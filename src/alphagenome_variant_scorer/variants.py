"""Prepare literal API requests using the pinned SDK normalizer, without RPCs."""

from dataclasses import dataclass
import re

from alphagenome.data import genome

from .errors import UserError

WINDOW = 2**20
CONTIGS = frozenset([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
PREPARATION_VERSION = 1


@dataclass
class RowIssue(Exception):
    state: str
    reason: str


def invalid(reason):
    raise RowIssue("invalid", reason)


def unsupported(reason):
    raise RowIssue("unsupported", reason)


def integer(value):
    if len(value) > 20 or not re.fullmatch(r"[+-]?[0-9]+", value):
        invalid("invalid_integer")
    return int(value)


def classify(ref, alt):
    if ref == alt:
        invalid("identical_alleles")
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    # Ordinary anchored indels retain one entire allele at one end.
    if alt.startswith(ref) or alt.endswith(ref):
        return "INS"
    if ref.startswith(alt) or ref.endswith(alt):
        return "DEL"
    return "COMPLEX_REPLACEMENT"


def variant_dict(v):
    return {
        "chrom": v.chromosome, "pos": v.position,
        "ref": v.reference_bases, "alt": v.alternate_bases,
    }


def verified_normalize(v, reference, contig_length):
    observed = reference.extract(v.reference_interval)
    if set(observed) - set("ACGT"):
        unsupported("ambiguous_reference")
    if observed != v.reference_bases:
        invalid("ref_mismatch")
    if v.is_snv:
        return v
    if v.start + max(len(v.reference_bases), len(v.alternate_bases)) > contig_length:
        unsupported("normalization_context_out_of_bounds")
    try:
        normalized = genome.normalize_variant(v, reference)
    except UserError:
        # SDK 0.9.0 can request negative coordinates while left shifting.
        unsupported("normalization_reference_bounds_or_bases")
    if reference.extract(normalized.reference_interval) != normalized.reference_bases:
        raise UserError("Normalizer REF verification failed; preparation stopped.")
    if set(normalized.reference_bases + normalized.alternate_bases) - set("ACGT"):
        unsupported("ambiguous_normalized_allele")
    start = max(0, min(v.start, normalized.start) - 32)
    end = min(contig_length, max(v.end, normalized.end) + 32)
    sequence = reference.sequence(v.chromosome, start, end)

    def apply(variant):
        return sequence[:variant.start - start] + variant.alternate_bases + sequence[variant.end - start:]

    if apply(v) != apply(normalized):
        raise UserError("Normalized variant changes the resulting sequence; preparation stopped.")
    return normalized


def prepare_variant(row, input_format, reference, *, dup_model, snv_source):
    if input_format == "literal-tsv":
        chrom = row["#CHROM"]
    else:
        chrom = row["SV_chrom"]
    # Deliberately no chr-prefix or MT alias inference.
    if chrom not in reference.contigs:
        invalid("unknown_reference_contig")
    if chrom not in CONTIGS:
        unsupported("contig_not_in_initial_api_scope")
    length = reference.contigs[chrom]
    focus = None
    source_bounds = None
    assumption = None
    if input_format == "literal-tsv":
        pos = integer(row["POS"])
        ref, alt = row["REF"].upper(), row["ALT"].upper()
        if "," in alt:
            unsupported("multiallelic_input")
        if "<" in alt or ">" in alt or "[" in alt or "]" in alt or alt == "*":
            unsupported("symbolic_allele")
        if not ref or not alt or set(ref + alt) - set("ACGT"):
            invalid("invalid_literal_allele")
        kind = classify(ref, alt)
        if not 1 <= pos <= length or pos - 1 + len(ref) > length:
            invalid("variant_out_of_bounds")
        if max(len(ref), len(alt)) >= WINDOW:
            unsupported("allele_exceeds_initial_size_limit")
        original = genome.Variant(chrom, pos, ref, alt)
    else:
        kind = "SV_" + row["SV_type"]
        if kind not in {"SV_DEL", "SV_DUP"}:
            unsupported("sv_type_not_supported")
        start, end, declared = (integer(row[k]) for k in ("SV_start", "SV_end", "SV_length"))
        if not 1 <= start <= end <= length:
            invalid("variant_out_of_bounds")
        size = end - start + 1
        if abs(declared) != size:
            invalid("sv_length_mismatch")
        if start == 1:
            unsupported("left_anchor_unavailable")
        if size + 1 >= WINDOW or (kind == "SV_DUP" and 2 * size + 1 >= WINDOW):
            unsupported("sv_exceeds_initial_size_limit")
        if kind == "SV_DUP" and dup_model != "tandem-extra-copy":
            unsupported("dup_structure_unspecified")
        seq = reference.sequence(chrom, start - 2, end)
        if set(seq) - set("ACGT"):
            unsupported("ambiguous_reference")
        if kind == "SV_DEL":
            original = genome.Variant(chrom, start - 1, seq, seq[0])
            focus = genome.Interval(chrom, start - 2, end)
            assumption = "simple_deletion"
        else:
            original = genome.Variant(chrom, end, seq[-1], seq[-1] + seq[1:])
            focus = genome.Interval(chrom, start - 1, end + size)
            assumption = "same_orientation_tandem_duplication_one_extra_copy"
        source_bounds = (start - 2, end)
    normalized = verified_normalize(original, reference, length)
    api = "Atlas" if kind == "SNV" and snv_source == "atlas" else "Prediction"
    request = {
        "api": api, "assembly": "GRCh38", "organism": "HOMO_SAPIENS",
        "variant": variant_dict(normalized), "requested_scorer": "RNA_SEQ",
    }
    if api == "Prediction":
        if length < WINDOW:
            unsupported("contig_shorter_than_prediction_window")
        proposed = (focus if focus is not None else normalized.reference_interval).resize(WINDOW)
        start = max(0, min(proposed.start, length - WINDOW))
        end = start + WINDOW
        if not start <= normalized.start <= normalized.end <= end:
            unsupported("normalized_variant_outside_prediction_window")
        if source_bounds and not (start <= source_bounds[0] and source_bounds[1] <= end):
            unsupported("sv_outside_prediction_window")
        request.update(
            interval={"start0": start, "end0": end},
            model_version="ALL_FOLDS",
            scorer={"type": "GeneMaskLFCScorer", "requested_output": "RNA_SEQ"},
            merge_stranded_gene_tracks=True,
        )
    return kind, request, {
        "original_variant": variant_dict(original),
        "normalized_variant": variant_dict(normalized),
        "reference_verified": True,
        "haplotype_equality_verified": True,
        "structure_assumption": assumption,
        "coordinate_convention": "1-based-closed" if input_format == "sv-tsv" else "1-based-POS",
    }
