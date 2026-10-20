"""Pinned SDK adapter. Secrets and raw exception messages are never persisted."""

from collections import namedtuple
import importlib.metadata
from threading import Lock

import grpc

from .errors import UserError
from .scores import encode_scores


class FetchError(Exception):
    def __init__(self, code):
        self.code = code
        super().__init__(code)


def error_code(exc):
    cause = exc.__cause__
    if isinstance(cause, grpc.RpcError):
        exc = cause
    if isinstance(exc, grpc.RpcError):
        return {
            grpc.StatusCode.UNAUTHENTICATED: "authentication",
            grpc.StatusCode.PERMISSION_DENIED: "authentication",
            grpc.StatusCode.RESOURCE_EXHAUSTED: "rate_limit",
            grpc.StatusCode.DEADLINE_EXCEEDED: "timeout",
            grpc.StatusCode.UNAVAILABLE: "unavailable",
            grpc.StatusCode.NOT_FOUND: "not_found",
            grpc.StatusCode.INVALID_ARGUMENT: "invalid_request",
            grpc.StatusCode.OUT_OF_RANGE: "invalid_request",
        }.get(exc.code(), "api_error")
    if isinstance(exc, (TimeoutError, grpc.FutureTimeoutError)):
        return "timeout"
    return "api_error"


_Details = namedtuple("_Details", "method timeout metadata credentials wait_for_ready compression")


class Deadline(grpc.UnaryUnaryClientInterceptor, grpc.StreamStreamClientInterceptor):
    def __init__(self, seconds):
        self.seconds = seconds

    def details(self, details):
        return _Details(details.method, min(details.timeout, self.seconds)
                        if details.timeout is not None else self.seconds,
                        details.metadata, details.credentials, details.wait_for_ready,
                        getattr(details, "compression", None))

    def intercept_unary_unary(self, continuation, details, request):
        return continuation(self.details(details), request)

    def intercept_stream_stream(self, continuation, details, iterator):
        return continuation(self.details(details), iterator)


def validate_request(request):
    common = {"api", "assembly", "organism", "variant", "requested_scorer"}
    extra = {"interval", "model_version", "scorer", "merge_stranded_gene_tracks"}
    route = request.get("api")
    if route not in {"Atlas", "Prediction"} or set(request) != (
        common if route == "Atlas" else common | extra
    ):
        raise UserError("Unsupported prepared API request.")
    if (request["assembly"] != "GRCh38" or request["organism"] != "HOMO_SAPIENS"
            or request["requested_scorer"] != "RNA_SEQ"):
        raise UserError("Unsupported assembly, organism or scorer.")
    if route == "Prediction" and (
        request["model_version"] != "ALL_FOLDS"
        or request["scorer"] != {"type": "GeneMaskLFCScorer", "requested_output": "RNA_SEQ"}
        or request["merge_stranded_gene_tracks"] is not True
    ):
        raise UserError("Unsupported prediction conditions.")


class ScoreClient:
    def __init__(self, *, key_provider=None, timeout=120):
        self.key_provider = key_provider
        self.timeout = timeout
        self.channel = None
        self.atlas = self.model = None
        self._initialization_lock = Lock()

    def close(self):
        if self.channel is not None:
            self.channel.close()
            self.channel = None
            self.atlas = self.model = None

    def _initialize(self):
        with self._initialization_lock:
            self._initialize_locked()

    def _initialize_locked(self):
        if self.channel is not None:
            return
        from alphagenome.atlas import atlas
        from alphagenome.models import dna_client
        from alphagenome.protos import atlas_service_pb2_grpc

        try:
            if self.key_provider is None:
                raise UserError("No config credential provider supplied.")
            key = self.key_provider()
        except UserError:
            raise FetchError("credential_unavailable") from None
        channel = grpc.secure_channel(
            "dns:///gdmscience.googleapis.com:443", grpc.ssl_channel_credentials(),
            options=[("grpc.max_receive_message_length", -1), ("grpc.max_send_message_length", -1),
                     ("grpc.enable_retries", 0)])
        try:
            grpc.channel_ready_future(channel).result(timeout=self.timeout)
        except BaseException:
            channel.close()
            raise
        self.channel = channel
        limited = grpc.intercept_channel(channel, Deadline(self.timeout))
        metadata = [("x-goog-api-key", key)]
        self.atlas = atlas.AtlasClient(
            atlas_service_pb2_grpc.AtlasServiceStub(limited), metadata=metadata)
        self.model = dna_client.DnaClient(
            channel=limited, model_version=dna_client.ModelVersion.ALL_FOLDS, metadata=metadata)

    def __call__(self, request):
        validate_request(request)
        if importlib.metadata.version("alphagenome") != "0.9.0":
            raise UserError("API adapter requires SDK 0.9.0.")
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers

        v = request["variant"]
        variant = genome.Variant(v["chrom"], v["pos"], v["ref"], v["alt"])
        try:
            self._initialize()
            if request["api"] == "Atlas":
                result = self.atlas.query_variant(variant, requested_scorers=["RNA_SEQ"])
                if "RNA_SEQ" not in result:
                    raise FetchError("not_returned")
                adata = result["RNA_SEQ"]
            else:
                interval = request["interval"]
                # Pinned SDK wrapper retries internally; the runner alone owns retries.
                adata, = dna_client.DnaClient.score_variant.__wrapped__(
                    self.model, genome.Interval(v["chrom"], interval["start0"], interval["end0"]),
                    variant, [variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]],
                    organism=dna_client.Organism.HOMO_SAPIENS,
                    merge_stranded_gene_tracks=True)
            return encode_scores(adata), {
                "sdk_version": "0.9.0", "api": request["api"], "scorer": "RNA_SEQ",
                "model_revision": None,
            }
        except FetchError:
            raise
        except Exception as exc:
            raise FetchError(error_code(exc)) from None
