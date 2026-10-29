# alphagenome-variant-scorer

AlphaGenomeを利用してSNV・indel・単純DEL・tandem DUPのRNA_SEQスコアを取得する、Windows/Linux向けCLIです。

準備・進捗の保存、取得・再開、スコアの集計・TSV/JSON出力に対応しています。
`run/resume` は未取得分のAPIアクセスを行います。他のコマンドはAPIキー不要・通信なしで動作します。

## インストール

Python 3.12以上が必要です。各OSで個別に仮想環境を作成してください。
正規化・API取得にはAlphaGenome SDK 0.9.0を使用します。
以下の手順はリポジトリのルートで実行してください。

Windows PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\alphagenome-variant-scorer.exe --help
```

Linux:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/alphagenome-variant-scorer --help
```

以下はCLIがPATHに入った環境での例です。仮想環境を有効化するか、上記の実行ファイルを指定してください。

## 初期設定

### FASTA

利用者が保有する非圧縮のGRCh38/hg38 FASTAを指定します。参照配列は同梱・自動ダウンロードしません。
REF照合、indelの正規化、DEL/tandem DUPの配列構築に使用し、API側の参照ゲノムは変更しません。

```bash
alphagenome-variant-scorer config set reference.fasta "/data/reference/hg38.fa"
alphagenome-variant-scorer reference check --create-index
```

Windowsでは例えば `"D:/reference/hg38.fa"` と指定できます。設定は次回以降も使用されます。
既存の `.fai` を使用し、なければ作成します。既存索引は上書きしません。
`reference check` は索引と配列の部分読出しを確認します。GRCh38やAPI側参照との一致は検証しません。

### APIキー

APIキーは[公式サイト](https://alphagenome.google/api)で取得します。
次のコマンドの `config_file` に表示された `config.toml` を開き、`[auth]` の `api_key` を追加・更新してください。

```bash
alphagenome-variant-scorer config show
```

```toml
[auth]
api_key = "ここにAPIキーを記入"
```

キーはconfigからのみ読み込みます。キー入りconfigを共有・Git登録しないでください。
Linuxでは `chmod 600` でconfigを所有者のみ読み書きできるようにします。

### 共通キャッシュ

取得済みスコアを解析間で再利用します。未作成の専用フォルダを指定して初期化してください。

```bash
alphagenome-variant-scorer config set cache.directory "/data/cache/alphagenome"
alphagenome-variant-scorer cache init
alphagenome-variant-scorer cache status
```

初期化済みなら状態を確認し、内容を保持します。
同一ホストのローカルディスクを使用してください。ネットワーク共有での運用は対象外です。
スコア本体は共通キャッシュに保存します。再開・再出力には実行フォルダと対応する共通キャッシュの両方を保管してください。

## TSVを準備する

入力はUTF-8のタブ区切りTSVです。prepareでREF照合・正規化とAPI送信の準備を行います。
未作成の実行フォルダを指定してください。失敗時に作成されたフォルダも再利用できません。

### REF/ALTを指定する形式（literal-tsv）

必須列は `#CHROM`, `POS`, `REF`, `ALT`。POSは1-basedです。
REF/ALTはアンカー付きの非空ACGT文字列とし、小文字は準備時に大文字化します。

```bash
alphagenome-variant-scorer prepare --input variants.tsv --input-format literal-tsv --run-dir /data/runs/example
alphagenome-variant-scorer status --run-dir /data/runs/example
```

SNVは既定でAtlas（事前計算スコア）、indel/MNV/複合置換はPrediction（予測API）を使用します。
SNVもPredictionで評価する場合は `--snv-source prediction` を指定します。MNVをSNVへ分割しません。

### 座標とSV種別を指定する形式（sv-tsv）

必須列は `SV_chrom`, `SV_start`, `SV_end`, `SV_type`, `SV_length`。
start/endは1-based closed、長さは `end - start + 1`（SV_lengthの絶対値）と一致させます。
単純DELと、明示的に仮定した同方向・追加1コピーのtandem DUPを準備できます。

```bash
alphagenome-variant-scorer prepare --input structural.tsv --input-format sv-tsv --dup-model tandem-extra-copy --run-dir /data/runs/structural
```

`--dup-model` がなければDUPは `dup_structure_unspecified` として保存し、要求を作りません。

### 入力条件と準備結果

対応するcontig表記は `chr1`〜`chr22`, `chrX`, `chrY` に限定し、別名を自動変換しません。
予測windowは1,048,576 bpです。windowより短いcontig、入力長制限、染色体端でSDKの正規化に必要な配列が得られない場合は未対応として保存します。
DUPの長さは `2 * 長さ + 1 < 1,048,576` を満たす必要があります。

各行はAPI送信準備済み `ready`、不正 `invalid`、未対応 `unsupported` の状態で保存します。
必須列欠落等は処理全体のエラーです。prepare完了時は不正・未対応行があっても終了コード0なので、状態別の件数と理由を確認してください。

元の列値・追加列・行順を保持し、実行フォルダには次を保存します。

- `input.private.tsv`: 明示的に指定した入力のコピー
- `run.private.sqlite`: 元行、行別状態、重複除去したジョブと準備済み要求
- `manifest.private.json`: 入力・参照のハッシュ、使用版、構造仮定と解析条件

prepareごとにFASTA全体のSHA-256を計算します。

## 取得と再開

prepare済みの実行フォルダと初期化済みの共通キャッシュを使用します。

```bash
alphagenome-variant-scorer run --run-dir /data/runs/example
alphagenome-variant-scorer status --run-dir /data/runs/example
```

中断した処理はresumeで再開します。

```bash
alphagenome-variant-scorer resume --run-dir /data/runs/example
```

run/resumeの主なオプション:

| オプション | 既定値 | 内容 |
| --- | --- | --- |
| `--max-jobs` | 残りすべて | 今回処理するジョブ数の上限。キャッシュ利用も数える |
| `--workers` | 1 | 同時処理数（1〜128） |
| `--min-interval` | 1秒 | API取得の試行間隔 |
| `--timeout` | 120秒 | 接続・RPCの制限時間 |
| `--max-attempts` | 3 | ジョブ当たりの試行上限 |

一時的なタイムアウト・接続障害のみ再試行します。認証・クォータ等の失敗は停止し、原因解消後にresumeします。
Atlas未収載（`not_found`）や要求したscorerが返らない場合（`not_returned`）は状態を保存して他のジョブを続け、同じ実行のresumeでは再照会しません。
送信間隔は全worker・Atlas/Prediction・再試行に共通で、キャッシュ利用には適用しません。複数実行では送信頻度の合計に注意してください。

異常やCtrl+Cでは新規送信を止め、通信中の処理の完了またはタイムアウトを待ち、取得できた応答を保存して終了します。
再開には最初と同じ共通キャッシュを使います。FASTAは不要で、並列数・送信間隔は変更できます。

## 保存結果の出力

```bash
alphagenome-variant-scorer export --run-dir /data/runs/example --output results.tsv --threshold 0.5
alphagenome-variant-scorer export --run-dir /data/runs/example --output results.json --format json --threshold 0.5
# gene×trackごとのraw値をすべて出す場合（大きなファイルになります）
alphagenome-variant-scorer export --run-dir /data/runs/example --output scores.tsv --format scores-tsv
```

exportは保存済みのrawスコアから再集計します。閾値を変更してもAPIの再取得は不要です。
TSVには元列に `AGVS_` 列を追加し、JSONには元入力値・行番号・取得状態・gene別要約を保存します。
不正行・重複行も順番を保ちます。列数が不正な入力行の全値は `AGVS_original_values_json` で保持します。
TSVは元の文字列を保持するため、Excel等では文字列として取り込んでください。
`AGVS_` の出力列名と元列名が衝突する場合はJSONを使ってください。
詳細TSVには集計対象外も含む全raw値と、gene/track情報・非有限値の種類を出力します。

### 集計と評価状態

Predictionでは `GeneMaskLFCScorer` を使用します。
遺伝子ごとに、返却された全trackのraw値が有限値（NaN/Infを含まない）のときだけ集計対象とします。
`|raw| >= threshold` を丸め前の値で判定し、正・負両方向の該当を保持します。既定閾値は0.5です。
これは保存スコアの抽出条件であり、病原性判定ではありません。

| 評価状態 / コード | 意味 |
| --- | --- |
| complete | 返却された全遺伝子の全trackが集計可能 |
| partial / P01 | 一部遺伝子が非有限値を含み、集計対象外。残りの有効遺伝子を集計 |
| unavailable / E02 | 返却遺伝子を1つも集計できない |
| no_genes / E02 | 遺伝子が返っていない |
| no_tracks / E02 | trackが返っていない |

TSVでは評価状態を `AGVS_evaluation`、P01/E02を `AGVS_code`、閾値判定を `AGVS_threshold_result` に記録します。
「該当なし」は集計できた範囲の判定で、未算出部分は判定対象外です。

未完了の解析は `--allow-partial` で出力でき、未取得・失敗行も要約TSV/JSONに含めます。
各出力の `<出力名>.manifest.json` に解析条件・ハッシュ・途中結果かどうかを記録します。
既存出力やmanifestは上書きしないため、新しい出力名を指定してください。
出力先はハードリンクに対応したローカルファイルシステムを使用してください。

## 詳細設定

### 設定ファイル

configの既定の保存先:

- Windows: `%LOCALAPPDATA%/alphagenome-variant-scorer/config.toml`
- Linux: `$XDG_CONFIG_HOME/alphagenome-variant-scorer/config.toml`
  （未設定時は `~/.config/alphagenome-variant-scorer/config.toml`）

設定例は [examples/config.example.toml](examples/config.example.toml) を参照してください。
別のconfigを使う場合は、コマンドの前に `--config FILE` を置きます。

```bash
alphagenome-variant-scorer --config /path/to/settings.toml config set reference.fasta /data/reference/hg38.fa
```

指定したconfigとユーザーconfigはマージしません。手編集したconfig内の相対パスはconfigの親ディレクトリを基準に解決します。
`config set` は他の設定値を保持しますが、TOMLのコメント・書式は保持しません。
キャッシュは `cache init/status --cache-dir DIR` でconfigを使わず指定することもできます。

### 参照配列と索引

FASTAの配置先が読取専用なら、書込み可能な別の索引パスを設定します。

```bash
alphagenome-variant-scorer config set reference.index "/data/indices/hg38.fa.fai"
alphagenome-variant-scorer reference check --create-index
```

古い・不正な既存索引がある場合も、新しい索引パスを指定してください。
索引保存先はハードリンク対応のローカルファイルシステムを使用します。

一時的に別のFASTAを使う場合は、configを変更せず指定できます。以前のFASTA用索引は引き継ぎません。

```bash
alphagenome-variant-scorer reference check --reference-fasta /data/other/reference.fa --create-index
```

`config set reference.fasta` で保存するFASTAを変更すると、以前の `reference.index` 設定を解除します。FASTA・索引ファイル自体は保持します。
FASTA全体のSHA-256を確認する場合は `reference check --sha256` を使用します。

## 開発用の確認

```bash
python -m unittest discover -s tests -v
python tools/check_public_files.py
python tools/check_public_files.py --staged
```

テストは合成データを使用し、実データ・既存キャッシュ・APIキーは不要です。

## ライセンス

本プログラムは [Apache License 2.0](LICENSE) で提供します。
AlphaGenome APIの利用には、提供元の利用条件が別途適用されます。
