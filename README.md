# alphagenome-variant-scorer

AlphaGenomeを利用してSNV・indel・単純DEL・tandem DUPのスコアを取得する、Windows/Linux向けCLIです。
取得対象はRNA_SEQスコアです。Predictionでは `GeneMaskLFCScorer` を使用します。

準備結果と進捗を実行単位のSQLiteへ保存し、スコア応答は解析間で共有します。保存済みスコアの集計・TSV/JSON出力に対応しています。
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

## FASTAを一度設定する

利用者が保有する、非圧縮のGRCh38/hg38 FASTAを指定します。
参照配列は同梱・自動ダウンロードしません。
ローカルFASTAはREF照合、indelの正規化、DEL/tandem DUPの配列構築に使用します。
FASTAを指定してもAlphaGenome API側の参照ゲノムは変更されません。

```bash
alphagenome-variant-scorer config set reference.fasta "/data/reference/hg38.fa"
alphagenome-variant-scorer config show
```

Windowsでは例えば `"D:/reference/hg38.fa"` と指定できます。
設定は保存され、次回以降も使用されます。索引の確認・作成は次節で行います。

保存先:

- Windows: `%LOCALAPPDATA%/alphagenome-variant-scorer/config.toml`
- Linux: `$XDG_CONFIG_HOME/alphagenome-variant-scorer/config.toml`
  （未設定時は `~/.config/alphagenome-variant-scorer/config.toml`）

別のconfigを使う場合は、コマンドの前に `--config FILE` を置きます。
指定したファイルとユーザーconfigを暗黙にマージしません。

```bash
alphagenome-variant-scorer --config /path/to/settings.toml config set reference.fasta /data/reference/hg38.fa
```

手編集したconfig内の相対パスは、configの親ディレクトリを基準に解決します。
設定例は [examples/config.example.toml](examples/config.example.toml) を参照してください。
`config set` は他の対応済み設定値を保持しますが、TOMLのコメント・書式は保持しません。

## 参照配列を確認する

```bash
alphagenome-variant-scorer reference check
```

既存の `.fai` を使用します。索引がなければ、明示的に作成できます。

```bash
alphagenome-variant-scorer reference check --create-index
```

FASTAの配置先が読取専用なら、書込み可能な別の索引パスを設定します。

```bash
alphagenome-variant-scorer config set reference.index "/data/indices/hg38.fa.fai"
alphagenome-variant-scorer reference check --create-index
```

既存索引は上書きしません。古い・不正な索引がある場合は、新しい索引パスを指定してください。
索引保存先はハードリンク対応のローカルファイルシステムを使用してください。

一時的なFASTA変更は次のように指定できます。configは変更されず、以前のFASTA用索引も引き継ぎません。

```bash
alphagenome-variant-scorer reference check --reference-fasta /data/other/reference.fa --create-index
```

保存するFASTAを `config set reference.fasta` で変更した場合は、以前の `reference.index` 設定を解除します。
以前のFASTAや索引ファイル自体は保持します。

FASTA全体のSHA-256も確認する場合は `reference check --sha256` を使用します。

`reference check` は索引構造と各contigの先頭・末尾の部分読出しを確認します。
全塩基の妥当性、GRCh38であること、API側参照との完全一致を検証した結果ではありません。
参照と索引を対応づけて管理してください。

## TSVを準備する

FASTAの設定と索引確認を済ませてから、新しい実行フォルダを指定します。
以下はAPIへ送信せず、REF照合、正規化、重複要求の共有、入力行とジョブの対応付けを行います。

```bash
alphagenome-variant-scorer prepare --input variants.tsv --input-format literal-tsv --run-dir /data/runs/example
alphagenome-variant-scorer status --run-dir /data/runs/example
```

literal TSVの必須列は `#CHROM`, `POS`, `REF`, `ALT`。POSは1-basedです。
REF/ALTはアンカー付きの非空ACGT文字列とし、小文字は準備時に大文字化します。
元の列値・追加列・行順は保持します。
SNVは既定でAtlas、indel/MNV/複合置換はPrediction用のリクエストにします。
SNVもPredictionで評価する計画なら `--snv-source prediction` を指定します。
MNVをSNVへ分割しません。

SV TSVの必須列は `SV_chrom`, `SV_start`, `SV_end`, `SV_type`, `SV_length`。
start/endは1-based closed、長さは `end - start + 1`（SV_lengthの絶対値）と一致させます。
単純DELと、明示的に仮定した同方向・追加1コピーのtandem DUPを準備できます。

```bash
alphagenome-variant-scorer prepare --input structural.tsv --input-format sv-tsv --dup-model tandem-extra-copy --run-dir /data/runs/structural
```

`--dup-model` がなければDUPは `dup_structure_unspecified` として保存し、要求を作りません。
対応するcontig表記は `chr1`〜`chr22`, `chrX`, `chrY` に限定し、別名を自動変換しません。
予測windowは2^20 bpです。windowより短いcontig、入力長制限、染色体端でSDKの正規化に必要な配列が得られない場合は未対応として保存します。
DUPの長さは `2 * 長さ + 1 < 2^20` を満たす必要があります。

REF不一致・不正座標等は `invalid`、未対応形式・構造未指定等は `unsupported`、
準備できた行は `ready` としてDBに残します。`ready` はAPI取得成功を意味しません。
必須列欠落等は処理全体のエラーですが、行ごとの不正・未対応は他の行の準備を妨げません。
検証処理が完了した場合の終了コードは0なので、ready件数と理由の内訳も確認してください。

実行フォルダには次を保存します。これらは非公開の解析データとして扱ってください。

- `input.private.tsv`: 明示的に指定した入力のコピー
- `run.private.sqlite`: 元行、行別状態、重複除去したジョブと準備済み要求
- `manifest.private.json`: 入力・FASTA・索引のハッシュ、使用版、構造仮定、条件と件数

prepareごとにFASTA全体のSHA-256を計算します。
既存の実行フォルダは上書きしません。prepareに失敗した場合も、修正後は新しいフォルダを指定してください。

## 共通キャッシュ

取得済みスコアを解析間で共有し、同じ要求のAPI再送を避けます。
まだ存在しない専用ディレクトリを指定して初期化してください。

```bash
alphagenome-variant-scorer config set cache.directory /path/to/shared-cache
alphagenome-variant-scorer cache init
alphagenome-variant-scorer cache status
```

`cache init/status --cache-dir DIR` で設定を使わず指定することもできます。
初期化済みのキャッシュに対するinitは状態を確認し、内容をリセットしません。
同一ホストのローカルディスクを使用してください。ネットワーク共有での運用は対象外です。

## 取得と再開

prepare済みの実行フォルダと初期化済みの共通キャッシュを使用します。
APIキーは非公開のconfig.tomlの `[auth]` に設定します。
`config show` に表示されるconfig_fileを開き、次を追加してください（既にauthがあればその項目を編集）。
Windowsのパッケージ版アプリでは設定フォルダが転送される場合があるため、表示された実体のパスを使用してください。

```toml
[auth]
api_key = "ここにAPIキーを記入"
```

環境変数・envファイルからのキー読込みには対応しません。`config show` ではキーを伏せます。
キー入りconfigはGitへ追加しないでください。
configは既定でリポジトリ外に保存します。Linuxでは所有者のみ読み書きできるよう `chmod 600` を設定してください。

```bash
alphagenome-variant-scorer run --run-dir /path/to/run --max-jobs 10
# 並列数を指定する場合
alphagenome-variant-scorer run --run-dir /path/to/run --workers 4
alphagenome-variant-scorer status --run-dir /path/to/run
alphagenome-variant-scorer resume --run-dir /path/to/run
```

`--max-jobs` は今回処理するジョブ数の上限で、キャッシュ利用も数えます。省略すると残りすべてを処理します。
`--workers` は同時処理数（1〜128、既定1）です。128はプログラムが受け付ける上限であり、推奨値ではありません。
`--min-interval` は実行内の全workerを合わせた取得試行間隔（既定1秒）、
`--timeout` は接続・RPCの制限時間（既定120秒）、`--max-attempts` はジョブ当たりの試行上限（既定3回）です。
一時的なタイムアウト・接続障害のみ再試行します。認証・クォータ等の失敗は停止し、原因解消後にresumeします。
Atlas未収載（not_found）や要求したscorerが返らない場合（not_returned）は行の取得状態として保存し、
他のジョブを続けます。同じ実行のresumeでは再送しません。これらの負の結果は共通キャッシュへ保存せず、別解析では再照会します。
送信間隔はAtlas/Predictionと再試行に共通で、キャッシュ利用には適用しません。
複数実行の同時起動では、実行ごとに管理される送信頻度の合計に注意してください。
異常やCtrl+Cでは新規送信を止め、通信中の応答を保存して終了します。
通信の完了またはタイムアウトを待つ場合があります。中断後は保存済み応答を再利用して再開します。
runとresumeは同じ再開可能な処理です。再開には最初と同じ共通キャッシュ保存先を使用します。
再開時にworkers・送信間隔を変更できます。FASTAは再開時に不要です。
`status` のfetch_attemptsはアダプターの試行数です。接続前失敗も含み、サーバーが受信したRPC数とは異なります。

## 保存結果の出力

```bash
alphagenome-variant-scorer export --run-dir /path/to/run --output results.tsv --threshold 0.5
alphagenome-variant-scorer export --run-dir /path/to/run --output results.json --format json --threshold 0.5
# gene×trackごとのraw値をすべて出す場合（大きなファイルになります）
alphagenome-variant-scorer export --run-dir /path/to/run --output scores.tsv --format scores-tsv
```

exportは保存結果から再集計するため、閾値を変えてもAPIへの再送は不要です。
TSVには元列に `AGVS_` 列を追加し、JSONには元入力値・行番号・取得状態・gene別要約を保存します。
不正行・重複行も順番を保ちます。列数が不正な入力行の全値は `AGVS_original_values_json` で保持します。
TSVは元の文字列を保持します。第三者由来のデータをExcel等で開く場合は、数式として解釈されないよう文字列として取り込んでください。
`AGVS_` の出力列名と元列名が衝突する場合はJSONを使ってください。
詳細TSVは取得済み応答のrawセルのみを出し、gene/trackメタデータと非有限値の種類を付けます。
quantile等のlayerはキャッシュに保持しますが、exportの判定・詳細値はrawのみです。

遺伝子ごとに、返却された全trackが有限値のときだけ集計対象とします。
`|raw| >= threshold` を丸め前の値で判定し、正・負両方向の該当を保持します。既定閾値は0.5です。
これは保存スコアの抽出条件であり、病原性判定ではありません。

| 評価状態 | 意味 |
| --- | --- |
| complete | 返却された全遺伝子の全trackが集計可能 |
| partial / P01 | 一部遺伝子が非有限値を含み、集計対象外。残りの有効遺伝子を集計 |
| unavailable / E02 | 返却遺伝子を1つも集計できない |
| no_genes / E02 | 遺伝子が返っていない |
| no_tracks / E02 | trackが返っていない |

P01/E02は本プログラムの評価状態で、APIの通信エラーではありません。P01でも閾値該当の有無を別項目で残します。
「該当なし」は集計できた範囲内の結果で、未算出部分の陰性を意味しません。
一部trackのみ非有限の遺伝子も全体を除外しますが、元の有限値を含む全raw値はキャッシュと詳細TSVに残します。
未算出の原因はこの判定では特定できません。completeも返却集合内での完全性を意味します。

未完了の解析は `--allow-partial` を指定すると出力でき、未取得・失敗行も要約TSV/JSONに含めます。各出力に `<出力名>.manifest.json` を付け、
閾値・集計方式の版・入力/要求ハッシュ・出力SHA-256・途中結果かどうかを記録します。
新しい出力名を使ってください。既存出力やmanifestを上書きしません。
出力先はハードリンクに対応したローカルファイルシステムを使用してください。

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
