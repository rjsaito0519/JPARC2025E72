# JPARC2025E72

J-PARC E72 の解析で使う較正・解析ツール群。analyzer（k18-analyzer の E72 ブランチ）で decode した ROOT を読み、較正パラメータを analyzer の `param/` に書き戻す。

使い方の詳細（較正の流れ: TDC window → Hodoscope → DC → D5）は [Wiki](https://github.com/rjsaito0519/JPARC2025E72/wiki) を参照。

## セットアップ

```bash
cp common/lib/paths_example.py common/lib/paths.py
./build.sh                             # C++ のツールを bin/ にビルド
```

Python 共通は `from common.lib import …`（トップに `lib` は置かない）。トップの symlink は人が回す入口だけ（`./calib.py` / `./ops.py` / `run_*.py` / `create_runlist.py` / `dst_create_runlist.py`）。自動化の日常起動は `./calib.py <run>`。状態操作は `./ops.py` または Slack。`iterate_*` / `adopt_*` / `extra_*` / 旧 session・dispatch は `calibration/scripts/…` 直叩き。

ローカル用 symlink（gitignore）: `group_root` → グループデータ領域、`param` → analyzer の `submodule/param`（コード上の正は `PARAM_DIR`）。

## 構成

| パス | 中身 |
|---|---|
| `calibration/scripts/` | 較正フロント本体（トップ symlink の先。`calib` / `ops` / `run_*` / `iterate_*` / 補助の `adopt_*`・`extra_*`・`post_*` など） |
| `calibration/config/` | 例: `calib_dispatch_example.yml` / `calib_session_example.yml` |
| `calibration/tdc_window_web/` | TDC window をブラウザで決めるツール（Wiki の TDC-Window） |
| `calibration/src/` | 較正の fit 本体（C++） |
| `analysis/` | 較正以外の解析・確認用（C++ / Python）。補助スクリプトは `analysis/scripts/<theme>/`（tpcpos, tpc, kpsc, hodo_kvc 等） |
| `common/` | 共通部品。C++ は `include/`・`src/`、Python は `common.lib`（`common/lib/`） |
| `macro/` | 表示・手作業用。`macro/evdisp/`（helix 3D 等）、`macro/calib/`（検査用 `.C`）。本番数値は C++/`bin` 優先。**ACLiC しない** |
| `archive/` | 使っていないソースの退避 |
| `wiki/` | 使い方（GitHub wiki）。手順の正本 |
| `bin/`, `results/` | ビルドした実行ファイル、出力（PDF・ROOT） |

## 較正の基本の形

```bash
# 自動化（推奨）: 起動は calib、判断は Slack / ops
./calib.py 2977
./calib.py 2977 --kaon --local
./ops.py status 2977                 # 別ターミナル / 別 cw から可

# 単体（従来どおり）
./create_runlist.py --hodo <run> Pi_hdprm     # conf・runlist を作る
# analyzer の runmanager/myrun.py で decode（既定 bsub、--local 可）
./run_hodo.py <run> hdprm                     # fit して param を更新
```

ドキュメント（物理的な背景・検討の経過）は analyzer 側の `docs/` に置いている。
