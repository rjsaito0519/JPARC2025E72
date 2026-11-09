# JPARC2025E72

J-PARC E72 の解析で使う較正・解析ツール群。analyzer（k18-analyzer の E72 ブランチ）で decode した ROOT を読み、較正パラメータを analyzer の `param/` に書き戻す。

使い方の詳細（較正の流れ: TDC window → Hodoscope → DC → D5）は [Wiki](https://github.com/rjsaito0519/JPARC2025E72/wiki) を参照。

## セットアップ

```bash
cp common/lib/paths_example.py common/lib/paths.py   # または cp lib/paths_example.py lib/paths.py（symlink）
./build.sh                             # C++ のツールを bin/ にビルド
./link.sh                              # 主なスクリプトを最上位にリンク
```

## 構成

| パス | 中身 |
|---|---|
| `calibration/scripts/` | 較正のフロントエンド（`create_runlist.py`、`run_*.py`、`iterate_*.py`、`calib_dispatch.py`、`update_param.py` など） |
| `calibration/config/` | 例: `calib_dispatch_example.yml` |
| `calibration/tdc_window_web/` | TDC window をブラウザで決めるツール（Wiki の TDC-Window） |
| `calibration/src/` | 較正の fit 本体（C++） |
| `analysis/` | 較正以外の解析・確認用（C++ / Python）。補助スクリプトは `analysis/scripts/<theme>/`（tpcpos, tpc, kpsc, hodo_kvc 等） |
| `common/` | 共通部品。C++ は `include/`・`src/`、Python は `lib/`（トップの `lib` は symlink） |
| `macro/` | 表示・手作業用。`macro/evdisp/`（helix 3D 等）、`macro/calib/`（検査用 `.C`）。本番数値は C++/`bin` 優先。**ACLiC しない** |
| `archive/` | 使っていないソースの退避 |
| `wiki/` | 使い方（GitHub wiki）。手順の正本 |
| `docs/` | analyzer `docs/` へのポインタのみ（文書は書かない） |
| `bin/`, `results/` | ビルドした実行ファイル、出力（PDF・ROOT） |

## 較正の基本の形

```bash
./create_runlist.py --hodo <run> Pi_hdprm     # conf・runlist を作る
# analyzer の runmanager/myrun.py で decode
./run_hodo.py <run> hdprm                     # fit して param を更新
```

ドキュメント（物理的な背景・検討の経過）は analyzer 側の `docs/` に置いている。
