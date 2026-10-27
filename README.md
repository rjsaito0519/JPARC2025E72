# JPARC2025E72

J-PARC E72 の解析で使う較正・解析ツール群。analyzer（k18-analyzer の E72 ブランチ）で decode した ROOT を読み、較正パラメータを analyzer の `param/` に書き戻す。

使い方の詳細（較正の流れ: TDC window → Hodoscope → DC → D5）は [Wiki](https://github.com/rjsaito0519/JPARC2025E72/wiki) を参照。

## セットアップ

```bash
cp lib/paths_example.py lib/paths.py   # analyzer・param・データ置き場のパスを設定
./build.sh                             # C++ のツールを bin/ にビルド
./link.sh                              # 主なスクリプトを最上位にリンク
```

## 構成

| パス | 中身 |
|---|---|
| `calibration/scripts/` | 較正のフロントエンド（`create_runlist.py`、`set_tdc_window.py`、`run_hodo.py`、`run_dcprm.py`、`run_d5.py`、`run_tpc.py`、`update_param.py` など） |
| `calibration/src/` | 較正の fit 本体（C++） |
| `analysis/` | 較正以外の解析・確認用（C++ / Python） |
| `common/` | 共通部品。`common/include/params.h` に fit 範囲のヒントや run 固有の設定 |
| `lib/` | Python の共通部品（パス、param ファイルの読み書き、入力チェック） |
| `macro/` | ROOT マクロ・表示用スクリプト |
| `bin/`, `results/` | ビルドした実行ファイル、出力（PDF・ROOT） |

## 較正の基本の形

```bash
./create_runlist.py --hodo <run> Pi_hdprm     # conf・runlist を作る
# analyzer の runmanager/myrun.py で decode
./run_hodo.py <run> hdprm                     # fit して param を更新
```

ドキュメント（物理的な背景・検討の経過）は analyzer 側の `docs/` に置いている。
