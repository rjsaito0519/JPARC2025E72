# archive

使っていない・置き場を移したファイルの退避。現行の本番コードはここを import / ビルドしない。

| パス | 由来 |
|---|---|
| `unused_cpp/` | 旧 `src/unused/` と `analysis/src/unused/` |
| `param_backup/` | 旧トップ `scripts/param_backup/` |
| `scripts_debug/` | 試験用スクリプト（例: `tpcpos_write_test_maps.py`） |
| `macro_legacy/` | 非 helix の旧 TPC viewer（`check_tpc_track*`） |
| `macro_manuals/` | evdisp 用 manual PDF/tex/fig |
| `scripts_compare_bt/` | 旧 beam-through 比較 PDF（`compare_*`；run 表は `common/lib/tpcpos_runs.py`） |

生きている補助スクリプトは `analysis/scripts/`（例: `lh2eloss/`、`beam_summary/`）。較正フロントは `calibration/scripts/`。表示は `macro/{evdisp,calib}/`。
