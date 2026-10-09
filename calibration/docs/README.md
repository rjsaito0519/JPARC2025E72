# calibration/docs — 較正関連のドキュメント

myanalysis のツール（`HTOF_Calib`、`run_hodo.py`、較正スクリプト、`tpc_phase_plot` など）と結びついた較正の文書だけを置く。解析全般の文書は analyzer 側 `/gpfs/home/had/sryuta/analyzer/JPARC2025E72/docs/` に集約している。テーマごとにディレクトリを分けている。

| ディレクトリ / ファイル | 内容 |
|---|---|
| [htof/](htof/) | HTOF の較正（TDC / ADC / PHC / 絶対 TOF オフセット、`run_hodo.py ... htof`）と、HTOF・散乱粒子側 TOF 関連コードのレビュー |
| [dc/](dc/) | DC（BLC 等）のグローバルオフセット調整の記録（DCGEO の意味の説明は analyzer の `docs/geometry/dcgeo_tracking_geometry.md`） |
| [tpc/](tpc/) | TPC の位相プロットの仕様 |
| [conf_templates_restore.md](conf_templates_restore.md) | runlist 用 conf テンプレートの探索順・既存 conf 保護・TPCPOS 既定 |

- HTOF: [htof/htof_calibration.md](htof/htof_calibration.md)（手順・現状・既知の問題）、[htof/htof_code_review_20260926.md](htof/htof_code_review_20260926.md)（コードレビューと整理の提案）
- DC: [dc/blc_global_ofs_tuning_run03773.md](dc/blc_global_ofs_tuning_run03773.md)
- TPC: [tpc/tpc_phase_plot_contract.md](tpc/tpc_phase_plot_contract.md)
- conf: [conf_templates_restore.md](conf_templates_restore.md)
