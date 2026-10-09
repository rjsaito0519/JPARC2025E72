# conf テンプレート配置と復元メモ

最終更新: 2026-10-09

`create_runlist.py` / `dst_create_runlist.py` は param の analyzer conf テンプレートを読む。

## 探索順

1. `param/conf/example/<name>`（文書・実環境で想定している正）
2. 無ければ `param/conf/<name>`（一部の param バックアップ配置）

どちらも無い場合はエラー終了する。

## 既存 run conf

- 既に `param/conf/<SUB_DIR>/analyzer_runNNNNN_*.conf` がある場合、スクリプトは **上書きしない**（手追加の optional key や run 固有 TPCPOS を消さない）。
- **新規生成**の DST conf だけ `DEFAULTS["TPCPOS:"]` = `TPCPositionCorrectionMap_0` を使う。既存 conf が KinFit map 等を指していてもそのまま残る。

## 復元時の推奨

param をバックアップから戻すときは、可能なら `conf/example/` にテンプレート（era conf・DST/TPCHit example）を揃える。無くても上記 fallback で動くが、example 配置の方が文書と一致する。
