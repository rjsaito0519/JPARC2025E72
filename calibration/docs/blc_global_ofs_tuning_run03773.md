# BLC グローバル位置（Ofs）調整の流れ — run03773 (K)

基準ラン: **run03773（Kaon）**  
対象パラメータ: `param/DCGEO/e72/DCGeomParam_run03773_K`  
（グローバル平行移動は `X/Y` を直接書き換えず、`apply_dcgeo_offset.py` で **Ofs に射影**する）

関連ツール:

- `myanalysis/calibration/scripts/apply_dcgeo_offset.py`
- `myanalysis/analysis/src/blc_bh2_correlation_review.cpp`
- `myanalysis/analysis/src/d5_wire_resi_review.cpp`

PDF 出力例: `/home/had/sryuta/JPARC2025E72/results/img/run03773/`

---

## 方針（最初に決めたこと）

1. **BLC2（BcOut）は基本固定**し、まず **BLC1（BcIn）の Y** を合わせる。
2. **X は指標が弱いので触らない**（Y のみ）。
3. チャンバー平行移動 `(dx, dy)` → 各面の局所 s 軸へ投影して **Ofs のみ更新**（`X/Y/Z/TA/RA` は変更しない）。
4. 変更のたびにバックアップを取り、`--skip-xy` なしの再適用で Ofs を二重計上しない。

---

## 調整順序

### 0. 準備

- DCGEO のバックアップ（例: `DCGeomParam_run03773_K.bk`）。
- サーベイ由来の `X/Y` がまだ Ofs に載っていなければ、先に  
  `apply_dcgeo_offset.py <DCGEO> --write`  
  （ファイルの `X/Y` → Ofs）。以降の残差補正は必ず `--skip-xy`。

### 1. 現状確認（調整前）

| 確認 | コマンド／見どころ |
|------|-------------------|
| BLC1–BLC2 相対 Y | `blc_bh2_correlation_review 3773 --kaon` |
| | chamber-center Δy、外挿 Δy（BLC1@z ≈ 2900、BLC2@z ≈ −1299.815） |
| D5 wire residual | `d5_wire_resi_review`（decode ROOT） |
| | LocalFit（面内）vs D5Fit（ジョイント） |

調整前の目安（run03773）:

- 外挿 Δy @ 2900: Gauss μ ≈ **+3.6 mm** → 「BLC1 y を +μ（BLC2 固定）」
- chamber-center Δy: μ ≈ **+4.6 mm**
- X は専用指標が弱く、今回は固定

### 2. BLC1 の Y を合わせる（主調整）

指標: **BLC2 固定**のまま、BLC1 外挿先での

\[
\Delta y = y_{\mathrm{BLC2}}(z_{\mathrm{BLC2}}) - y_{\mathrm{BLC1}}(z)
\]

を 0 に近づける（コード推奨どおり `resiy` は Gauss μ と同符号）。

例（初回は外挿 @ 2900 基準）:

```bash
python3 myanalysis/calibration/scripts/apply_dcgeo_offset.py \
  param/DCGEO/e72/DCGeomParam_run03773_K \
  --skip-xy --bcin --resiy <μ>   # 例: 3.57
  # 確認後 --write
```

その後:

1. BcIn（必要なら D5）再 decode  
2. `blc_bh2_correlation_review` / `d5_wire_resi_review` 再生成  
3. 符号が逆なら `--scale -1` または `resiy` の符号反転

結果（このキャンペーン）: 外挿 Δy @ 2900 の μ ≈ **0**。

### 3. 外挿 z 依存の確認（2900 が適当か）

`blc_bh2_correlation_review` に追加した hist:

- BLC2 は **z = −1299.815（チャンバー中心）固定**
- BLC1 の外挿 z を **2900 ± 300** でスキャン
- 縦軸: 上記と同じ Δy

数値チェック（調整後・同一イベント選択）:

| z_BLC1 | Gauss μ [mm] |
|--------|--------------|
| 2600   | ≈ −0.12      |
| 2800   | ≈ −0.04      |
| 2900   | ≈ 0.00       |
| 3000   | ≈ +0.04      |
| 3200   | ≈ +0.12      |

→ ±300 mm 動かしても **0.1 mm 級**。2900 一点合わせでこの窓は十分。  
一方 **chamber-center（z=0）Δy ≈ +1 mm** は別定義の残り（外挿窓の再チューニングでは消えない）。

### 4. BLC2 の微調（D5Fit の U/V ジグザグ）

BLC1 合わせ後も `d5_wire_resi_review` で:

- LocalFit: きれい（面内 OK）
- D5Fit: BLC2 で U/V が ±0.2 mm 程度ジグザグ
- **BLC2a と BLC2b で符号が逆** → 共通 `--bcout --resiy` は不可

面平均の D5Fit mean から |dy| ≈ 0.20 × √2 ≈ **0.28 mm**、a/b 逆符号で試す:

```bash
python3 myanalysis/calibration/scripts/apply_dcgeo_offset.py \
  param/DCGEO/e72/DCGeomParam_run03773_K \
  --skip-xy \
  --residual-offset BLC2a:0,0.28 \
  --residual-offset BLC2b:0,-0.28
  # dry-run で ΔOfs 確認後 --write
```

符号が逆に動いたら a/b の dy を入れ替え（または `--scale -1`）。

再 decode → `d5_wire_resi_review` でジグザグ減少、`blc_bh2` で外挿 Δy が壊れないことを確認。

→ このステップ後、全体として **いい感じに収束**。

---

## 順序の要約（チェックリスト）

```text
[0] バックアップ / （必要なら）X,Y→Ofs のサーベイ投影
[1] blc_bh2 + d5_wire_resi で現状確認
[2] BLC1 のみ --skip-xy --bcin --resiy  （X は触らない、BLC2 は触らない）
[3] 再 decode → 外挿 Δy≈0 を確認
[4] Δy vs z（2900±300）で z 依存が小さいことを確認
[5] D5Fit で BLC2 U/V ジグザグが残る場合のみ
    BLC2a / BLC2b を別 dy で微調（共通 bcout は使わない）
[6] 再 decode → D5Fit + 外挿 Δy の両立を確認
```

---

## 触らなかった／残したもの

- **X オフセット**: 指標が弱いため未調整。
- **chamber-center Δy ~1 mm**: 外挿（2900 付近）とは別指標。今回の主目標外として様子見可。
- **面内 LocalFit 用の細かい Ofs**: Local が良いので今回は触らない（必要なら後段で `update_param.py residual`）。

---

## 座標・定義の覚え書き

- BLC2 チャンバー中心: **z = −1299.815 mm**（`blc_bh2_correlation_review` の `kBlc2BeamProfileZ`）
- 一点外挿の従来値: BLC1 **z = 2900 mm**（`kBlc1ExtrapZ`）
- Δy の符号: `y_BLC2 − y_BLC1`。Gauss μ が正 → BLC1 の y を **+μ**（BLC2 固定）がコード上の推奨。
- `apply_dcgeo_offset.py`: `delta_ofs = scale * (dsdx*dx + dsdy*dy + dsdz*dz)`（TA/RA から s 軸）。典型 BLC は RA=0 で dz は Ofs に効かない。

詳細な幾何・Ofs の意味は `calibration/docs/dcgeo_tracking_geometry.md` を参照。
