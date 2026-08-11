# TPC phase: fit 出力と plot の契約

`tpc_phase_from_tpcbcout` が書く `TpcPhase_*.root` と、`tpc_phase_plot` が読む前提を固定する。  
ここを破ると close-up の fit 赤線が消えるなど、修正のたびに時間が溶ける。

関連ソース:

- `myanalysis/calibration/src/tpc_phase_from_tpcbcout.cpp`
- `myanalysis/calibration/src/tpc_phase_plot.cpp`

## 運用前提

- **正準は 1-step**（`--fit-step 1`）。多段 `N>1` は CLI 互換のみ（非推奨・想定外）。
- デフォルト（baseline ON）は **2 段階**:
  1. pre-fit（`--base` 相当）で `t`, `w` を決める
  2. baseline 本 fit で `t`,`w` 固定・`c0`/`A` を当てる（採用は 2）
- **ステップ幅**: fit は `--step-width`（既定 **0.1 ns**）で安定化。保存時は `A,t`（と `c0`）を取り、**`w=0.001` の新しい TF1 を作り直して** TGraph / Fit / Fallback `p2` に書く（Fix 済みパラへの `SetParameter` は使わない）。TGraph は全域サンプルに加え **ステップ近傍を密サンプリング**する。
- **MINUIT `fit_status != 0`（非収束）でも最終パラメータを採用**して TGraph / Fallback を書く（Warning 付き）。従来の `--base` / StepSum 互換。非収束＝「結果を捨てる」ではない。

## TPCPRM

- aty=3（kCobo PhaseShift）は常に **`(p0,p1,p2)=(A,t,w)` の 3 パラメータ**。個数を変えない。
- `c0` は aty=3 に入れない。TGraph（baseline 込み）と `--ofs-update`（aty=2 p0）のみ。

## phase ファイルに書くもの

| オブジェクト | 必須 | 用途 |
|---|---|---|
| `TpcPhase_Cobo%d` (TGraph) | **必須** | 補正本体 / plot の **赤実線** |
| `TpcPhase_Profile_Cobo%d` | 推奨 | オレンジ profile 点 |
| `TpcPhase_ProfileFitUsed_Cobo%d` | 推奨 | fit 窓・青縦線 |
| `TpcPhase_CoboFallback` (TTree) | **必須** | `p0,p1,p2,nstep_fit,t0_init,c0,...` |
| `TpcPhase_Fit_Cobo%d` (TF1, 3par A,t,w) | 推奨 | plot フォールバック |
| `TpcPhase_FitRaw_Cobo%d` (TF1, c0+A,t,w) | 推奨 | baseline 時の raw |

**`TpcPhase_Cobo*` を書かないと close-up の fit 赤線は出ない。**  
Profile / `t0_init` だけ残ってもオレンジ点とマゼンタ縦線だけ見える、という壊れ方をする。

### Fallback フィールド（1-step）

- `nstep_fit == 1`: `p0,p1,p2` = `(A,t,w)` 有効。plot の **赤点線縦線 = p1**
- `nstep_fit > 1`: 多段（非推奨）。`p0–p2` は無効扱い
- `nstep_fit == -1`: 失敗（TGraph も無い想定）
- `t0_init`: 採用ステップ時刻（正準 2 段階では pre-fit の `t`）。マゼンタ縦線

## plot の線の意味

| 見た目 | 意味 |
|---|---|
| 赤実線 | fit 曲線（主に `TpcPhase_Cobo*`） |
| 赤点線縦線 | ステップ位置 **p1** |
| マゼンタ点線縦線 | `t0_init` |
| 青点線縦線 | fit 窓（ProfileFitUsed） |
| オレンジ丸 | profile 点 |

close-up は **p1 を付け直さない**（Fallback / Fit の `t` をそのままズーム中心・縦線に使う）。

## 修正時チェックリスト

1. `TpcPhase_*.root` に `TpcPhase_Cobo0`… があるか
2. `TpcPhase_CoboFallback` で `nstep_fit==1` かつ有限の `p1` か
3. ログに `non-converged; using last params` が出ても TGraph が書かれているか
4. `tpc_phase_plot --phase ...` で page1 / close-up に **赤実線 + 赤点線(p1)** が出るか
5. aty=3 更新が `(A,t,w)` の 3 列のままか
