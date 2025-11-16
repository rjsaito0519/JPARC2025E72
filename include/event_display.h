#ifndef EVENTDISPLAY_H
#define EVENTDISPLAY_H

#include <TROOT.h>
#include <TString.h> // 修正: TString のために追加
#include <map>

// ROOTクラスの前方宣言
class TH2Poly;
class TCanvas;

/**
 * @brief E72セットアップのイベントディスプレイを管理するヘルパークラス
*/
class EventDisplay {
public:
    /**
     * @brief コンストラクタ
     * @param canvasName キャンバスのタイトル
     */
    EventDisplay(const TString& canvasName = "E72 Event Display");
    
    /**
     * @brief デストラクタ
     */
    virtual ~EventDisplay();

    /**
     * @brief 検出器のジオメトリ（余白）を描画します (最初に1回だけ呼び出す)
     */
    void DrawSetup();

    // --- イベントループ内で使用する関数 ---

    /**
     * @brief 検出器名とセグメントID (0始まり) を指定してADC値でFillする
     */
    void FillByDetector(const TString& detectorName, Int_t segmentId, Double_t adc_value); // 修正: const char* -> const TString&

    /**
     * @brief 全てのビンのADC値（内容）をリセットします (イベント毎に呼ぶ)
     */
    void ResetDisplay();

    void ClearEvent();

    /**
     * @brief キャンバスの描画を更新します
     */
    void UpdateDisplay();

    /**
     * @brief TCanvasオブジェクトへのポインタを取得します
     */
    TCanvas* GetCanvas() { return fCanvas; }

private:
    // --- プライベートヘルパー関数 (ジオメトリ定義用) ---
    void addDetectorPoly(std::map<Int_t, Int_t>& map, Double_t x_pos, Double_t width, Double_t y_min, Double_t y_max, Int_t n_segments);
    void rotate_point(Double_t x, Double_t y, Double_t center_x, Double_t center_y, Double_t theta, Double_t &out_x, Double_t &out_y);
    void addTrapezoidalHTOF(Double_t x_center, Double_t y_center);

    // --- メンバー変数 (グローバル変数の代わり) ---
    TCanvas* fCanvas;       // このディスプレイが描画するキャンバス
    TH2Poly* fDetectorPoly; // ジオメトリ本体
    
    // IDマップ
    std::map<Int_t, Int_t> fHtofSegmentMap;
    std::map<Int_t, Int_t> fBh2SegmentMap;
    std::map<Int_t, Int_t> fBacSegmentMap;
    std::map<Int_t, Int_t> fKvcSegmentMap;
};

#endif // EVENTDISPLAY_H