#include "event_display.h" // 対応するヘッダーをインクルード

#include <TCanvas.h>
#include <TH2Poly.h>      
#include <TLine.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TStyle.h>
#include <TPad.h>         
#include <TMath.h>        
#include <TString.h>      
#include <iostream>       
#include <cstring> // (TStringへの移行により不要になったが、念のため残す)

EventDisplay::EventDisplay(const TString& canvasName) 
    : fCanvas(nullptr), fDetectorPoly(nullptr) 
{
    // キャンバスを作成
    fCanvas = new TCanvas("c1", canvasName.Data(), 1200, 900); 
    
    // ヒストグラム（ジオメトリ）を作成
    fDetectorPoly = new TH2Poly("gDetectorPoly", "Experimental Setup;X [cm];Y [cm]", -80, 80, -80, 80); 
    fDetectorPoly->SetFillColor(kWhite); 
    fDetectorPoly->SetLineColor(kBlack); 
    fDetectorPoly->SetMinimum(0.1); // 0のビンを白く描画
}

EventDisplay::~EventDisplay() {
    // fDetectorPoly (TH2Poly) は gDirectory が所有権を持つため、deleteしない
    // fCanvas は TApplication 終了時に自動で閉じられる
}

/**
 * @brief 検出器のジオメトリ（余白）を描画します
 */
void EventDisplay::DrawSetup() {
    if (!fCanvas || !fDetectorPoly) return;
    
    fCanvas->cd(); 
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(kTemperatureMap); 

    addTrapezoidalHTOF(0.0, 0.0);
    addDetectorPoly(fKvcSegmentMap, 61.0, 3.0, 11.3, 32.1, 8);
    addDetectorPoly(fBh2SegmentMap, -59.57, 3.0, -10.5, 10.5, 15);
    addDetectorPoly(fBacSegmentMap, -50.36, 3.0, -7.45, 4.05, 1);

    fDetectorPoly->Draw("COLZ L"); 
    fCanvas->Modified();
    fCanvas->Update();
}

void EventDisplay::FillByDetector(const TString& detectorName,
                                  Int_t segmentId,
                                  Double_t adc_value)
{
    if (!fDetectorPoly) return;

    // どのマップを使うか選ぶ
    std::map<Int_t, Int_t>* segMap = nullptr;

    if (detectorName == "BH2") {
        segMap = &fBh2SegmentMap;
    } else if (detectorName == "BAC") {
        segMap = &fBacSegmentMap;
    } else if (detectorName == "HTOF") {
        segMap = &fHtofSegmentMap;
    } else if (detectorName == "KVC") {
        segMap = &fKvcSegmentMap;
    } else {
        std::cerr << "EventDisplay::FillByDetector: unknown detector '"
                  << detectorName << "'" << std::endl;
        return;
    }

    // segmentId → bin 番号 を探す
    auto it = segMap->find(segmentId);
    if (it == segMap->end()) {
        std::cerr << "EventDisplay::FillByDetector: segmentId "
                  << segmentId << " not found for detector '"
                  << detectorName << "'" << std::endl;
        return;
    }

    Int_t bin = it->second;

    // bin の妥当性チェック（お守り）
    if (bin <= 0 || bin > fDetectorPoly->GetNumberOfBins()) {
        std::cerr << "EventDisplay::FillByDetector: invalid bin "
                  << bin << " for detector '" << detectorName << "'" << std::endl;
        return;
    }

    // ★ここで「このイベントの値」を上書きセット
    fDetectorPoly->SetBinContent(bin, adc_value);
}


void EventDisplay::ResetDisplay() {
    if (!fDetectorPoly) return;
    fDetectorPoly->Reset("");
}

void EventDisplay::ClearEvent()
{
    // ジオメトリ（形）はそのまま、中身（Bin Content）だけゼロにする
    if (fDetectorPoly) {
        // TH2Poly は TH1 の派生なので Reset が使える
        fDetectorPoly->Reset("ICES");
    }
}


void EventDisplay::UpdateDisplay() {
    if (!fCanvas || !fDetectorPoly) return;

    fCanvas->cd();
    // ジオメトリ (fDetectorPoly) 自体が "bin content = ADC" を持っている想定
    // 最初の DrawSetup() で軸や枠は描いてあるはずですが、
    // 毎回上書きしたいならここでもう一度 Draw して OK
    fDetectorPoly->Draw("COLZ");
    
    gPad->Modified();
    gPad->Update();
}

void EventDisplay::addDetectorPoly(std::map<Int_t, Int_t>& map, Double_t x_pos, Double_t width, Double_t y_min, Double_t y_max, Int_t n_segments) {
    if (!fDetectorPoly) return;
    Double_t total_height = y_max - y_min;
    Double_t seg_height = total_height / n_segments;
    for (Int_t i = 0; i < n_segments; ++i) { 
        Double_t y1 = y_min + i * seg_height;
        Double_t y2 = y1 + seg_height;
        Double_t px[4] = {x_pos, x_pos + width, x_pos + width, x_pos};
        Double_t py[4] = {y1, y1, y2, y2};
        
        fDetectorPoly->AddBin(4, px, py);
        Int_t bin_index = fDetectorPoly->GetNumberOfBins();
        map[i] = bin_index; 
    }
}

void EventDisplay::rotate_point(Double_t x, Double_t y, Double_t center_x, Double_t center_y, Double_t theta, Double_t &out_x, Double_t &out_y) {
    Double_t translated_x = x - center_x;
    Double_t translated_y = y - center_y;
    Double_t rotated_x = translated_x * TMath::Cos(theta) - translated_y * TMath::Sin(theta);
    Double_t rotated_y = translated_x * TMath::Sin(theta) + translated_y * TMath::Cos(theta);
    out_x = rotated_x + center_x;
    out_y = rotated_y + center_y;
}

void EventDisplay::addTrapezoidalHTOF(Double_t x_center, Double_t y_center) {
    if (!fDetectorPoly) return;
    
    const Double_t l = 687.264 / 2.0 / 10.0; 
    const Double_t htof_width = 70.0 / 10.0;   
    const Double_t htof_thick = 2.0; 
    const Double_t htof_gap = 1.0 / 10.0;      
    const Double_t tan3pi8 = TMath::Tan(3.0 * TMath::Pi() / 8.0); 
    const Double_t l_inner = l; 
    const Double_t l_outer = l + htof_thick; 
    const Double_t l_split = l + htof_thick / 2.0; 

    Double_t px_world[4];
    Double_t py_world[4];
    
    for (Int_t i = 0; i < 8; ++i) { 
        Double_t rotation_angle = TMath::Pi() / 4.0 * (i + 2.0); 

        for (Int_t j = 0; j < 4; ++j) { 
            
            Int_t segment_id_base = -1;
            if (i == 0) { // Left (Seg 0-5)
                if (j == 0) segment_id_base = 5;
                else if (j == 1) segment_id_base = 3; // [3, 4]
                else if (j == 2) segment_id_base = 1; // [1, 2]
                else if (j == 3) segment_id_base = 0;
            }
            else if (i == 1) segment_id_base =  9 - j;  // (Seg 6-9)
            else if (i == 2) segment_id_base = 13 - j; // Bottom (Seg 10-13)
            else if (i == 3) segment_id_base = 17 - j; // (Seg 14-17)
            else if (i == 4) segment_id_base = 21 - j; // Right (Seg 18-21)
            else if (i == 5) segment_id_base = 25 - j; // (Seg 22-25)
            else if (i == 6) segment_id_base = 29 - j; // Top (Seg 26-29)
            else if (i == 7) segment_id_base = 33 - j; // (Seg 30-33)

            Double_t y1 = l_inner, y2 = l_outer; 
            Double_t x_start = (htof_width + htof_gap) * (j - 2);
            Double_t x_end = x_start + htof_width;
            
            Double_t local_x_verts[4];
            Double_t local_y_verts[4];

            if (j == 0) { 
                local_x_verts[0] = x_start;                       local_y_verts[0] = y1;
                local_x_verts[1] = x_end;                         local_y_verts[1] = y1;
                local_x_verts[2] = x_end;                         local_y_verts[2] = y2;
                local_x_verts[3] = x_start - htof_thick / tan3pi8; local_y_verts[3] = y2;
            } else if (j == 3) { 
                local_x_verts[0] = x_start;                       local_y_verts[0] = y1;
                local_x_verts[1] = x_end;                         local_y_verts[1] = y1;
                local_x_verts[2] = x_end + htof_thick / tan3pi8;   local_y_verts[2] = y2;
                local_x_verts[3] = x_start;                       local_y_verts[3] = y2;
            } else { 
                local_x_verts[0] = x_start; local_y_verts[0] = y1;
                local_x_verts[1] = x_end;   local_y_verts[1] = y1;
                local_x_verts[2] = x_end;   local_y_verts[2] = y2;
                local_x_verts[3] = x_start; local_y_verts[3] = y2;
            }

            Bool_t isSplit = (i == 0 && (j == 1 || j == 2)); 

            if (isSplit) {
                // ビンA (内側)
                local_y_verts[0] = l_inner; local_y_verts[1] = l_inner;
                local_y_verts[2] = l_split; local_y_verts[3] = l_split;
                for(Int_t k=0; k<4; ++k) {
                    rotate_point(local_x_verts[k], local_y_verts[k], 0, 0, rotation_angle, px_world[k], py_world[k]);
                    px_world[k] += x_center; py_world[k] += y_center;
                }
                fDetectorPoly->AddBin(4, px_world, py_world);
                Int_t bin_index_A = fDetectorPoly->GetNumberOfBins(); 
                
                // ビンB (外側)
                local_y_verts[0] = l_split; local_y_verts[1] = l_split;
                local_y_verts[2] = l_outer; local_y_verts[3] = l_outer;
                for(Int_t k=0; k<4; ++k) {
                    rotate_point(local_x_verts[k], local_y_verts[k], 0, 0, rotation_angle, px_world[k], py_world[k]);
                    px_world[k] += x_center; py_world[k] += y_center;
                }
                fDetectorPoly->AddBin(4, px_world, py_world);
                Int_t bin_index_B = fDetectorPoly->GetNumberOfBins();

                if (segment_id_base == 1) { 
                    fHtofSegmentMap[1] = bin_index_A; 
                    fHtofSegmentMap[2] = bin_index_B; 
                } else if (segment_id_base == 3) { 
                    fHtofSegmentMap[3] = bin_index_A; 
                    fHtofSegmentMap[4] = bin_index_B; 
                }
            
            } else {
                for(Int_t k=0; k<4; ++k) {
                    rotate_point(local_x_verts[k], local_y_verts[k], 0, 0, rotation_angle, px_world[k], py_world[k]);
                    px_world[k] += x_center; py_world[k] += y_center;
                }
                fDetectorPoly->AddBin(4, px_world, py_world);
                Int_t bin_index = fDetectorPoly->GetNumberOfBins(); 
                fHtofSegmentMap[segment_id_base] = bin_index;
            }
        }
    }
}