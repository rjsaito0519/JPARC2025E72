#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH2Poly.h"
#include "TRandom.h"

#include "../include/event_display.h"

static TString path;

void draw_track(Int_t n_rand)
{
    // +------------------------------------+
    // | Load ROOT file and set TTreeReader |
    // +------------------------------------+
    TFile *f = new TFile(path.Data());
    TTreeReader reader("hodo", f);
    Int_t tot_num = reader.GetEntries();

    TTreeReaderValue<std::vector<Double_t>> bh2_raw_seg(reader, "bh2_raw_seg");
    TTreeReaderValue<std::vector<Double_t>> bh2_adc_u(reader, "bh2_adc_u");
    TTreeReaderValue<std::vector<Double_t>> bh2_adc_d(reader, "bh2_adc_d");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> bh2_tdc(reader, "bh2_tdc_s");

    TTreeReaderValue<std::vector<std::vector<Double_t>>> bac_tdc(reader, "bac_tdc_u");
    TTreeReaderValue<std::vector<Double_t>> bac_adc(reader, "bac_adc_u");

    TTreeReaderValue<std::vector<Double_t>> htof_raw_seg(reader, "htof_raw_seg");
    TTreeReaderValue<std::vector<Double_t>> htof_adc_u(reader, "htof_adc_u");
    TTreeReaderValue<std::vector<Double_t>> htof_adc_d(reader, "htof_adc_d");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> htof_tdc(reader, "htof_tdc_s");

    TTreeReaderValue<std::vector<Double_t>> kvc_raw_seg(reader, "kvc_raw_seg");
    TTreeReaderValue<std::vector<Double_t>> kvc_adc(reader, "kvc_adc_s");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> kvc_tdc(reader, "kvc_tdc_s");

    Int_t n = n_rand % tot_num;
    reader.Restart();
    reader.SetEntry(n);

    // +---------------------+
    // | Prepare Time window |
    // +---------------------+
    Double_t bh2_tdc_min  = 713417.0, bh2_tdc_max  =  734174.0;
    Double_t bac_tdc_min  = 692328.0, bac_tdc_max  =  723194.0;
    Double_t htof_tdc_min = 682178.0, htof_tdc_max =  730195.0;
    Double_t kvc_tdc_min  = 500000.0, kvc_tdc_max  = 1500000.0;
    
    // +-------------------+
    // | Prepare histogram |
    // +-------------------+
    // 1. セットアップ（余白）を描画
    EventDisplay* display = new EventDisplay("E72 Event Display");
    display->DrawSetup(); // これでジオメトリが描画される

    // +---------------------+
    // | Fill event and draw |
    // +---------------------+
    for (Int_t i = 0, n = (*bh2_raw_seg).size(); i < n; i++) {
        Bool_t has_hit = false;
        for (Int_t ihit = 0, nhit = (*bh2_tdc)[i].size(); ihit < nhit; ihit++) {
            if (bh2_tdc_min < (*bh2_tdc)[i][ihit] && (*bh2_tdc)[i][ihit] < bh2_tdc_max) has_hit = true;
        }
        if (has_hit) display->FillByDetector("BH2", i, TMath::Sqrt((*bh2_adc_u)[i]*(*bh2_adc_d)[i]));
    }
    {
        Bool_t has_hit = false;
        for (Int_t ihit = 0, nhit = (*bac_tdc)[4].size(); ihit < nhit; ihit++) {
            if (bac_tdc_min < (*bac_tdc)[4][ihit] && (*bac_tdc)[4][ihit] < bac_tdc_max) has_hit = true;
        }
        if (has_hit) display->FillByDetector("BAC", 0, (*bac_adc)[4]);
    }
    for (Int_t i = 0, n = (*htof_raw_seg).size(); i < n; i++) {
        Bool_t has_hit = false;
        for (Int_t ihit = 0, nhit = (*htof_tdc)[i].size(); ihit < nhit; ihit++) {
            if (htof_tdc_min < (*htof_tdc)[i][ihit] && (*htof_tdc)[i][ihit] < htof_tdc_max) has_hit = true;
        }
        if (has_hit) display->FillByDetector("HTOF", i, TMath::Sqrt((*htof_adc_u)[i]*(*htof_adc_d)[i]));
    }
    for (Int_t i = 0, n = (*kvc_raw_seg).size(); i < n; i++) {
        Bool_t has_hit = false;
        for (Int_t ihit = 0, nhit = (*kvc_tdc)[i].size(); ihit < nhit; ihit++) {
            if (kvc_tdc_min < (*kvc_tdc)[i][ihit] && (*kvc_tdc)[i][ihit] < kvc_tdc_max) has_hit = true;
        }
        if (has_hit) display->FillByDetector("KVC", i, (*kvc_adc)[i]);
    }
    display->UpdateDisplay();

    // Close the ROOT file to release memory
    f->Close();
    delete f;
}

void event(Int_t n = -1)
{
    if (n == -1) {
        // Use ROOT global random generator
        // 大きい範囲の乱数が欲しければ適当な上限を指定
        Int_t n_rand = gRandom->Integer(1000000000); // 0〜1e9-1 の整数

        draw_track(n_rand);
    } else {
        draw_track(n);
    }
}

void set_path(TString rootfile_path)
{
    path = rootfile_path;
}