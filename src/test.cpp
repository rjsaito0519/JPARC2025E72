// src/test.cpp
//
// Standalone event viewer using TTreeReader and EventDisplay.
//
// Usage:
//   ./test /path/to/hodo.root
//
// Controls (type on terminal, then Enter):
//   (empty line) : next event
//   b            : previous event
//   r            : random event
//   q            : quit
//

#include <iostream>
#include <string>
#include <vector>
#include <random>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TSystem.h"
#include "TMath.h"

#include "event_display.h"  // EventDisplay class

// 時間窓パラメータ
struct TimeWindows {
    Double_t bh2_tdc_min  = 713417.0;
    Double_t bh2_tdc_max  = 734174.0;
    Double_t bac_tdc_min  = 692328.0;
    Double_t bac_tdc_max  = 723194.0;
    Double_t htof_tdc_min = 682178.0;
    Double_t htof_tdc_max = 730195.0;
    Double_t kvc_tdc_min  = 500000.0;
    Double_t kvc_tdc_max  = 1500000.0;
};

// ------------------------------------------------------
// Helper: fill EventDisplay for the current tree entry
// ------------------------------------------------------
static void FillDisplay(EventDisplay& display,
                        const TimeWindows& win,
                        TTreeReaderValue<std::vector<Double_t>>&              bh2_raw_seg,
                        TTreeReaderValue<std::vector<Double_t>>&              bh2_adc_u,
                        TTreeReaderValue<std::vector<Double_t>>&              bh2_adc_d,
                        TTreeReaderValue<std::vector<std::vector<Double_t>>>& bh2_tdc,
                        TTreeReaderValue<std::vector<std::vector<Double_t>>>& bac_tdc,
                        TTreeReaderValue<std::vector<Double_t>>&              bac_adc,
                        TTreeReaderValue<std::vector<Double_t>>&              htof_raw_seg,
                        TTreeReaderValue<std::vector<Double_t>>&              htof_adc_u,
                        TTreeReaderValue<std::vector<Double_t>>&              htof_adc_d,
                        TTreeReaderValue<std::vector<std::vector<Double_t>>>& htof_tdc,
                        TTreeReaderValue<std::vector<Double_t>>&              kvc_raw_seg,
                        TTreeReaderValue<std::vector<Double_t>>&              kvc_adc,
                        TTreeReaderValue<std::vector<std::vector<Double_t>>>& kvc_tdc)
{
    // ------------ BH2 ------------
    {
        const auto& raw_seg = *bh2_raw_seg;
        const auto& adc_u   = *bh2_adc_u;
        const auto& adc_d   = *bh2_adc_d;
        const auto& tdc     = *bh2_tdc;

        const Int_t nseg = static_cast<Int_t>(raw_seg.size());
        for (Int_t i = 0; i < nseg; ++i) {
            Bool_t has_hit = kFALSE;

            if (i < static_cast<Int_t>(tdc.size())) {
                const auto& tvec = tdc[i];
                for (size_t ihit = 0; ihit < tvec.size(); ++ihit) {
                    Double_t t = tvec[ihit];
                    if (win.bh2_tdc_min < t && t < win.bh2_tdc_max) {
                        has_hit = kTRUE;
                        break;
                    }
                }
            }

            if (has_hit &&
                i < static_cast<Int_t>(adc_u.size()) &&
                i < static_cast<Int_t>(adc_d.size())) {
                Double_t au = adc_u[i];
                Double_t ad = adc_d[i];
                Double_t adc = TMath::Sqrt(au * ad);
                display.FillByDetector("BH2", i, adc);
            }
        }
    }

    // ------------ BAC ------------
    // もとのマクロと同じく、ch=4 を BAC セグメント 0 にマッピング
    {
        const auto& tdc = *bac_tdc;
        const auto& adc = *bac_adc;

        if (tdc.size() > 4 && adc.size() > 4) {
            Bool_t has_hit = kFALSE;
            const auto& tvec = tdc[4];

            for (size_t ihit = 0; ihit < tvec.size(); ++ihit) {
                Double_t t = tvec[ihit];
                if (win.bac_tdc_min < t && t < win.bac_tdc_max) {
                    has_hit = kTRUE;
                    break;
                }
            }

            if (has_hit) {
                display.FillByDetector("BAC", 0, adc[4]);
            }
        }
    }

    // ------------ HTOF ------------
    {
        const auto& raw_seg = *htof_raw_seg;
        const auto& adc_u   = *htof_adc_u;
        const auto& adc_d   = *htof_adc_d;
        const auto& tdc     = *htof_tdc;

        const Int_t nseg = static_cast<Int_t>(raw_seg.size());
        for (Int_t i = 0; i < nseg; ++i) {
            Bool_t has_hit = kFALSE;

            if (i < static_cast<Int_t>(tdc.size())) {
                const auto& tvec = tdc[i];
                for (size_t ihit = 0; ihit < tvec.size(); ++ihit) {
                    Double_t t = tvec[ihit];
                    if (win.htof_tdc_min < t && t < win.htof_tdc_max) {
                        has_hit = kTRUE;
                        break;
                    }
                }
            }

            if (has_hit &&
                i < static_cast<Int_t>(adc_u.size()) &&
                i < static_cast<Int_t>(adc_d.size())) {
                Double_t au = adc_u[i];
                Double_t ad = adc_d[i];
                Double_t adc = TMath::Sqrt(au * ad);
                display.FillByDetector("HTOF", i, adc);
            }
        }
    }

    // ------------ KVC ------------
    {
        const auto& raw_seg = *kvc_raw_seg;
        const auto& adc     = *kvc_adc;
        const auto& tdc     = *kvc_tdc;

        const Int_t nseg = static_cast<Int_t>(raw_seg.size());
        for (Int_t i = 0; i < nseg; ++i) {
            Bool_t has_hit = kFALSE;

            if (i < static_cast<Int_t>(tdc.size())) {
                const auto& tvec = tdc[i];
                for (size_t ihit = 0; ihit < tvec.size(); ++ihit) {
                    Double_t t = tvec[ihit];
                    if (win.kvc_tdc_min < t && t < win.kvc_tdc_max) {
                        has_hit = kTRUE;
                        break;
                    }
                }
            }

            if (has_hit && i < static_cast<Int_t>(adc.size())) {
                display.FillByDetector("KVC", i, adc[i]);
            }
        }
    }
}

// -------------------------
// main
// -------------------------
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " /path/to/hodo.root" << std::endl;
        return 1;
    }

    const char* filename = argv[1];

    // ROOT application (GUI 用)
    TApplication app("event_viewer", &argc, argv);

    // ファイルを開く
    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return 1;
    }

    // TTree を取得
    TTree* tree = nullptr;
    f->GetObject("hodo", tree);
    if (!tree) {
        std::cerr << "Error: cannot find TTree 'hodo' in file " << filename << std::endl;
        f->Close();
        delete f;
        return 1;
    }

    Long64_t nEntries = tree->GetEntries();
    if (nEntries <= 0) {
        std::cerr << "Error: tree has no entries." << std::endl;
        f->Close();
        delete f;
        return 1;
    }

    // TTreeReader とブランチをセット
    TTreeReader reader(tree);

    TTreeReaderValue<std::vector<Double_t>>              bh2_raw_seg(reader, "bh2_raw_seg");
    TTreeReaderValue<std::vector<Double_t>>              bh2_adc_u(reader,   "bh2_adc_u");
    TTreeReaderValue<std::vector<Double_t>>              bh2_adc_d(reader,   "bh2_adc_d");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> bh2_tdc(reader,     "bh2_tdc_s");

    TTreeReaderValue<std::vector<std::vector<Double_t>>> bac_tdc(reader,     "bac_tdc_u");
    TTreeReaderValue<std::vector<Double_t>>              bac_adc(reader,     "bac_adc_u");

    TTreeReaderValue<std::vector<Double_t>>              htof_raw_seg(reader, "htof_raw_seg");
    TTreeReaderValue<std::vector<Double_t>>              htof_adc_u(reader,   "htof_adc_u");
    TTreeReaderValue<std::vector<Double_t>>              htof_adc_d(reader,   "htof_adc_d");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> htof_tdc(reader,     "htof_tdc_s");

    TTreeReaderValue<std::vector<Double_t>>              kvc_raw_seg(reader,  "kvc_raw_seg");
    TTreeReaderValue<std::vector<Double_t>>              kvc_adc(reader,      "kvc_adc_s");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> kvc_tdc(reader,      "kvc_tdc_s");

    TimeWindows win;

    // EventDisplay を 1 回だけ生成（キャンバス使い回し）
    EventDisplay display("E72 Event Display");
    display.DrawSetup();
    gSystem->ProcessEvents();

    // イベント番号制御
    Long64_t entry = 0;
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<Long64_t> dist(0, nEntries - 1);

    std::cout << "Loaded file: " << filename << std::endl;
    std::cout << "Total entries: " << nEntries << std::endl;
    std::cout << "Controls: [Enter: next, b: prev, r: random, q: quit]" << std::endl;

    while (true) {
        if (entry < 0) entry = 0;
        if (entry >= nEntries) entry = nEntries - 1;

        auto status = reader.SetEntry(entry);
        if (status != TTreeReader::kEntryValid) {
            std::cerr << "Warning: entry " << entry
                      << " is not valid (status=" << status << ")" << std::endl;
        } else {
            // 前イベントの内容をクリア
            display.ResetDisplay();

            // 今のエントリで塗り直し
            FillDisplay(display, win,
                        bh2_raw_seg, bh2_adc_u, bh2_adc_d, bh2_tdc,
                        bac_tdc, bac_adc,
                        htof_raw_seg, htof_adc_u, htof_adc_d, htof_tdc,
                        kvc_raw_seg, kvc_adc, kvc_tdc);

            // キャンバス更新（fCanvas には触らない）
            display.UpdateDisplay();
            gSystem->ProcessEvents();

            std::cout << "Showing entry " << entry << " / " << (nEntries - 1) << std::endl;
        }

        // ターミナルからコマンド入力
        std::cout << "[Enter/b/r/q] > " << std::flush;
        std::string line;
        if (!std::getline(std::cin, line)) {
            break; // EOF
        }

        if (line.empty()) {
            // Enter → 次
            ++entry;
        } else if (line == "b") {
            --entry;
        } else if (line == "r") {
            entry = dist(gen);
        } else if (line == "q") {
            break;
        } else {
            std::cout << "Unknown command: " << line << std::endl;
        }
    }

    f->Close();
    delete f;
    return 0;
}
