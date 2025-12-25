#ifndef ANA_HELPER_
#define ANA_HELPER_

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <tuple>
#include <algorithm>

// ROOTヘッダー
#include <Rtypes.h>
#include <TCanvas.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
#include <TLine.h>
#include <TBox.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TKDE.h>
#include <TGraph.h> 

#include "config.h"
#include "params.h"

// raw と trig をペアで管理する構造体
struct HistPair {
    TH1D* raw;  // raw histogram
    TH1D* trig; // trigged histogram

    HistPair(const TString& object_name, const TString& title, Int_t bins, Double_t range_min, Double_t range_max) {
        raw  = new TH1D(object_name+"_raw",  title, bins, range_min, range_max);
        trig = new TH1D(object_name+"_trig", title, bins, range_min, range_max);
    }

    // // デストラクタでメモリを解放 (これをするとグラフが消えてしまう。ちゃんと動作はするので、グラフ化しないんだったらOK)
    // ~HistPair() {
    //     delete raw;
    //     delete trig;
    // }
};

// raw と trig をペアで管理する構造体
struct HistPair2D {
    TH2D* raw;  // raw histogram
    TH2D* trig; // trigged histogram

    HistPair2D(const TString& object_name, const TString& title, Int_t bins1, Double_t range_min1, Double_t range_max1, Int_t bins2, Double_t range_min2, Double_t range_max2) {
        raw  = new TH2D(object_name+"_raw",  title, bins1, range_min1, range_max1, bins2, range_min2, range_max2);
        trig = new TH2D(object_name+"_trig", title, bins1, range_min1, range_max1, bins2, range_min2, range_max2);
    }

    // // デストラクタでメモリを解放 (これをするとグラフが消えてしまう。ちゃんと動作はするので、グラフ化しないんだったらOK)
    // ~HistPair() {
    //     delete raw;
    //     delete trig;
    // }
};

namespace ana_helper {
    TCanvas* add_tab(TGTab *tab, const char* tabName);

    // -- hodo -----
    void set_tdc_search_range(TH1D *h);
    FitResult tdc_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult adc_fit(TH1D *h, TCanvas *c, Int_t n_c, Int_t n_rebin = 1);
    FitResult htof_adc_fit(TH1D *h, TCanvas *c, Int_t& n_c, TString key, Int_t n_rebin = 1);
    FitResult bht_tot_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult t0_offset_fit(TH1D *h, TCanvas *c, Int_t n_c);
    std::pair<Double_t, Double_t> find_phc_range(TH1D* h, Double_t ratio);
    FitResult phc_fit(TH2D *h, TCanvas *c, Int_t n_c);

    // -- dcana -----
    FitResult dc_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c);
    TGraph* make_drift_function(TH1D *h, TCanvas *c, Int_t n_c, Int_t plane);
    FitResult residual_fit(TH1D *h, TCanvas *c, Int_t n_c); // same as t0 offset fit

    // -- opg -----
    std::vector<Int_t> get_kvc_info(Int_t run_number);
    Int_t get_kvc_segment(Int_t run_number);
}

#endif  // ANA_HELPER_