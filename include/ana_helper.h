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

#include "config.h"

// -- fitting result container -----
struct FitResult {
    std::vector<Double_t> par;
    std::vector<Double_t> err;
    Double_t chi_square;
    Int_t ndf;
    Int_t migrad_stats;
    std::vector<Double_t> additional; // 何か追加で値を返したいときのコンテナ
};

namespace ana_helper {
    TCanvas* add_tab(TGTab *tab, const char* tabName);

    // -- hodo -----
    void set_tdc_search_range(TH1D *h);
    FitResult tdc_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult adc_fit(TH1D *h, TCanvas *c, Int_t n_c, Int_t n_rebin = 1);
    FitResult bht_tot_fit(TH1D *h, TCanvas *c, Int_t n_c);
    FitResult t0_offset_fit(TH1D *h, TCanvas *c, Int_t n_c);
    std::pair<Double_t, Double_t> find_phc_range(TH1D* h, Double_t ratio);
    FitResult phc_fit(TH2D *h, TCanvas *c, Int_t n_c);
}

#endif  // ANA_HELPER_