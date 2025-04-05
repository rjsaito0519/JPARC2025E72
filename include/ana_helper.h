#ifndef ANA_HELPER_
#define ANA_HELPER_

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <tuple>

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

// -- fitting result container -----
struct PwaOption {
    Bool_t debug_cout;
    Bool_t fix_s_wave;
    Int_t min_print_level;
    Int_t min_strategy;
    TString anl_osaka_model;
};


namespace ana_helper {
    TCanvas* add_tab(TGTab *tab, const char* tabName);
    std::vector<std::vector<Double_t>> load_data(TString path);
}

#endif  // ANA_HELPER_