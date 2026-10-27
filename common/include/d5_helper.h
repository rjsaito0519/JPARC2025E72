#ifndef D5_HELPER_
#define D5_HELPER_

// Helpers shared by the D5 / BLC1-position review tools (d5_pos_review, d5_plane_resi_qa).

#include <cmath>
#include <map>

#include <Rtypes.h>
#include <TString.h>

class TH1;

namespace d5_helper {

    // beam selection of the D5 tools: --pion (default), --kaon, --all
    constexpr Int_t kBeamFlagPion = 1;
    constexpr Int_t kBeamFlagKaon = 2;
    constexpr Int_t kBeamFlagAll = -1;

    struct BeamSel {
        Int_t flag = kBeamFlagPion;   // d5 tree beam_flag to select (kBeamFlagAll: no selection)
        TString param_suffix = "Pi";  // suffix of the params / decode files (Pi or K; K for --all)
        TString hist_suffix = "_Pi";  // suffix of the per-beam-flag histograms ("" for --all)
        TString out_tag = "Pi";       // tag in the output file names (Pi, K, All)
        TString label = "pi (beam_flag=1)";
    };

    // true if arg is one of --pion / --kaon / --all (sel is set), false otherwise
    Bool_t parse_beam_option(const TString& arg, BeamSel& sel);

    // Gauss fit in mean +- 1.5 sigma, iterated n_iter times from (mean0, sigma0).
    // The function stays attached to h as "g_<name>" (kNotDraw) for drawing.
    // Fewer than min_entries entries: NaN.
    struct GaussCore {
        Double_t mean = NAN;
        Double_t mean_err = NAN;
        Double_t sigma = NAN;
    };
    GaussCore fit_core(TH1* h, Double_t mean0, Double_t sigma0, Int_t n_iter,
                       Int_t min_entries = 100);

    // DCGEO layer name -> (Z, TA); BLC layers only
    struct DCLayer {
        Double_t z = 0.0;
        Double_t ta = 0.0;
    };
    std::map<TString, DCLayer> load_dcgeo_blc(const TString& path);

    // mean Z of the layers named "<chamber>-*" (e.g. BLC1a); false if none
    Bool_t chamber_mean_z(const std::map<TString, DCLayer>& geo, const TString& chamber,
                          Double_t& z);

} // namespace d5_helper

#endif
