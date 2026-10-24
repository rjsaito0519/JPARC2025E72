#include "d5_helper.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include <TF1.h>
#include <TH1.h>

namespace d5_helper {

    // ____________________________________________________________________________________________
    Bool_t parse_beam_option(const TString& arg, BeamSel& sel) {
        if (arg == "--pion") {
            sel = BeamSel();
        } else if (arg == "--kaon") {
            sel.flag = kBeamFlagKaon;
            sel.param_suffix = "K";
            sel.hist_suffix = "_K";
            sel.label = "K window (beam_flag=2)";
        } else if (arg == "--all") {
            sel.flag = kBeamFlagAll;
            sel.param_suffix = "K";
            sel.hist_suffix = "";
            sel.label = "all beam_flag";
        } else {
            return kFALSE;
        }
        return kTRUE;
    }

    // ____________________________________________________________________________________________
    GaussCore fit_core(TH1* h, Double_t mean0, Double_t sigma0, Int_t n_iter, Int_t min_entries) {
        GaussCore g;
        if (!h || h->GetEntries() < min_entries || !(sigma0 > 0.)) return g;
        Double_t m = mean0, s = sigma0;
        auto* f = new TF1(Form("g_%s", h->GetName()), "gaus", m - 1.5 * s, m + 1.5 * s);
        f->SetLineColor(kRed);
        for (Int_t it = 0; it < n_iter; ++it) {
            f->SetRange(m - 1.5 * s, m + 1.5 * s);
            f->SetParameters(h->GetMaximum(), m, s);
            h->Fit(f, "QR0");
            m = f->GetParameter(1);
            s = std::fabs(f->GetParameter(2));
        }
        g.mean = m;
        g.mean_err = f->GetParError(1);
        g.sigma = s;
        return g;
    }

    // ____________________________________________________________________________________________
    std::map<TString, DCLayer> load_dcgeo_blc(const TString& path) {
        std::map<TString, DCLayer> geo;
        std::ifstream ifs(path.Data());
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            std::string id, name;
            Double_t x, y, z, ta;
            if (!(iss >> id >> name >> x >> y >> z >> ta)) continue;
            if (name.rfind("BLC", 0) != 0 || name.find('-') == std::string::npos) continue;
            geo[name.c_str()] = DCLayer{z, ta};
        }
        return geo;
    }

    // ____________________________________________________________________________________________
    Bool_t chamber_mean_z(const std::map<TString, DCLayer>& geo, const TString& chamber, Double_t& z) {
        Double_t sum = 0.0;
        Int_t n = 0;
        for (const auto& [name, layer] : geo) {
            if (!name.BeginsWith(chamber + "-")) continue;
            sum += layer.z;
            ++n;
        }
        if (n == 0) return kFALSE;
        z = sum / n;
        return kTRUE;
    }

} // namespace d5_helper
