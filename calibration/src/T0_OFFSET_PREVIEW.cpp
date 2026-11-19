#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

#include <TFile.h>
#include <TROOT.h>

#include "ana_helper.h"
#include "hodo_input.h"

//_____________________________________________________________________________
int
main(
    int argc,
    char** argv
)
{
    try {
        gROOT->SetBatch(kTRUE);
        TH1::AddDirectory(kFALSE);
        const bool curve_only = argc == 5 && std::string(argv[1]) == "--curve";
        ana_helper::MipFitSettings settings;
        ana_helper::MipFitState state;
        std::string output_path;
        if (curve_only) {
            if (!ana_helper::read_mip_settings(argv[2], argv[3], settings)
                || settings.function != "gaussian" || settings.parameters[3] != 0.)
                throw std::runtime_error("invalid T0 Gaussian settings");
            state.parameters = settings.parameters;
            state.lower = settings.lower;
            state.upper = settings.upper;
            state.status = 0;
            output_path = argv[4];
        } else {
            if (argc != 6 && argc != 7) throw std::runtime_error("invalid T0 preview arguments");
            TFile input(argv[1], "READ");
            const std::string particle = argv[2];
        const std::string input_suffix = hodo_input::histogram_suffix(particle.c_str()).Data();
            const int segment = std::stoi(argv[3]);
            if (input.IsZombie() || (particle != "Pi" && particle != "K") || segment < 0 || segment >= 15
                || std::string(argv[4]) != "T0") throw std::runtime_error("invalid T0 preview input");
            const std::string name = "T0_seg"+std::to_string(segment)+"_TimeOffset"+input_suffix;
            auto* source = dynamic_cast<TH1D*>(input.Get(name.c_str()));
            if (!source) throw std::runtime_error("missing histogram: "+name);
            std::unique_ptr<TH1D> histogram(static_cast<TH1D*>(source->Clone("t0_preview")));
            auto& conf = Config::getInstance();
            conf.detector = "bh2";
            const std::string key = "bh2-"+std::to_string(segment)+"-t0:"+particle;
            if (argc == 7 && !ana_helper::read_mip_settings(argv[6], key, settings))
                throw std::runtime_error("missing T0 preview settings");
            std::unique_ptr<TH1D> sum(static_cast<TH1D*>(source->Clone("t0_sum")));
            sum->Reset();
            for (int channel = 0; channel < 15; ++channel) {
                auto* hist = dynamic_cast<TH1D*>(input.Get(Form("T0_seg%d_TimeOffset%s", channel, input_suffix.c_str())));
                if (!hist) throw std::runtime_error("missing T0 input segment");
                sum->Add(hist);
            }
            ana_helper::set_tdc_search_range(sum.get());
            TCanvas canvas("t0_preview_canvas", "", 800, 600);
            ana_helper::t0_offset_fit(histogram.get(), &canvas, 1, argc == 7 ? &settings : nullptr, &state);
            output_path = argv[5];
        }
        std::ofstream output(output_path);
        if (!output) throw std::runtime_error("cannot write T0 preview");
        output << std::setprecision(17);
        output << "{\"function\":\"gaussian\",\"success\":" << (state.success ? "true" : "false")
            << ",\"status\":" << state.status << ",\"covariance\":" << state.covariance
            << ",\"ndf\":" << state.ndf << ",\"parameters\":[";
        for (int parameter = 0; parameter < 3; ++parameter) {
            if (parameter) output << ',';
            output << (std::isfinite(state.parameters[parameter]) ? state.parameters[parameter] : 0.);
        }
        output << "],\"range\":[" << state.lower << ',' << state.upper << "],\"curve\":[";
        TF1 curve("t0_curve", "gausn", state.lower, state.upper);
        curve.SetParameters(state.parameters.data());
        for (int point = 0; state.status >= 0 && point <= 300; ++point) {
            const double x = state.lower+(state.upper-state.lower)*point/300.;
            const double y = curve.Eval(x);
            if (point) output << ',';
            output << '[' << x << ',' << (std::isfinite(y) ? y : 0.) << ']';
        }
        output << "]}";
        return output ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
