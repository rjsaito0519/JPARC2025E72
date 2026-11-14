#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

#include <TFile.h>
#include <TROOT.h>

#include "ana_helper.h"

//_____________________________________________________________________________
int
main(
    int argc,
    char** argv
)
{
    if (argc == 5 && std::string(argv[1]) == "--curve") {
        try {
            ana_helper::MipFitSettings settings;
            if (!ana_helper::read_mip_settings(argv[2], argv[3], settings)) return 1;
            TF1 curve("initial_curve", settings.function == "landau" ? "landaun(0)+pol0(3)" : "gausn(0)+pol0(3)", settings.lower, settings.upper);
            curve.SetParameters(settings.parameters.data());
            std::ofstream output(argv[4]);
            output << std::setprecision(17) << "{\"curve\":[";
            for (int i = 0; i <= 300; ++i) {
                double x = settings.lower + (settings.upper-settings.lower)*i/300.;
                double y = curve.Eval(x);
                if (i) output << ',';
                output << '[' << x << ',' << (std::isfinite(y) ? y : 0.) << ']';
            }
            output << "]}";
            return output ? 0 : 1;
        } catch (const std::exception& error) {
            std::cerr << error.what() << std::endl;
            return 1;
        }
    }
    if (argc < 6 || argc > 7) return 1;
    try {
        gROOT->SetBatch(kTRUE);
        TH1::AddDirectory(kFALSE);
        TFile input(argv[1], "READ");
        const std::string particle = argv[2];
        const int segment = std::stoi(argv[3]);
        const std::string side = argv[4];
        if (input.IsZombie() || (particle != "Pi" && particle != "K")
            || segment < 0 || segment >= 15 || (side != "U" && side != "D"))
            throw std::runtime_error("invalid preview input");
        const std::string name = "BH2_ADC_seg" + std::to_string(segment) + side + "_" + particle;
        auto* source = dynamic_cast<TH1D*>(input.Get(name.c_str()));
        if (!source) throw std::runtime_error("missing histogram: " + name);
        std::unique_ptr<TH1D> histogram(static_cast<TH1D*>(source->Clone("mip_preview")));
        auto& conf = Config::getInstance();
        conf.detector = "bh2";
        const std::string key = "bh2-" + std::to_string(segment) + "-" + (side == "U" ? "u" : "d");
        conf.hdprm_mip_range_left = ana_helper::hdprm_mip_left_override(key,
            param::hdprm_params.count(key) ? param::hdprm_params.at(key)[0] : -1.);
        ana_helper::MipFitSettings settings;
        if (argc == 7 && !ana_helper::read_mip_settings(argv[6], key + ":" + particle, settings))
            throw std::runtime_error("missing preview settings");
        ana_helper::MipFitState state;
        TCanvas canvas("mip_canvas", "", 800, 600);
        if (argc == 6) {
            const std::string typical_name = "BH2_ADC_seg6U_" + particle;
            auto* typical_source = dynamic_cast<TH1D*>(input.Get(typical_name.c_str()));
            if (!typical_source) throw std::runtime_error("missing typical BH2 histogram");
            std::unique_ptr<TH1D> typical(static_cast<TH1D*>(typical_source->Clone("mip_typical")));
            const double channel_left = conf.hdprm_mip_range_left;
            const std::string typical_key = "bh2-6-u";
            conf.hdprm_mip_range_left = ana_helper::hdprm_mip_left_override(typical_key,
                param::hdprm_params.count(typical_key) ? param::hdprm_params.at(typical_key)[0] : -1.);
            conf.hdprm_typical_value = ana_helper::adc_fit(typical.get(), &canvas, 1);
            conf.hdprm_mip_range_left = channel_left;
        }
        ana_helper::adc_fit(histogram.get(), &canvas, 1, 2, argc == 7 ? &settings : nullptr, &state);
        std::ofstream output(argv[5]);
        if (!output) throw std::runtime_error("cannot write preview result");
        output << std::setprecision(17);
        output << "{\"success\":" << (state.success ? "true" : "false")
            << ",\"status\":" << state.status << ",\"covariance\":" << state.covariance
            << ",\"ndf\":" << state.ndf << ",\"function\":\"" << state.function << "\",\"parameters\":[";
        for (int i = 0; i < 4; ++i) {
            if (i) output << ',';
            output << (std::isfinite(state.parameters[i]) ? state.parameters[i] : 0.);
        }
        output << "],\"range\":[" << state.lower << ',' << state.upper << "],\"curve\":[";
        TF1 curve("preview_curve", state.function == "landau" ? "landaun(0)+pol0(3)" : "gausn(0)+pol0(3)", state.lower, state.upper);
        curve.SetParameters(state.parameters.data());
        for (int i = 0; i <= 300; ++i) {
            const double x = state.lower + (state.upper-state.lower)*i/300.;
            const double y = curve.Eval(x);
            if (i) output << ',';
            output << '[' << x << ',' << (std::isfinite(y) ? y : 0.) << ']';
        }
        output << "]}";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
