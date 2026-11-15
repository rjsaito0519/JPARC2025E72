#ifndef MIP_FIT_H
#define MIP_FIT_H

#include <array>
#include <string>

#include <TH1.h>
#include <TF1.h>

namespace ana_helper {
struct MipFitSettings {
    std::array<double, 4> parameters{{0., 0., 1., 0.}};
    double lower = 0.;
    double upper = 0.;
    std::string function = "auto";
    std::string range_mode = "auto";
};

struct MipStatistics {
    double window_lower = 0.;
    double window_upper = 0.;
    double lower = 0.;
    double upper = 0.;
    double mean = 0.;
    double median = 0.;
    double window_median = 0.;
    double entries = 0.;
    double selected_entries = 0.;
    double error = 0.;
    bool valid = false;
};

struct MipFitState {
    std::array<double, 4> parameters{{0., 0., 1., 0.}};
    double lower = 0.;
    double upper = 0.;
    double chi_square = 0.;
    int ndf = 0;
    int status = -1;
    int covariance = -1;
    std::string function = "gaussian";
    bool success = false;
    MipStatistics statistics;
};

MipStatistics summarize_mip_statistics(TH1D* histogram, double lower, double upper, double qlo = 0., double qhi = 1.);
bool read_mip_settings(const std::string& path, const std::string& key, MipFitSettings& settings);
bool load_mip_settings(const std::string& key, MipFitSettings& settings);
MipFitState describe_mip_fit(TF1* function, int status, int covariance);
TF1* fit_manual_mip(TH1D* histogram, const MipFitSettings& settings, MipFitState& state);
}
#endif
