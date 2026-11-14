#include <cmath>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "mip_fit.h"

namespace ana_helper {
//_____________________________________________________________________________
bool
read_mip_settings(
    const std::string& path,
    const std::string& key,
    MipFitSettings& settings
)
{
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("MIP settings file is not readable: " + path);
    std::string line;
    while (std::getline(stream, line)) {
        std::istringstream row(line);
        std::string candidate;
        row >> candidate;
        if (candidate != key) continue;
        row >> settings.function >> settings.range_mode >> settings.lower >> settings.upper;
        for (auto& parameter : settings.parameters) row >> parameter;
        if (!row || settings.lower >= settings.upper || settings.parameters[0] <= 0.
            || settings.parameters[2] <= 0. || settings.parameters[3] < 0.
            || (settings.function != "auto" && settings.function != "gaussian" && settings.function != "landau")
            || (settings.range_mode != "auto" && settings.range_mode != "fixed"))
            throw std::runtime_error("invalid MIP fit override: " + key);
        for (auto parameter : settings.parameters)
            if (!std::isfinite(parameter)) throw std::runtime_error("nonfinite MIP parameter");
        if (!std::isfinite(settings.lower) || !std::isfinite(settings.upper))
            throw std::runtime_error("nonfinite MIP range");
        return true;
    }
    return false;
}

//_____________________________________________________________________________
bool
load_mip_settings(
    const std::string& key,
    MipFitSettings& settings
)
{
    const char* path = std::getenv("E72_MIP_FIT_OVERRIDE");
    return path && read_mip_settings(path, key, settings);
}

//_____________________________________________________________________________
MipFitState
describe_mip_fit(
    TF1* function,
    int status,
    int covariance
)
{
    MipFitState state;
    for (int i = 0; i < 4; ++i) state.parameters[i] = function->GetParameter(i);
    state.lower = function->GetXmin();
    state.upper = function->GetXmax();
    state.chi_square = function->GetChisquare();
    state.ndf = function->GetNDF();
    state.status = status;
    state.covariance = covariance;
    state.function = std::string(function->GetName()).find("landau") != std::string::npos ? "landau" : "gaussian";
    // 4000 reports only the optional IMPROVE failure.
    state.success = (status == 0 || status == 4000) && covariance >= 2 && state.ndf > 0
        && std::isfinite(state.chi_square) && state.parameters[0] > 0.
        && state.parameters[2] > 0. && state.parameters[1] >= state.lower
        && state.parameters[1] <= state.upper;
    for (auto parameter : state.parameters) state.success = state.success && std::isfinite(parameter);
    return state;
}

//_____________________________________________________________________________
TF1*
fit_manual_mip(
    TH1D* histogram,
    const MipFitSettings& settings,
    MipFitState& state
)
{
    auto gaussian = std::make_unique<TF1>(Form("manual_gauss_%s", histogram->GetName()), "gausn(0)+pol0(3)", settings.lower, settings.upper);
    auto landau = std::make_unique<TF1>(Form("manual_landau_%s", histogram->GetName()), "landaun(0)+pol0(3)", settings.lower, settings.upper);
    for (auto* function : {gaussian.get(), landau.get()}) {
        function->SetParameters(settings.parameters.data());
        function->SetParLimits(1, settings.lower, histogram->GetXaxis()->GetXmax());
        function->SetParLimits(3, 0., std::max(1., histogram->GetMaximum()));
    }
    double lower = settings.lower, upper = settings.upper;
    TF1* chosen = gaussian.get();
    MipFitState gaussian_state, landau_state;
    for (int round = 0; round < 2; ++round) {
        if (lower >= upper) throw std::runtime_error("invalid adaptive MIP range");
        for (auto* function : {gaussian.get(), landau.get()}) {
            if (settings.function == "gaussian" && function == landau.get()) continue;
            if (settings.function == "landau" && function == gaussian.get()) continue;
            function->SetRange(lower, upper);
            TFitResultPtr fitted = histogram->Fit(function, "0QEMRSN", "", lower, upper);
            auto described = describe_mip_fit(function, int(fitted), fitted.Get() ? fitted->CovMatrixStatus() : -1);
            if (function == gaussian.get()) gaussian_state = described;
            else landau_state = described;
        }
        if (settings.function == "landau") chosen = landau.get();
        else if (settings.function == "gaussian") chosen = gaussian.get();
        else if (gaussian_state.success != landau_state.success)
            chosen = gaussian_state.success ? gaussian.get() : landau.get();
        else chosen = gaussian->GetChisquare() <= landau->GetChisquare() ? gaussian.get() : landau.get();
        state = chosen == gaussian.get() ? gaussian_state : landau_state;
        if (settings.range_mode == "auto" && round == 0) {
            lower = std::max(settings.lower, state.parameters[1] - 1.8*state.parameters[2]);
            upper = std::min(histogram->GetXaxis()->GetXmax(), state.parameters[1] + 2.2*state.parameters[2]);
        }
    }
    chosen->SetLineColor(kOrange);
    chosen->SetLineWidth(2);
    chosen->SetNpx(1000);
    return chosen == gaussian.get() ? gaussian.release() : landau.release();
}
}
