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
            || (settings.function != "auto" && settings.function != "gaussian" && settings.function != "landau"
                && settings.function != "mean")
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
MipStatistics
summarize_mip_statistics(
    TH1D* histogram,
    double lower,
    double upper,
    double qlo,
    double qhi
)
{
    MipStatistics summary;
    summary.window_lower = lower;
    summary.window_upper = upper;
    if (!(std::isfinite(lower) && std::isfinite(upper) && lower < upper
          && 0. <= qlo && qlo < qhi && qhi <= 1.)) return summary;
    const int first = std::max(1, histogram->FindBin(lower));
    const int last = std::min(histogram->GetNbinsX(), histogram->FindBin(upper));
    if (first > last) return summary;
    summary.entries = histogram->Integral(first, last);
    if (!(summary.entries > 0.) || !std::isfinite(summary.entries)) return summary;
    auto quantile = [&](double probability) {
        double cumulative = 0.;
        for (int bin = first; bin <= last; ++bin) {
            const double count = histogram->GetBinContent(bin);
            if (count > 0. && cumulative + count >= probability*summary.entries)
                return histogram->GetBinLowEdge(bin)
                    + histogram->GetBinWidth(bin)*(probability*summary.entries-cumulative)/count;
            cumulative += count;
        }
        return histogram->GetBinLowEdge(last+1);
    };
    summary.lower = qlo == 0. ? lower : quantile(qlo);
    summary.upper = qhi == 1. ? upper : quantile(qhi);
    summary.window_median = quantile(.5);
    double cumulative = 0., sum = 0., squares = 0.;
    summary.median = quantile(.5*(qlo+qhi));
    for (int bin = first; bin <= last; ++bin) {
        const double count = histogram->GetBinContent(bin);
        const double q_in = cumulative/summary.entries;
        const double q_out = (cumulative+count)/summary.entries;
        cumulative += count;
        const double fraction = std::min(q_out, qhi)-std::max(q_in, qlo);
        if (fraction <= 0.) continue;
        const double weight = fraction*summary.entries;
        const double center = histogram->GetBinCenter(bin);
        sum += weight*center;
        squares += weight*center*center;
        summary.selected_entries += weight;
    }
    if (summary.selected_entries > 0.) {
        summary.mean = sum/summary.selected_entries;
        const double variance = std::max(0., squares/summary.selected_entries-summary.mean*summary.mean);
        summary.error = std::sqrt(variance/summary.selected_entries);
        summary.valid = std::isfinite(summary.mean) && std::isfinite(summary.median);
    }
    return summary;
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
    for (int i = 0; i < std::min(4, function->GetNpar()); ++i) state.parameters[i] = function->GetParameter(i);
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
fit_manual_t0(
    TH1D* histogram,
    const MipFitSettings& settings,
    MipFitState& state
)
{
    if (settings.function != "gaussian" || settings.parameters[3] != 0.)
        throw std::runtime_error("T0 uses a three-parameter Gaussian");
    const double minimum = histogram->GetXaxis()->GetXmin();
    const double maximum = histogram->GetXaxis()->GetXmax();
    if (settings.lower < minimum || settings.upper > maximum)
        throw std::runtime_error("T0 fit range is outside the histogram");
    auto function = std::make_unique<TF1>(Form("manual_t0_gauss_%s", histogram->GetName()),
        "gausn", settings.lower, settings.upper);
    function->SetParameters(settings.parameters.data());
    function->SetParLimits(1, settings.lower, settings.upper);
    function->SetParLimits(2, histogram->GetBinWidth(1)*.01, maximum-minimum);
    double lower = settings.lower, upper = settings.upper;
    histogram->GetXaxis()->UnZoom();
    for (int round = 0; round < 2; ++round) {
        if (lower >= upper) throw std::runtime_error("invalid T0 adaptive fit range");
        function->SetRange(lower, upper);
        const TString option = histogram->GetMaximum() > 50. ? "0QEMRSN" : "0QEMRLSN";
        TFitResultPtr fitted = histogram->Fit(function.get(), option, "", lower, upper);
        state = describe_mip_fit(function.get(), int(fitted), fitted.Get() ? fitted->CovMatrixStatus() : -1);
        if (round == 0 && settings.range_mode == "auto") {
            lower = std::max(settings.lower, state.parameters[1]-2.*state.parameters[2]);
            upper = std::min(settings.upper, state.parameters[1]+2.*state.parameters[2]);
        }
    }
    function->SetLineColor(kOrange);
    function->SetLineWidth(2);
    function->SetNpx(1000);
    return function.release();
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
