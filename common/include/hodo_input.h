#ifndef HODO_INPUT_H
#define HODO_INPUT_H

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <TString.h>
#include <TFile.h>
#include <TObjString.h>

#include "config.h"

namespace hodo_input {
//_____________________________________________________________________________
inline bool
all_histograms()
{
    const char* selection = std::getenv("E72_CALIB_HIST_SELECTION");
    if (!selection || !*selection) return false;
    if (std::string(selection) != "All")
        throw std::runtime_error("Unsupported calibration histogram selection");
    return true;
}

//_____________________________________________________________________________
inline TString
histogram_suffix(
    const TString& parameter_suffix
)
{
    return all_histograms() ? TString("") : "_" + parameter_suffix;
}

//_____________________________________________________________________________
inline void
configure_phc(
    Config& settings
)
{
    const char* range = std::getenv("E72_CALIB_PHC_TIME_WINDOW");
    if (!all_histograms() || !range || !*range) return;
    std::istringstream values(range);
    double lower, upper;
    std::string extra;
    if (!(values >> lower >> upper) || (values >> extra)
        || !std::isfinite(lower) || !std::isfinite(upper) || lower >= upper)
        throw std::runtime_error("Invalid calibration PHC time window");
    settings.phc_time_window[settings.detector.Data()] = {lower, upper};
}
//_____________________________________________________________________________
inline std::pair<double, double>
t0_window()
{
    const char* range = std::getenv("E72_CALIB_T0_TIME_WINDOW");
    if (!range || !*range)
        throw std::runtime_error("All T0 requires an explicit raw-time fit window");
    std::istringstream values(range);
    double lower, upper;
    std::string extra;
    if (!(values >> lower >> upper) || (values >> extra)
        || !std::isfinite(lower) || !std::isfinite(upper) || lower >= upper)
        throw std::runtime_error("Invalid calibration T0 time window");
    return {lower, upper};
}

//_____________________________________________________________________________
inline void
write_metadata(
    TFile* output,
    const TString& parameter_suffix
)
{
    if (!all_histograms()) return;
    output->cd();
    TObjString("All").Write("calib_input_selection");
    TObjString(parameter_suffix).Write("calib_parameter_particle");
}

}
#endif
