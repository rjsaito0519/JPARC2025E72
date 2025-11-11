#ifndef PARAMS_H
#define PARAMS_H

#include <unordered_map>
#include <vector>
#include <string>
#include <Rtypes.h>

namespace param
{
    static const std::unordered_map<std::string, std::vector<Double_t>> hdprm_params{
        // key        mip_range_left
        // name-seg-ud
        { "bh2-0-u", { 160.0 } },
        { "bh2-0-d", { 200.0 } },
        { "bh2-1-u", { 200.0 } },
        { "bh2-1-d", { 190.0 } },
        { "bh2-2-u", { 180.0 } },
        { "bh2-2-d", { 200.0 } },
        { "bh2-3-u", { 220.0 } },
        { "bh2-3-d", { 190.0 } },
        { "bh2-4-u", { 180.0 } },
        { "bh2-4-d", { 180.0 } },
        { "bh2-5-u", { 200.0 } },
        { "bh2-5-d", { 180.0 } },
        { "bh2-6-u", { 160.0 } },
        { "bh2-6-d", { 170.0 } },
        { "bh2-7-u", { 200.0 } },
        { "bh2-7-d", { 180.0 } },
        { "bh2-8-u", { 200.0 } },
        { "bh2-8-d", { 170.0 } },
        { "bh2-9-u", { 200.0 } },
        { "bh2-9-d", { 180.0 } },
        { "bh2-10-u", { 220.0 } },
        { "bh2-10-d", { 180.0 } },
        { "bh2-11-u", { 170.0 } },
        { "bh2-11-d", { 180.0 } },
        { "bh2-12-u", { 200.0 } },
        { "bh2-12-d", { 190.0 } },
        { "bh2-13-u", { 200.0 } },
        { "bh2-13-d", { 200.0 } },
        { "bh2-14-u", { 210.0 } },
        { "bh2-14-d", { 210.0 } },

        { "htof-18-u", { 350.0 } },
        { "htof-18-d", { 460.0 } },
        { "htof-19-u", { 280.0 } },
        { "htof-19-d", { 280.0 } },
        { "htof-20-u", { 280.0 } },
        { "htof-20-d", { 280.0 } },
        { "htof-21-u", { 280.0 } },
        { "htof-21-d", { 280.0 } },
        
    };
}

#endif // PARAMS_H
