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
        { "bh2-0-u", { 250.0 } },
        { "bh2-0-d", { 250.0 } },
        { "bh2-1-u", { 250.0 } },
        { "bh2-1-d", { 220.0 } },
        { "bh2-2-u", { 250.0 } },
        { "bh2-2-d", { 250.0 } },
        { "bh2-3-u", { 220.0 } },
        { "bh2-3-d", { 220.0 } },
        { "bh2-4-u", { 220.0 } },
        { "bh2-4-d", { 190.0 } },
        { "bh2-5-u", { 220.0 } },
        { "bh2-5-d", { 190.0 } },
        { "bh2-6-u", { 190.0 } },
        { "bh2-6-d", { 190.0 } },
        { "bh2-7-u", { 240.0 } },
        { "bh2-7-d", { 230.0 } },
        { "bh2-8-u", { 250.0 } },
        { "bh2-8-d", { 260.0 } },
        { "bh2-9-u", { 220.0 } },
        { "bh2-9-d", { 220.0 } },
        { "bh2-10-u", { 230.0 } },
        { "bh2-10-d", { 220.0 } },
        
    };
}

#endif // PARAMS_H
