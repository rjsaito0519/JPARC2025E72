// config.h
#ifndef CONFIG_H
#define CONFIG_H

// c++
#include <vector>
#include <string>
#include <unordered_map>
#include <utility>

// ROOT
#include "TString.h" 
#include "Rtypes.h"

// -- fitting result container -----
struct FitResult {
    std::vector<Double_t> par;
    std::vector<Double_t> err;
    Double_t chi_square;
    Int_t ndf;
    Int_t migrad_stats;
    std::vector<Double_t> additional; // 何か追加で値を返したいときのコンテナ
};

class Config {
public:
    static Config& getInstance() {
        static Config instance;
        return instance;
    }

    const Int_t max_tdc_hit = 16;

    TString detector = "none";
    Double_t adc_ped_remove_nsigma = 15.0;
    Double_t hdprm_mip_range_left  = -1.0;
    FitResult hdprm_typical_value;

    const std::pair<Int_t, Int_t> htof_adc_exist_seg = {0, 33};

    const std::unordered_map<std::string, Int_t> num_of_ch{
        { "bht", 63  },
        {  "t0",  5  },
        { "bac",  4  },
        { "kvc",  8  },
        { "bh2",  15 },
        { "htof", 34 },
        { "cvc",   8 },
        { "sac3",  1 },
        { "sfv",   1 },

        { "blc",  8 },
    };

    // 検出器ごとの MaxDriftTime [ns]
    std::unordered_map<std::string, Double_t> max_drift_time = {
        {"BLC1a", 200.0},
        {"BLC1b", 200.0},
        {"BLC2a", 200.0},
        {"BLC2b", 200.0},
    };

    // 検出器ごとの MaxDriftLength [mm]
    std::unordered_map<std::string, Double_t> max_drift_length = {
        {"BLC1a", 4.0},
        {"BLC1b", 4.0},
        {"BLC2a", 2.5},
        {"BLC2b", 2.5},
    };

    std::unordered_map<std::string, std::pair<Double_t, Double_t>> tdc_search_range{
        { "bht", {0.0, 0.0} },
        {  "t0", {0.0, 0.0} },
        { "bh2", {0.0, 0.0} },
        { "htof", {0.0, 0.0} },
        {  "cvc", {0.0, 0.0} },
        { "sac3", {0.0, 0.0} },
        {  "sfv", {0.0, 0.0} },

        // { "BAC",  4 },
        // { "SAC",  8 },
        // { "KVC",  4 },
        // { "BH2",  1 },
    };

    std::unordered_map<std::string, std::pair<Double_t, Double_t>> phc_time_window{
        { "bht",  {-4.0, 2.0} },
        {  "t0",  {-2.5, 1.5} },
        { "bh2",  {-2.0, 1.0} },
        { "htof", {-2.0, 1.0} },
        { "cvc",  {-4.0, 2.0} },
        
        // { "BAC",  4 },
        // { "SAC",  8 },
        // { "KVC",  4 },
        // { "BH2",  1 },
    };

    Double_t phc_de_range_min = 0.0;
    Double_t phc_de_range_max = 5.0;
    std::unordered_map<std::string, Double_t> phc_de_range_ratio{
        { "bht", 0.001 },
        {  "t0", 0.01 },
        { "bh2", 0.001 },
        { "htof", 0.01 },
        { "cvc", 0.005 },
        
        // { "BAC",  4 },
        // { "SAC",  8 },
        // { "KVC",  4 },
        // { "BH2",  1 },
    };

    

    // void beam_initialize() {
    //     beam_generator = -1;
    // }
    
private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
