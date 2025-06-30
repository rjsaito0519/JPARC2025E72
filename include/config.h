// config.h
#ifndef CONFIG_H
#define CONFIG_H

class Config {
public:
    static Config& getInstance() {
        static Config instance;
        return instance;
    }

    const Double_t HRTDC_factor = -0.000939002;
    const Int_t max_tdc_hit = 16;

    TString detector = "none";
    Double_t adc_ped_remove_nsigma = 15.0;
    Double_t hdprm_mip_range_left = -1.0;

    const std::pair<Int_t, Int_t> htof_adc_exist_seg = {18, 21};

    const std::unordered_map<std::string, Int_t> num_of_ch{
        { "bht", 63 },
        {  "t0",  5 },
        { "BAC",  4 },
        { "SAC",  8 },
        { "KVC",  4 },
        { "bh2",  11 },
        { "htof",  34 },
    };

    std::unordered_map<std::string, std::pair<Double_t, Double_t>> tdc_search_range{
        { "bht", {0.0, 0.0} },
        {  "t0", {0.0, 0.0} },
        { "bh2", {0.0, 0.0} },
        { "htof", {0.0, 0.0} },
        
        // { "BAC",  4 },
        // { "SAC",  8 },
        // { "KVC",  4 },
        // { "BH2",  1 },
    };

    std::unordered_map<std::string, std::pair<Double_t, Double_t>> phc_time_window{
        { "bht", {-4.0, 2.0} },
        {  "t0", {-2.5, 1.5} },
        { "bh2", {-2.0, 1.0} },
        { "htof", {-2.0, 1.0} },
        
        // { "BAC",  4 },
        // { "SAC",  8 },
        // { "KVC",  4 },
        // { "BH2",  1 },
    };

    Double_t phc_de_range_min = 0.0;
    std::unordered_map<std::string, Double_t> phc_de_range_ratio{
        { "bht", 0.001 },
        {  "t0", 0.01 },
        { "bh2", 0.001 },
        { "htof", 0.01 },
        
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
