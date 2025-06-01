// config.h
#ifndef CONFIG_H
#define CONFIG_H

class Config {
public:
    static Config& getInstance() {
        static Config instance;
        return instance;
    }

    Double_t HRTDC_factor = -0.000939002;

    TString detector = "none";

    Double_t t0_adc_ped_remove_nsigma = 15.0;

    const Int_t max_tdc_hit = 16;

    const std::unordered_map<std::string, Int_t> num_of_ch{
        { "bht", 63 },
        {  "t0",  5 },
        { "BAC",  4 },
        { "SAC",  8 },
        { "KVC",  4 },
        { "bh2",  11 },
    };

    std::unordered_map<std::string, std::pair<Double_t, Double_t>> tdc_search_range{
        { "bht", {0.0, 0.0} },
        {  "t0", {0.0, 0.0} },
        { "bh2", {0.0, 0.0} },
        
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
