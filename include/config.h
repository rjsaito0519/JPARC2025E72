// config.h
#ifndef CONFIG_H
#define CONFIG_H

class Config {
public:
    static Config& getInstance() {
        static Config instance;
        return instance;
    }

    const Int_t max_tdc_hit = 16;

    const std::unordered_map<std::string, Int_t> num_of_ch{
        { "bht", 63 },
        {  "t0",  5 },
        { "BAC",  4 },
        { "SAC",  8 },
        { "KVC",  4 },
        { "BH2",  1 },
    };

    const std::vector<std::vector<Double_t>> t0_ped_mip_distance{
        {77.0, 52.0, 75.0, 50.0, 71.0}, // UP
        {55.0, 35.0, 55.0, 53.0, 45.0}  // DOWN
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
