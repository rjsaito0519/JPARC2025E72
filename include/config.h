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

    // void beam_initialize() {
    //     beam_generator = -1;
    // }
    
private:
    Config() = default; // コンストラクタをプライベートにして外部からのインスタンス生成を禁止
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
};

#endif
