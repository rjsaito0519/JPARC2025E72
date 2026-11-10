#include "ana_helper.h"

#include <cstdlib>
#include <sstream>
#include <sys/stat.h>
#include <iostream>
#include "paths.h"
#include <TObjArray.h>
#include <TObjString.h>

namespace ana_helper {

    namespace {
        bool g_fit_override_loaded = false;
        std::unordered_map<std::string, Double_t> g_mip_left_override;
    }

    // ____________________________________________________________________________________________
    void load_fit_override_file()
    {
        if (g_fit_override_loaded) {
            return;
        }
        g_fit_override_loaded = true;
        g_mip_left_override.clear();
        const char* path = std::getenv("E72_FIT_OVERRIDE");
        if (!path || !*path) {
            return;
        }
        std::ifstream in(path);
        if (!in) {
            std::cerr << "[ana_helper] E72_FIT_OVERRIDE not readable: "
                      << path << std::endl;
            return;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            std::istringstream iss(line);
            std::string key;
            Double_t value = 0.0;
            if (!(iss >> key >> value)) {
                continue;
            }
            g_mip_left_override[key] = value;
        }
        std::cerr << "[ana_helper] loaded " << g_mip_left_override.size()
                  << " fit override(s) from " << path << std::endl;
    }

    // ____________________________________________________________________________________________
    Double_t hdprm_mip_left_override(const std::string& key, Double_t fallback)
    {
        load_fit_override_file();
        auto it = g_mip_left_override.find(key);
        if (it == g_mip_left_override.end()) {
            return fallback;
        }
        return it->second;
    }

    // ____________________________________________________________________________________________
    TString get_run_dir(const TString& base_dir, Int_t run_num) {
        TString dir = Form("%s/run%05d", base_dir.Data(), run_num);
        struct stat st = {0};
        if (stat(dir.Data(), &st) == -1) {
            mkdir(dir.Data(), 0755);
        }
        return dir;
    }

    TString get_img_dir(const TString& base_output_dir, Int_t run_num) {
        TString img_base = Form("%s/img", base_output_dir.Data());
        struct stat st = {0};
        if (stat(img_base.Data(), &st) == -1) {
            mkdir(img_base.Data(), 0755);
        }
        return get_run_dir(img_base, run_num);
    }

    TString get_root_dir(const TString& base_output_dir, Int_t run_num) {
        TString root_base = Form("%s/root", base_output_dir.Data());
        struct stat st = {0};
        if (stat(root_base.Data(), &st) == -1) {
            mkdir(root_base.Data(), 0755);
        }
        return get_run_dir(root_base, run_num);
    }
    
    TString get_scratch_dir(const TString& base_output_dir, Int_t run_num) {
        TString scratch_base = Form("%s/scratch_root", base_output_dir.Data());
        struct stat st = {0};
        if (stat(scratch_base.Data(), &st) == -1) {
            mkdir(scratch_base.Data(), 0755);
        }
        return get_run_dir(scratch_base, run_num);
    }

    // ____________________________________________________________________________________________
    TCanvas* add_tab(TGTab *tab, const char* tabName) {
        // タブを作成し、キャンバスを埋め込む
        TGCompositeFrame *tf = tab->AddTab(tabName);
        TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas(tabName, tf, 1000, 800);
       tf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        return embeddedCanvas->GetCanvas();
    }
    // ____________________________________________________________________________________________
    std::pair<Double_t, Double_t> get_user_param(const TString& key, Int_t run_num) {
        // Try to locate UserParam file
         TString param_file = Form("param/USER/e72/UserParam_run%05d", run_num);
         struct stat st;
         // Try path relative to Project Root (CWD)
         if (stat(param_file.Data(), &st) != 0) {
              // Try relative to bin/
              param_file = Form("../param/USER/e72/UserParam_run%05d", run_num);
              if (stat(param_file.Data(), &st) != 0) {
                   // Try the analyzer tree (works from any CWD)
                   param_file = Form("%s/param/USER/e72/UserParam_run%05d", ANALYZER_DIR.Data(), run_num);
                   if (stat(param_file.Data(), &st) != 0) {
                        std::cerr << "Warning: no UserParam_run" << Form("%05d", run_num)
                                  << " found; " << key << " not taken from USER" << std::endl;
                        return {0.0, 0.0};
                   }
              }
         }

         std::ifstream ifs(param_file.Data());
         if (!ifs.is_open()) return {0.0, 0.0};

         std::string line;
         while (std::getline(ifs, line)) {
             TString tline = line;
             if (tline.BeginsWith(key)) {
                 // Found key
                 // Format: KEY min max
                 // e.g. BHT_TDC 1400 1500
                 TObjArray *tokens = tline.Tokenize(" \t"); // split by space or tab
                 if (tokens->GetEntries() >= 3) {
                     TString s_min = ((TObjString*)tokens->At(1))->GetString();
                     TString s_max = ((TObjString*)tokens->At(2))->GetString();
                     delete tokens;
                     return {s_min.Atof(), s_max.Atof()};
                 }
                 delete tokens;
             }
         }
         return {0.0, 0.0};
    }

}