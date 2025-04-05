#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    TCanvas* add_tab(TGTab *tab, const char* tabName) {
        // タブを作成し、キャンバスを埋め込む
        TGCompositeFrame *tf = tab->AddTab(tabName);
        TRootEmbeddedCanvas *embeddedCanvas = new TRootEmbeddedCanvas(tabName, tf, 1000, 800);
       tf->AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        return embeddedCanvas->GetCanvas();
    }

    // ____________________________________________________________________________________________
    std::vector<std::vector<Double_t>> load_data(TString path) {
        std::ifstream file(path.Data());
        std::string line;

        std::vector<std::vector<Double_t>> data;

        // Skip the header line
        std::getline(file, line);

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string value;

            std::vector<Double_t> buf;
            while (std::getline(ss, value, ',')) {
                buf.push_back(std::stod(value));
            }
            data.push_back(buf);
        }
        return data;
    }

}