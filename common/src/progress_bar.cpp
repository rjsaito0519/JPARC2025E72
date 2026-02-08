#include <iostream>
#include <iomanip>
#include <Rtypes.h>

void displayProgressBar(Int_t current, Int_t total) {
    static Int_t lastPercent = -1;
    Double_t progress = (Double_t)current / total;
    Int_t percent = static_cast<Int_t>(progress * 100);

    // Update only when the percentage changes
    if (percent != lastPercent) {
        lastPercent = percent;
        Int_t barWidth = 50;  // Width of the progress bar

        std::cout << "[";
        Int_t pos = barWidth * progress;
        for (Int_t i = 0; i < barWidth; ++i) {
            if (i < pos) {
                std::cout << "=";
            } else if (i == pos) {
                std::cout << ">";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << " %\r";
        std::cout.flush();  // Flush the output to display immediately
        if (percent == 100) {
            std::cout << std::endl;
        }
    }
}
