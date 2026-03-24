// Inspect TpcPhase_Fit_Cobo%d functions in a TpcPhase.root file
// and print basic step parameters (N, C, A0, t0, w0) for each CoBo.
//
// Usage (from ROOT):
//   .L inspect_tpc_phase_steps.C+
//   inspect_tpc_phase_steps("param/TPCPHASE/TpcPhase_02601.root");
//
// 出力は 1 行/CoBo で、run ごとの代表振幅や位置分布の確認に使うことを想定。

#include <TFile.h>
#include <TF1.h>
#include <TString.h>
#include <iostream>

void inspect_tpc_phase_steps(const char* phase_path,
                             Int_t max_cobo = 8)
{
  TFile* f = TFile::Open(phase_path, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << phase_path << std::endl;
    return;
  }

  std::cout << "# CoBo  N  C[ns]  A0[ns]  t0[ns]  w0[ns]" << std::endl;

  for (Int_t cobo = 0; cobo < max_cobo; ++cobo) {
    TString name = Form("TpcPhase_Fit_Cobo%d", cobo);
    TF1* f_fit = dynamic_cast<TF1*>(f->Get(name));
    if (!f_fit) {
      std::cout << cobo << "  -  (no function)" << std::endl;
      continue;
    }

    const Double_t N   = f_fit->GetParameter(0);
    const Double_t C   = f_fit->GetParameter(1);
    Double_t A0 = 0.0, t0 = 0.0, w0 = 0.0;
    if (N >= 1.0) {
      A0 = f_fit->GetParameter(2);
      t0 = f_fit->GetParameter(3);
      w0 = f_fit->GetParameter(4);
    }

    std::cout << cobo
              << "  " << N
              << "  " << C
              << "  " << A0
              << "  " << t0
              << "  " << w0
              << std::endl;
  }

  f->Close();
  delete f;
}

