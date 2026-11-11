// Draw unsuffixed CBtof0 (all events) as hdphc QA book front matter (1 page).
#include <iostream>
#include <sys/stat.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "ana_helper.h"
#include "paths.h"

//_____________________________________________________________________________
static void
ApplyPhcStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.06, "XY");
  gStyle->SetTitleSize(0.06, "x");
  gStyle->SetTitleSize(0.06, "y");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gROOT->GetColor(0)->SetAlpha(0.01);
}

//_____________________________________________________________________________
static Int_t
ReadRunNumber(TFile* f)
{
  if (!f) {
    return -1;
  }
  TTreeReader reader("hodo", f);
  if (!reader.GetTree()) {
    return -1;
  }
  TTreeReaderValue<unsigned int> run_number(reader, "run_number");
  if (!reader.Next()) {
    return -1;
  }
  return static_cast<Int_t>(*run_number);
}

//_____________________________________________________________________________
static Bool_t
DrawPage(
  TCanvas* c,
  const TString& pdf_path,
  TObject* obj,
  const TString& title
)
{
  if (!obj) {
    return kFALSE;
  }
  c->cd();
  c->Clear();
  if (obj->InheritsFrom(TH1::Class())) {
    auto* h1 = static_cast<TH1*>(obj);
    h1->SetTitle(title);
    h1->SetLineWidth(2);
    h1->SetFillStyle(0);
    h1->Draw("HIST");
  } else {
    return kFALSE;
  }
  c->Print(pdf_path);
  return kTRUE;
}

//_____________________________________________________________________________
static Int_t
Analyze(TString path, Int_t run_override)
{
  ApplyPhcStyle();

  auto* f = TFile::Open(path.Data(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: Could not open file : " << path << std::endl;
    return 1;
  }

  Int_t run_num = run_override > 0 ? run_override : ReadRunNumber(f);
  if (run_num < 1) {
    std::cerr << "Error: cannot determine run number" << std::endl;
    f->Close();
    return 1;
  }

  // All-events hist (no _Pi / _K suffix).
  TH1* h_cbtof0 = dynamic_cast<TH1*>(f->Get("CBtof0"));
  if (!h_cbtof0) {
    std::cerr << "Error: missing CBtof0 in " << path << std::endl;
    f->Close();
    return 1;
  }

  TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
  TString pdf_path = Form(
    "%s/run%05d_hodo_btof0.pdf",
    img_base_dir.Data(), run_num
  );

  auto* c = new TCanvas("cbtof0_qa", "", 1500, 1200);
  c->Print(pdf_path + "[");
  const Bool_t drawn = DrawPage(
    c, pdf_path, h_cbtof0,
    Form("run%05d CBtof0 (all)", run_num)
  );
  c->Print(pdf_path + "]");
  delete c;
  f->Close();

  if (!drawn) {
    std::cerr << "Error: no CBtof0 pages drawn" << std::endl;
    return 1;
  }
  std::cout << "wrote " << pdf_path << " (1 pages)" << std::endl;
  return 0;
}

//_____________________________________________________________________________
Int_t
main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <Hodo.root> [run]" << std::endl;
    return 1;
  }
  TString path = argv[1];
  Int_t run_override = -1;
  if (argc >= 3) {
    TString arg2 = argv[2];
    // Legacy: <Hodo.root> <Pi|K> [run] — ignore particle, take run if present.
    if (arg2 == "Pi" || arg2 == "K") {
      if (argc >= 4) {
        run_override = TString(argv[3]).Atoi();
      }
    } else {
      run_override = arg2.Atoi();
    }
  }
  return Analyze(path, run_override);
}
