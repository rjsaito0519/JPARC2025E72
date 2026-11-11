// Draw CBtof0 (+ related 2D) pages for hdphc QA book front matter.
#include <iostream>
#include <sys/stat.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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
  if (obj->InheritsFrom(TH2::Class())) {
    auto* h2 = static_cast<TH2*>(obj);
    h2->SetTitle(title);
    h2->Draw("colz");
  } else if (obj->InheritsFrom(TH1::Class())) {
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
Analyze(TString path, TString particle, Int_t run_override)
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

  const TString tag = "_" + particle;
  TH1* h_cbtof0 = dynamic_cast<TH1*>(f->Get("CBtof0" + tag));
  if (!h_cbtof0) {
    std::cerr << "Error: missing CBtof0" << tag << " in " << path << std::endl;
    f->Close();
    return 1;
  }
  TH2* h_vs_t0 = dynamic_cast<TH2*>(f->Get("CBtof0_vs_deT0Seg" + tag));
  TH2* h_vs_btof = dynamic_cast<TH2*>(f->Get("CBtof0_vs_deBtof0Seg" + tag));

  TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
  TString pdf_path = Form(
    "%s/run%05d_hodo_btof0_%s.pdf",
    img_base_dir.Data(), run_num, particle.Data()
  );

  auto* c = new TCanvas("cbtof0_qa", "", 1500, 1200);
  c->Print(pdf_path + "[");
  Int_t pages = 0;
  if (DrawPage(
        c, pdf_path, h_cbtof0,
        Form("run%05d CBtof0 (%s)", run_num, particle.Data())
      )) {
    ++pages;
  }
  if (DrawPage(
        c, pdf_path, h_vs_t0,
        Form("run%05d CBtof0 vs dE(T0) (%s)", run_num, particle.Data())
      )) {
    ++pages;
  }
  if (DrawPage(
        c, pdf_path, h_vs_btof,
        Form("run%05d CBtof0 vs dE(Btof0) (%s)", run_num, particle.Data())
      )) {
    ++pages;
  }
  c->Print(pdf_path + "]");
  delete c;
  f->Close();

  if (pages < 1) {
    std::cerr << "Error: no CBtof0 pages drawn" << std::endl;
    return 1;
  }
  std::cout << "wrote " << pdf_path << " (" << pages << " pages)" << std::endl;
  return 0;
}

//_____________________________________________________________________________
Int_t
main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " <Hodo.root> <Pi|K> [run]" << std::endl;
    return 1;
  }
  TString path = argv[1];
  TString particle = argv[2];
  if (particle != "Pi" && particle != "K") {
    std::cerr << "Error: Unexpected particle name: " << particle << std::endl;
    return 1;
  }
  Int_t run_override = -1;
  if (argc >= 4) {
    run_override = TString(argv[3]).Atoi();
  }
  return Analyze(path, particle, run_override);
}
