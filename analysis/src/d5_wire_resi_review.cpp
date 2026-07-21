// D5 wire residual review (UserD5Tracking histograms) -> multi-page PDF
// Usage: d5_wire_resi_review [-o <out.pdf>] <decode.root>
// Default PDF: OUTPUT_DIR/img/runXXXXX/runXXXXX_D5_wire_resi_review.pdf

#include "ana_helper.h"
#include "config.h"
#include "paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <cstdio>
#include <cstring>
#include <iostream>

namespace
{

Config& conf = Config::getInstance();

const char* kDetectors[] = {"BLC1a", "BLC1b", "BLC2a", "BLC2b"};
const Int_t kNumDetectors = 4;

const char* kPlaneLabels[] = {
  "U1", "UP1", "V1", "VP1", "U2", "UP2", "V2", "VP2"
};

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " [-o <out.pdf>] <decode.root>\n"
            << "  Reads D5WireResi_* histograms from D5 Tracking decode ROOT.\n"
            << "  Default PDF: " << OUTPUT_DIR << "/img/runXXXXX/"
            << "runXXXXX_D5_wire_resi_review.pdf\n";
}

TObject*
get_object(TFile* f, const char* name)
{
  if (!f || !name) return nullptr;
  return f->Get(name);
}

Long64_t
hist_entries(TFile* f, const char* name)
{
  auto* h = dynamic_cast<TH1*>(get_object(f, name));
  if (!h) return 0;
  return static_cast<Long64_t>(h->GetEntries());
}

void
setup_style()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleSize(0.045, "XYZ");
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
}

Int_t
parse_run_from_path(const TString& path)
{
  const std::string s(path.Data());
  const std::size_t pos = s.find("run");
  if (pos == std::string::npos || pos + 8 > s.size()) return -1;
  Int_t run = -1;
  if (std::sscanf(s.c_str() + pos + 3, "%5d", &run) == 1) return run;
  return -1;
}

TString
default_pdf_path(Int_t run_num)
{
  const Int_t run = run_num >= 0 ? run_num : 0;
  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run);
  return Form("%s/run%05d_D5_wire_resi_review.pdf", img_dir.Data(), run);
}

void
print_page(TCanvas* c, const TString& pdf_path)
{
  c->Print(pdf_path);
}

void
draw_placeholder(const char* message)
{
  auto* tex = new TLatex(0.5, 0.5, message);
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->Draw();
}

void
draw_plane_label(const char* label)
{
  auto* tex = new TLatex(0.18, 0.88, label);
  tex->SetNDC();
  tex->SetTextAlign(11);
  tex->SetTextSize(0.045);
  tex->Draw();
}

void
draw_title_page(TCanvas* c, const TString& pdf_path, TFile* f, Int_t run_num)
{
  c->Clear();
  c->cd();
  TTimeStamp now;
  now.Set();

  Long64_t total = 0;
  for (Int_t id = 0; id < kNumDetectors; ++id) {
    total += hist_entries(f, Form("D5WireResi_%s_LocalFit_vs_plane", kDetectors[id]));
  }

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->DrawLatex(0.5, 0.84, "D5 wire residual review (rank0 D5 pairs)");
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.5, 0.74, f->GetName());
  tex->DrawLatex(0.5, 0.68, Form("run_number = %d", run_num));
  tex->DrawLatex(0.5, 0.62,
                 Form("D5WireResi LocalFit_vs_plane entries (sum) = %lld",
                      static_cast<long long>(total)));
  tex->DrawLatex(0.5, 0.54, "LocalFit = chamber DoFit | D5Fit = joint-fit prediction");
  tex->DrawLatex(0.5, 0.46, Form("created: %s", now.AsString("s")));
  print_page(c, pdf_path);
}

void
draw_vs_plane_page(TCanvas* c, const TString& pdf_path, TFile* f,
                   const char* detector)
{
  c->Clear();
  c->Divide(2, 1);
  const char* names[] = {"LocalFit", "D5Fit"};
  for (Int_t i = 0; i < 2; ++i) {
    c->cd(i + 1);
    const TString hname =
      Form("D5WireResi_%s_%s_vs_plane", detector, names[i]);
    auto* h = get_object(f, hname);
    if (!h) {
      draw_placeholder(Form("NOT FOUND:\n%s", hname.Data()));
      continue;
    }
    gPad->SetLogz();
    h->Draw("colz");
    auto* tex = new TLatex(0.18, 0.88, Form("%s %s", detector, names[i]));
    tex->SetNDC();
    tex->SetTextAlign(11);
    tex->SetTextSize(0.045);
    tex->Draw();
  }
  print_page(c, pdf_path);
}

void
draw_plane_overlay_page(TCanvas* c, const TString& pdf_path, TFile* f,
                        const char* detector)
{
  const Int_t nplane = conf.num_of_ch.at("blc");
  c->Clear();
  c->Divide(4, 2);
  for (Int_t i = 0; i < nplane; ++i) {
    c->cd(i + 1);
    const TString h_local =
      Form("D5WireResi_%s_LocalFit_plane%d", detector, i);
    const TString h_d5 =
      Form("D5WireResi_%s_D5Fit_plane%d", detector, i);
    auto* hl = dynamic_cast<TH1*>(get_object(f, h_local));
    auto* hd = dynamic_cast<TH1*>(get_object(f, h_d5));
    if (!hl && !hd) {
      draw_placeholder(Form("NOT FOUND:\nplane %d", i));
      continue;
    }
    if (hl) {
      hl->SetLineColor(kBlack);
      hl->Draw("hist");
    }
    if (hd) {
      hd->SetLineColor(kRed + 1);
      if (hl) hd->Draw("hist same");
      else hd->Draw("hist");
    }
    if (hl || hd) {
      auto* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
      if (hl) leg->AddEntry(hl, "LocalFit", "l");
      if (hd) leg->AddEntry(hd, "D5Fit", "l");
      leg->Draw();
    }
    if (i < 8) draw_plane_label(kPlaneLabels[i]);
  }
  print_page(c, pdf_path);
}

Int_t
analyze(const TString& input_path, const TString& output_path_in)
{
  setup_style();

  auto* f = TFile::Open(input_path);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open: " << input_path << std::endl;
    return 1;
  }

  Int_t run_num = -1;
  auto* d5 = dynamic_cast<TTree*>(f->Get("d5"));
  if (d5 && d5->GetEntries() > 0) {
    UInt_t run_number = 0;
    if (d5->GetBranch("run_number"))
      d5->SetBranchAddress("run_number", &run_number);
    d5->GetEntry(0);
    run_num = static_cast<Int_t>(run_number);
  }
  if (run_num < 0) run_num = parse_run_from_path(input_path);

  TString output_path = output_path_in;
  if (output_path.IsNull()) {
    output_path = default_pdf_path(run_num);
    std::cout << "Default output: " << output_path << std::endl;
  }

  if (hist_entries(f, "D5WireResi_BLC1a_LocalFit_vs_plane") <= 0) {
    std::cerr << "Warning: D5WireResi_* histograms not found or empty.\n"
              << "  Re-decode with updated UserD5Tracking.\n";
  }

  auto* canvas = new TCanvas("d5_wire_resi_review", output_path, 1600, 1200);
  canvas->Print(output_path + "[");
  draw_title_page(canvas, output_path, f, run_num);
  for (Int_t id = 0; id < kNumDetectors; ++id) {
    draw_vs_plane_page(canvas, output_path, f, kDetectors[id]);
    draw_plane_overlay_page(canvas, output_path, f, kDetectors[id]);
  }
  canvas->Print(output_path + "]");
  delete canvas;
  f->Close();

  std::cout << "Wrote " << output_path << " (run=" << run_num << ")" << std::endl;
  return 0;
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  TString input_path;
  TString output_path;
  bool has_output = false;

  for (Int_t i = 1; i < argc; ++i) {
    const char* arg = argv[i];
    if (strcmp(arg, "-o") == 0) {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      output_path = argv[++i];
      has_output = true;
      continue;
    }
    if (arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << std::endl;
      usage(argv[0]);
      return 1;
    }
    if (!input_path.IsNull()) {
      std::cerr << "Multiple input files specified." << std::endl;
      usage(argv[0]);
      return 1;
    }
    input_path = arg;
  }

  if (input_path.IsNull()) {
    usage(argv[0]);
    return 1;
  }

  gSystem->ExpandPathName(input_path);
  return analyze(input_path, has_output ? output_path : TString());
}
