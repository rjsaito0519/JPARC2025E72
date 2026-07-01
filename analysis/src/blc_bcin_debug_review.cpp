// UserBcInTracking ROOT -> BLC1 geometry debug multi-page PDF
// Usage: blc_bcin_debug_review [-o <out.pdf>] <input.root>
// Default PDF: OUTPUT_DIR/img/runXXXXX/runXXXXX_BLC1_bcin_debug_review.pdf

#include "ana_helper.h"
#include "config.h"
#include "paths.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TTree.h>

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <vector>

namespace
{

Config& conf = Config::getInstance();

const char* kPlaneLabels[] = {
  "U1", "UP1", "V1", "VP1", "U2", "UP2", "V2", "VP2"
};

const char* kSuffixCandidates[] = {"", "_Pi", "_K"};
const Int_t kNumSuffix = 3;

struct SuffixInfo
{
  TString suffix;
  Long64_t entries[kNumSuffix];
};

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " [-o <out.pdf>] <input.root>\n"
            << "  Reads UserBcInTracking histograms (bcin tree required).\n"
            << "  Default PDF: " << OUTPUT_DIR << "/img/runXXXXX/"
            << "runXXXXX_BLC1_bcin_debug_review.pdf\n";
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

SuffixInfo
resolve_suffix(TFile* f)
{
  SuffixInfo info;
  Long64_t best_entries = -1;
  for (Int_t i = 0; i < kNumSuffix; ++i) {
    const char* sfx = kSuffixCandidates[i];
    info.entries[i] = hist_entries(f, Form("BcInTrack_Y0%s", sfx));
    if (info.entries[i] > best_entries) {
      best_entries = info.entries[i];
      info.suffix = sfx;
    } else if (info.entries[i] == best_entries && best_entries > 0) {
      // tie: prefer earlier candidate ("" > _Pi > _K)
      info.suffix = sfx;
    }
  }
  if (best_entries <= 0) info.suffix = "";
  return info;
}

void
draw_placeholder(const char* message)
{
  auto* tex = new TLatex(0.5, 0.5, message);
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.06);
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
draw_hist_or_placeholder(TFile* f, const char* hname, const char* opt,
                         const char* plane_label, Bool_t logz = false)
{
  auto* obj = get_object(f, hname);
  if (!obj) {
    draw_placeholder(Form("NOT FOUND:\n%s", hname));
    return;
  }
  if (logz) gPad->SetLogz();
  obj->Draw(opt);
  if (plane_label) draw_plane_label(plane_label);
}

void
setup_style()
{
  gROOT->SetBatch();
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
  return Form("%s/run%05d_BLC1_bcin_debug_review.pdf", img_dir.Data(), run);
}

void
print_page(TCanvas* c, const TString& pdf_path)
{
  c->Print(pdf_path);
}

void
draw_title_page(TCanvas* c, const TString& pdf_path, TFile* f,
                Int_t run_num, const SuffixInfo& sfx)
{
  c->Clear();
  c->cd();

  const TString input_name = f->GetName();
  TTimeStamp now;
  now.Set();

  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(22);
  tex->SetTextSize(0.05);
  tex->DrawLatex(0.5, 0.82, "BLC1 debug review (UserBcInTracking)");
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.5, 0.72, input_name);
  tex->DrawLatex(0.5, 0.66, Form("run_number = %d", run_num));
  tex->DrawLatex(0.5, 0.60, Form("selected suffix = \"%s\"", sfx.suffix.Data()));
  tex->DrawLatex(0.5, 0.54,
                 Form("BcInTrack_Y0 entries: all=%lld  Pi=%lld  K=%lld",
                      sfx.entries[0], sfx.entries[1], sfx.entries[2]));
  tex->DrawLatex(0.5, 0.46, Form("created: %s", now.AsString("s")));
  print_page(c, pdf_path);
}

void
draw_status_beam_page(TCanvas* c, const TString& pdf_path, TFile* f)
{
  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  auto* status = get_object(f, "Status");
  if (status) status->Draw();
  else draw_placeholder("NOT FOUND: Status");
  c->cd(2);
  auto* beam = get_object(f, "BeamFlag");
  if (beam) beam->Draw();
  else draw_placeholder("NOT FOUND: BeamFlag");
  print_page(c, pdf_path);
}

void
draw_track_summary_page(TCanvas* c, const TString& pdf_path, TFile* f,
                        const TString& suffix)
{
  const char* names[] = {
    "BcInTrack_NHit", "BcInTrack_ChiSquare", "BcInTrack_X0",
    "BcInTrack_Y0", "BcInTrack_U0", "BcInTrack_V0"
  };
  c->Clear();
  c->Divide(3, 2);
  for (Int_t i = 0; i < 6; ++i) {
    c->cd(i + 1);
    const TString hname = Form("%s%s", names[i], suffix.Data());
    draw_hist_or_placeholder(f, hname, "", names[i], false);
  }
  print_page(c, pdf_path);
}

void
draw_residual_page(TCanvas* c, const TString& pdf_path, TFile* f,
                   const char* chamber, const TString& suffix)
{
  const Int_t nplane = conf.num_of_ch.at("blc");
  c->Clear();
  c->Divide(4, 2);
  for (Int_t i = 0; i < nplane; ++i) {
  c->cd(i + 1);
    const TString hname =
      Form("%s_Track_Residual_plane%d%s", chamber, i, suffix.Data());
    auto* h = dynamic_cast<TH1D*>(get_object(f, hname));
    if (!h || h->GetEntries() <= 0) {
      draw_placeholder(h ? "empty histogram" : Form("NOT FOUND:\n%s", hname.Data()));
      continue;
    }
    ana_helper::residual_fit(h, c, i + 1);
    draw_plane_label(kPlaneLabels[i]);
  }
  print_page(c, pdf_path);
}

void
draw_colz_page(TCanvas* c, const TString& pdf_path, TFile* f,
               const char* chamber, const TString& suffix,
               const char* hist_middle)
{
  const Int_t nplane = conf.num_of_ch.at("blc");
  c->Clear();
  c->Divide(4, 2);
  for (Int_t i = 0; i < nplane; ++i) {
    c->cd(i + 1);
    const TString hname =
      Form("%s_%s_plane%d%s", chamber, hist_middle, i, suffix.Data());
    draw_hist_or_placeholder(f, hname, "colz", kPlaneLabels[i], true);
  }
  print_page(c, pdf_path);
}

Int_t
analyze(const TString& input_path, const TString& output_path_in)
{
  setup_style();

  auto* f = TFile::Open(input_path);
  if (!f || f->IsZombie()) {
    std::cerr << "Error: could not open " << input_path << std::endl;
    return 1;
  }
  auto* bcin = dynamic_cast<TTree*>(f->Get("bcin"));
  if (!bcin) {
    std::cerr << "Error: 'bcin' tree not found (UserBcInTracking ROOT required)."
              << std::endl;
    f->Close();
    return 1;
  }

  Int_t run_num = -1;
  if (bcin->GetEntries() > 0) {
    UInt_t run_number = 0;
    if (bcin->GetBranch("run_number"))
      bcin->SetBranchAddress("run_number", &run_number);
    bcin->GetEntry(0);
    run_num = static_cast<Int_t>(run_number);
  }
  if (run_num < 0) run_num = parse_run_from_path(input_path);

  TString output_path = output_path_in;
  if (output_path.IsNull()) {
    output_path = default_pdf_path(run_num);
    std::cout << "Default output (calib-style img dir): " << output_path << std::endl;
  }

  const SuffixInfo sfx = resolve_suffix(f);
  if (sfx.entries[0] + sfx.entries[1] + sfx.entries[2] <= 0) {
    std::cerr << "Warning: BcInTrack_Y0* histograms are empty; PDF may be sparse."
              << std::endl;
  }

  auto* canvas = new TCanvas("blc_bcin_debug_review", output_path, 1600, 1200);
  canvas->Print(output_path + "[");

  draw_title_page(canvas, output_path, f, run_num, sfx);
  draw_status_beam_page(canvas, output_path, f);
  draw_track_summary_page(canvas, output_path, f, sfx.suffix);

  const char* chambers[] = {"BLC1a", "BLC1b"};
  for (const char* chamber : chambers) {
    draw_residual_page(canvas, output_path, f, chamber, sfx.suffix);
    draw_colz_page(canvas, output_path, f, chamber, sfx.suffix,
                   "Track_Residual_vs_DriftLength");
    draw_colz_page(canvas, output_path, f, chamber, sfx.suffix,
                   "Track_DriftLength_vs_HitPat");
    draw_colz_page(canvas, output_path, f, chamber, sfx.suffix,
                   "Hit_DriftLength_vs_HitPat");
  }

  canvas->Print(output_path + "]");
  delete canvas;
  f->Close();

  std::cout << "Wrote " << output_path << " (run=" << run_num
            << ", suffix=\"" << sfx.suffix << "\")" << std::endl;
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
