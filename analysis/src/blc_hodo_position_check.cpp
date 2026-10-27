// Independent check of the BLC global x positions against the hodoscope segment geometry
// Usage: blc_hodo_position_check <run> [--pion | --kaon | --all] [-o <out.pdf>]
// Inputs (S = Pi for --pion (default), K for --kaon / --all):
//           BcIn / BcOut / Hodo  run%05d_{BcIn,BcOut,Hodo}.root under DATA_DIR,
//                                else DECODE_DIR/runXXXXX/; a Hodo decode of the same run is required
//           DCGEO  ANALYZER_DIR/param/DCGEO/e72/DCGeomParam_runXXXXX_S  (z of BH2, BHT, BLC layers)
//           pdf    OUTPUT_DIR/img/runXXXXX/blc_hodo_position_check_runXXXXX_<Pi|K|All>.pdf
// Events: Hodo beam_flag selection; BcIn / BcOut entries joined by event_number; exactly one
// track in the chamber used. BH2 needs exactly one BH2 hit segment; BHT uses every hit segment.
//
// BH2 (15 seg along x, nominal centre x_seg = (seg - 7) * 14 mm, as Bh2XPos of the analyzer):
// the BLC2 (BcOut) track is extrapolated to z_BH2. The per-segment core mean of x is fitted with
// <x> = offset + slope * x_seg. slope ~ 1 checks the scale; offset = x_BLC2 - x_BH2 at the BH2
// centre, i.e. a common x shift of BLC2 relative to BH2 (or of BH2). BcOut tracking does not use
// the BH2 segment filter, so the two measurements are independent. y vs segment is flat unless x
// and y are swapped. BHT (63 seg) has no segment centres here: BLC1 (BcIn) at z_BHT vs segment
// is shown as a correlation only.
// Pages: 1 summary; 2 BLC2 x, y at z_BH2 vs BH2 seg; 3 BH2 per-segment mean x vs x_seg with the
//        linear fit and residual; 4 BLC1 x, y at z_BHT vs BHT seg; 5 BLC1 / BLC2 beam profiles at
//        the chamber centres.
// Fit numbers are printed to stdout.

#include "ana_helper.h"
#include "d5_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{

constexpr Int_t kNumBh2Seg = 15;
constexpr Double_t kBh2Width = 14.0; // [mm]
constexpr Int_t kNumBhtSeg = 63;

constexpr Double_t kPosRange = 120.0;      // [mm] drawn |x|, |y| at BH2 / BLC2
constexpr Double_t kBlc1Range = 200.0;     // [mm] drawn |x|, |y| at BHT / BLC1 (wider beam)
constexpr Double_t kCoreHalfWindow = 14.0; // [mm] window around the median for the core mean
constexpr Int_t kMinSegEntries = 50;       // entries needed for a segment to enter the fit

Double_t
bh2_seg_x(Int_t seg)
{
  return (seg - (kNumBh2Seg - 1) / 2.0) * kBh2Width;
}

// first existing of DATA_DIR/run%05d_<tag>.root, DECODE_DIR/runXXXXX/run%05d_<tag>.root
TString
resolve_input_path(Int_t run, const char* tag)
{
  const TString candidates[] = {
    Form("%s/run%05d_%s.root", DATA_DIR.Data(), run, tag),
    Form("%s/run%05d/run%05d_%s.root", DECODE_DIR.Data(), run, run, tag),
  };
  for (const auto& path : candidates) {
    if (gSystem->AccessPathName(path) == 0) return path;
  }
  return candidates[0];
}

// Z column of the DCGEO row named name (same format as d5_helper::load_dcgeo_blc)
Bool_t
load_dcgeo_z(const TString& path, const TString& name, Double_t& z)
{
  std::ifstream ifs(path.Data());
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string id, nm;
    Double_t x, y, zz;
    if (!(iss >> id >> nm >> x >> y >> zz)) continue;
    if (name == nm.c_str()) { z = zz; return kTRUE; }
  }
  return kFALSE;
}

std::unordered_map<UInt_t, Long64_t>
build_event_index(TTree* tree)
{
  std::unordered_map<UInt_t, Long64_t> index;
  UInt_t event_number = 0;
  tree->SetBranchAddress("event_number", &event_number);
  const Long64_t nent = tree->GetEntries();
  index.reserve(static_cast<std::size_t>(nent));
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tree->GetEntry(ie);
    index[event_number] = ie;
  }
  tree->ResetBranchAddress(tree->GetBranch("event_number"));
  return index;
}

// branches of a bcin / bcout local-track tree used here
struct TrackTree
{
  TFile* file = nullptr;
  TTree* tree = nullptr;
  Int_t ntrack = 0;
  std::vector<Double_t> *x0 = nullptr, *y0 = nullptr, *u0 = nullptr, *v0 = nullptr;

  Bool_t open(const TString& path, const char* name)
  {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return kFALSE;
    tree = dynamic_cast<TTree*>(file->Get(name));
    return tree != nullptr;
  }
  void set_branches()
  {
    tree->SetBranchAddress("ntrack", &ntrack);
    tree->SetBranchAddress("x0", &x0); tree->SetBranchAddress("y0", &y0);
    tree->SetBranchAddress("u0", &u0); tree->SetBranchAddress("v0", &v0);
  }
  Bool_t single() const { return ntrack == 1 && x0 && x0->size() == 1; }
  Double_t x(Double_t z) const { return (*x0)[0] + (*u0)[0] * z; }
  Double_t y(Double_t z) const { return (*y0)[0] + (*v0)[0] * z; }
};

// mean of h within median +- kCoreHalfWindow, re-centred 3 times; false if too few entries
Bool_t
core_mean(TH1D* h, Double_t& mean, Double_t& err)
{
  if (!h || h->GetEntries() < kMinSegEntries) return kFALSE;
  Double_t q = 0.5, med = 0.0;
  h->GetQuantiles(1, &med, &q);
  mean = med;
  for (Int_t it = 0; it < 3; ++it) {
    h->GetXaxis()->SetRangeUser(mean - kCoreHalfWindow, mean + kCoreHalfWindow);
    mean = h->GetMean();
    err = h->GetMeanError();
  }
  h->GetXaxis()->SetRange();
  return kTRUE;
}

// per-segment core mean of the y axis of h (x axis = segment); points at x_of(seg), or seg
TGraphErrors*
seg_means(TH2D* h, Int_t nseg, Double_t (*x_of)(Int_t))
{
  auto* g = new TGraphErrors();
  for (Int_t s = 0; s < nseg; ++s) {
    TH1D* p = h->ProjectionY(Form("%s_s%d", h->GetName(), s), s + 1, s + 1);
    Double_t m = 0.0, e = 0.0;
    if (core_mean(p, m, e)) {
      const Int_t n = g->GetN();
      g->SetPoint(n, x_of ? x_of(s) : s, m);
      g->SetPointError(n, 0.0, e);
    }
    delete p;
  }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  return g;
}

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " <run> [--pion | --kaon | --all] [-o <out.pdf>]\n"
            << "  needs run%05d_{BcIn,BcOut,Hodo}.root (Hodo decode of the same run) under\n"
            << "  " << DATA_DIR << " or " << DECODE_DIR << "/runXXXXX/\n";
}

} // namespace

Int_t
main(Int_t argc, char** argv)
{
  if (argc < 2) { usage(argv[0]); return 1; }
  const Int_t run = std::atoi(argv[1]);
  if (run <= 0) { usage(argv[0]); return 1; }
  d5_helper::BeamSel beam;
  TString path_pdf;
  for (Int_t i = 2; i < argc; ++i) {
    const TString a = argv[i];
    if (d5_helper::parse_beam_option(a, beam)) continue;
    if (a == "-o" && i + 1 < argc) path_pdf = argv[++i];
    else { usage(argv[0]); return 1; }
  }
  if (path_pdf.IsNull())
    path_pdf = Form("%s/blc_hodo_position_check_run%05d_%s.pdf",
                    ana_helper::get_img_dir(OUTPUT_DIR, run).Data(), run, beam.out_tag.Data());

  // --- geometry from the run's DCGEO ---
  const TString path_geo = Form("%s/param/DCGEO/e72/DCGeomParam_run%05d_%s",
                                ANALYZER_DIR.Data(), run, beam.param_suffix.Data());
  Double_t z_bh2 = 0.0, z_bht = 0.0;
  if (!load_dcgeo_z(path_geo, "BH2", z_bh2) || !load_dcgeo_z(path_geo, "BHT", z_bht)) {
    std::cerr << "cannot read BH2 / BHT z from DCGEO: " << path_geo << std::endl;
    return 1;
  }
  const auto geo = d5_helper::load_dcgeo_blc(path_geo);
  Double_t z1a = 0.0, z1b = 0.0, z2a = 0.0, z2b = 0.0;
  if (!d5_helper::chamber_mean_z(geo, "BLC1a", z1a) || !d5_helper::chamber_mean_z(geo, "BLC1b", z1b) ||
      !d5_helper::chamber_mean_z(geo, "BLC2a", z2a) || !d5_helper::chamber_mean_z(geo, "BLC2b", z2b)) {
    std::cerr << "cannot read BLC layer z from DCGEO: " << path_geo << std::endl;
    return 1;
  }
  const Double_t z_blc1 = 0.5 * (z1a + z1b); // chamber centre in the BcIn frame
  const Double_t z_blc2 = 0.5 * (z2a + z2b); // chamber centre in the BcOut frame

  // --- inputs ---
  const TString path_bcin = resolve_input_path(run, "BcIn");
  const TString path_bcout = resolve_input_path(run, "BcOut");
  const TString path_hodo = resolve_input_path(run, "Hodo");
  TrackTree bcin, bcout;
  if (!bcin.open(path_bcin, "bcin")) { std::cerr << "cannot read bcin tree of " << path_bcin << std::endl; return 1; }
  if (!bcout.open(path_bcout, "bcout")) { std::cerr << "cannot read bcout tree of " << path_bcout << std::endl; return 1; }
  auto* f_hodo = TFile::Open(path_hodo);
  auto* tree_hodo = (f_hodo && !f_hodo->IsZombie()) ? dynamic_cast<TTree*>(f_hodo->Get("hodo")) : nullptr;
  if (!tree_hodo) {
    std::cerr << "cannot read hodo tree of " << path_hodo
              << " (a Hodo decode of the same run is required)" << std::endl;
    return 1;
  }
  std::cout << "BcIn : " << path_bcin << "\nBcOut: " << path_bcout << "\nHodo : " << path_hodo
            << "\nDCGEO: " << path_geo << "\nbeam : " << beam.label << std::endl;

  const auto bcin_idx = build_event_index(bcin.tree);
  const auto bcout_idx = build_event_index(bcout.tree);
  bcin.set_branches();
  bcout.set_branches();

  UInt_t event_number = 0;
  Int_t beam_flag = 0;
  std::vector<Double_t>* bh2_seg = nullptr;
  std::vector<Double_t>* bht_seg = nullptr;
  tree_hodo->SetBranchAddress("event_number", &event_number);
  tree_hodo->SetBranchAddress("beam_flag", &beam_flag);
  tree_hodo->SetBranchAddress("bh2_hit_seg", &bh2_seg);
  tree_hodo->SetBranchAddress("bht_hit_seg", &bht_seg);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(29.0, 20.5);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetTitleOffset(1.4, "Y");

  const Int_t npos = static_cast<Int_t>(2 * kPosRange / 0.5);
  const Int_t npos1 = static_cast<Int_t>(2 * kBlc1Range / 0.5);
  auto* h_bh2_x = new TH2D("h_bh2_x", Form("BLC2 x @ z_{BH2} = %.1f mm vs BH2 seg; BH2 seg; x_{BLC2} [mm]", z_bh2),
                           kNumBh2Seg, -0.5, kNumBh2Seg - 0.5, npos, -kPosRange, kPosRange);
  auto* h_bh2_y = new TH2D("h_bh2_y", Form("BLC2 y @ z_{BH2} = %.1f mm vs BH2 seg; BH2 seg; y_{BLC2} [mm]", z_bh2),
                           kNumBh2Seg, -0.5, kNumBh2Seg - 0.5, npos, -kPosRange, kPosRange);
  auto* h_bht_x = new TH2D("h_bht_x", Form("BLC1 x @ z_{BHT} = %.1f mm vs BHT seg; BHT seg; x_{BLC1} [mm]", z_bht),
                           kNumBhtSeg, -0.5, kNumBhtSeg - 0.5, npos1, -kBlc1Range, kBlc1Range);
  auto* h_bht_y = new TH2D("h_bht_y", Form("BLC1 y @ z_{BHT} = %.1f mm vs BHT seg; BHT seg; y_{BLC1} [mm]", z_bht),
                           kNumBhtSeg, -0.5, kNumBhtSeg - 0.5, npos1, -kBlc1Range, kBlc1Range);
  auto* h_blc1_xy = new TH2D("h_blc1_xy", Form("BLC1 beam profile @ z = %.1f mm (BcIn frame); x [mm]; y [mm]", z_blc1),
                             200, -kBlc1Range, kBlc1Range, 200, -kBlc1Range, kBlc1Range);
  auto* h_blc2_xy = new TH2D("h_blc2_xy", Form("BLC2 beam profile @ z = %.1f mm (BcOut frame); x [mm]; y [mm]", z_blc2),
                             200, -kPosRange, kPosRange, 200, -kPosRange, kPosRange);

  Long64_t n_sel = 0, n_matched = 0, n_bh2 = 0, n_bht = 0;
  const Long64_t n_hodo = tree_hodo->GetEntries();
  for (Long64_t ie = 0; ie < n_hodo; ++ie) {
    tree_hodo->GetEntry(ie);
    if (beam.flag >= 0 && beam_flag != beam.flag) continue;
    ++n_sel;
    const auto it_in = bcin_idx.find(event_number);
    const auto it_out = bcout_idx.find(event_number);
    if (it_in == bcin_idx.end() || it_out == bcout_idx.end()) continue;
    ++n_matched;
    bcin.tree->GetEntry(it_in->second);
    bcout.tree->GetEntry(it_out->second);

    if (bcout.single()) {
      h_blc2_xy->Fill(bcout.x(z_blc2), bcout.y(z_blc2));
      if (bh2_seg && bh2_seg->size() == 1) {
        const Int_t seg = static_cast<Int_t>(std::lround((*bh2_seg)[0]));
        h_bh2_x->Fill(seg, bcout.x(z_bh2));
        h_bh2_y->Fill(seg, bcout.y(z_bh2));
        ++n_bh2;
      }
    }
    if (bcin.single()) {
      h_blc1_xy->Fill(bcin.x(z_blc1), bcin.y(z_blc1));
      if (bht_seg && !bht_seg->empty()) {
        for (const Double_t s : *bht_seg) {
          const Int_t seg = static_cast<Int_t>(std::lround(s));
          h_bht_x->Fill(seg, bcin.x(z_bht));
          h_bht_y->Fill(seg, bcin.y(z_bht));
        }
        ++n_bht;
      }
    }
  }

  // --- BH2: per-segment mean x vs nominal segment centre ---
  TGraphErrors* g_bh2 = seg_means(h_bh2_x, kNumBh2Seg, bh2_seg_x);
  g_bh2->SetTitle("BH2: core mean of x_{BLC2} vs segment centre; x_{seg} = (seg - 7) #times 14 [mm]; <x_{BLC2}> [mm]");
  auto* f_bh2 = new TF1("f_bh2", "pol1", -kPosRange, kPosRange);
  f_bh2->SetLineColor(kRed);
  Double_t offset = NAN, offset_err = NAN, slope = NAN, slope_err = NAN;
  if (g_bh2->GetN() >= 3) {
    g_bh2->Fit(f_bh2, "Q");
    offset = f_bh2->GetParameter(0); offset_err = f_bh2->GetParError(0);
    slope = f_bh2->GetParameter(1); slope_err = f_bh2->GetParError(1);
  }
  auto* g_res = new TGraphErrors();
  g_res->SetTitle("BH2: <x_{BLC2}> - x_{seg}; x_{seg} [mm]; <x_{BLC2}> - x_{seg} [mm]");
  for (Int_t i = 0; i < g_bh2->GetN(); ++i) {
    g_res->SetPoint(i, g_bh2->GetX()[i], g_bh2->GetY()[i] - g_bh2->GetX()[i]);
    g_res->SetPointError(i, 0.0, g_bh2->GetEY()[i]);
  }
  g_res->SetMarkerStyle(20);
  g_res->SetMarkerSize(0.8);
  TGraphErrors* g_bh2x_seg = seg_means(h_bh2_x, kNumBh2Seg, nullptr);
  TGraphErrors* g_bh2y = seg_means(h_bh2_y, kNumBh2Seg, nullptr);
  TGraphErrors* g_bhtx = seg_means(h_bht_x, kNumBhtSeg, nullptr);
  TGraphErrors* g_bhty = seg_means(h_bht_y, kNumBhtSeg, nullptr);
  for (auto* g : { g_bh2x_seg, g_bh2y, g_bhtx, g_bhty }) { g->SetMarkerColor(kRed); g->SetLineColor(kRed); }

  std::cout << Form("BH2 fit (%d segments, %lld events): <x_BLC2> = offset + slope * x_seg\n"
                    "  slope  = %.5f +- %.5f\n  offset = %+.3f +- %.3f mm  (x_BLC2 - x_BH2 at x = 0)",
                    g_bh2->GetN(), n_bh2, slope, slope_err, offset, offset_err) << std::endl;
  for (Int_t i = 0; i < g_res->GetN(); ++i)
    std::cout << Form("  x_seg %+6.1f  <x> - x_seg = %+.3f +- %.3f mm", g_res->GetX()[i], g_res->GetY()[i],
                      g_res->GetEY()[i]) << std::endl;
  std::cout << Form("BLC1 beam @ z=%.1f: x %.2f (rms %.2f), y %.2f (rms %.2f) mm", z_blc1, h_blc1_xy->GetMean(1),
                    h_blc1_xy->GetRMS(1), h_blc1_xy->GetMean(2), h_blc1_xy->GetRMS(2)) << std::endl;
  std::cout << Form("BLC2 beam @ z=%.1f: x %.2f (rms %.2f), y %.2f (rms %.2f) mm", z_blc2, h_blc2_xy->GetMean(1),
                    h_blc2_xy->GetRMS(1), h_blc2_xy->GetMean(2), h_blc2_xy->GetRMS(2)) << std::endl;

  // --- PDF ---
  auto* c = new TCanvas("c", "c", 1600, 1130);
  c->Print(path_pdf + "[");
  c->cd();
  auto* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.024);
  Double_t yy = 0.94;
  auto line = [&](const TString& s) { tex->DrawLatex(0.05, yy, s); yy -= 0.040; };
  line(Form("BLC position vs hodoscope segments: run %05d, %s", run, beam.label.Data()));
  line(Form("BcIn: %s", path_bcin.Data()));
  line(Form("BcOut: %s", path_bcout.Data()));
  line(Form("Hodo: %s", path_hodo.Data()));
  line(Form("DCGEO: %s", path_geo.Data()));
  line(Form("z: BH2 %.1f, BLC2 centre %.1f mm (BcOut frame);  BHT %.1f, BLC1 centre %.1f mm (BcIn frame)",
            z_bh2, z_blc2, z_bht, z_blc1));
  line(Form("events: hodo %lld, beam selected %lld, in BcIn and BcOut %lld", n_hodo, n_sel, n_matched));
  line(Form("  BH2: 1 BcOut track and 1 BH2 hit segment: %lld;  BHT: 1 BcIn track, every BHT hit: %lld", n_bh2, n_bht));
  line("BH2 fit: <x_{BLC2}> (core mean, median #pm 14 mm, per segment) = offset + slope #times x_{seg}");
  line(Form("  slope = %.4f #pm %.4f,  offset = %+.3f #pm %.3f mm  (%d segments with #geq %d entries)",
            slope, slope_err, offset, offset_err, g_bh2->GetN(), kMinSegEntries));
  line("  slope #approx 1 expected; offset = x_{BLC2} - x_{BH2} at the BH2 centre (common x shift of BLC2 or BH2)");
  line("  the per-segment mean follows the beam profile inside a segment: edge segments pull the slope");
  line("BH2 / BHT y vs segment should be flat (a band means an x/y swap)");
  line("BHT: no segment centres used here; correlation only");
  c->Print(path_pdf, "Title:summary");

  auto draw2 = [&](TH2D* h, TGraphErrors* g) {
    gPad->SetLogz();
    h->Draw("COLZ");
    if (g && g->GetN() > 0) g->Draw("P same");
  };

  c->Clear();
  c->Divide(2, 1);
  c->cd(1); draw2(h_bh2_x, g_bh2x_seg);
  c->cd(2); draw2(h_bh2_y, g_bh2y);
  c->Print(path_pdf, "Title:BH2 correlation");

  c->Clear();
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetGrid();
  g_bh2->Draw("AP");
  auto* t_fit = new TLatex();
  t_fit->SetNDC();
  t_fit->SetTextSize(0.035);
  t_fit->DrawLatex(0.16, 0.84, Form("slope = %.4f #pm %.4f", slope, slope_err));
  t_fit->DrawLatex(0.16, 0.79, Form("offset = %+.3f #pm %.3f mm", offset, offset_err));
  c->cd(2);
  gPad->SetGrid();
  g_res->Draw("AP");
  if (std::isfinite(slope)) {
    auto* f_res = new TF1("f_res", "[0] + ([1] - 1) * x", -kPosRange, kPosRange);
    f_res->SetParameters(offset, slope);
    f_res->SetLineColor(kRed);
    f_res->Draw("same");
  }
  c->Print(path_pdf, "Title:BH2 fit");

  c->Clear();
  c->Divide(2, 1);
  c->cd(1); draw2(h_bht_x, g_bhtx);
  c->cd(2); draw2(h_bht_y, g_bhty);
  c->Print(path_pdf, "Title:BHT correlation");

  c->Clear();
  c->Divide(2, 1);
  c->cd(1); h_blc1_xy->Draw("COLZ");
  c->cd(2); h_blc2_xy->Draw("COLZ");
  c->Print(path_pdf, "Title:beam profiles");
  c->Print(path_pdf + "]");

  std::cout << "wrote " << path_pdf << std::endl;
  return 0;
}
