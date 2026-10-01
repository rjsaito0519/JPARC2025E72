// -*- C++ -*-
// Kp Geant4 true momentum-loss survey (p / K / beam).
// Check whether energy-loss correction looks necessary before designing a table.
//
// Usage:
//   kpsc_g4_eloss_trend <g4.root> [-o out.pdf]
//
// Input tree: g4hyptpc
// Output default:
//   OUTPUT_DIR/img/run00000/kpsc_misc/kpsc_g4_eloss_trend_run00000.pdf
//
// Loss definitions (G4 truth only; no helix / kpsc):
//   scattered p : PRM proton |p|  → first primary TPC proton hit |p|
//   scattered K : PRM K+/-  |p|   → first primary TPC K hit |p|
//   beam        : BEAM |p|        → TGT hit |p|  (generator==7201 beam-only)
// Units:
 //   Absolute momenta and #Deltap are shown in MeV/c.
 //   TParticle momenta in tree are MeV/c; TPC hit momenta are GeV/c
 //   (converted to MeV/c when used).
 //   delta = (p_end - p_start) / p_start   [dimensionless; 0.01 = 1%]
 //   Delta p = p_end - p_start             [MeV/c]

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TInterpreter.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{

constexpr Int_t kBeamGenerator = 7201;
constexpr Double_t kGeVtoMeV = 1.e3;
constexpr Int_t kRunOut = 0;

constexpr Double_t kDeltaLo = -0.40;
constexpr Double_t kDeltaHi = 0.10;
constexpr Double_t kDpPLo = -250.;  // MeV/c absolute for proton
constexpr Double_t kDpPHi = 50.;
constexpr Double_t kDpKLo = -100.;
constexpr Double_t kDpKHi = 50.;
constexpr Double_t kDpBLo = -30.;
constexpr Double_t kDpBHi = 10.;
constexpr Double_t kThetaLo = 0.0;
constexpr Double_t kThetaHi = 1.8;
constexpr Double_t kVtxR = 120.;   // mm display
constexpr Double_t kVtxZ = 220.;   // mm display
constexpr Double_t kP0PLo = 50.;   // MeV/c
constexpr Double_t kP0PHi = 900.;
constexpr Double_t kP0KLo = 100.;
constexpr Double_t kP0KHi = 900.;
constexpr Double_t kP0BLo = 600.;
constexpr Double_t kP0BHi = 850.;

void
usage(const char* a0)
{
  std::cerr
    << "Usage: " << a0 << " <g4.root> [-o out.pdf]\n"
    << "  tree: g4hyptpc\n"
    << "  Default PDF: OUTPUT_DIR/img/run00000/kpsc_misc/kpsc_g4_eloss_trend_run00000.pdf\n"
    << "  delta = (p_end - p_start) / p_start   [dimensionless; also shown in %]\n"
    << "  Delta p = p_end - p_start            [MeV/c]\n"
    << "  p/K  : PRM → TPC 1st primary hit (reaction events)\n"
    << "  beam : BEAM → TGT (beam-only generator==7201)\n";
}

class PdfWriter
{
public:
  explicit PdfWriter(TString path) : path_(std::move(path)) {}

  void Print(TCanvas& c)
  {
    ++page_;
    if (page_ == 1)
      c.Print(path_ + "(");
    else
      c.Print(path_);
  }

  void Close(TCanvas& c)
  {
    if (page_ <= 0)
      return;
    c.Print(path_ + "]");
  }

private:
  TString path_;
  Int_t page_ = 0;
};

Bool_t
MatchPdg(Int_t hit_pdg, Int_t want_pdg)
{
  if (want_pdg == 2212)
    return hit_pdg == 2212;
  if (std::abs(want_pdg) == 321)
    return std::abs(hit_pdg) == 321;
  return hit_pdg == want_pdg;
}

// First hit of the longest primary TPC track with the requested PDG.
Bool_t
FirstHitMom(const std::vector<int>* pid, const std::vector<int>* parentid,
            const std::vector<int>* trackid, const std::vector<double>* px,
            const std::vector<double>* py, const std::vector<double>* pz,
            Int_t want_pdg, TVector3& mom, TVector3& pos,
            const std::vector<double>* x0, const std::vector<double>* y0,
            const std::vector<double>* z0)
{
  mom.SetXYZ(0, 0, 0);
  pos.SetXYZ(0, 0, 0);
  if (!pid || !parentid || !trackid || !px || !py || !pz)
    return kFALSE;
  if (pid->size() != px->size())
    return kFALSE;

  std::map<int, int> cnt;
  std::map<int, size_t> first;
  for (size_t j = 0; j < pid->size(); ++j) {
    if ((*parentid)[j] != 0)
      continue;
    if (!MatchPdg((*pid)[j], want_pdg))
      continue;
    const int tid = (*trackid)[j];
    if (cnt.find(tid) == cnt.end())
      first[tid] = j;
    cnt[tid]++;
  }
  int best = -1;
  int bn = 0;
  for (const auto& kv : cnt) {
    if (kv.second > bn) {
      bn = kv.second;
      best = kv.first;
    }
  }
  if (best < 0)
    return kFALSE;
  const size_t j = first[best];
  mom.SetXYZ((*px)[j], (*py)[j], (*pz)[j]);
  if (x0 && y0 && z0 && j < x0->size())
    pos.SetXYZ((*x0)[j], (*y0)[j], (*z0)[j]);
  return mom.Mag() > 1e-9;
}

Double_t
MedianSorted(std::vector<Double_t> v)
{
  if (v.empty())
    return TMath::QuietNaN();
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  if (n % 2 == 1)
    return v[n / 2];
  return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

struct LossSample
{
  std::vector<Double_t> delta;
  std::vector<Double_t> dabs;  // MeV/c
  std::vector<Double_t> p0;    // MeV/c
  std::vector<Double_t> theta;
  std::vector<Double_t> vx;
  std::vector<Double_t> vy;
  std::vector<Double_t> vz;
  std::vector<Double_t> path;  // |r_end - r_start| mm

  void Add(Double_t d, Double_t da, Double_t p_start, Double_t th, Double_t x,
           Double_t y, Double_t z, Double_t L)
  {
    if (!std::isfinite(d) || !std::isfinite(da) || !(p_start > 0.))
      return;
    delta.push_back(d);
    dabs.push_back(da);
    p0.push_back(p_start);
    theta.push_back(th);
    vx.push_back(x);
    vy.push_back(y);
    vz.push_back(z);
    path.push_back(L);
  }

  Long64_t N() const { return static_cast<Long64_t>(delta.size()); }

  Double_t MedDelta() const { return MedianSorted(delta); }
  Double_t MedDabs() const { return MedianSorted(dabs); }
};

void
FillScatter(const LossSample& s, TH1D& h_d, TH1D& h_da, TProfile& pr_th,
            TProfile& pr_p0, TProfile& pr_vz, TProfile& pr_vx, TProfile& pr_path,
            TH2D& h_th, TProfile& pr_th_a, TProfile& pr_p0_a, TProfile& pr_vz_a,
            TProfile& pr_vx_a, TProfile& pr_path_a, TH2D& h_th_a)
{
  for (size_t i = 0; i < s.delta.size(); ++i) {
    h_d.Fill(s.delta[i]);
    h_da.Fill(s.dabs[i]);
    pr_th.Fill(s.theta[i], s.delta[i]);
    pr_p0.Fill(s.p0[i], s.delta[i]);
    pr_vz.Fill(s.vz[i], s.delta[i]);
    pr_vx.Fill(s.vx[i], s.delta[i]);
    pr_path.Fill(s.path[i], s.delta[i]);
    h_th.Fill(s.theta[i], s.delta[i]);
    pr_th_a.Fill(s.theta[i], s.dabs[i]);
    pr_p0_a.Fill(s.p0[i], s.dabs[i]);
    pr_vz_a.Fill(s.vz[i], s.dabs[i]);
    pr_vx_a.Fill(s.vx[i], s.dabs[i]);
    pr_path_a.Fill(s.path[i], s.dabs[i]);
    h_th_a.Fill(s.theta[i], s.dabs[i]);
  }
}

void
FillBeam(const LossSample& s, TH1D& h_d, TH1D& h_da, TProfile& pr_p0,
         TProfile& pr_vz, TProfile& pr_vx, TProfile& pr_path, TProfile& pr_p0_a,
         TProfile& pr_vz_a, TProfile& pr_vx_a, TProfile& pr_path_a)
{
  for (size_t i = 0; i < s.delta.size(); ++i) {
    h_d.Fill(s.delta[i]);
    h_da.Fill(s.dabs[i]);
    pr_p0.Fill(s.p0[i], s.delta[i]);
    pr_vz.Fill(s.vz[i], s.delta[i]);
    pr_vx.Fill(s.vx[i], s.delta[i]);
    pr_path.Fill(s.path[i], s.delta[i]);
    pr_p0_a.Fill(s.p0[i], s.dabs[i]);
    pr_vz_a.Fill(s.vz[i], s.dabs[i]);
    pr_vx_a.Fill(s.vx[i], s.dabs[i]);
    pr_path_a.Fill(s.path[i], s.dabs[i]);
  }
}

void
SetupPad(Bool_t with_z = kFALSE)
{
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(with_z ? 0.16 : 0.08);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.08);
}

// unit_mode: 0 = dimensionless delta (also print %), 1 = MeV/c absolute Delta p
void
DrawH1(TCanvas& c, PdfWriter& w, TH1D& h, const char* title, Double_t med,
       Int_t unit_mode)
{
  c.Clear();
  c.cd();
  SetupPad();
  h.SetTitle(title);
  h.GetXaxis()->SetTitleOffset(1.15);
  h.GetYaxis()->SetTitleOffset(1.35);
  h.SetLineWidth(2);
  h.Draw("hist");
  TLine ln(med, 0., med, h.GetMaximum() * 1.05);
  ln.SetLineColor(kRed);
  ln.SetLineStyle(2);
  ln.Draw("same");
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.032);
  if (unit_mode == 0) {
    lat.DrawLatex(0.50, 0.86, Form("median = %+.4f  (%+.2f%%)", med, 100. * med));
    lat.DrawLatex(0.50, 0.81,
                  Form("mean   = %+.4f  (%+.2f%%)", h.GetMean(),
                       100. * h.GetMean()));
  } else {
    lat.DrawLatex(0.50, 0.86, Form("median = %+.1f MeV/c", med));
    lat.DrawLatex(0.50, 0.81, Form("mean   = %+.1f MeV/c", h.GetMean()));
  }
  lat.DrawLatex(0.50, 0.76, Form("N      = %.0f", h.GetEntries()));
  w.Print(c);
}

void
DrawProf(TCanvas& c, PdfWriter& w, TProfile& pr, const char* title,
         Double_t y_lo = kDeltaLo, Double_t y_hi = kDeltaHi)
{
  c.Clear();
  c.cd();
  SetupPad();
  pr.SetTitle(title);
  pr.GetXaxis()->SetTitleOffset(1.15);
  pr.GetYaxis()->SetTitleOffset(1.35);
  pr.SetLineWidth(2);
  pr.SetMinimum(y_lo);
  pr.SetMaximum(y_hi);
  pr.Draw("hist");
  TLine z(pr.GetXaxis()->GetXmin(), 0., pr.GetXaxis()->GetXmax(), 0.);
  z.SetLineStyle(2);
  z.Draw("same");
  w.Print(c);
}

void
DrawTH2(TCanvas& c, PdfWriter& w, TH2D& h, const char* title)
{
  c.Clear();
  c.cd();
  SetupPad(kTRUE);
  h.SetTitle(title);
  h.GetXaxis()->SetTitleOffset(1.15);
  h.GetYaxis()->SetTitleOffset(1.35);
  h.Draw("colz");
  w.Print(c);
}

TString
ShortPath(const TString& path, Int_t max_len = 55)
{
  if (path.Length() <= max_len)
    return path;
  return TString("...") + path(path.Length() - max_len + 3, max_len - 3);
}

TString
VerdictLine(const char* name, Double_t med_delta, Double_t med_dabs)
{
  // Rough thresholds for "correction looks needed" discussion only.
  const Bool_t large = std::isfinite(med_delta) && std::abs(med_delta) > 0.02;
  const Bool_t mild = std::isfinite(med_delta) && std::abs(med_delta) > 0.005;
  const char* tag = large ? "LIKELY NEED" : (mild ? "MAYBE" : "SMALL");
  return Form("%s: med #delta=%+.4f (%+.2f%%), med #Deltap=%+.1f MeV/c [%s]",
              name, med_delta, 100. * med_delta, med_dabs, tag);
}

}  // namespace

int
main(int argc, char** argv)
{
  TString inPath;
  TString outPdf;
  for (int i = 1; i < argc; ++i) {
    const TString a = argv[i];
    if (a == "-h" || a == "--help") {
      usage(argv[0]);
      return 0;
    }
    if (a == "-o") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      outPdf = argv[++i];
      continue;
    }
    if (a.BeginsWith("-")) {
      std::cerr << "Unknown option: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
    if (inPath.IsNull())
      inPath = a;
    else {
      std::cerr << "Extra argument: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
  }
  if (inPath.IsNull()) {
    usage(argv[0]);
    return 1;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gInterpreter->GenerateDictionary("vector<TParticle>", "vector;TParticle.h");

  TFile* fin = TFile::Open(inPath, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open " << inPath << "\n";
    return 1;
  }
  auto* tree = dynamic_cast<TTree*>(fin->Get("g4hyptpc"));
  if (!tree) {
    std::cerr << "Error: tree g4hyptpc not found in " << inPath << "\n";
    return 1;
  }

  Int_t generator = 0;
  std::vector<TParticle>* BEAM = nullptr;
  std::vector<TParticle>* PRM = nullptr;
  std::vector<TParticle>* TGT = nullptr;
  std::vector<int>* pidtpc = nullptr;
  std::vector<int>* parentidtpc = nullptr;
  std::vector<int>* trackidtpc = nullptr;
  std::vector<double>* pxtpc = nullptr;
  std::vector<double>* pytpc = nullptr;
  std::vector<double>* pztpc = nullptr;
  std::vector<double>* x0tpc = nullptr;
  std::vector<double>* y0tpc = nullptr;
  std::vector<double>* z0tpc = nullptr;

  tree->SetBranchAddress("generator", &generator);
  tree->SetBranchAddress("BEAM", &BEAM);
  tree->SetBranchAddress("PRM", &PRM);
  tree->SetBranchAddress("TGT", &TGT);
  tree->SetBranchAddress("pidtpc", &pidtpc);
  tree->SetBranchAddress("parentidtpc", &parentidtpc);
  tree->SetBranchAddress("trackidtpc", &trackidtpc);
  tree->SetBranchAddress("pxtpc", &pxtpc);
  tree->SetBranchAddress("pytpc", &pytpc);
  tree->SetBranchAddress("pztpc", &pztpc);
  tree->SetBranchAddress("x0tpc", &x0tpc);
  tree->SetBranchAddress("y0tpc", &y0tpc);
  tree->SetBranchAddress("z0tpc", &z0tpc);

  LossSample samp_p, samp_k, samp_beam;

  const Long64_t nent = tree->GetEntries();
  Long64_t n_reac = 0, n_beam = 0;
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tree->GetEntry(ie);

    if (generator == kBeamGenerator) {
      ++n_beam;
      if (!BEAM || BEAM->empty() || !TGT || TGT->empty())
        continue;
      const TParticle& b0 = (*BEAM)[0];
      const TParticle& tg = (*TGT)[0];
      // TParticle::P() is MeV/c
      const Double_t p_start = b0.P();
      const Double_t p_end = tg.P();
      if (!(p_start > 1e-3) || !(p_end > 1e-3))
        continue;
      const Double_t d = (p_end - p_start) / p_start;
      const Double_t da = p_end - p_start;  // MeV/c
      const TVector3 v0(b0.Vx(), b0.Vy(), b0.Vz());
      const TVector3 v1(tg.Vx(), tg.Vy(), tg.Vz());
      const TVector3 pb(b0.Px(), b0.Py(), b0.Pz());
      samp_beam.Add(d, da, p_start, pb.Theta(), tg.Vx(), tg.Vy(), tg.Vz(),
                    (v1 - v0).Mag());
      continue;
    }

    // reaction
    ++n_reac;
    if (!PRM || PRM->empty())
      continue;

    TVector3 p_prm(0, 0, 0), v_prm(0, 0, 0);
    TVector3 k_prm(0, 0, 0);
    Bool_t has_p = kFALSE, has_k = kFALSE;
    for (const auto& part : *PRM) {
      if (part.GetPdgCode() == 2212 && !has_p) {
        // PRM momenta are MeV/c
        p_prm.SetXYZ(part.Px(), part.Py(), part.Pz());
        v_prm.SetXYZ(part.Vx(), part.Vy(), part.Vz());
        has_p = kTRUE;
      }
      if (std::abs(part.GetPdgCode()) == 321 && !has_k) {
        k_prm.SetXYZ(part.Px(), part.Py(), part.Pz());
        if (!has_p)
          v_prm.SetXYZ(part.Vx(), part.Vy(), part.Vz());
        has_k = kTRUE;
      }
    }

    if (has_p && p_prm.Mag() > 1e-3) {
      TVector3 p1, r1;
      if (FirstHitMom(pidtpc, parentidtpc, trackidtpc, pxtpc, pytpc, pztpc, 2212,
                      p1, r1, x0tpc, y0tpc, z0tpc)) {
        // TPC hit momenta are GeV/c -> MeV/c
        p1 *= kGeVtoMeV;
        const Double_t p0 = p_prm.Mag();
        const Double_t d = (p1.Mag() - p0) / p0;
        samp_p.Add(d, p1.Mag() - p0, p0, p_prm.Theta(), v_prm.X(), v_prm.Y(),
                   v_prm.Z(), (r1 - v_prm).Mag());
      }
    }
    if (has_k && k_prm.Mag() > 1e-3) {
      TVector3 p1, r1;
      if (FirstHitMom(pidtpc, parentidtpc, trackidtpc, pxtpc, pytpc, pztpc, -321,
                      p1, r1, x0tpc, y0tpc, z0tpc)) {
        p1 *= kGeVtoMeV;
        const Double_t p0 = k_prm.Mag();
        const Double_t d = (p1.Mag() - p0) / p0;
        samp_k.Add(d, p1.Mag() - p0, p0, k_prm.Theta(), v_prm.X(), v_prm.Y(),
                   v_prm.Z(), (r1 - v_prm).Mag());
      }
    }
  }

  // default outputs go to the kpsc_misc/ sub-directory of img/runNNNNN
  const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, kRunOut) + "/kpsc_misc";
  gSystem->mkdir(imgDir.Data(), kTRUE);
  if (outPdf.IsNull())
    outPdf = Form("%s/kpsc_g4_eloss_trend_run%05d.pdf", imgDir.Data(), kRunOut);
  TString outTxt = outPdf;
  if (outTxt.EndsWith(".pdf") || outTxt.EndsWith(".PDF"))
    outTxt.Replace(outTxt.Length() - 4, 4, ".txt");
  else
    outTxt += ".txt";

  // Axis convention:
  //   #delta = (p_end-p_start)/p_start   [dimensionless]
  //   #Deltap = p_end-p_start            [MeV/c]
  TH1D h_dp("h_dp",
            "scattered p;#delta=(p_{1st}-p_{PRM})/p_{PRM} [dimensionless];"
            "entries",
            80, kDeltaLo, kDeltaHi);
  TH1D h_dk("h_dk",
            "scattered K;#delta=(p_{1st}-p_{PRM})/p_{PRM} [dimensionless];"
            "entries",
            80, kDeltaLo, kDeltaHi);
  TH1D h_db("h_db",
            "beam;#delta=(p_{TGT}-p_{BEAM})/p_{BEAM} [dimensionless];entries",
            80, -0.05, 0.02);
  TH1D h_dap("h_dap", "scattered p;#Deltap = p_{1st}-p_{PRM} [MeV/c];entries",
            80, kDpPLo, kDpPHi);
  TH1D h_dak("h_dak", "scattered K;#Deltap = p_{1st}-p_{PRM} [MeV/c];entries",
            80, kDpKLo, kDpKHi);
  TH1D h_dab("h_dab", "beam;#Deltap = p_{TGT}-p_{BEAM} [MeV/c];entries", 80,
            kDpBLo, kDpBHi);

  TProfile pr_p_th("pr_p_th",
                   "p;#theta_{lab} [rad];#delta = #Deltap/p [dimensionless]",
                   18, kThetaLo, kThetaHi);
  TProfile pr_k_th("pr_k_th",
                   "K;#theta_{lab} [rad];#delta = #Deltap/p [dimensionless]",
                   18, kThetaLo, kThetaHi);
  TProfile pr_p_p0("pr_p_p0",
                   "p;p_{PRM} [MeV/c];#delta = #Deltap/p [dimensionless]", 20,
                   kP0PLo, kP0PHi);
  TProfile pr_k_p0("pr_k_p0",
                   "K;p_{PRM} [MeV/c];#delta = #Deltap/p [dimensionless]", 20,
                   kP0KLo, kP0KHi);
  TProfile pr_b_p0("pr_b_p0",
                   "beam;p_{BEAM} [MeV/c];#delta = #Deltap/p [dimensionless]",
                   20, kP0BLo, kP0BHi);
  TProfile pr_p_vz("pr_p_vz",
                   "p;vtx z [mm];#delta = #Deltap/p [dimensionless]", 22,
                   -kVtxZ, kVtxZ);
  TProfile pr_k_vz("pr_k_vz",
                   "K;vtx z [mm];#delta = #Deltap/p [dimensionless]", 22,
                   -kVtxZ, kVtxZ);
  TProfile pr_b_vz("pr_b_vz",
                   "beam;TGT z [mm];#delta = #Deltap/p [dimensionless]", 22,
                   -kVtxZ, kVtxZ);
  TProfile pr_p_vx("pr_p_vx",
                   "p;vtx x [mm];#delta = #Deltap/p [dimensionless]", 24,
                   -kVtxR, kVtxR);
  TProfile pr_k_vx("pr_k_vx",
                   "K;vtx x [mm];#delta = #Deltap/p [dimensionless]", 24,
                   -kVtxR, kVtxR);
  TProfile pr_b_vx("pr_b_vx",
                   "beam;TGT x [mm];#delta = #Deltap/p [dimensionless]", 24,
                   -kVtxR, kVtxR);
  TProfile pr_p_path("pr_p_path",
                     "p;|r_{1st}-v_{PRM}| [mm];#delta = #Deltap/p "
                     "[dimensionless]",
                     20, 0, 250);
  TProfile pr_k_path("pr_k_path",
                     "K;|r_{1st}-v_{PRM}| [mm];#delta = #Deltap/p "
                     "[dimensionless]",
                     20, 0, 250);
  TProfile pr_b_path("pr_b_path",
                     "beam;|r_{TGT}-r_{BEAM}| [mm];#delta = #Deltap/p "
                     "[dimensionless]",
                     20, 0, 1200);

  // Absolute #Deltap [MeV/c] profiles (same x axes)
  TProfile pr_p_th_a("pr_p_th_a", "p;#theta_{lab} [rad];#Deltap [MeV/c]", 18,
                     kThetaLo, kThetaHi);
  TProfile pr_k_th_a("pr_k_th_a", "K;#theta_{lab} [rad];#Deltap [MeV/c]", 18,
                     kThetaLo, kThetaHi);
  TProfile pr_p_p0_a("pr_p_p0_a", "p;p_{PRM} [MeV/c];#Deltap [MeV/c]", 20,
                     kP0PLo, kP0PHi);
  TProfile pr_k_p0_a("pr_k_p0_a", "K;p_{PRM} [MeV/c];#Deltap [MeV/c]", 20,
                     kP0KLo, kP0KHi);
  TProfile pr_b_p0_a("pr_b_p0_a", "beam;p_{BEAM} [MeV/c];#Deltap [MeV/c]", 20,
                     kP0BLo, kP0BHi);
  TProfile pr_p_vz_a("pr_p_vz_a", "p;vtx z [mm];#Deltap [MeV/c]", 22, -kVtxZ,
                     kVtxZ);
  TProfile pr_k_vz_a("pr_k_vz_a", "K;vtx z [mm];#Deltap [MeV/c]", 22, -kVtxZ,
                     kVtxZ);
  TProfile pr_b_vz_a("pr_b_vz_a", "beam;TGT z [mm];#Deltap [MeV/c]", 22, -kVtxZ,
                     kVtxZ);
  TProfile pr_p_vx_a("pr_p_vx_a", "p;vtx x [mm];#Deltap [MeV/c]", 24, -kVtxR,
                     kVtxR);
  TProfile pr_k_vx_a("pr_k_vx_a", "K;vtx x [mm];#Deltap [MeV/c]", 24, -kVtxR,
                     kVtxR);
  TProfile pr_b_vx_a("pr_b_vx_a", "beam;TGT x [mm];#Deltap [MeV/c]", 24, -kVtxR,
                     kVtxR);
  TProfile pr_p_path_a("pr_p_path_a",
                       "p;|r_{1st}-v_{PRM}| [mm];#Deltap [MeV/c]", 20, 0, 250);
  TProfile pr_k_path_a("pr_k_path_a",
                       "K;|r_{1st}-v_{PRM}| [mm];#Deltap [MeV/c]", 20, 0, 250);
  TProfile pr_b_path_a("pr_b_path_a",
                       "beam;|r_{TGT}-r_{BEAM}| [mm];#Deltap [MeV/c]", 20, 0,
                       1200);

  TH2D h2_p_th("h2_p_th",
               "p;#theta_{lab} [rad];#delta = #Deltap/p [dimensionless]", 36,
               kThetaLo, kThetaHi, 60, kDeltaLo, kDeltaHi);
  TH2D h2_k_th("h2_k_th",
               "K;#theta_{lab} [rad];#delta = #Deltap/p [dimensionless]", 36,
               kThetaLo, kThetaHi, 60, kDeltaLo, kDeltaHi);
  TH2D h2_p_th_a("h2_p_th_a", "p;#theta_{lab} [rad];#Deltap [MeV/c]", 36,
                 kThetaLo, kThetaHi, 60, kDpPLo, kDpPHi);
  TH2D h2_k_th_a("h2_k_th_a", "K;#theta_{lab} [rad];#Deltap [MeV/c]", 36,
                 kThetaLo, kThetaHi, 60, kDpKLo, kDpKHi);

  FillScatter(samp_p, h_dp, h_dap, pr_p_th, pr_p_p0, pr_p_vz, pr_p_vx, pr_p_path,
              h2_p_th, pr_p_th_a, pr_p_p0_a, pr_p_vz_a, pr_p_vx_a, pr_p_path_a,
              h2_p_th_a);
  FillScatter(samp_k, h_dk, h_dak, pr_k_th, pr_k_p0, pr_k_vz, pr_k_vx, pr_k_path,
              h2_k_th, pr_k_th_a, pr_k_p0_a, pr_k_vz_a, pr_k_vx_a, pr_k_path_a,
              h2_k_th_a);
  FillBeam(samp_beam, h_db, h_dab, pr_b_p0, pr_b_vz, pr_b_vx, pr_b_path,
           pr_b_p0_a, pr_b_vz_a, pr_b_vx_a, pr_b_path_a);

  const Double_t med_p = samp_p.MedDelta();
  const Double_t med_k = samp_k.MedDelta();
  const Double_t med_b = samp_beam.MedDelta();
  const Double_t med_pa = samp_p.MedDabs();
  const Double_t med_ka = samp_k.MedDabs();
  const Double_t med_ba = samp_beam.MedDabs();

  {
    std::ofstream ofs(outTxt.Data());
    ofs << "# kpsc_g4_eloss_trend\n";
    ofs << "# input: " << inPath << "\n";
    ofs << "# entries: " << nent << "  reaction_loops: " << n_reac
        << "  beam_loops: " << n_beam << "\n";
    ofs << "# definitions:\n";
    ofs << "#   p/K: PRM |p| -> TPC 1st primary hit |p|\n";
    ofs << "#   beam: BEAM |p| -> TGT |p| (generator==7201)\n";
    ofs << "# delta=(p_end-p_start)/p_start  [dimensionless]\n";
    ofs << "# Delta_p=p_end-p_start          [MeV/c]\n";
    ofs << "N_p " << samp_p.N() << " median_delta " << med_p
        << " median_delta_percent " << (100. * med_p) << " median_Dp_MeV "
        << med_pa << "\n";
    ofs << "N_K " << samp_k.N() << " median_delta " << med_k
        << " median_delta_percent " << (100. * med_k) << " median_Dp_MeV "
        << med_ka << "\n";
    ofs << "N_beam " << samp_beam.N() << " median_delta " << med_b
        << " median_delta_percent " << (100. * med_b) << " median_Dp_MeV "
        << med_ba << "\n";
    ofs << VerdictLine("proton", med_p, med_pa) << "\n";
    ofs << VerdictLine("kaon", med_k, med_ka) << "\n";
    ofs << VerdictLine("beam", med_b, med_ba) << "\n";
    ofs << "# note: LIKELY NEED if |median delta|>2%; MAYBE if >0.5%; "
           "thresholds are discussion aids only\n";
  }

  TCanvas c("c_eloss", "kpsc g4 eloss trend", 900, 700);
  PdfWriter writer(outPdf);

  // title / summary
  {
    c.Clear();
    c.cd();
    SetupPad();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.10, 0.90, "Kp Geant4 true momentum loss survey");
    lat.SetTextSize(0.026);
    lat.DrawLatex(0.10, 0.84, Form("input: %s", ShortPath(inPath).Data()));
    lat.DrawLatex(0.10, 0.79,
                  Form("tree g4hyptpc  entries=%lld  reac~%lld  beam~%lld", nent,
                       n_reac, n_beam));
    lat.DrawLatex(0.10, 0.72, "Units:");
    lat.DrawLatex(0.12, 0.67,
                  "#delta = (p_{end}-p_{start})/p_{start}   "
                  "[dimensionless; 0.01 = 1%]");
    lat.DrawLatex(0.12, 0.62, "#Deltap = p_{end}-p_{start}   [MeV/c]");
    lat.DrawLatex(0.10, 0.55, "Definitions:");
    lat.DrawLatex(0.12, 0.50, "p / K: PRM |p|  ->  TPC 1st primary hit |p|");
    lat.DrawLatex(0.12, 0.45, "beam: BEAM |p|  ->  TGT |p|  (generator==7201)");
    lat.DrawLatex(0.10, 0.38, "Results (median):");
    lat.SetTextSize(0.024);
    lat.DrawLatex(0.12, 0.33, VerdictLine("proton", med_p, med_pa));
    lat.DrawLatex(0.12, 0.28, VerdictLine("kaon  ", med_k, med_ka));
    lat.DrawLatex(0.12, 0.23, VerdictLine("beam  ", med_b, med_ba));
    lat.SetTextSize(0.022);
    lat.DrawLatex(0.10, 0.16,
                  "LIKELY NEED: |med #delta|>2%   MAYBE: >0.5%   (rough guide)");
    lat.DrawLatex(0.10, 0.11,
                  Form("PDF: %s", ShortPath(outPdf, 60).Data()));
    writer.Print(c);
  }

  DrawH1(c, writer, h_dp, "scattered proton #delta", med_p, 0);
  DrawH1(c, writer, h_dk, "scattered kaon #delta", med_k, 0);
  DrawH1(c, writer, h_db, "beam #delta (BEAM #rightarrow TGT)", med_b, 0);

  DrawH1(c, writer, h_dap, "scattered proton #Deltap", med_pa, 1);
  DrawH1(c, writer, h_dak, "scattered kaon #Deltap", med_ka, 1);
  DrawH1(c, writer, h_dab, "beam #Deltap (BEAM #rightarrow TGT)", med_ba, 1);

  DrawProf(c, writer, pr_p_th, "proton <#delta> vs #theta_{lab}");
  DrawProf(c, writer, pr_k_th, "kaon <#delta> vs #theta_{lab}");
  DrawTH2(c, writer, h2_p_th, "proton #delta vs #theta_{lab}");
  DrawTH2(c, writer, h2_k_th, "kaon #delta vs #theta_{lab}");

  DrawProf(c, writer, pr_p_th_a, "proton <#Deltap> vs #theta_{lab} [MeV/c]",
           kDpPLo, kDpPHi);
  DrawProf(c, writer, pr_k_th_a, "kaon <#Deltap> vs #theta_{lab} [MeV/c]",
           kDpKLo, kDpKHi);
  DrawTH2(c, writer, h2_p_th_a, "proton #Deltap vs #theta_{lab} [MeV/c]");
  DrawTH2(c, writer, h2_k_th_a, "kaon #Deltap vs #theta_{lab} [MeV/c]");

  DrawProf(c, writer, pr_p_p0, "proton <#delta> vs p_{PRM}");
  DrawProf(c, writer, pr_k_p0, "kaon <#delta> vs p_{PRM}");
  DrawProf(c, writer, pr_b_p0, "beam <#delta> vs p_{BEAM}", -0.05, 0.02);
  DrawProf(c, writer, pr_p_p0_a, "proton <#Deltap> vs p_{PRM} [MeV/c]", kDpPLo,
           kDpPHi);
  DrawProf(c, writer, pr_k_p0_a, "kaon <#Deltap> vs p_{PRM} [MeV/c]", kDpKLo,
           kDpKHi);
  DrawProf(c, writer, pr_b_p0_a, "beam <#Deltap> vs p_{BEAM} [MeV/c]", kDpBLo,
           kDpBHi);

  DrawProf(c, writer, pr_p_vz, "proton <#delta> vs vtx z");
  DrawProf(c, writer, pr_k_vz, "kaon <#delta> vs vtx z");
  DrawProf(c, writer, pr_b_vz, "beam <#delta> vs TGT z", -0.05, 0.02);
  DrawProf(c, writer, pr_p_vz_a, "proton <#Deltap> vs vtx z [MeV/c]", kDpPLo,
           kDpPHi);
  DrawProf(c, writer, pr_k_vz_a, "kaon <#Deltap> vs vtx z [MeV/c]", kDpKLo,
           kDpKHi);
  DrawProf(c, writer, pr_b_vz_a, "beam <#Deltap> vs TGT z [MeV/c]", kDpBLo,
           kDpBHi);

  DrawProf(c, writer, pr_p_vx, "proton <#delta> vs vtx x");
  DrawProf(c, writer, pr_k_vx, "kaon <#delta> vs vtx x");
  DrawProf(c, writer, pr_b_vx, "beam <#delta> vs TGT x", -0.05, 0.02);
  DrawProf(c, writer, pr_p_vx_a, "proton <#Deltap> vs vtx x [MeV/c]", kDpPLo,
           kDpPHi);
  DrawProf(c, writer, pr_k_vx_a, "kaon <#Deltap> vs vtx x [MeV/c]", kDpKLo,
           kDpKHi);
  DrawProf(c, writer, pr_b_vx_a, "beam <#Deltap> vs TGT x [MeV/c]", kDpBLo,
           kDpBHi);

  DrawProf(c, writer, pr_p_path, "proton <#delta> vs |r_{1st}-v|");
  DrawProf(c, writer, pr_k_path, "kaon <#delta> vs |r_{1st}-v|");
  DrawProf(c, writer, pr_b_path, "beam <#delta> vs path BEAM#rightarrowTGT",
           -0.05, 0.02);
  DrawProf(c, writer, pr_p_path_a, "proton <#Deltap> vs |r_{1st}-v| [MeV/c]",
           kDpPLo, kDpPHi);
  DrawProf(c, writer, pr_k_path_a, "kaon <#Deltap> vs |r_{1st}-v| [MeV/c]",
           kDpKLo, kDpKHi);
  DrawProf(c, writer, pr_b_path_a,
           "beam <#Deltap> vs path BEAM#rightarrowTGT [MeV/c]", kDpBLo, kDpBHi);

  // closing summary page (correction need)
  {
    c.Clear();
    c.cd();
    SetupPad();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.10, 0.90, "Correction-need checklist");
    lat.SetTextSize(0.024);
    lat.DrawLatex(0.10, 0.82, VerdictLine("proton", med_p, med_pa));
    lat.DrawLatex(0.10, 0.76, VerdictLine("kaon  ", med_k, med_ka));
    lat.DrawLatex(0.10, 0.70, VerdictLine("beam  ", med_b, med_ba));
    lat.SetTextSize(0.026);
    lat.DrawLatex(0.10, 0.60, "How to read profiles:");
    lat.SetTextSize(0.024);
    lat.DrawLatex(0.12, 0.54,
                  "- Strong #theta or vtx dependence #rightarrow path / table");
    lat.DrawLatex(0.12, 0.48,
                  "- Small global offset only #rightarrow single scale, or ignore");
    lat.DrawLatex(0.12, 0.42,
                  "- Beam BEAM#rightarrowTGT is upstream-of-vertex loss only");
    lat.DrawLatex(0.10, 0.32, Form("N: p=%lld  K=%lld  beam=%lld", samp_p.N(),
                                  samp_k.N(), samp_beam.N()));
    lat.DrawLatex(0.10, 0.26,
                  Form("summary txt: %s", ShortPath(outTxt, 60).Data()));
    writer.Print(c);
  }

  writer.Close(c);

  std::cout << "N_p=" << samp_p.N() << " median_delta=" << med_p
            << " median_Dp=" << med_pa << "\n";
  std::cout << "N_K=" << samp_k.N() << " median_delta=" << med_k
            << " median_Dp=" << med_ka << "\n";
  std::cout << "N_beam=" << samp_beam.N() << " median_delta=" << med_b
            << " median_Dp=" << med_ba << "\n";
  std::cout << "wrote " << outPdf << "\n";
  std::cout << "wrote " << outTxt << "\n";

  fin->Close();
  return 0;
}
