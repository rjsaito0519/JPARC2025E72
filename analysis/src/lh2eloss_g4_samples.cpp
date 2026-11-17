// -*- C++ -*-
// Write the LH2 momentum-loss samples of a G4 single-particle gun (FirstGenerator 7216) for
// analysis/scripts/lh2eloss/build_lh2_eloss_map.py.
//
// Usage:
//   lh2eloss_g4_samples <g4.root> [g4.root ...] -o samples.root --pdg c1,c2,...
//                       [--holdout-frac f] [--min-hits N]
//
// One sample per primary of the given PDG codes that reaches the TPC (TTree "samples"):
//   p0    |p| at the vertex (PRM) [MeV/c]
//   p1    |p| at the first TPC hit of that primary [MeV/c]
//   loss  p0 - p1 [MeV/c] (samples with loss < 0 are dropped)
//   vx vy vz [mm], ux uy uz: vertex and unit direction at the vertex
//   L, cap  chord from the vertex through the liquid cylinder [mm] (0 if the vertex is outside it)
//           and its exit surface (0 side, 1 end cap)
//   pdg, nhit (TPC hits of the primary), file (index of the input file)
//   holdout 1 for a fixed fraction of samples (hash of the rounded kinematics), 0 otherwise
// Samples with identical rounded kinematics (p0, p1, vertex, L, loss, cap at 0.01) are written once.

#include <TFile.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <unordered_set>
#include <vector>

namespace
{

constexpr Int_t kBeamGenerator = 7201;
constexpr Double_t kGeVtoMeV = 1.e3;
constexpr Double_t kTargetRadius = 40.0;      // mm (XZ)
constexpr Double_t kTargetHalfHeight = 50.0;  // mm (along Y)
constexpr Double_t kTargetCenterZ = -143.0;   // mm; centre x = y = 0
constexpr Double_t kEps = 1.e-12;

struct Options
{
  std::vector<TString> inputs;
  TString output;
  std::vector<Int_t> pdg;
  Double_t holdout_frac = 0.3;
  Int_t min_hits = 0;
};

struct Sample
{
  Float_t p0 = 0, p1 = 0, L = 0, loss = 0;
  Float_t vx = 0, vy = 0, vz = 0, ux = 0, uy = 0, uz = 0;
  Int_t pdg = 0, cap = 0, holdout = 0, file = 0, nhit = 0;
};

struct TpcHits
{
  std::vector<int>* pid = nullptr;
  std::vector<int>* parentid = nullptr;
  std::vector<int>* trackid = nullptr;
  std::vector<double>* px = nullptr;
  std::vector<double>* py = nullptr;
  std::vector<double>* pz = nullptr;
};

void
Usage(const char* a0)
{
  std::cerr << "Usage: " << a0 << " <g4.root> [g4.root ...] -o samples.root --pdg c1,c2,..."
            << " [--holdout-frac f (0.3)] [--min-hits N (0)]\n";
}

std::vector<TString>
SplitComma(const TString& s)
{
  std::vector<TString> out;
  std::unique_ptr<TObjArray> tok(s.Tokenize(","));
  for (Int_t k = 0; k < tok->GetEntries(); ++k)
    out.push_back(static_cast<TObjString*>(tok->At(k))->GetString());
  return out;
}

Bool_t
ParseArgs(int argc, char** argv, Options& o)
{
  for (int i = 1; i < argc; ++i) {
    const TString a = argv[i];
    const Bool_t has_value = i + 1 < argc;
    if (a == "-o" && has_value)
      o.output = argv[++i];
    else if (a == "--pdg" && has_value)
      for (const auto& t : SplitComma(argv[++i]))
        o.pdg.push_back(t.Atoi());
    else if (a == "--holdout-frac" && has_value)
      o.holdout_frac = TString(argv[++i]).Atof();
    else if (a == "--min-hits" && has_value)
      o.min_hits = TString(argv[++i]).Atoi();
    else if (a.BeginsWith("-"))
      return kFALSE;
    else
      o.inputs.push_back(a);
  }
  if (o.inputs.empty() || o.output.IsNull() || o.pdg.empty())
    return kFALSE;
  if (!(o.holdout_frac >= 0.) || o.holdout_frac >= 1.) {
    std::cerr << "--holdout-frac must be in [0, 1)\n";
    return kFALSE;
  }
  return kTRUE;
}

// Momentum [GeV/c] at the first TPC hit of the primary (parent 0) with the most hits of this PDG code.
Bool_t
FirstHitMom(const TpcHits& h, Int_t pdg, TVector3& mom, Int_t& nhit)
{
  if (!h.pid || h.pid->size() != h.px->size())
    return kFALSE;
  std::map<int, int> count;
  std::map<int, size_t> first;
  for (size_t j = 0; j < h.pid->size(); ++j) {
    if ((*h.parentid)[j] != 0 || (*h.pid)[j] != pdg)
      continue;
    const int tid = (*h.trackid)[j];
    if (count.find(tid) == count.end())
      first[tid] = j;
    ++count[tid];
  }
  int best = -1;
  nhit = 0;
  for (const auto& kv : count) {
    if (kv.second > nhit) {
      nhit = kv.second;
      best = kv.first;
    }
  }
  if (best < 0)
    return kFALSE;
  const size_t j = first[best];
  mom.SetXYZ((*h.px)[j], (*h.py)[j], (*h.pz)[j]);
  return mom.Mag() > 1e-9;
}

// Chord from vtx along the unit vector u to where it leaves the liquid cylinder (axis // y);
// 0 if vtx is outside the cylinder.
Double_t
ChordInLiquid(const TVector3& vtx, const TVector3& u, Bool_t& exit_cap)
{
  exit_cap = kFALSE;
  const Double_t R = kTargetRadius;
  const Double_t H = kTargetHalfHeight;
  const Double_t x0 = vtx.X();
  const Double_t y0 = vtx.Y();
  const Double_t z0 = vtx.Z() - kTargetCenterZ;
  const Double_t rho2 = x0 * x0 + z0 * z0;
  if (!(rho2 <= R * R + 1.e-9 && std::abs(y0) <= H + 1.e-9))
    return 0.;
  Double_t t_exit = 1.e300;
  Bool_t cap = kFALSE;
  const Double_t a = u.X() * u.X() + u.Z() * u.Z();
  const Double_t b = 2. * (x0 * u.X() + z0 * u.Z());
  const Double_t c = rho2 - R * R;
  const Double_t disc = b * b - 4. * a * c;
  if (a > kEps && disc >= 0.) {
    const Double_t s = TMath::Sqrt(disc);
    for (Double_t t : {(-b - s) / (2. * a), (-b + s) / (2. * a)}) {
      if (t > kEps && std::abs(y0 + t * u.Y()) <= H + 1.e-6 && t < t_exit)
        t_exit = t;
    }
  }
  if (std::abs(u.Y()) > kEps) {
    for (Double_t ycap : {H, -H}) {
      const Double_t t = (ycap - y0) / u.Y();
      const Double_t x = x0 + t * u.X();
      const Double_t z = z0 + t * u.Z();
      if (t > kEps && x * x + z * z <= R * R + 1.e-6 && t < t_exit) {
        t_exit = t;
        cap = kTRUE;
      }
    }
  }
  if (!(t_exit < 1.e299))
    return 0.;
  exit_cap = cap;
  return t_exit;
}

// Kinematics of one sample in double precision (the fingerprint rounds these, not the stored floats).
struct Kinematics
{
  Double_t p0 = 0, p1 = 0, vx = 0, vy = 0, vz = 0, L = 0, loss = 0;
  Int_t cap = 0;
};

// FNV-1a over the kinematics rounded to 0.01 (mm, MeV/c): duplicates and the holdout flag.
std::size_t
Fingerprint(const Kinematics& k)
{
  std::size_t h = 1469598103934665603ull;
  const auto mix = [&h](Int_t v) {
    h ^= static_cast<std::size_t>(static_cast<std::uint32_t>(v));
    h *= 1099511628211ull;
  };
  for (Double_t x : {k.p0, k.p1, k.vx, k.vy, k.vz, k.L, k.loss})
    mix(static_cast<Int_t>(std::llround(x * 100.)));
  mix(k.cap);
  return h;
}

// Fill s from one primary; false if it does not give a sample.
Bool_t
MakeSample(const TParticle& part, const TpcHits& hits, Int_t min_hits, Sample& s, Kinematics& k)
{
  const TVector3 p_prm(part.Px(), part.Py(), part.Pz());
  if (!(p_prm.Mag() > 1e-3))
    return kFALSE;
  TVector3 p1;
  Int_t nhit = 0;
  if (!FirstHitMom(hits, part.GetPdgCode(), p1, nhit) || nhit < min_hits)
    return kFALSE;
  p1 *= kGeVtoMeV;
  const Double_t p0 = p_prm.Mag();
  const Double_t mag1 = p1.Mag();
  const Double_t loss = p0 - mag1;
  if (!(loss >= 0.) || !std::isfinite(loss))
    return kFALSE;
  const TVector3 v(part.Vx(), part.Vy(), part.Vz());
  const TVector3 u = p_prm.Unit();
  Bool_t cap = kFALSE;
  const Double_t L = ChordInLiquid(v, u, cap);
  k = {p0, mag1, v.X(), v.Y(), v.Z(), L, loss, cap ? 1 : 0};
  s.p0 = p0; s.p1 = mag1; s.L = L; s.loss = loss;
  s.vx = v.X(); s.vy = v.Y(); s.vz = v.Z();
  s.ux = u.X(); s.uy = u.Y(); s.uz = u.Z();
  s.pdg = part.GetPdgCode();
  s.cap = cap ? 1 : 0;
  s.nhit = nhit;
  return kTRUE;
}

}  // namespace

int
main(int argc, char** argv)
{
  Options opt;
  if (!ParseArgs(argc, argv, opt)) {
    Usage(argv[0]);
    return 1;
  }
  std::unique_ptr<TFile> fout(TFile::Open(opt.output, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    std::cerr << "cannot write " << opt.output << "\n";
    return 1;
  }
  auto* out = new TTree("samples", "LH2 momentum-loss samples of a G4 particle gun");
  Sample s;
  out->Branch("p0", &s.p0);
  out->Branch("p1", &s.p1);
  out->Branch("L", &s.L);
  out->Branch("loss", &s.loss);
  out->Branch("vx", &s.vx);
  out->Branch("vy", &s.vy);
  out->Branch("vz", &s.vz);
  out->Branch("ux", &s.ux);
  out->Branch("uy", &s.uy);
  out->Branch("uz", &s.uz);
  out->Branch("pdg", &s.pdg);
  out->Branch("cap", &s.cap);
  out->Branch("holdout", &s.holdout);
  out->Branch("file", &s.file);
  out->Branch("nhit", &s.nhit);

  std::unordered_set<std::size_t> seen;
  const auto n_holdout = static_cast<std::size_t>(opt.holdout_frac * 10000. + 0.5);
  Long64_t n_dup = 0;
  for (size_t fi = 0; fi < opt.inputs.size(); ++fi) {
    std::unique_ptr<TFile> fin(TFile::Open(opt.inputs[fi], "READ"));
    auto* tree = fin && !fin->IsZombie() ? dynamic_cast<TTree*>(fin->Get("g4hyptpc")) : nullptr;
    if (!tree) {
      std::cerr << "cannot read tree g4hyptpc from " << opt.inputs[fi] << "\n";
      return 1;
    }
    Int_t generator = 0;
    std::vector<TParticle>* prm = nullptr;
    TpcHits hits;
    tree->SetBranchStatus("*", 0);
    for (const char* bn : {"generator", "PRM", "pidtpc", "parentidtpc", "trackidtpc", "pxtpc", "pytpc", "pztpc"})
      tree->SetBranchStatus(bn, 1);
    tree->SetBranchAddress("generator", &generator);
    tree->SetBranchAddress("PRM", &prm);
    tree->SetBranchAddress("pidtpc", &hits.pid);
    tree->SetBranchAddress("parentidtpc", &hits.parentid);
    tree->SetBranchAddress("trackidtpc", &hits.trackid);
    tree->SetBranchAddress("pxtpc", &hits.px);
    tree->SetBranchAddress("pytpc", &hits.py);
    tree->SetBranchAddress("pztpc", &hits.pz);
    for (Long64_t ie = 0; ie < tree->GetEntries(); ++ie) {
      tree->GetEntry(ie);
      if (generator == kBeamGenerator || !prm)
        continue;
      for (const auto& part : *prm) {
        if (std::find(opt.pdg.begin(), opt.pdg.end(), part.GetPdgCode()) == opt.pdg.end())
          continue;
        Kinematics k;
        if (!MakeSample(part, hits, opt.min_hits, s, k))
          continue;
        const std::size_t fp = Fingerprint(k);
        if (!seen.insert(fp).second) {
          ++n_dup;
          continue;
        }
        s.holdout = (fp % 10000ull) < n_holdout ? 1 : 0;
        s.file = static_cast<Int_t>(fi);
        out->Fill();
      }
    }
  }
  fout->cd();
  out->Write();
  std::cout << "wrote " << out->GetEntries() << " samples (" << n_dup << " duplicates skipped) to "
            << opt.output << "\n";
  fout->Close();
  return 0;
}
