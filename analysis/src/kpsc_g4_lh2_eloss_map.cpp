// -*- C++ -*-
// Build an LH2 energy-loss map for one particle species from Geant4 truth
// (primary particles of the given PDG codes; default: protons).
//
// Usage:
//   kpsc_g4_lh2_eloss_map <g4.root> [g4.root ...] [-o mapfile] [--pdf out.pdf]
//                         [--paxis prm|1st] [--pgrid coarse|fine] [--min-hits N]
//                         [--holdout-frac f] [--split-exit]
//                         [--pdg c1,c2,...] [--particle proton|pion|kaon]
//                         [--pedges e1,e2,...] [--vessel-radius R]
//
// Output: the Lh2ElossMan 5-column format
//   particle_id exit_id p[MeV/c] L[mm] Loss[MeV/c]
//   (particle_id 0=proton 1=pion 2=kaon, exit_id 0=side 1=cap)
// with the Target* and VesselRadius header keys. Without --split-exit the
// cap rows repeat the side grid. Several single-particle files are combined
// into one map with scripts/lh2eloss/merge_lh2_eloss_maps.py.
// --pdg: primary PDG codes used for the map (default 2212).
// --particle: species written as particle_id (default proton).
// --pedges: p grid edges [MeV/c] (overrides --pgrid).
// --ledges: L grid edges [mm] (default 0,10,...,80,100).
// --vessel-radius: VesselRadius header value [mm] (default 48.5).
// --dump-samples f.root: write every unique sample (train and holdout) to a
//             TTree "samples" for offline studies of the map.
// --min-hits: require at least N TPC hits on the primary track
//             (mimics the reconstructed sample; default 0).
// --split-exit: side grid (exit_id 0) from particles leaving the cylinder
//             through the side wall, cap grid (exit_id 1) from particles
//             leaving through a Y end cap.
// --holdout-frac f: fraction of unique events held out of the map median
//             (0..1, default 0). Holdout residuals are summarized after the
//             map is built. Duplicate kinematics never enter either set.
//
// Loss = p_PRM - p_TPC1st  [MeV/c] (>=0; negative samples dropped)
// L    = PathLength in finite LH2 cylinder from vtx along the track dir [mm]
// Grid: median Loss in (p, L) bins over unique kinematics; written as flat text map.
// --paxis selects the p used for the grid axis:
//   prm (default): p_PRM (vertex momentum)
//   1st          : p_TPC1st (momentum at the first TPC hit, close to what the
//                  helix measures; Lh2ElossMan::Loss is called with |p|_TPC)

#include "ana_helper.h"
#include "paths.h"

#include <TCanvas.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH2D.h>
#include <TInterpreter.h>
#include <TLatex.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParticle.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

namespace
{

constexpr Int_t kBeamGenerator = 7201;
constexpr Double_t kGeVtoMeV = 1.e3;
constexpr Int_t kRunOut = 0;
constexpr Double_t kTargetRadius = 40.0;       // mm (XZ)
constexpr Double_t kTargetHalfHeight = 50.0;  // mm (along Y)
constexpr Double_t kTargetCenterX = 0.0;      // mm
constexpr Double_t kTargetCenterY = 0.0;      // mm
constexpr Double_t kTargetCenterZ = -143.0;   // mm (= tpc::Z_TARGET)
constexpr Double_t kEps = 1.e-12;

const std::vector<Double_t> kPEdgesCoarse = {
  100., 150., 200., 250., 300., 350., 400., 450.,
  500., 550., 600., 650., 700., 750., 800., 900.
};
// Finer steps at low p, where Loss(p) is steep and strongly convex.
const std::vector<Double_t> kPEdgesFine = {
  100., 110., 120., 130., 140., 150., 165., 180., 200., 225., 250., 275.,
  300., 325., 350., 375., 400., 450., 500., 550., 600., 650., 700., 750.,
  800., 900.
};
const std::vector<Double_t> kLEdgesDefault = {
  0., 10., 20., 30., 40., 50., 60., 70., 80., 100.
};

void
usage(const char* a0)
{
  std::cerr
    << "Usage: " << a0
    << " <g4.root> [g4.root ...] [-o mapfile] [--pdf out.pdf]"
    << " [--paxis prm|1st] [--pgrid coarse|fine] [--min-hits N]"
    << " [--holdout-frac f] [--split-exit]"
    << " [--pdg c1,c2,...] [--particle proton|pion|kaon]"
    << " [--pedges e1,e2,...] [--ledges e1,e2,...] [--vessel-radius R]"
    << " [--dump-samples f.root]\n"
    << "  output: Lh2ElossMan 5-column map (particle_id exit_id p L Loss)\n"
    << "  tree: g4hyptpc\n"
    << "  Loss = p_PRM - p_TPC1st [MeV/c], L = cylinder PathLength [mm]\n"
    << "  Unique kinematics only (rounded fingerprint); optional holdout\n"
    << "  Default map/PDF under OUTPUT_DIR/img/run00000/kpsc_lh2eloss/g4/\n";
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
FirstHitMom(const std::vector<int>* pid, const std::vector<int>* parentid,
            const std::vector<int>* trackid, const std::vector<double>* px,
            const std::vector<double>* py, const std::vector<double>* pz,
            TVector3& mom, TVector3& pos,
            const std::vector<double>* x0, const std::vector<double>* y0,
            const std::vector<double>* z0, Int_t* nhit = nullptr,
            Int_t pdg = 2212)
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
    if ((*pid)[j] != pdg)
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
  if (nhit)
    *nhit = bn;
  const size_t j = first[best];
  mom.SetXYZ((*px)[j], (*py)[j], (*pz)[j]);
  if (x0 && y0 && z0 && j < x0->size())
    pos.SetXYZ((*x0)[j], (*y0)[j], (*z0)[j]);
  return mom.Mag() > 1e-9;
}

// Finite Y-axis cylinder PathLength (same convention as Lh2ElossMan / tpc::Z_TARGET).
Double_t
PathLengthCylinder(const TVector3& vtx, const TVector3& dir_in,
                   Double_t R, Double_t H,
                   Double_t cx, Double_t cy, Double_t cz,
                   Bool_t* exit_cap = nullptr)
{
  if (exit_cap)
    *exit_cap = kFALSE;
  TVector3 dir = dir_in;
  const Double_t n2 = dir.Mag2();
  if (!(n2 > kEps))
    return 0.;
  dir *= 1. / TMath::Sqrt(n2);

  const Double_t x0 = vtx.X() - cx;
  const Double_t y0 = vtx.Y() - cy;
  const Double_t z0 = vtx.Z() - cz;
  const Double_t rho2 = x0 * x0 + z0 * z0;
  if (!(rho2 <= R * R + 1.e-9 && std::abs(y0) <= H + 1.e-9))
    return 0.;

  const Double_t dx = dir.X();
  const Double_t dy = dir.Y();
  const Double_t dz = dir.Z();
  Double_t t_exit = 1.e300;
  Bool_t cap = kFALSE;

  const Double_t a = dx * dx + dz * dz;
  const Double_t b = 2. * (x0 * dx + z0 * dz);
  const Double_t c = rho2 - R * R;
  if (a > kEps) {
    const Double_t disc = b * b - 4. * a * c;
    if (disc >= 0.) {
      const Double_t s = TMath::Sqrt(disc);
      for (Double_t t : {(-b - s) / (2. * a), (-b + s) / (2. * a)}) {
        if (t <= kEps)
          continue;
        const Double_t y = y0 + t * dy;
        if (std::abs(y) <= H + 1.e-6 && t < t_exit) {
          t_exit = t;
          cap = kFALSE;
        }
      }
    }
  }
  if (std::abs(dy) > kEps) {
    for (Double_t ycap : {H, -H}) {
      const Double_t t = (ycap - y0) / dy;
      if (t <= kEps)
        continue;
      const Double_t x = x0 + t * dx;
      const Double_t z = z0 + t * dz;
      if (x * x + z * z <= R * R + 1.e-6 && t < t_exit) {
        t_exit = t;
        cap = kTRUE;
      }
    }
  }
  if (!(t_exit < 1.e299) || t_exit <= 0.)
    return 0.;
  if (exit_cap)
    *exit_cap = cap;
  return t_exit;
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

Int_t
BinIndex(const std::vector<Double_t>& edges, Double_t x)
{
  if (x < edges.front() || x > edges.back())
    return -1;
  // nearest edge index for map grid (store at edge centers? Plan stores at edge values)
  // Assign to nearest edge.
  Int_t best = 0;
  Double_t best_d = std::abs(x - edges[0]);
  for (Int_t i = 1; i < static_cast<Int_t>(edges.size()); ++i) {
    const Double_t d = std::abs(x - edges[static_cast<size_t>(i)]);
    if (d < best_d) {
      best_d = d;
      best = i;
    }
  }
  return best;
}

//_____________________________________________________________________________
Int_t
Quantize(
  Double_t x,
  Double_t scale
)
{
  return static_cast<Int_t>(std::llround(x * scale));
}

// Rounded kinematics key so duplicated G4 blocks (same reaction, many copies)
// contribute once to the median. Scale keeps mm / MeV/c at ~0.01 resolution.
struct EventFingerprint
{
  Int_t p0_q = 0;
  Int_t p1_q = 0;
  Int_t vx_q = 0;
  Int_t vy_q = 0;
  Int_t vz_q = 0;
  Int_t L_q = 0;
  Int_t loss_q = 0;
  Int_t cap = 0;

  Bool_t
  operator==(
    const EventFingerprint& o
  ) const
  {
    return p0_q == o.p0_q && p1_q == o.p1_q && vx_q == o.vx_q &&
           vy_q == o.vy_q && vz_q == o.vz_q && L_q == o.L_q &&
           loss_q == o.loss_q && cap == o.cap;
  }
};

struct EventFingerprintHash
{
  std::size_t
  operator()(
    const EventFingerprint& f
  ) const
  {
    std::size_t h = 1469598103934665603ull;
    const auto mix = [&](Int_t v) {
      h ^= static_cast<std::size_t>(static_cast<std::uint32_t>(v));
      h *= 1099511628211ull;
    };
    mix(f.p0_q);
    mix(f.p1_q);
    mix(f.vx_q);
    mix(f.vy_q);
    mix(f.vz_q);
    mix(f.L_q);
    mix(f.loss_q);
    mix(f.cap);
    return h;
  }
};

struct HoldoutSample
{
  Double_t p = 0.;    // same p axis as the map (--paxis)
  Double_t L = 0.;
  Double_t loss = 0.;
  Bool_t exit_cap = kFALSE;
};

//_____________________________________________________________________________
Int_t
FindLowerIndex(
  const std::vector<Double_t>& edges,
  Double_t x
)
{
  // Same clamp convention as Lh2ElossMan::FindLowerIndex.
  if (x <= edges.front())
    return 0;
  if (x >= edges.back())
    return static_cast<Int_t>(edges.size()) - 2;
  auto it = std::upper_bound(edges.begin(), edges.end(), x);
  Int_t i = static_cast<Int_t>(std::distance(edges.begin(), it)) - 1;
  if (i < 0)
    i = 0;
  if (i > static_cast<Int_t>(edges.size()) - 2)
    i = static_cast<Int_t>(edges.size()) - 2;
  return i;
}

//_____________________________________________________________________________
Double_t
BilinearLoss(
  const std::vector<Double_t>& p_edges,
  const std::vector<Double_t>& l_edges,
  const std::vector<Double_t>& grid,
  Double_t p_MeV,
  Double_t L_mm
)
{
  // Mirror of Lh2ElossMan::Loss bilinear clamp (without Cap/main selection).
  if (p_edges.size() < 2 || l_edges.size() < 2)
    return 0.;
  if (!(p_MeV > 0.) || !(L_mm >= 0.) || !std::isfinite(p_MeV) || !std::isfinite(L_mm))
    return 0.;
  const Int_t nl = static_cast<Int_t>(l_edges.size());
  const Int_t ip = FindLowerIndex(p_edges, p_MeV);
  const Int_t il = FindLowerIndex(l_edges, L_mm);
  const Double_t p0 = p_edges[static_cast<size_t>(ip)];
  const Double_t p1 = p_edges[static_cast<size_t>(ip + 1)];
  const Double_t l0 = l_edges[static_cast<size_t>(il)];
  const Double_t l1 = l_edges[static_cast<size_t>(il + 1)];
  auto at = [&](Int_t i, Int_t j) -> Double_t {
    return grid[static_cast<size_t>(i * nl + j)];
  };
  const Double_t f00 = at(ip, il);
  const Double_t f10 = at(ip + 1, il);
  const Double_t f01 = at(ip, il + 1);
  const Double_t f11 = at(ip + 1, il + 1);
  Double_t tp = 0.;
  if (std::abs(p1 - p0) > kEps) {
    tp = (p_MeV - p0) / (p1 - p0);
    if (tp < 0.)
      tp = 0.;
    if (tp > 1.)
      tp = 1.;
  }
  Double_t tl = 0.;
  if (std::abs(l1 - l0) > kEps) {
    tl = (L_mm - l0) / (l1 - l0);
    if (tl < 0.)
      tl = 0.;
    if (tl > 1.)
      tl = 1.;
  }
  const Double_t loss
    = (1. - tp) * (1. - tl) * f00
    + tp * (1. - tl) * f10
    + (1. - tp) * tl * f01
    + tp * tl * f11;
  if (!(loss > 0.) || !std::isfinite(loss))
    return 0.;
  return loss;
}

//_____________________________________________________________________________
EventFingerprint
MakeFingerprint(
  Double_t p0,
  Double_t p1,
  const TVector3& vtx,
  Double_t L,
  Double_t loss,
  Bool_t exit_cap
)
{
  EventFingerprint f;
  f.p0_q = Quantize(p0, 100.);
  f.p1_q = Quantize(p1, 100.);
  f.vx_q = Quantize(vtx.X(), 100.);
  f.vy_q = Quantize(vtx.Y(), 100.);
  f.vz_q = Quantize(vtx.Z(), 100.);
  f.L_q = Quantize(L, 100.);
  f.loss_q = Quantize(loss, 100.);
  f.cap = exit_cap ? 1 : 0;
  return f;
}

}  // namespace

int
main(int argc, char** argv)
{
  std::vector<TString> in_paths;
  TString out_map;
  TString out_pdf;
  Bool_t paxis_1st = kFALSE;
  Bool_t pgrid_fine = kFALSE;
  Int_t min_hits = 0;
  Bool_t split_exit = kFALSE;
  Double_t holdout_frac = 0.;
  std::vector<Int_t> pdg_list = {2212};
  TString particle_name = "proton";
  std::vector<Double_t> p_edges_user;
  std::vector<Double_t> l_edges_user;
  Double_t vessel_radius = 48.5;  // mm, vacuum-vessel window (G4 TargetAl)
  TString dump_path;
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
      out_map = argv[++i];
      continue;
    }
    if (a == "--pdf") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      out_pdf = argv[++i];
      continue;
    }
    if (a == "--paxis") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      const TString v = argv[++i];
      if (v == "1st")
        paxis_1st = kTRUE;
      else if (v != "prm") {
        std::cerr << "Unknown --paxis: " << v << "\n";
        return 1;
      }
      continue;
    }
    if (a == "--split-exit") {
      split_exit = kTRUE;
      continue;
    }
    if (a == "--min-hits") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      min_hits = TString(argv[++i]).Atoi();
      continue;
    }
    if (a == "--holdout-frac") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      holdout_frac = TString(argv[++i]).Atof();
      if (!(holdout_frac >= 0.) || holdout_frac >= 1.) {
        std::cerr << "Error: --holdout-frac must be in [0, 1)\n";
        return 1;
      }
      continue;
    }
    if (a == "--pdg" && i + 1 < argc) {
      // comma separated PDG codes of the primary particles to use
      pdg_list.clear();
      TString v = argv[++i];
      std::unique_ptr<TObjArray> tok(v.Tokenize(","));
      for (Int_t k = 0; k < tok->GetEntries(); ++k)
        pdg_list.push_back(static_cast<TObjString*>(tok->At(k))->GetString().Atoi());
      continue;
    }
    if (a == "--particle" && i + 1 < argc) {
      particle_name = argv[++i];
      continue;
    }
    if (a == "--ledges" && i + 1 < argc) {
      // comma separated L grid edges [mm]
      TString v = argv[++i];
      std::unique_ptr<TObjArray> tok(v.Tokenize(","));
      for (Int_t k = 0; k < tok->GetEntries(); ++k)
        l_edges_user.push_back(static_cast<TObjString*>(tok->At(k))->GetString().Atof());
      continue;
    }
    if (a == "--dump-samples" && i + 1 < argc) {
      dump_path = argv[++i];
      continue;
    }
    if (a == "--vessel-radius" && i + 1 < argc) {
      vessel_radius = TString(argv[++i]).Atof();
      continue;
    }
    if (a == "--pedges" && i + 1 < argc) {
      // comma separated p grid edges [MeV/c]
      TString v = argv[++i];
      std::unique_ptr<TObjArray> tok(v.Tokenize(","));
      for (Int_t k = 0; k < tok->GetEntries(); ++k)
        p_edges_user.push_back(static_cast<TObjString*>(tok->At(k))->GetString().Atof());
      continue;
    }
    if (a == "--pgrid") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        return 1;
      }
      const TString v = argv[++i];
      if (v == "fine")
        pgrid_fine = kTRUE;
      else if (v != "coarse") {
        std::cerr << "Unknown --pgrid: " << v << "\n";
        return 1;
      }
      continue;
    }
    if (a.BeginsWith("-")) {
      std::cerr << "Unknown option: " << a << "\n";
      usage(argv[0]);
      return 1;
    }
    in_paths.push_back(a);
  }
  if (in_paths.empty()) {
    usage(argv[0]);
    return 1;
  }
  // particle_id of the Lh2ElossMan 5-column map
  Int_t particle_id = -1;
  if (particle_name == "proton")
    particle_id = 0;
  else if (particle_name == "pion")
    particle_id = 1;
  else if (particle_name == "kaon")
    particle_id = 2;
  if (particle_id < 0) {
    std::cerr << "Unknown --particle: " << particle_name
              << " (proton|pion|kaon)\n";
    return 1;
  }
  if (!(vessel_radius > kTargetRadius)) {
    std::cerr << "--vessel-radius must be larger than TargetRadius ("
              << kTargetRadius << " mm)\n";
    return 1;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gInterpreter->GenerateDictionary("vector<TParticle>", "vector;TParticle.h");

  const std::vector<Double_t> kPEdges = !p_edges_user.empty() ? p_edges_user
    : (pgrid_fine ? kPEdgesFine : kPEdgesCoarse);
  const Int_t np = static_cast<Int_t>(kPEdges.size());
  const std::vector<Double_t> kLEdges = !l_edges_user.empty() ? l_edges_user
                                                              : kLEdgesDefault;
  const Int_t nl = static_cast<Int_t>(kLEdges.size());
  std::vector<std::vector<Double_t>> cells(
    static_cast<size_t>(np * nl));
  std::vector<std::vector<Double_t>> cells_cap(
    static_cast<size_t>(np * nl));
  Long64_t n_used_cap = 0;

  Long64_t n_raw = 0;
  Long64_t n_used = 0;
  Long64_t n_dup_skipped = 0;
  Long64_t n_holdout = 0;
  Long64_t n_skip_neg = 0;
  Long64_t n_skip_L0 = 0;
  Long64_t n_skip_bin = 0;
  std::unordered_set<EventFingerprint, EventFingerprintHash> seen;
  std::vector<HoldoutSample> holdout;
  const EventFingerprintHash fp_hash;

  // optional per-sample dump (MeV/c, mm)
  std::unique_ptr<TFile> dump_file;
  TTree* dump_tree = nullptr;
  Float_t d_p0 = 0, d_p1 = 0, d_L = 0, d_loss = 0;
  Float_t d_vx = 0, d_vy = 0, d_vz = 0, d_ux = 0, d_uy = 0, d_uz = 0;
  Int_t d_pdg = 0, d_cap = 0, d_holdout = 0, d_file = 0, d_nhit = 0;
  if (!dump_path.IsNull()) {
    dump_file.reset(TFile::Open(dump_path, "RECREATE"));
    if (!dump_file || dump_file->IsZombie()) {
      std::cerr << "Error: cannot write " << dump_path << "\n";
      return 1;
    }
    dump_tree = new TTree("samples", "unique LH2 eloss samples");
    dump_tree->Branch("p0", &d_p0);
    dump_tree->Branch("p1", &d_p1);
    dump_tree->Branch("L", &d_L);
    dump_tree->Branch("loss", &d_loss);
    dump_tree->Branch("vx", &d_vx);
    dump_tree->Branch("vy", &d_vy);
    dump_tree->Branch("vz", &d_vz);
    dump_tree->Branch("ux", &d_ux);
    dump_tree->Branch("uy", &d_uy);
    dump_tree->Branch("uz", &d_uz);
    dump_tree->Branch("pdg", &d_pdg);
    dump_tree->Branch("cap", &d_cap);
    dump_tree->Branch("holdout", &d_holdout);
    dump_tree->Branch("file", &d_file);
    dump_tree->Branch("nhit", &d_nhit);
  }
  Int_t file_index = -1;

  for (const auto& in_path : in_paths) {
    ++file_index;
    TFile* fin = TFile::Open(in_path, "READ");
    if (!fin || fin->IsZombie()) {
      std::cerr << "Error: cannot open " << in_path << "\n";
      return 1;
    }
    auto* tree = dynamic_cast<TTree*>(fin->Get("g4hyptpc"));
    if (!tree) {
      std::cerr << "Error: tree g4hyptpc not found in " << in_path << "\n";
      return 1;
    }

    Int_t generator = 0;
    std::vector<TParticle>* PRM = nullptr;
    std::vector<int>* pidtpc = nullptr;
    std::vector<int>* parentidtpc = nullptr;
    std::vector<int>* trackidtpc = nullptr;
    std::vector<double>* pxtpc = nullptr;
    std::vector<double>* pytpc = nullptr;
    std::vector<double>* pztpc = nullptr;
    std::vector<double>* x0tpc = nullptr;
    std::vector<double>* y0tpc = nullptr;
    std::vector<double>* z0tpc = nullptr;

    // read only the branches used below
    tree->SetBranchStatus("*", 0);
    for (const char* bn : {"generator", "PRM", "pidtpc", "parentidtpc", "trackidtpc",
                           "pxtpc", "pytpc", "pztpc", "x0tpc", "y0tpc", "z0tpc"})
      tree->SetBranchStatus(bn, 1);
    tree->SetBranchAddress("generator", &generator);
    tree->SetBranchAddress("PRM", &PRM);
    tree->SetBranchAddress("pidtpc", &pidtpc);
    tree->SetBranchAddress("parentidtpc", &parentidtpc);
    tree->SetBranchAddress("trackidtpc", &trackidtpc);
    tree->SetBranchAddress("pxtpc", &pxtpc);
    tree->SetBranchAddress("pytpc", &pytpc);
    tree->SetBranchAddress("pztpc", &pztpc);
    tree->SetBranchAddress("x0tpc", &x0tpc);
    tree->SetBranchAddress("y0tpc", &y0tpc);
    tree->SetBranchAddress("z0tpc", &z0tpc);

    const Long64_t nent = tree->GetEntries();
    for (Long64_t ie = 0; ie < nent; ++ie) {
      tree->GetEntry(ie);
      if (generator == kBeamGenerator)
        continue;
      if (!PRM || PRM->empty())
        continue;

      for (const auto& part : *PRM) {
      const Int_t pdg = part.GetPdgCode();
      if (std::find(pdg_list.begin(), pdg_list.end(), pdg) == pdg_list.end())
        continue;
      TVector3 p_prm(part.Px(), part.Py(), part.Pz());
      TVector3 v_prm(part.Vx(), part.Vy(), part.Vz());
      const Bool_t has_p = kTRUE;
      if (!has_p || !(p_prm.Mag() > 1e-3))
        continue;

      TVector3 p1, r1;
      Int_t nhit_p = 0;
      if (!FirstHitMom(pidtpc, parentidtpc, trackidtpc, pxtpc, pytpc, pztpc,
                       p1, r1, x0tpc, y0tpc, z0tpc, &nhit_p, pdg))
        continue;
      if (nhit_p < min_hits)
        continue;
      p1 *= kGeVtoMeV;
      const Double_t p0 = p_prm.Mag();
      const Double_t loss = p0 - p1.Mag();  // MeV/c, want >=0
      if (!(loss >= 0.) || !std::isfinite(loss)) {
        ++n_skip_neg;
        continue;
      }
      Bool_t exit_cap = kFALSE;
      const Double_t L = PathLengthCylinder(v_prm, p_prm.Unit(),
                                           kTargetRadius, kTargetHalfHeight,
                                           kTargetCenterX, kTargetCenterY,
                                           kTargetCenterZ, &exit_cap);
      if (!(L > 0.)) {
        ++n_skip_L0;
        continue;
      }
      const Int_t ip = BinIndex(kPEdges, paxis_1st ? p1.Mag() : p0);
      const Int_t il = BinIndex(kLEdges, L);
      if (ip < 0 || il < 0) {
        ++n_skip_bin;
        continue;
      }
      ++n_raw;
      const EventFingerprint fp = MakeFingerprint(
        p0, p1.Mag(), v_prm, L, loss, exit_cap);
      if (!seen.insert(fp).second) {
        ++n_dup_skipped;
        continue;
      }
      const Bool_t is_holdout =
        (holdout_frac > 0.) &&
        ((fp_hash(fp) % 10000ull) <
         static_cast<std::size_t>(holdout_frac * 10000. + 0.5));
      if (dump_tree) {
        const TVector3 u = p_prm.Unit();
        d_p0 = p0;
        d_p1 = p1.Mag();
        d_L = L;
        d_loss = loss;
        d_vx = v_prm.X();
        d_vy = v_prm.Y();
        d_vz = v_prm.Z();
        d_ux = u.X();
        d_uy = u.Y();
        d_uz = u.Z();
        d_pdg = pdg;
        d_cap = exit_cap ? 1 : 0;
        d_holdout = is_holdout ? 1 : 0;
        d_file = file_index;
        d_nhit = nhit_p;
        dump_tree->Fill();
      }
      if (is_holdout) {
        HoldoutSample s;
        s.p = paxis_1st ? p1.Mag() : p0;
        s.L = L;
        s.loss = loss;
        s.exit_cap = exit_cap;
        holdout.push_back(s);
        ++n_holdout;
        continue;
      }
      if (split_exit && exit_cap) {
        cells_cap[static_cast<size_t>(ip * nl + il)].push_back(loss);
        ++n_used_cap;
      } else {
        cells[static_cast<size_t>(ip * nl + il)].push_back(loss);
      }
      ++n_used;
      }
    }
    fin->Close();
  }

  if (dump_tree) {
    dump_file->cd();
    dump_tree->Write();
    dump_file->Close();
    std::cout << "Wrote samples: " << dump_path << "\n";
  }

  // Fill empty cells by nearest non-empty along L at same p, else 0.
  std::vector<Double_t> loss_grid(static_cast<size_t>(np * nl), 0.);
  std::vector<Long64_t> n_grid(static_cast<size_t>(np * nl), 0);
  for (Int_t ip = 0; ip < np; ++ip) {
    for (Int_t il = 0; il < nl; ++il) {
      const size_t idx = static_cast<size_t>(ip * nl + il);
      n_grid[idx] = static_cast<Long64_t>(cells[idx].size());
      if (!cells[idx].empty())
        loss_grid[idx] = MedianSorted(cells[idx]);
    }
  }
  for (Int_t ip = 0; ip < np; ++ip) {
    for (Int_t il = 0; il < nl; ++il) {
      const size_t idx = static_cast<size_t>(ip * nl + il);
      if (n_grid[idx] > 0)
        continue;
      // prefer same-p nearest L with stats
      Int_t best_il = -1;
      Int_t best_d = 999;
      for (Int_t jl = 0; jl < nl; ++jl) {
        const size_t jdx = static_cast<size_t>(ip * nl + jl);
        if (n_grid[jdx] == 0)
          continue;
        const Int_t d = std::abs(jl - il);
        if (d < best_d) {
          best_d = d;
          best_il = jl;
        }
      }
      if (best_il >= 0)
        loss_grid[idx] = loss_grid[static_cast<size_t>(ip * nl + best_il)];
      else if (il == 0)
        loss_grid[idx] = 0.;
    }
  }

  // Cap grid: median per cell; empty cells take the nearest filled L at the
  // same p, and rows without any cap sample fall back to the main grid.
  std::vector<Double_t> loss_cap(static_cast<size_t>(np * nl), 0.);
  if (split_exit) {
    for (Int_t ip = 0; ip < np; ++ip) {
      Bool_t any = kFALSE;
      for (Int_t il = 0; il < nl; ++il)
        any = any || !cells_cap[static_cast<size_t>(ip * nl + il)].empty();
      for (Int_t il = 0; il < nl; ++il) {
        const size_t idx = static_cast<size_t>(ip * nl + il);
        if (!any) {
          loss_cap[idx] = loss_grid[idx];
          continue;
        }
        if (!cells_cap[idx].empty()) {
          loss_cap[idx] = MedianSorted(cells_cap[idx]);
          continue;
        }
        Int_t best_il = -1, best_d = 999;
        for (Int_t jl = 0; jl < nl; ++jl) {
          if (cells_cap[static_cast<size_t>(ip * nl + jl)].empty())
            continue;
          if (std::abs(jl - il) < best_d) {
            best_d = std::abs(jl - il);
            best_il = jl;
          }
        }
        loss_cap[idx] = MedianSorted(cells_cap[static_cast<size_t>(ip * nl + best_il)]);
      }
    }
  }

  // default outputs go to the kpsc_lh2eloss/g4/ sub-directory of img/runNNNNN
  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, kRunOut) + "/kpsc_lh2eloss/g4";
  gSystem->mkdir(img_dir.Data(), kTRUE);
  if (out_map.IsNull())
    out_map = Form("%s/Lh2ElossMap_%s", img_dir.Data(), particle_name.Data());
  if (out_pdf.IsNull())
    out_pdf = Form("%s/kpsc_g4_lh2_eloss_map_run%05d.pdf", img_dir.Data(),
                   kRunOut);

  {
    std::ofstream ofs(out_map.Data());
    if (!ofs) {
      std::cerr << "Error: cannot write " << out_map << "\n";
      return 1;
    }
    ofs << "#\n";
    ofs << "#  Lh2ElossMap_" << particle_name << "\n";
    ofs << "#  primary PDG codes:";
    for (const auto c : pdg_list) ofs << " " << c;
    ofs << "\n";
    ofs << "#  Loss = p_start - p_end  [MeV/c]  (>=0)\n";
    ofs << "#  apply: |p|_vtx [MeV/c] = |p|_TPC [MeV/c] + Loss\n";
    ofs << "#  L from vtx+dir vs cylinder (TargetRadius, TargetHalfHeight)\n";
    ofs << "#  built by kpsc_g4_lh2_eloss_map  N_train=" << n_used
        << "  N_raw=" << n_raw
        << "  N_unique=" << static_cast<Long64_t>(seen.size())
        << "  N_dup_skipped=" << n_dup_skipped
        << "  N_holdout=" << n_holdout << "\n";
    ofs << "#  dedupe: rounded kinematics fingerprint (cross-file)\n";
    ofs << "#  holdout_frac=" << holdout_frac << "\n";
    ofs << "#  min TPC hits on primary track: " << min_hits << "\n";
    if (split_exit)
      ofs << "#  split by exit surface: exit_id 0 = side exits,"
          << " exit_id 1 = end-cap exits (N_cap=" << n_used_cap << ")\n";
    else
      ofs << "#  no exit split: exit_id 1 rows repeat the exit_id 0 grid\n";
    ofs << "#  p axis: "
        << (paxis_1st ? "p_TPC1st (momentum at first TPC hit)"
                      : "p_PRM (vertex momentum)") << "\n";
    for (const auto& p : in_paths)
      ofs << "#  input: " << p << "\n";
    ofs << "#\n";
    ofs << "#______________________________________________________________________________\n";
    ofs << "TargetRadius       " << kTargetRadius << "    # mm\n";
    ofs << "TargetHalfHeight   " << kTargetHalfHeight << "    # mm\n";
    ofs << "TargetCenterX      " << kTargetCenterX << "    # mm\n";
    ofs << "TargetCenterY      " << kTargetCenterY << "    # mm\n";
    ofs << "TargetCenterZ      " << kTargetCenterZ << "    # mm\n";
    ofs << "VesselRadius       " << vessel_radius << "    # mm\n";
    ofs << "\n";
    ofs << "#______________________________________________________________________________\n";
    ofs << "# particle  exit  p[MeV/c]  L[mm]  Loss[MeV/c]"
        << "  (particle 0=proton 1=pion 2=kaon, exit 0=side 1=cap)\n";
    ofs << std::fixed;
    for (Int_t exit_id = 0; exit_id < 2; ++exit_id) {
      const std::vector<Double_t>& grid
        = (exit_id == 1 && split_exit) ? loss_cap : loss_grid;
      for (Int_t ip = 0; ip < np; ++ip) {
        for (Int_t il = 0; il < nl; ++il) {
          const size_t idx = static_cast<size_t>(ip * nl + il);
          ofs << particle_id << "  " << exit_id << "  "
              << kPEdges[static_cast<size_t>(ip)] << "  "
              << kLEdges[static_cast<size_t>(il)] << "  "
              << grid[idx] << "\n";
        }
      }
    }
  }

  const TString p_title = paxis_1st ? "p_{TPC1st} [MeV/c]" : "p_{PRM} [MeV/c]";
  TH2D h_loss("h_loss", "median Loss;" + p_title + ";L_{LH2} [mm]",
              np - 1, kPEdges.front(), kPEdges.back(),
              nl - 1, kLEdges.front(), kLEdges.back());
  TH2D h_n("h_n", "entries;" + p_title + ";L_{LH2} [mm]",
           np - 1, kPEdges.front(), kPEdges.back(),
           nl - 1, kLEdges.front(), kLEdges.back());
  for (Int_t ip = 0; ip < np; ++ip) {
    for (Int_t il = 0; il < nl; ++il) {
      const size_t idx = static_cast<size_t>(ip * nl + il);
      h_loss.Fill(kPEdges[static_cast<size_t>(ip)],
                  kLEdges[static_cast<size_t>(il)], loss_grid[idx]);
      h_n.Fill(kPEdges[static_cast<size_t>(ip)],
               kLEdges[static_cast<size_t>(il)],
               static_cast<Double_t>(n_grid[idx]));
    }
  }

  TCanvas c("c", "c", 900, 700);
  PdfWriter pdf(out_pdf);
  c.Clear();
  c.cd();
  gPad->SetRightMargin(0.14);
  h_loss.SetTitle(Form("LH2 %s Loss median [MeV/c]", particle_name.Data()));
  h_loss.Draw("colz");
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.03);
  lat.DrawLatex(
    0.12, 0.92,
    Form("N_train=%lld  N_raw=%lld  dup=%lld  holdout=%lld",
         n_used, n_raw, n_dup_skipped, n_holdout));
  pdf.Print(c);
  c.Clear();
  c.cd();
  gPad->SetRightMargin(0.14);
  h_n.SetTitle("LH2 map unique train entries per (p,L) cell");
  h_n.Draw("colz");
  pdf.Print(c);
  pdf.Close(c);

  Double_t holdout_mean_abs = 0.;
  Double_t holdout_rms = 0.;
  Long64_t n_holdout_eval = 0;
  if (!holdout.empty()) {
    Double_t sum = 0.;
    Double_t sum2 = 0.;
    for (const auto& s : holdout) {
      // Same bilinear + Cap selection as Lh2ElossMan::Loss for this map.
      const std::vector<Double_t>& grid =
        (split_exit && s.exit_cap) ? loss_cap : loss_grid;
      const Double_t pred = BilinearLoss(kPEdges, kLEdges, grid, s.p, s.L);
      const Double_t resid = s.loss - pred;
      sum += std::abs(resid);
      sum2 += resid * resid;
      ++n_holdout_eval;
    }
    if (n_holdout_eval > 0) {
      holdout_mean_abs = sum / static_cast<Double_t>(n_holdout_eval);
      holdout_rms = std::sqrt(sum2 / static_cast<Double_t>(n_holdout_eval));
    }
  }

  std::cout << "Wrote map: " << out_map << "\n"
            << "Wrote PDF: " << out_pdf << "\n"
            << "N_raw=" << n_raw
            << " N_unique=" << seen.size()
            << " N_train=" << n_used
            << " N_dup_skipped=" << n_dup_skipped
            << " N_holdout=" << n_holdout
            << " skip_neg=" << n_skip_neg
            << " L0=" << n_skip_L0
            << " bin=" << n_skip_bin << "\n";
  if (n_holdout_eval > 0) {
    std::cout << "holdout residual vs Lh2ElossMan-style bilinear Loss: "
              << "N=" << n_holdout_eval
              << " mean|d|=" << holdout_mean_abs
              << " rms=" << holdout_rms << " MeV/c\n";
  }
  return 0;
}
