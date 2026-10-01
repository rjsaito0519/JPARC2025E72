// -*- C++ -*-
// Build Lh2ElossMap_proton from Geant4 truth (scattered protons).
//
// Usage:
//   kpsc_g4_lh2_eloss_map <g4.root> [g4.root ...] [-o mapfile] [--pdf out.pdf]
//                         [--paxis prm|1st] [--pgrid coarse|fine] [--min-hits N]
// --min-hits: require at least N TPC hits on the primary proton track
//             (mimics the reconstructed sample; default 0).
// --split-exit: main grid from protons leaving the cylinder through the side
//             wall, plus "Cap p L Loss" lines from protons leaving through a
//             Y end cap (read by Lh2ElossMan as the cap grid).
//
// Loss = p_PRM - p_TPC1st  [MeV/c] (>=0; negative samples dropped)
// L    = PathLength in finite LH2 cylinder from vtx along proton dir [mm]
// Grid: median Loss in (p, L) bins; written as flat text map.
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
#include <TParticle.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>

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
const std::vector<Double_t> kLEdges = {
  0., 10., 20., 30., 40., 50., 60., 70., 80., 100.
};

void
usage(const char* a0)
{
  std::cerr
    << "Usage: " << a0
    << " <g4.root> [g4.root ...] [-o mapfile] [--pdf out.pdf]"
    << " [--paxis prm|1st] [--pgrid coarse|fine] [--min-hits N] [--split-exit]\n"
    << "  tree: g4hyptpc\n"
    << "  Loss = p_PRM - p_TPC1st [MeV/c], L = cylinder PathLength [mm]\n"
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
            const std::vector<double>* z0, Int_t* nhit = nullptr)
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
    if ((*pid)[j] != 2212)
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

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gInterpreter->GenerateDictionary("vector<TParticle>", "vector;TParticle.h");

  const std::vector<Double_t>& kPEdges = pgrid_fine ? kPEdgesFine : kPEdgesCoarse;
  const Int_t np = static_cast<Int_t>(kPEdges.size());
  const Int_t nl = static_cast<Int_t>(kLEdges.size());
  std::vector<std::vector<Double_t>> cells(
    static_cast<size_t>(np * nl));
  std::vector<std::vector<Double_t>> cells_cap(
    static_cast<size_t>(np * nl));
  Long64_t n_used_cap = 0;

  Long64_t n_used = 0;
  Long64_t n_skip_neg = 0;
  Long64_t n_skip_L0 = 0;
  Long64_t n_skip_bin = 0;

  for (const auto& in_path : in_paths) {
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

      TVector3 p_prm(0, 0, 0), v_prm(0, 0, 0);
      Bool_t has_p = kFALSE;
      for (const auto& part : *PRM) {
        if (part.GetPdgCode() == 2212) {
          p_prm.SetXYZ(part.Px(), part.Py(), part.Pz());
          v_prm.SetXYZ(part.Vx(), part.Vy(), part.Vz());
          has_p = kTRUE;
          break;
        }
      }
      if (!has_p || !(p_prm.Mag() > 1e-3))
        continue;

      TVector3 p1, r1;
      Int_t nhit_p = 0;
      if (!FirstHitMom(pidtpc, parentidtpc, trackidtpc, pxtpc, pytpc, pztpc,
                       p1, r1, x0tpc, y0tpc, z0tpc, &nhit_p))
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
      if (split_exit && exit_cap) {
        cells_cap[static_cast<size_t>(ip * nl + il)].push_back(loss);
        ++n_used_cap;
      } else {
        cells[static_cast<size_t>(ip * nl + il)].push_back(loss);
      }
      ++n_used;
    }
    fin->Close();
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
    out_map = Form("%s/Lh2ElossMap_proton", img_dir.Data());
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
    ofs << "#  Lh2ElossMap_proton\n";
    ofs << "#  Loss = p_start - p_end  [MeV/c]  (>=0)\n";
    ofs << "#  apply: |p|_vtx [MeV/c] = |p|_TPC [MeV/c] + Loss\n";
    ofs << "#  L from vtx+dir vs cylinder (TargetRadius, TargetHalfHeight)\n";
    ofs << "#  built by kpsc_g4_lh2_eloss_map  N_used=" << n_used << "\n";
    ofs << "#  min TPC hits on proton track: " << min_hits << "\n";
    if (split_exit)
      ofs << "#  split by exit surface: main grid = side exits,"
          << " Cap lines = end-cap exits (N_cap=" << n_used_cap << ")\n";
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
    ofs << "\n";
    ofs << "#______________________________________________________________________________\n";
    ofs << "# p[MeV/c]   L[mm]   Loss[MeV/c]\n";
    ofs << std::fixed;
    for (Int_t ip = 0; ip < np; ++ip) {
      for (Int_t il = 0; il < nl; ++il) {
        const size_t idx = static_cast<size_t>(ip * nl + il);
        ofs << kPEdges[static_cast<size_t>(ip)] << "  "
            << kLEdges[static_cast<size_t>(il)] << "  "
            << loss_grid[idx] << "\n";
      }
    }
    if (split_exit) {
      ofs << "\n";
      ofs << "#______________________________________________________________________________\n";
      ofs << "# end-cap exits:  Cap  p[MeV/c]   L[mm]   Loss[MeV/c]\n";
      for (Int_t ip = 0; ip < np; ++ip) {
        for (Int_t il = 0; il < nl; ++il) {
          const size_t idx = static_cast<size_t>(ip * nl + il);
          ofs << "Cap  " << kPEdges[static_cast<size_t>(ip)] << "  "
              << kLEdges[static_cast<size_t>(il)] << "  "
              << loss_cap[idx] << "\n";
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
  h_loss.SetTitle("LH2 proton Loss median [MeV/c]");
  h_loss.Draw("colz");
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.03);
  lat.DrawLatex(0.12, 0.92, Form("N_used=%lld  skip_neg=%lld  L0=%lld  bin=%lld",
                                 n_used, n_skip_neg, n_skip_L0, n_skip_bin));
  pdf.Print(c);
  c.Clear();
  c.cd();
  gPad->SetRightMargin(0.14);
  h_n.SetTitle("LH2 map entries per (p,L) cell");
  h_n.Draw("colz");
  pdf.Print(c);
  pdf.Close(c);

  std::cout << "Wrote map: " << out_map << "\n"
            << "Wrote PDF: " << out_pdf << "\n"
            << "N_used=" << n_used
            << " skip_neg=" << n_skip_neg
            << " L0=" << n_skip_L0
            << " bin=" << n_skip_bin << "\n";
  return 0;
}
