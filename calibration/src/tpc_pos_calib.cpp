// TPCPOS pad-level fit from TPCHelix inclusive residuals.
// Spec: docs/tpc/TPCPOS_ITER_PREP.ja.md
// formula_id: helix_inclusive_tangent_phi_v1
//
// Usage:
//   tpc_pos_calib <helix.root> [--model tan|xz] [--threshold N] [--ave] [--debug]
//                 [--parent <map>] [--overwrite]
//
// The residuals are measured with the TPCPOS map that was used to reconstruct
// <helix.root>, so the update is incremental: new = parent + (-alpha * <r>).
// --parent must name that map (default: the shared e72 map).
// --overwrite replaces the touched pads by -alpha * <r> instead (valid only
// when the helix was reconstructed with a zero map).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <TArrow.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TPRegexp.h>
#include <TString.h>
#include <TVector3.h>

#include "TPCPadHelper.hh"
#include "ana_helper.h"
#include "paths.h"

namespace {

constexpr Double_t kAlpha = 1.0;
constexpr Int_t kMinHitsTrack = 15;
constexpr Double_t kMaxAbsDz = 0.05;
constexpr const char* kFormulaId = "helix_inclusive_tangent_phi_v1";

enum class Model {
  kTan,
  kXz
};

struct PadAgg {
  Long64_t n_hit = 0;
  Double_t sum_rx = 0.0;
  Double_t sum_rz = 0.0;
  Double_t sum_rt = 0.0;
};

struct PadCorr {
  Bool_t fitted = kFALSE;
  Bool_t from_ave = kFALSE;
  Bool_t center_frame = kFALSE;
  Double_t vx = 0.0;
  Double_t vz = 0.0;
  Double_t ct = 0.0;
  Long64_t n_hit = 0;
};

using PadKey = std::pair<Int_t, Int_t>;

Bool_t
ParseModel(
  const TString& s,
  Model& out
)
{
  const TString v = s.Strip(TString::kBoth);
  if (v == "tan" || v == "tangent" || v == "rt") {
    out = Model::kTan;
    return kTRUE;
  }
  if (v == "xz" || v == "cartesian" || v == "cart") {
    out = Model::kXz;
    return kTRUE;
  }
  return kFALSE;
}

const char*
ModelName(
  Model m
)
{
  return (m == Model::kTan) ? "tan" : "xz";
}

void
PadUnitTangent(
  Int_t layer,
  Int_t row,
  Double_t& e_tx,
  Double_t& e_tz
)
{
  const Double_t th = tpc::getTheta(layer, static_cast<Double_t>(row)) * TMath::DegToRad();
  e_tx = std::cos(th);
  e_tz = -std::sin(th);
}

Int_t
DetectRunNumber(
  TFile* f,
  const TString& input_root
)
{
  Int_t run_number = -1;
  TTree* tree = dynamic_cast<TTree*>(f->Get("tpc"));
  if (tree) {
    try {
      TTreeReader reader(tree);
      TTreeReaderValue<UInt_t> rv(reader, "run_number");
      if (reader.Next()) {
        run_number = static_cast<Int_t>(*rv);
      }
    } catch (...) {
    }
  }
  if (run_number < 0) {
    TPRegexp run_regex("run([0-9]+)");
    TString run_str = input_root(run_regex);
    if (!run_str.IsNull()) {
      run_str.ReplaceAll("run", "");
      run_number = run_str.Atoi();
    }
  }
  return run_number;
}

Bool_t
ReadParentMap(
  const TString& path,
  Int_t& n_y,
  Double_t& y0,
  Double_t& dy,
  std::map<PadKey, std::vector<std::pair<Double_t, std::pair<Double_t, Double_t>>>>& grid
)
{
  std::ifstream ifs(path.Data());
  if (!ifs.is_open()) {
    return kFALSE;
  }
  std::string line;
  if (!std::getline(ifs, line)) {
    return kFALSE;
  }
  {
    std::istringstream iss(line);
    if (!(iss >> n_y >> y0 >> dy) || n_y <= 0) {
      return kFALSE;
    }
  }
  grid.clear();
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    Int_t layer = 0;
    Int_t row = 0;
    Double_t y = 0.0;
    Double_t vx = 0.0;
    Double_t vz = 0.0;
    if (!(iss >> layer >> row >> y >> vx >> vz)) {
      continue;
    }
    grid[{layer, row}].push_back({y, {vx, vz}});
  }
  return kTRUE;
}

Bool_t
WriteMap(
  const TString& path,
  Int_t n_y,
  Double_t y0,
  Double_t dy,
  const std::map<PadKey, std::vector<std::pair<Double_t, std::pair<Double_t, Double_t>>>>& grid
)
{
  std::ofstream ofs(path.Data());
  if (!ofs.is_open()) {
    return kFALSE;
  }
  ofs << n_y << "  " << y0 << "  " << dy << "\n";
  ofs << std::setprecision(10);
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < n_pad; ++row) {
      const auto it = grid.find({layer, row});
      if (it == grid.end()) {
        continue;
      }
      for (const auto& node : it->second) {
        ofs << layer << "   " << row << "   " << node.first << "   "
            << node.second.first << "   " << node.second.second << "\n";
      }
    }
  }
  return kTRUE;
}

void
PrintProgress(
  Long64_t i,
  Long64_t n
)
{
  if (n <= 0) {
    return;
  }
  const Long64_t step = std::max<Long64_t>(1, n / 200);
  if ((i % step) != 0 && i + 1 != n) {
    return;
  }
  const Double_t pct = 100.0 * static_cast<Double_t>(i + 1) / static_cast<Double_t>(n);
  std::fprintf(stderr, "\r[pos] events %lld / %lld (%.1f%%)", static_cast<long long>(i + 1),
               static_cast<long long>(n), pct);
  std::fflush(stderr);
}

}  // namespace

//_____________________________________________________________________________
int
main(
  Int_t argc,
  char** argv
)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <helix.root> [--model tan|xz] [--threshold N] [--ave] [--debug]"
              << " [--parent <map>] [--overwrite]"
              << std::endl;
    return EXIT_FAILURE;
  }

  const TString input_root = argv[1];
  Model model = Model::kTan;
  Int_t threshold = 100;
  Bool_t use_ave = kFALSE;
  Bool_t debug_mode = kFALSE;
  TString parent_arg;
  Bool_t overwrite = kFALSE;

  for (Int_t i = 2; i < argc; ++i) {
    const TString arg = argv[i];
    if (arg == "--model" && i + 1 < argc) {
      Model parsed;
      if (!ParseModel(argv[++i], parsed)) {
        std::cerr << "Error: unknown --model (use tan|xz)" << std::endl;
        return EXIT_FAILURE;
      }
      model = parsed;
    } else if (arg == "--threshold" && i + 1 < argc) {
      threshold = std::atoi(argv[++i]);
    } else if (arg == "--ave" || arg == "--asad-avg") {
      use_ave = kTRUE;
    } else if (arg == "--debug") {
      debug_mode = kTRUE;
    } else if (arg == "--overwrite") {
      overwrite = kTRUE;
    } else if (arg == "--parent" && i + 1 < argc) {
      parent_arg = argv[++i];
    } else {
      std::cerr << "Error: unknown arg: " << arg << std::endl;
      return EXIT_FAILURE;
    }
  }

  TFile* f = TFile::Open(input_root.Data(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << input_root << std::endl;
    return EXIT_FAILURE;
  }

  Int_t run_number = DetectRunNumber(f, input_root);
  if (run_number < 0) {
    run_number = 0;
    std::cout << "[WARN] run number not detected; using 0 for output paths" << std::endl;
  } else {
    std::cout << "Run number: " << run_number << std::endl;
  }

  const TString map_path = ANALYZER_DIR + "/param/TPCPOS/e72/TPCPositionCorrectionMap";
  const TString parent_path = parent_arg.IsNull() ? map_path : parent_arg;
  if (parent_arg.IsNull()) {
    std::cout << "[WARN] --parent not given; assuming the helix was reconstructed with "
              << map_path << std::endl;
  }

  Int_t n_y = 0;
  Double_t y0 = 0.0;
  Double_t dy = 0.0;
  std::map<PadKey, std::vector<std::pair<Double_t, std::pair<Double_t, Double_t>>>> grid;
  if (!ReadParentMap(parent_path, n_y, y0, dy, grid)) {
    std::cerr << "Error: cannot read parent map " << parent_path << std::endl;
    f->Close();
    return EXIT_FAILURE;
  }
  std::cout << "Parent map: " << parent_path << " (n_y=" << n_y << ")" << std::endl;
  std::cout << "model=" << ModelName(model) << " threshold=" << threshold
            << " ave=" << (use_ave ? 1 : 0) << " debug=" << (debug_mode ? 1 : 0)
            << " alpha=" << kAlpha << " update=" << (overwrite ? "overwrite" : "add")
            << std::endl;

  TTree* tree = dynamic_cast<TTree*>(f->Get("tpc"));
  if (!tree) {
    std::cerr << "Error: tree 'tpc' not found" << std::endl;
    f->Close();
    return EXIT_FAILURE;
  }

  std::map<PadKey, PadAgg> aggs;
  std::map<PadKey, TH1D*> h_rt;
  std::map<PadKey, TH1D*> h_rx;
  std::map<PadKey, TH1D*> h_rz;
  Long64_t n_tracks_sel = 0;
  Long64_t n_hits_sel = 0;
  const Long64_t nentries = tree->GetEntries();

  auto EnsureHist = [](
    std::map<PadKey, TH1D*>& store,
    const PadKey& key,
    const char* tag,
    const char* xlabel
  ) -> TH1D* {
    auto it = store.find(key);
    if (it != store.end()) {
      return it->second;
    }
    TH1D* h = new TH1D(
      Form("hpos_%s_L%02d_R%03d", tag, key.first, key.second),
      Form("L%02d R%03d %s;%s [mm];entries", key.first, key.second, tag, xlabel),
      100, -10.0, 10.0
    );
    h->SetDirectory(nullptr);
    store[key] = h;
    return h;
  };

  {
    TTreeReader reader(tree);
    TTreeReaderValue<std::vector<Int_t>> nhtrack(reader, "nhtrack");
    TTreeReaderValue<std::vector<Double_t>> helix_dz(reader, "helix_dz");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> hitlayer(reader, "hitlayer");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> residual_x(reader, "residual_x");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> residual_z(reader, "residual_z");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_x(reader, "hitpos_x");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_z(reader, "hitpos_z");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> row_center(
      reader, "track_cluster_row_center"
    );

    Long64_t iev = 0;
    while (reader.Next()) {
      PrintProgress(iev, nentries);
      ++iev;
      const auto& nh = *nhtrack;
      const auto& dz = *helix_dz;
      const auto& layers = *hitlayer;
      const auto& rxv = *residual_x;
      const auto& rzv = *residual_z;
      const auto& hx = *hitpos_x;
      const auto& hz = *hitpos_z;
      const auto& rowc = *row_center;
      const Int_t n_trk = static_cast<Int_t>(nh.size());
      for (Int_t it = 0; it < n_trk; ++it) {
        if (nh[it] < kMinHitsTrack) {
          continue;
        }
        if (it >= static_cast<Int_t>(dz.size()) || std::abs(dz[it]) >= kMaxAbsDz) {
          continue;
        }
        if (it >= static_cast<Int_t>(layers.size())) {
          continue;
        }
        ++n_tracks_sel;
        const auto& ly = layers[it];
        const auto& rx = rxv[it];
        const auto& rz = rzv[it];
        const auto& px = hx[it];
        const auto& pz = hz[it];
        const auto& rc = rowc[it];
        const Int_t n_hit = static_cast<Int_t>(ly.size());
        for (Int_t ih = 0; ih < n_hit; ++ih) {
          if (ih >= static_cast<Int_t>(rc.size()) || ih >= static_cast<Int_t>(rx.size()) ||
              ih >= static_cast<Int_t>(rz.size()) || ih >= static_cast<Int_t>(px.size()) ||
              ih >= static_cast<Int_t>(pz.size())) {
            continue;
          }
          const Int_t layer = static_cast<Int_t>(std::lround(ly[ih]));
          const Int_t row = static_cast<Int_t>(std::lround(rc[ih]));
          if (layer < 0 || layer >= tpc::NumOfLayersTPC || row < 0) {
            continue;
          }
          const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
          if (row >= n_pad) {
            continue;
          }
          const PadKey key{layer, row};
          const Double_t phi = std::atan2(px[ih], pz[ih] - tpc::ZTarget);
          const Double_t r_t = rx[ih] * std::cos(phi) - rz[ih] * std::sin(phi);
          auto& a = aggs[key];
          a.n_hit += 1;
          a.sum_rx += rx[ih];
          a.sum_rz += rz[ih];
          a.sum_rt += r_t;
          ++n_hits_sel;
          if (model == Model::kTan) {
            EnsureHist(h_rt, key, "rt", "r_{t}")->Fill(r_t);
          } else {
            EnsureHist(h_rx, key, "rx", "r_{x}")->Fill(rx[ih]);
            EnsureHist(h_rz, key, "rz", "r_{z}")->Fill(rz[ih]);
          }
        }
      }
    }
    std::fprintf(stderr, "\n");
  }
  f->Close();

  std::cout << "Selected tracks=" << n_tracks_sel << " hits=" << n_hits_sel
            << " occupied_pads=" << aggs.size() << std::endl;

  std::map<PadKey, PadCorr> corrs;
  Int_t n_fitted = 0;
  for (const auto& kv : aggs) {
    const Int_t layer = kv.first.first;
    const Int_t row = kv.first.second;
    const PadAgg& a = kv.second;
    if (a.n_hit < threshold) {
      continue;
    }
    if (tpc::IsPadOnCenterFrame(layer, row)) {
      continue;
    }
    PadCorr c;
    c.fitted = kTRUE;
    c.n_hit = a.n_hit;
    const Double_t mean_rx = a.sum_rx / static_cast<Double_t>(a.n_hit);
    const Double_t mean_rz = a.sum_rz / static_cast<Double_t>(a.n_hit);
    const Double_t mean_rt = a.sum_rt / static_cast<Double_t>(a.n_hit);
    if (model == Model::kTan) {
      c.ct = -kAlpha * mean_rt;
      Double_t e_tx = 0.0;
      Double_t e_tz = 0.0;
      PadUnitTangent(layer, row, e_tx, e_tz);
      c.vx = c.ct * e_tx;
      c.vz = c.ct * e_tz;
    } else {
      c.vx = -kAlpha * mean_rx;
      c.vz = -kAlpha * mean_rz;
      c.ct = 0.0;
    }
    corrs[kv.first] = c;
    ++n_fitted;
  }
  std::cout << "Fitted pads (n_hit>=" << threshold << "): " << n_fitted << std::endl;

  Int_t n_ave = 0;
  if (use_ave) {
    std::map<PadKey, std::pair<Double_t, Int_t>> asad_ct;   // tan
    std::map<PadKey, std::pair<Double_t, Int_t>> asad_vx;   // xz
    std::map<PadKey, std::pair<Double_t, Int_t>> asad_vz;
    for (const auto& kv : corrs) {
      if (!kv.second.fitted) {
        continue;
      }
      const Int_t layer = kv.first.first;
      const Int_t row = kv.first.second;
      const Int_t asad = tpc::GetASADId(layer, row);
      const PadKey ak{layer, asad};
      if (model == Model::kTan) {
        asad_ct[ak].first += kv.second.ct;
        asad_ct[ak].second += 1;
      } else {
        asad_vx[ak].first += kv.second.vx;
        asad_vx[ak].second += 1;
        asad_vz[ak].first += kv.second.vz;
        asad_vz[ak].second += 1;
      }
    }
    for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
      const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
      for (Int_t row = 0; row < n_pad; ++row) {
        if (tpc::IsPadOnCenterFrame(layer, row)) {
          continue;
        }
        const PadKey key{layer, row};
        if (corrs.count(key) && corrs[key].fitted) {
          continue;
        }
        const Int_t asad = tpc::GetASADId(layer, row);
        const PadKey ak{layer, asad};
        PadCorr c;
        c.from_ave = kTRUE;
        c.n_hit = aggs.count(key) ? aggs[key].n_hit : 0;
        if (model == Model::kTan) {
          const auto it = asad_ct.find(ak);
          if (it == asad_ct.end() || it->second.second <= 0) {
            continue;
          }
          c.ct = it->second.first / static_cast<Double_t>(it->second.second);
          Double_t e_tx = 0.0;
          Double_t e_tz = 0.0;
          PadUnitTangent(layer, row, e_tx, e_tz);
          c.vx = c.ct * e_tx;
          c.vz = c.ct * e_tz;
        } else {
          const auto itx = asad_vx.find(ak);
          const auto itz = asad_vz.find(ak);
          if (itx == asad_vx.end() || itx->second.second <= 0 || itz == asad_vz.end()) {
            continue;
          }
          c.vx = itx->second.first / static_cast<Double_t>(itx->second.second);
          c.vz = itz->second.first / static_cast<Double_t>(itz->second.second);
        }
        corrs[key] = c;
        ++n_ave;
      }
    }
    std::cout << "ASAD-mean filled pads: " << n_ave << std::endl;
  }

  Int_t n_center = 0;
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < n_pad; ++row) {
      if (!tpc::IsPadOnCenterFrame(layer, row)) {
        continue;
      }
      PadCorr c;
      c.center_frame = kTRUE;
      c.vx = 0.0;
      c.vz = 0.0;
      c.ct = 0.0;
      corrs[{layer, row}] = c;
      ++n_center;
    }
  }
  std::cout << "Center-frame pads forced to (0,0): " << n_center << std::endl;

  // Add the (Y-flat) increments onto the parent grid (or replace with --overwrite);
  // center-frame pads are reset to (0,0).
  Int_t n_map_touch = 0;
  for (const auto& kv : corrs) {
    auto it = grid.find(kv.first);
    if (it == grid.end()) {
      // Missing pad in parent: treat parent as zero on the header lattice.
      std::vector<std::pair<Double_t, std::pair<Double_t, Double_t>>> nodes;
      for (Int_t iy = 0; iy < n_y; ++iy) {
        nodes.push_back({y0 + iy * dy, {kv.second.vx, kv.second.vz}});
      }
      grid[kv.first] = nodes;
      ++n_map_touch;
      continue;
    }
    for (auto& node : it->second) {
      if (kv.second.center_frame) {
        node.second = {0.0, 0.0};
      } else if (overwrite) {
        node.second = {kv.second.vx, kv.second.vz};
      } else {
        node.second = {node.second.first + kv.second.vx, node.second.second + kv.second.vz};
      }
    }
    ++n_map_touch;
  }

  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
  const TString pdf_path =
    img_dir + Form("/tpc_pos_calib_%s.pdf", ModelName(model));

  std::cout << "formula_id=" << kFormulaId
            << " n_fitted=" << n_fitted << " n_ave=" << n_ave
            << " n_center=" << n_center << " n_map_touch=" << n_map_touch
            << " pdf=" << pdf_path << std::endl;

  // Sort fitted pads for stable PDF page order (layer, row).
  std::vector<PadKey> fitted_keys;
  fitted_keys.reserve(static_cast<std::size_t>(n_fitted));
  for (const auto& kv : corrs) {
    if (kv.second.fitted) {
      fitted_keys.push_back(kv.first);
    }
  }
  std::sort(fitted_keys.begin(), fitted_keys.end());

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // Per-pad residual pages (offset-style 4x4).
  TCanvas canvas("c_pos", "c_pos", 1200, 1000);
  canvas.Print(pdf_path + "[");
  const Int_t nx = 4;
  const Int_t ny = 4;
  const Int_t n_per_page = nx * ny;
  Int_t ipad = 0;
  std::vector<TLine*> page_lines;
  Int_t n_pad_pages = 0;

  auto FlushPadPage = [&]() {
    if (ipad <= 0) {
      return;
    }
    canvas.Print(pdf_path);
    for (auto* ln : page_lines) {
      delete ln;
    }
    page_lines.clear();
    ipad = 0;
    ++n_pad_pages;
  };

  auto DrawMeanLine = [&](TH1D* h, Double_t mean) {
    if (!h) {
      return;
    }
    const Double_t y0h = h->GetMinimum();
    const Double_t y1h = std::max(h->GetMaximum() * 1.05, 1.0);
    auto* vLine = new TLine(mean, y0h, mean, y1h);
    vLine->SetLineStyle(2);
    vLine->SetLineColor(kBlue + 1);
    vLine->SetLineWidth(2);
    vLine->Draw("same");
    page_lines.push_back(vLine);
  };

  auto DrawOneHist = [&](TH1D* h, Double_t mean, const TString& title) {
    if (ipad == 0) {
      canvas.Clear();
      canvas.Divide(nx, ny);
    }
    canvas.cd(++ipad);
    if (h) {
      h->SetTitle(title);
      h->Draw("HIST");
      DrawMeanLine(h, mean);
    } else {
      TLatex lab;
      lab.SetNDC();
      lab.SetTextSize(0.06);
      lab.DrawLatex(0.15, 0.5, "no hist");
    }
    if (ipad >= n_per_page) {
      FlushPadPage();
    }
  };

  for (const PadKey& key : fitted_keys) {
    const auto& c = corrs[key];
    const auto& a = aggs[key];
    if (model == Model::kTan) {
      const Double_t mean_rt = a.sum_rt / static_cast<Double_t>(a.n_hit);
      TH1D* h = h_rt.count(key) ? h_rt[key] : nullptr;
      DrawOneHist(
        h, mean_rt,
        Form("L%02d R%03d r_{t}  n=%lld  #LT r_{t}#GT=%.3f  C_{t}=%.3f",
             key.first, key.second, static_cast<long long>(a.n_hit), mean_rt, c.ct)
      );
    } else {
      const Double_t mean_rx = a.sum_rx / static_cast<Double_t>(a.n_hit);
      const Double_t mean_rz = a.sum_rz / static_cast<Double_t>(a.n_hit);
      TH1D* hx = h_rx.count(key) ? h_rx[key] : nullptr;
      TH1D* hz = h_rz.count(key) ? h_rz[key] : nullptr;
      DrawOneHist(
        hx, mean_rx,
        Form("L%02d R%03d r_{x}  n=%lld  #LT r_{x}#GT=%.3f  v_{x}=%.3f",
             key.first, key.second, static_cast<long long>(a.n_hit), mean_rx, c.vx)
      );
      DrawOneHist(
        hz, mean_rz,
        Form("L%02d R%03d r_{z}  n=%lld  #LT r_{z}#GT=%.3f  v_{z}=%.3f",
             key.first, key.second, static_cast<long long>(a.n_hit), mean_rz, c.vz)
      );
    }
  }
  FlushPadPage();
  std::cout << "Per-pad PDF pages: " << n_pad_pages << " (fitted pads=" << n_fitted << ")"
            << std::endl;

  // Summary maps: square canvas, no stat box.
  TCanvas csum("c_pos_sum", "c_pos_sum", 1000, 1000);
  TH2Poly* hCt = new TH2Poly(
    "hTPC_Pos_Ct",
    Form("TPCPOS %s (Run %05d);Z [mm];X [mm]", ModelName(model), run_number),
    -300, 300, -300, 300
  );
  hCt->SetDirectory(nullptr);
  tpc::InitializeHistograms(hCt);

  TH1D* h1 = new TH1D(
    "hTPC_Pos_1D",
    Form("TPCPOS update magnitude (%s);value [mm];pads", ModelName(model)),
    100, -5.0, 5.0
  );
  h1->SetDirectory(nullptr);

  for (const auto& kv : corrs) {
    if (!kv.second.fitted && !kv.second.from_ave) {
      continue;
    }
    const Int_t bin = tpc::GetPadId(kv.first.first, kv.first.second) + 1;
    const Double_t val =
      (model == Model::kTan) ? kv.second.ct : std::hypot(kv.second.vx, kv.second.vz);
    hCt->SetBinContent(bin, val);
    h1->Fill(val);
  }

  csum.cd();
  csum.SetFixedAspectRatio();
  hCt->Draw("COLZ");
  gPad->SetFixedAspectRatio();
  csum.Print(pdf_path);

  csum.Clear();
  h1->Draw("HIST");
  csum.Print(pdf_path);

  // Arrow field: draw shaft (TLine) + small head so the vector body is visible.
  csum.Clear();
  TH2Poly* hArrowBase = new TH2Poly(
    "hTPC_Pos_ArrowBase",
    Form("TPCPOS (vx,vz) arrows (%s, Run %05d);Z [mm];X [mm]", ModelName(model), run_number),
    -300, 300, -300, 300
  );
  hArrowBase->SetDirectory(nullptr);
  tpc::InitializeHistograms(hArrowBase);
  hArrowBase->Draw("COLZ");
  gPad->SetFixedAspectRatio();
  std::vector<TObject*> arrow_objs;
  constexpr Double_t kArrowScale = 15.0;
  Int_t iarrow = 0;
  for (const auto& kv : corrs) {
    if (!kv.second.fitted && !kv.second.from_ave) {
      continue;
    }
    const Int_t pad_id = tpc::GetPadId(kv.first.first, kv.first.second);
    if ((pad_id % 4) != 0) {
      continue;
    }
    const TVector3 pos = tpc::getPosition(kv.first.first, static_cast<Double_t>(kv.first.second));
    const Double_t z0 = pos.Z();
    const Double_t x0 = pos.X();
    const Double_t z1 = z0 + kArrowScale * kv.second.vz;
    const Double_t x1 = x0 + kArrowScale * kv.second.vx;
    auto* ln = new TLine(z0, x0, z1, x1);
    ln->SetLineColor(kRed + 1);
    ln->SetLineWidth(2);
    ln->Draw();
    arrow_objs.push_back(ln);
    auto* ar = new TArrow(z0, x0, z1, x1, 0.008, "|>");
    ar->SetLineColor(kRed + 1);
    ar->SetFillColor(kRed + 1);
    ar->SetLineWidth(2);
    ar->Draw();
    arrow_objs.push_back(ar);
    ++iarrow;
  }
  TLatex note;
  note.SetNDC();
  note.SetTextSize(0.03);
  note.DrawLatex(0.12, 0.92, Form("arrow scale = %.0f #times (vx,vz) [mm]", kArrowScale));
  csum.Print(pdf_path);
  csum.Print(pdf_path + "]");

  for (auto* o : arrow_objs) {
    delete o;
  }
  delete hArrowBase;
  delete hCt;
  delete h1;
  for (auto& kv : h_rt) {
    delete kv.second;
  }
  for (auto& kv : h_rx) {
    delete kv.second;
  }
  for (auto& kv : h_rz) {
    delete kv.second;
  }
  std::cout << "Wrote " << pdf_path << " (arrows=" << iarrow << ")" << std::endl;

  if (debug_mode) {
    std::cout << "--debug: skip writing " << map_path << std::endl;
  } else {
    if (!WriteMap(map_path, n_y, y0, dy, grid)) {
      std::cerr << "Error: failed to write map " << map_path << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Wrote map " << map_path << " (touched pads=" << n_map_touch << ")"
              << std::endl;
  }

  return EXIT_SUCCESS;
}
