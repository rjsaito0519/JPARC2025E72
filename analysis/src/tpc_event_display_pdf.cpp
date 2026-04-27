// -*- C++ -*-
// UserTPCHit `tpc` tree -> X-Z pad map filled by dE (COLZ), multi-page PDF
// Usage: tpc_event_display_pdf -n <N> [-o <out.pdf>] [--start <entry>] <input.root>
// If -o is omitted, writes under OUTPUT_DIR/img/runXXXXX/ (same layout as calibration tools; see paths.h).
// Draw condition: hit time cut 5 <= tTpc <= 155 and dead-pad masking via tpc::IsDead().

#include "paths.h"
#include "ana_helper.h"
#include "TPCPadHelper.hh"

#include <TCanvas.h>
#include <TFile.h>
#include <TH2Poly.h>
#include <TGraph.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{

void
TPCPadTemplateXZ(TH2Poly* h)
{
  // TH2Poly X = Z (horizontal), TH2Poly Y = X (vertical)
  Double_t X[5];
  Double_t Y[5];
  for (Int_t l = 0; l < tpc::NumOfLayersTPC; ++l) {
    Double_t pLength = tpc::padParameter[l][tpc::kLength];
    Double_t st = (180. - (360. / tpc::padParameter[l][tpc::kNumOfDivision]) *
                   tpc::padParameter[l][tpc::kNumOfPad] / 2.);
    Double_t sTheta = (-1 + st / 180.) * TMath::Pi();
    Double_t dTheta = (360. / tpc::padParameter[l][tpc::kNumOfDivision]) / 180. * TMath::Pi();
    Double_t cRad = tpc::padParameter[l][tpc::kRadius];
    Int_t nPad = static_cast<Int_t>(tpc::padParameter[l][tpc::kNumOfPad]);
    for (Int_t j = 0; j < nPad; ++j) {
      Double_t theta1 = j * dTheta + sTheta;
      Double_t theta2 = (j + 1) * dTheta + sTheta;
      Double_t r_in = cRad - (pLength / 2.);
      Double_t r_out = cRad + (pLength / 2.);

      X[1] = r_out * TMath::Cos(theta1) + tpc::ZTarget;
      X[2] = r_out * TMath::Cos(theta2) + tpc::ZTarget;
      X[3] = r_in * TMath::Cos(theta2) + tpc::ZTarget;
      X[4] = r_in * TMath::Cos(theta1) + tpc::ZTarget;
      X[0] = X[4];

      Y[1] = r_out * TMath::Sin(theta1);
      Y[2] = r_out * TMath::Sin(theta2);
      Y[3] = r_in * TMath::Sin(theta2);
      Y[4] = r_in * TMath::Sin(theta1);
      Y[0] = Y[4];

      h->AddBin(5, X, Y);
    }
  }
  h->SetMaximum(0x1000);
}

void
usage(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " -n <N> [-o <out.pdf>] [--start <entry>] [--asad|-asad <id[,id...]>] <input.root>\n"
            << "  -n     number of consecutive events to draw (from --start)\n"
            << "  -o     output PDF path (default: OUTPUT_DIR/img/runXXXXX/tpc_event_display_runXXXXX.pdf)\n"
            << "  --start  first tree entry (default 0)\n"
            << "  --asad|-asad  highlight specified ASAD ids (0-30), repeatable and comma-separated\n";
}

void
BuildPadPolygon(Int_t layer, Int_t row, Double_t z[5], Double_t x[5])
{
  const Double_t pLength = tpc::padParameter[layer][tpc::kLength];
  const Double_t st = (180. - (360. / tpc::padParameter[layer][tpc::kNumOfDivision]) *
                       tpc::padParameter[layer][tpc::kNumOfPad] / 2.);
  const Double_t sTheta = (-1 + st / 180.) * TMath::Pi();
  const Double_t dTheta = (360. / tpc::padParameter[layer][tpc::kNumOfDivision]) / 180. * TMath::Pi();
  const Double_t cRad = tpc::padParameter[layer][tpc::kRadius];
  const Double_t theta1 = row * dTheta + sTheta;
  const Double_t theta2 = (row + 1) * dTheta + sTheta;
  const Double_t r_in = cRad - (pLength / 2.);
  const Double_t r_out = cRad + (pLength / 2.);

  z[0] = r_out * TMath::Cos(theta1) + tpc::ZTarget;
  z[1] = r_out * TMath::Cos(theta2) + tpc::ZTarget;
  z[2] = r_in * TMath::Cos(theta2) + tpc::ZTarget;
  z[3] = r_in * TMath::Cos(theta1) + tpc::ZTarget;
  z[4] = z[0];

  x[0] = r_out * TMath::Sin(theta1);
  x[1] = r_out * TMath::Sin(theta2);
  x[2] = r_in * TMath::Sin(theta2);
  x[3] = r_in * TMath::Sin(theta1);
  x[4] = x[0];
}

void
DrawPadPolygon(Int_t layer, Int_t row, Int_t fillColor, Double_t alpha)
{
  Double_t z[5];
  Double_t x[5];
  BuildPadPolygon(layer, row, z, x);

  TGraph* poly = new TGraph(5, z, x);
  poly->SetFillColorAlpha(fillColor, alpha);
  poly->SetLineWidth(0);
  poly->Draw("F same");
}

std::set<Int_t>
BuildFramePadSet()
{
  std::set<Int_t> framePads;
  framePads.insert(std::begin(tpc::padOnCenterFrame), std::end(tpc::padOnCenterFrame));
  framePads.insert(std::begin(tpc::padOnSectionFrame), std::end(tpc::padOnSectionFrame));
  return framePads;
}

void
DrawFramePadOverlay(const std::set<Int_t>& framePads)
{
  if (framePads.empty()) return;

  const Int_t fillColor = kWhite;
  const Double_t fillAlpha = 0.25;
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    const Int_t nPad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < nPad; ++row) {
      const Int_t padID = tpc::GetPadId(layer, row);
      if (framePads.find(padID) == framePads.end()) continue;
      Double_t z[5];
      Double_t x[5];
      BuildPadPolygon(layer, row, z, x);

      TGraph* poly = new TGraph(5, z, x);
      poly->SetFillColorAlpha(fillColor, fillAlpha);
      poly->SetLineWidth(0);
      poly->Draw("F same");
    }
  }
}

void
DrawAsadPadOverlay(const std::set<Int_t>& asadFilter, const std::map<Int_t, bool>& asadHasHit)
{
  if (asadFilter.empty()) return;

  const Int_t fillColor = kPink + 1;
  const Double_t alpha = 0.2;
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    const Int_t nPad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < nPad; ++row) {
      const Int_t asad = tpc::GetASADId(layer, row);
      if (asadFilter.find(asad) == asadFilter.end()) continue;

      const bool hasHit = (asadHasHit.find(asad) != asadHasHit.end() && asadHasHit.at(asad));
      (void)hasHit;
      DrawPadPolygon(layer, row, fillColor, alpha);
    }
  }
}

} // namespace

int
main(int argc, char** argv)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  Long64_t nEvents = -1;
  Long64_t startEntry = 0;
  std::string outPdf;
  std::string inPath;
  std::set<Int_t> asadFilter;
  const std::set<Int_t> framePads = BuildFramePadSet();

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
      nEvents = std::atoll(argv[++i]);
    } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      outPdf = argv[++i];
    } else if (std::strcmp(argv[i], "--start") == 0 && i + 1 < argc) {
      startEntry = std::atoll(argv[++i]);
    } else if ((std::strcmp(argv[i], "--asad") == 0 || std::strcmp(argv[i], "-asad") == 0) && i + 1 < argc) {
      std::stringstream ss(argv[++i]);
      std::string token;
      while (std::getline(ss, token, ',')) {
        if (token.empty()) continue;
        const Int_t asad = static_cast<Int_t>(std::strtol(token.c_str(), nullptr, 10));
        if (asad < 0 || asad >= tpc::NumOfAsadTPC) {
          std::cerr << "Error: --asad out of range: " << asad << " (valid: 0-" << (tpc::NumOfAsadTPC - 1) << ")\n";
          return 1;
        }
        asadFilter.insert(asad);
      }
    } else if (argv[i][0] == '-') {
      std::cerr << "Unknown option: " << argv[i] << std::endl;
      usage(argv[0]);
      return 1;
    } else {
      inPath = argv[i];
    }
  }

  if (nEvents <= 0 || inPath.empty()) {
    usage(argv[0]);
    return 1;
  }

  TFile f(inPath.c_str(), "READ");
  if (!f.IsOpen() || f.IsZombie()) {
    std::cerr << "Error: cannot open " << inPath << std::endl;
    return 1;
  }

  TTree* tree = dynamic_cast<TTree*>(f.Get("tpc"));
  if (!tree) {
    std::cerr << "Error: tree 'tpc' not found in " << inPath << std::endl;
    return 1;
  }

  if (!tree->GetBranch("run_number") || !tree->GetBranch("event_number") ||
      !tree->GetBranch("nhTpc") || !tree->GetBranch("padTpc") ||
      !tree->GetBranch("deTpc") || !tree->GetBranch("tTpc")) {
    std::cerr << "Error: required branches missing (run_number, event_number, nhTpc, padTpc, deTpc, tTpc)\n";
    return 1;
  }

  UInt_t run_number = 0;
  UInt_t event_number = 0;
  Int_t nhTpc = 0;
  std::vector<Int_t>* padTpc = nullptr;
  std::vector<Double_t>* deTpc = nullptr;
  std::vector<Double_t>* tTpc = nullptr;

  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchAddress("event_number", &event_number);
  tree->SetBranchAddress("nhTpc", &nhTpc);
  tree->SetBranchAddress("padTpc", &padTpc);
  tree->SetBranchAddress("deTpc", &deTpc);
  tree->SetBranchAddress("tTpc", &tTpc);

  const Long64_t nentries = tree->GetEntries();
  if (startEntry < 0 || startEntry >= nentries) {
    std::cerr << "Error: --start " << startEntry << " out of range [0, " << nentries << ")\n";
    return 1;
  }

  if (outPdf.empty()) {
    tree->GetEntry(startEntry);
    const Int_t run = static_cast<Int_t>(run_number);
    const TString imgDir = ana_helper::get_img_dir(OUTPUT_DIR, run);
    outPdf = std::string(TString::Format("%s/tpc_event_display_run%05d.pdf", imgDir.Data(), run).Data());
    std::cout << "Default output (calib-style img dir): " << outPdf << std::endl;
  }

  const Double_t zmin = -300.0;
  const Double_t zmax = 300.0;
  const Double_t xmin = -300.0;
  const Double_t xmax = 300.0;
  const Double_t tMin = 5.0;
  const Double_t tMax = 155.0;

  std::cout << "Applying hit selection: " << tMin << " <= tTpc <= " << tMax
            << ", and skipping dead pads via tpc::IsDead()." << std::endl;
  std::cout << "Note: tpc::IsNoisy() is not available in this TPCPadHelper; noisy masking is not applied." << std::endl;
  std::cout << "Frame pad highlight enabled: " << framePads.size() << " pad(s)." << std::endl;
  if (!asadFilter.empty()) {
    std::cout << "ASAD highlight enabled for:";
    for (const auto& asad : asadFilter) std::cout << " " << asad;
    std::cout << std::endl;
  }

  TCanvas canvas("c_tpc_ed", "TPC X-Z", 900, 800);
  const TString outPath(outPdf.c_str());

  Long64_t nDrawn = 0;
  for (Long64_t entry = startEntry; entry < nentries && nDrawn < nEvents; ++entry) {
    tree->GetEntry(entry);

    if (!padTpc || !deTpc || !tTpc) {
      std::cerr << "Error: null branch pointer at entry " << entry << std::endl;
      return 1;
    }

    const std::size_t nPad = padTpc->size();
    const std::size_t nDe = deTpc->size();
    const std::size_t nT = tTpc->size();
    if (nPad != nDe || nPad != nT) {
      std::cerr << "Warning: entry " << entry
                << " size mismatch (padTpc=" << nPad
                << ", deTpc=" << nDe
                << ", tTpc=" << nT
                << ") — using min length\n";
    }
    const std::size_t n = std::min(nPad, std::min(nDe, nT));

    std::map<Int_t, Double_t> pad_de_map;
    std::map<Int_t, bool> asadHasHit;
    for (const auto& asad : asadFilter) {
      asadHasHit[asad] = false;
    }
    for (std::size_t i = 0; i < n; ++i) {
      const Int_t pad = (*padTpc)[i];
      if (pad < 0) {
        continue;
      }
      if (tpc::IsDead(pad)) {
        continue;
      }
      const Double_t t = (*tTpc)[i];
      if (t < tMin || t > tMax) {
        continue;
      }
      if (!asadFilter.empty()) {
        const Int_t layer = tpc::getLayerID(pad);
        const Int_t row = tpc::getRowID(pad);
        const Int_t asad = tpc::GetASADId(layer, row);
        if (asadFilter.find(asad) != asadFilter.end()) {
          asadHasHit[asad] = true;
        }
      }

      const Double_t de = (*deTpc)[i];
      const auto it = pad_de_map.find(pad);
      if (it == pad_de_map.end() || it->second < de) {
        pad_de_map[pad] = de;
      }
    }

    // Event-level cut: skip events with no hits after tTpc/dead-pad selection.
    if (pad_de_map.empty()) {
      continue;
    }

    TH2Poly* h =
      new TH2Poly(Form("h_tpc_%lld", static_cast<long long>(entry)),
                  Form("TPC X-Z (Run %u, Ev %u, nh=%d);Z [mm];X [mm]", run_number, event_number,
                       nhTpc),
                  zmin, zmax, xmin, xmax);
    h->SetStats(0);
    TPCPadTemplateXZ(h);
    h->SetMinimum(0);
    h->SetMaximum(1000);

    for (const auto& pair : pad_de_map) {
      const Int_t padid = pair.first;
      const Double_t de = pair.second;
      TVector3 padPos = tpc::getPosition(padid);
      if (!std::isfinite(padPos.Z()) || !std::isfinite(padPos.X())) {
        continue;
      }
      const Int_t bin = h->FindBin(padPos.Z(), padPos.X());
      if (bin > 0) {
        h->SetBinContent(bin, de);
      }
    }

    h->Draw("COLZ");
    DrawFramePadOverlay(framePads);
    DrawAsadPadOverlay(asadFilter, asadHasHit);

    canvas.Modified();
    canvas.Update();

    if (nDrawn == 0) {
      canvas.Print(outPath + "(");
    } else {
      canvas.Print(outPath);
    }
    ++nDrawn;

    // TH2Poly / graphs owned by current pad; clear pad for next event
    canvas.Clear();
  }

  if (nDrawn > 0) {
    canvas.Print(outPath + ")");
  } else {
    std::cerr << "Warning: no events passed selection; no PDF was written." << std::endl;
  }

  if (nDrawn < nEvents) {
    std::cerr << "Warning: requested " << nEvents
              << " events, but only " << nDrawn
              << " passed the selection in range [" << startEntry << ", " << (nentries - 1) << "].\n";
  }
  std::cout << "Wrote " << nDrawn << " page(s) to " << outPdf << std::endl;
  return 0;
}
