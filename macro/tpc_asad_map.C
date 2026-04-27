// -*- C++ -*-
// TPC X-Z pad map colored by ASAD id + ASAD/CoBo table.
// Also supports multi-page PDF with per-ASAD highlight pages.
//
// Usage:
//   root -l
//   .L macro/tpc_asad_map.C
//   tpc_asad_map();
//   tpc_asad_map_save(); // default: results/img/tpc_asad_map.pdf
//   tpc_asad_map_save("results/img/custom.pdf");

#include <TCanvas.h>
#include <TColor.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include "../common/include/TPCPadHelper.hh"

namespace
{

void TPCPadTemplateXZ(TH2Poly* h)
{
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
}

void SetDiscreteAsadStyle(TH2Poly* h)
{
  // Integer-step levels: prevents smooth interpolation and clarifies boundaries.
  static Double_t levels[tpc::NumOfAsadTPC + 1];
  for (Int_t i = 0; i <= tpc::NumOfAsadTPC; ++i) {
    levels[i] = static_cast<Double_t>(i) - 0.5;
  }
  h->SetContour(tpc::NumOfAsadTPC + 1, levels);
  h->SetMinimum(-0.5);
  h->SetMaximum(tpc::NumOfAsadTPC - 0.5);
}

void SetAsadPaletteVividFromBlue()
{
  // Base starts from blue, then intentionally vivid colors.
  static Int_t colors[] = {
    kBlue + 1, kCyan + 1, kSpring + 9, kYellow + 1,
    kOrange + 7, kRed + 1, kMagenta + 1, kViolet + 1
  };
  gStyle->SetPalette(sizeof(colors) / sizeof(Int_t), colors);
}

void SetSingleAsadPaletteBlueVsFlashy()
{
  // 0: base blue, 1: highlighted flashy color.
  static Int_t colors[] = {kAzure + 1, kPink + 7};
  gStyle->SetPalette(sizeof(colors) / sizeof(Int_t), colors);
}

void FillAsadContent(TH2Poly* h, std::set<std::pair<Int_t, Int_t>>* pairSet = nullptr)
{
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    Int_t nPad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < nPad; ++row) {
      const Int_t asad = tpc::GetASADId(layer, row);
      const Int_t cobo = tpc::GetCoBoId(layer, row);
      if (pairSet) pairSet->insert({asad, cobo});

      TVector3 pos = tpc::getPosition(layer, static_cast<Double_t>(row));
      if (!std::isfinite(pos.Z()) || !std::isfinite(pos.X())) continue;

      Int_t bin = h->FindBin(pos.Z(), pos.X());
      if (bin > 0) h->SetBinContent(bin, static_cast<Double_t>(asad));
    }
  }
}

void FillSingleAsad(TH2Poly* h, Int_t selectedAsad)
{
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    Int_t nPad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < nPad; ++row) {
      const Int_t asad = tpc::GetASADId(layer, row);
      TVector3 pos = tpc::getPosition(layer, static_cast<Double_t>(row));
      if (!std::isfinite(pos.Z()) || !std::isfinite(pos.X())) continue;
      Int_t bin = h->FindBin(pos.Z(), pos.X());
      if (bin > 0) h->SetBinContent(bin, asad == selectedAsad ? 1.0 : 0.0);
    }
  }
}

void FillDeadChannelMap(TH2Poly* h)
{
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    Int_t nPad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
    for (Int_t row = 0; row < nPad; ++row) {
      TVector3 pos = tpc::getPosition(layer, static_cast<Double_t>(row));
      if (!std::isfinite(pos.Z()) || !std::isfinite(pos.X())) continue;
      Int_t bin = h->FindBin(pos.Z(), pos.X());
      if (bin > 0) h->SetBinContent(bin, tpc::IsDead(layer, row) ? 1.0 : 0.0);
    }
  }
}

TString EnsureDefaultSavePath()
{
  const TString dir = "results/img";
  gSystem->mkdir("results", kTRUE);
  gSystem->mkdir(dir, kTRUE);
  return dir + "/tpc_asad_map.pdf";
}

} // namespace

//______________________________________________________________________________
void tpc_asad_map()
{
  gStyle->SetOptStat(0);
  SetAsadPaletteVividFromBlue();
  gStyle->SetNumberContours(tpc::NumOfAsadTPC + 1);
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetTitleSize(0.045, "XYZ");

  const Double_t zmin = -300.0, zmax = 300.0;
  const Double_t xmin = -300.0, xmax = 300.0;

  auto* c = new TCanvas("c_tpc_asad_map", "TPC ASAD / CoBo map", 1400, 900);
  auto* padMap = new TPad("padMap", "padMap", 0.0, 0.0, 0.70, 1.0);
  auto* padTab = new TPad("padTab", "padTab", 0.70, 0.0, 1.0, 1.0);
  padMap->SetRightMargin(0.16);
  padMap->SetLeftMargin(0.10);
  padMap->Draw();
  padTab->Draw();

  std::set<std::pair<Int_t, Int_t>> asadCoboPairs;

  padMap->cd();
  auto* h = new TH2Poly("h_tpc_asad", "TPC pad map: ASAD id;Z [mm];X [mm]", zmin, zmax, xmin, xmax);
  TPCPadTemplateXZ(h);
  FillAsadContent(h, &asadCoboPairs);
  SetDiscreteAsadStyle(h);
  h->Draw("COLZ");

  padTab->cd();
  padTab->SetMargin(0.05, 0.02, 0.03, 0.02);
  auto* pt = new TPaveText(0.0, 0.0, 1.0, 1.0, "NDC");
  pt->SetFillColor(kWhite);
  pt->SetBorderSize(1);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.035); // larger than before
  pt->AddText("ASAD - CoBo table");
  pt->AddText(" ");

  std::vector<std::pair<Int_t, Int_t>> pairs(asadCoboPairs.begin(), asadCoboPairs.end());
  std::sort(pairs.begin(), pairs.end());
  for (const auto& p : pairs) {
    pt->AddText(TString::Format("ASAD %2d  <->  CoBo %d", p.first, p.second));
  }
  pt->Draw();

  c->Update();
}

//______________________________________________________________________________
void tpc_asad_map_save(const char* pdfPath = "")
{
  gStyle->SetOptStat(0);

  TString out = (pdfPath && TString(pdfPath).Length() > 0) ? TString(pdfPath) : EnsureDefaultSavePath();

  const Double_t zmin = -300.0, zmax = 300.0;
  const Double_t xmin = -300.0, xmax = 300.0;

  // page 1: full map + table
  tpc_asad_map();
  TCanvas* cMain = dynamic_cast<TCanvas*>(gROOT->FindObject("c_tpc_asad_map"));
  if (!cMain) {
    std::cerr << "tpc_asad_map_save: failed to build main canvas" << std::endl;
    return;
  }
  cMain->Print(out + "(");

  // page 2..: each ASAD highlighted
  auto* cSingle = new TCanvas("c_tpc_asad_single", "TPC single ASAD pages", 1200, 850);
  cSingle->SetRightMargin(0.16);
  cSingle->SetLeftMargin(0.10);

  static Double_t lv2[3] = {-0.5, 0.5, 1.5};
  for (Int_t asad = 0; asad < tpc::NumOfAsadTPC; ++asad) {
    auto* h = new TH2Poly(Form("h_tpc_asad_%02d", asad),
                          Form("TPC pads: ASAD %02d highlighted;Z [mm];X [mm]", asad),
                          zmin, zmax, xmin, xmax);
    TPCPadTemplateXZ(h);
    FillSingleAsad(h, asad);
    h->SetContour(3, lv2);
    h->SetMinimum(-0.5);
    h->SetMaximum(1.5);
    SetSingleAsadPaletteBlueVsFlashy();
    gStyle->SetNumberContours(3);
    h->Draw("COLZ");

    cSingle->Print(out);
    cSingle->Clear();
  }

  // last page: dead channel map
  auto* hDead = new TH2Poly("h_tpc_dead", "TPC dead channels;Z [mm];X [mm]", zmin, zmax, xmin, xmax);
  TPCPadTemplateXZ(hDead);
  FillDeadChannelMap(hDead);
  hDead->SetContour(3, lv2);
  hDead->SetMinimum(-0.5);
  hDead->SetMaximum(1.5);
  SetSingleAsadPaletteBlueVsFlashy(); // 0: blue(normal), 1: flashy(dead)
  gStyle->SetNumberContours(3);
  hDead->Draw("COLZ");
  cSingle->Print(out);
  cSingle->Clear();

  cSingle->Print(out + ")");
  std::cout << "Wrote " << out << std::endl;
}

