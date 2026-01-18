// -*- C++ -*-
// Macro to check TPC tracking results: clustering and Hough transformation
// Usage:
//   root
//   .L check_tpc_track.C
//   set_path("path/to/rootfile.root")
//   event()        // random event
//   event(evnum)   // specific event

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TH3F.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TText.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TVector3.h>
#include <TSystem.h>
#include <vector>
#include <map>
#include <iostream>
#include "../include/TPCPadHelper.hh"

// Global variables
TFile* gMacroFile = nullptr;
TTreeReader* gMacroReader = nullptr;
TRandom3* gMacroRandom = nullptr;
TCanvas* gMacroCanvas = nullptr;

// Event data structure
struct EventData {
  UInt_t runnum;
  UInt_t evnum;
  Int_t nhTpc;
  Int_t nclTpc;
  Int_t ntTpc;
  
  std::vector<Double_t> raw_hitpos_x;
  std::vector<Double_t> raw_hitpos_y;
  std::vector<Double_t> raw_hitpos_z;
  std::vector<Double_t> raw_de;
  std::vector<Int_t> raw_padid;
  std::vector<Int_t> raw_layer;
  std::vector<Int_t> raw_row;
  
  std::vector<Double_t> cluster_x;
  std::vector<Double_t> cluster_y;
  std::vector<Double_t> cluster_z;
  std::vector<Double_t> cluster_de;
  std::vector<Int_t> cluster_layer;
  std::vector<Int_t> cluster_row_center;
  std::vector<Int_t> cluster_houghflag;
  
  std::vector<Double_t> x0Tpc;
  std::vector<Double_t> y0Tpc;
  std::vector<Double_t> u0Tpc;
  std::vector<Double_t> v0Tpc;
  std::vector<Double_t> theta;
  std::vector<Int_t> nhtrack;
  std::vector<std::vector<Double_t>> hitpos_x;
  std::vector<std::vector<Double_t>> hitpos_y;
  std::vector<std::vector<Double_t>> hitpos_z;
};

EventData gEvent;

//______________________________________________________________________________
void
TPC_pad_template(TH2Poly* h)
{
  // For X-Z view: TH2Poly X = actual Z (horizontal), TH2Poly Y = actual X (vertical)
  // According to tpc::getPosition: X = radius * sin(theta), Z = radius * cos(theta) + ZTarget
  Double_t X[5]; // TH2Poly X = actual Z coordinate (horizontal axis)
  Double_t Y[5]; // TH2Poly Y = actual X coordinate (vertical axis)
  for(Int_t l = 0; l < tpc::NumOfLayersTPC; ++l) {
    Double_t pLength = tpc::padParameter[l][tpc::kLength];
    Double_t st = (180. - (360. / tpc::padParameter[l][tpc::kNumOfDivision]) *
                   tpc::padParameter[l][tpc::kNumOfPad] / 2.);
    Double_t sTheta = (-1 + st / 180.) * TMath::Pi();
    Double_t dTheta = (360. / tpc::padParameter[l][tpc::kNumOfDivision]) / 180. * TMath::Pi();
    Double_t cRad = tpc::padParameter[l][tpc::kRadius];
    Int_t nPad = static_cast<Int_t>(tpc::padParameter[l][tpc::kNumOfPad]);
    for(Int_t j = 0; j < nPad; ++j) {
      Double_t theta1 = j * dTheta + sTheta;
      Double_t theta2 = (j + 1) * dTheta + sTheta;
      Double_t r_in = cRad - (pLength / 2.);
      Double_t r_out = cRad + (pLength / 2.);
      
      // X = radius * sin(theta), Z = radius * cos(theta) + ZTarget
      // TH2Poly X = actual Z (horizontal axis)
      X[1] = r_out * TMath::Cos(theta1) + tpc::ZTarget;
      X[2] = r_out * TMath::Cos(theta2) + tpc::ZTarget;
      X[3] = r_in * TMath::Cos(theta2) + tpc::ZTarget;
      X[4] = r_in * TMath::Cos(theta1) + tpc::ZTarget;
      X[0] = X[4];
      
      // TH2Poly Y = actual X (vertical axis)
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

//______________________________________________________________________________
void
set_path(const char* path)
{
  if(gMacroFile) {
    gMacroFile->Close();
    delete gMacroFile;
  }
  if(gMacroReader) {
    delete gMacroReader;
  }
  
  gMacroFile = TFile::Open(path);
  if(!gMacroFile || gMacroFile->IsZombie()) {
    std::cerr << "Error: Cannot open file " << path << std::endl;
    return;
  }
  
  TTree* tree = dynamic_cast<TTree*>(gMacroFile->Get("tpc"));
  if(!tree) {
    std::cerr << "Error: Cannot find tree 'tpc' in file" << std::endl;
    return;
  }
  
  gMacroReader = new TTreeReader(tree);
  
  std::cout << "File opened: " << path << std::endl;
  std::cout << "Total entries: " << tree->GetEntries() << std::endl;
  
  if(!gMacroRandom) gMacroRandom = new TRandom3();
}

//______________________________________________________________________________
void
load_event(Long64_t entry = -1)
{
  if(!gMacroReader) {
    std::cerr << "Error: No file opened. Use set_path() first." << std::endl;
    return;
  }
  
  Long64_t nentries = gMacroReader->GetEntries(false);
  if(nentries == 0) {
    std::cerr << "Error: No entries in tree" << std::endl;
    return;
  }
  
  if(entry < 0) {
    if(!gMacroRandom) gMacroRandom = new TRandom3();
    entry = gMacroRandom->Integer(nentries);
  }
  
  if(entry >= nentries) {
    std::cerr << "Error: Entry " << entry << " is out of range [0, " << nentries << ")" << std::endl;
    return;
  }
  
  // Restart the reader to allow creating new TTreeReaderValue objects
  gMacroReader->Restart();
  
  // Create TTreeReaderValue objects BEFORE calling SetEntry()
  TTreeReaderValue<UInt_t> runnum(*gMacroReader, "run_number");
  TTreeReaderValue<UInt_t> evnum(*gMacroReader, "event_number");
  TTreeReaderValue<Int_t> nhTpc(*gMacroReader, "nhTpc");
  TTreeReaderValue<Int_t> nclTpc(*gMacroReader, "nclTpc");
  TTreeReaderValue<Int_t> ntTpc(*gMacroReader, "ntTpc");
  
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_x(*gMacroReader, "raw_hitpos_x");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_y(*gMacroReader, "raw_hitpos_y");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_z(*gMacroReader, "raw_hitpos_z");
  TTreeReaderValue<std::vector<Double_t>> raw_de(*gMacroReader, "raw_de");
  TTreeReaderValue<std::vector<Int_t>> raw_padid(*gMacroReader, "raw_padid");
  TTreeReaderValue<std::vector<Int_t>> raw_layer(*gMacroReader, "raw_layer");
  TTreeReaderValue<std::vector<Int_t>> raw_row(*gMacroReader, "raw_row");
  
  TTreeReaderValue<std::vector<Double_t>> cluster_x(*gMacroReader, "cluster_x");
  TTreeReaderValue<std::vector<Double_t>> cluster_y(*gMacroReader, "cluster_y");
  TTreeReaderValue<std::vector<Double_t>> cluster_z(*gMacroReader, "cluster_z");
  TTreeReaderValue<std::vector<Double_t>> cluster_de(*gMacroReader, "cluster_de");
  TTreeReaderValue<std::vector<Int_t>> cluster_layer(*gMacroReader, "cluster_layer");
  TTreeReaderValue<std::vector<Int_t>> cluster_row_center(*gMacroReader, "cluster_row_center");
  TTreeReaderValue<std::vector<Int_t>> cluster_houghflag(*gMacroReader, "cluster_houghflag");
  
  TTreeReaderValue<std::vector<Double_t>> x0Tpc(*gMacroReader, "x0Tpc");
  TTreeReaderValue<std::vector<Double_t>> y0Tpc(*gMacroReader, "y0Tpc");
  TTreeReaderValue<std::vector<Double_t>> u0Tpc(*gMacroReader, "u0Tpc");
  TTreeReaderValue<std::vector<Double_t>> v0Tpc(*gMacroReader, "v0Tpc");
  TTreeReaderValue<std::vector<Double_t>> theta(*gMacroReader, "theta");
  TTreeReaderValue<std::vector<Int_t>> nhtrack(*gMacroReader, "nhtrack");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_x(*gMacroReader, "hitpos_x");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_y(*gMacroReader, "hitpos_y");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_z(*gMacroReader, "hitpos_z");
  
  // Now set the entry after all TTreeReaderValue objects are created
  gMacroReader->SetEntry(entry);
  
  gEvent.runnum = *runnum;
  gEvent.evnum = *evnum;
  gEvent.nhTpc = *nhTpc;
  gEvent.nclTpc = *nclTpc;
  gEvent.ntTpc = *ntTpc;
  gEvent.raw_hitpos_x = *raw_hitpos_x;
  gEvent.raw_hitpos_y = *raw_hitpos_y;
  gEvent.raw_hitpos_z = *raw_hitpos_z;
  gEvent.raw_de = *raw_de;
  gEvent.raw_padid = *raw_padid;
  gEvent.raw_layer = *raw_layer;
  gEvent.raw_row = *raw_row;
  gEvent.cluster_x = *cluster_x;
  gEvent.cluster_y = *cluster_y;
  gEvent.cluster_z = *cluster_z;
  gEvent.cluster_de = *cluster_de;
  gEvent.cluster_layer = *cluster_layer;
  gEvent.cluster_row_center = *cluster_row_center;
  gEvent.cluster_houghflag = *cluster_houghflag;
  gEvent.x0Tpc = *x0Tpc;
  gEvent.y0Tpc = *y0Tpc;
  gEvent.u0Tpc = *u0Tpc;
  gEvent.v0Tpc = *v0Tpc;
  gEvent.theta = *theta;
  gEvent.nhtrack = *nhtrack;
  gEvent.hitpos_x = *hitpos_x;
  gEvent.hitpos_y = *hitpos_y;
  gEvent.hitpos_z = *hitpos_z;
  
  std::cout << "Event loaded: Run=" << gEvent.runnum 
       << ", Event=" << gEvent.evnum 
       << ", Entry=" << entry << std::endl;
  std::cout << "  Raw hits: " << gEvent.nhTpc << std::endl;
  std::cout << "  Clusters: " << gEvent.nclTpc << std::endl;
  std::cout << "  Tracks: " << gEvent.ntTpc << std::endl;
}

//______________________________________________________________________________
void
event(Long64_t evnum = -1)
{
  load_event(evnum);
  
  // Create new canvas each time to avoid overwriting
  if(gMacroCanvas) {
    delete gMacroCanvas;
  }
  gMacroCanvas = new TCanvas(Form("c1_%u_%u", gEvent.runnum, gEvent.evnum), 
                              Form("TPC Track Check (Run %u, Event %u)", gEvent.runnum, gEvent.evnum), 
                              1400, 700);
  gMacroCanvas->Divide(2, 1);
  
  gMacroCanvas->cd(1);
  {
    // X-Z view: Pad geometry with hits colored by de
    // Fixed range for consistent viewing
    // Note: Z is horizontal axis, X is vertical axis
    const Double_t zmin = -300.0, zmax = 300.0;
    const Double_t xmin = -300.0, xmax = 300.0;
    
    TH2Poly* h1 = new TH2Poly("h1", Form("X-Z View (Run %u, Event %u);Z [mm];X [mm]", 
                                          gEvent.runnum, gEvent.evnum),
                               zmin, zmax, xmin, xmax);
    h1->SetStats(0);
    
    // Create pad template (all pads with minimum value)
    TPC_pad_template(h1);
    h1->SetMinimum(0);
    h1->SetMaximum(1000); // Adjust based on typical de values
    
    // Fill pads with hit de values
    std::map<Int_t, Double_t> pad_de_map; // padid -> max de
    for(Int_t i = 0; i < gEvent.nhTpc; i++) {
      if(i < static_cast<Int_t>(gEvent.raw_padid.size()) && gEvent.raw_padid[i] >= 0) {
        Int_t padid = gEvent.raw_padid[i];
        Double_t de = (i < static_cast<Int_t>(gEvent.raw_de.size())) ? gEvent.raw_de[i] : 0.0;
        if(pad_de_map.find(padid) == pad_de_map.end() || pad_de_map[padid] < de) {
          pad_de_map[padid] = de;
        }
      }
    }
    
    // Set pad values based on de
    for(const auto& pair : pad_de_map) {
      Int_t padid = pair.first;
      Double_t de = pair.second;
      // Find bin corresponding to this pad
      // Note: Z is horizontal (first arg), X is vertical (second arg)
      TVector3 padPos = tpc::getPosition(padid);
      Int_t bin = h1->FindBin(padPos.Z(), padPos.X());
      if(bin > 0) {
        h1->SetBinContent(bin, de);
      }
    }
    
    h1->Draw("COLZ");
    
    // Draw clusters on top
    // Note: Z is horizontal (x), X is vertical (y)
    TGraph* gCl = new TGraph();
    gCl->SetMarkerStyle(21);
    gCl->SetMarkerSize(1.0);
    gCl->SetMarkerColor(kBlue);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      gCl->SetPoint(gCl->GetN(), gEvent.cluster_z[i], gEvent.cluster_x[i]);
    }
    gCl->Draw("P same");
    
    // Draw Hough-selected clusters
    TGraph* gClHough = new TGraph();
    gClHough->SetMarkerStyle(22);
    gClHough->SetMarkerSize(1.2);
    gClHough->SetMarkerColor(kRed);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(gEvent.cluster_houghflag[i] > 0) {
        gClHough->SetPoint(gClHough->GetN(), gEvent.cluster_z[i], gEvent.cluster_x[i]);
      }
    }
    gClHough->Draw("P same");
    
    // Draw tracks
    // Note: Z is horizontal (x), X is vertical (y)
    for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
      Double_t x0 = gEvent.x0Tpc[itrack];
      Double_t u0 = gEvent.u0Tpc[itrack];
      Double_t x1 = x0 + u0 * zmin;
      Double_t x2 = x0 + u0 * zmax;
      
      TGraph* gTrack = new TGraph();
      gTrack->SetLineWidth(2);
      gTrack->SetLineColor(itrack % 9 + 1);
      gTrack->SetPoint(0, zmin, x1);
      gTrack->SetPoint(1, zmax, x2);
      gTrack->Draw("L same");
    }
    
    TLegend* leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry((TObject*)0, "Pads: de (color)", "");
    leg1->AddEntry(gCl, "Clusters", "p");
    leg1->AddEntry(gClHough, "Hough selected", "p");
    leg1->Draw();
  }
  
  gMacroCanvas->cd(2);
  {
    // Y-Z view: Pad geometry with hits colored by de
    // Fixed range for consistent viewing
    // Note: Z is horizontal axis, Y is vertical axis
    const Double_t zmin = -300.0, zmax = 300.0;
    const Double_t ymin = -200.0, ymax = 200.0;
    
    TH2Poly* h2 = new TH2Poly("h2", Form("Y-Z View (Run %u, Event %u);Z [mm];Y [mm]", 
                                          gEvent.runnum, gEvent.evnum),
                               zmin, zmax, ymin, ymax);
    h2->SetStats(0);
    
    // Create pad template (all pads with minimum value)
    // For Y-Z view, we need to project pads to Y-Z plane
    // Since pads are in X-Z plane, we'll use the pad's Y position (which is 0 for pad center)
    // and use raw_hitpos_y for actual Y positions
    Double_t X[5];
    Double_t Y[5];
    for(Int_t l = 0; l < tpc::NumOfLayersTPC; ++l) {
      Double_t pLength = tpc::padParameter[l][tpc::kLength];
      Double_t st = (180. - (360. / tpc::padParameter[l][tpc::kNumOfDivision]) *
                     tpc::padParameter[l][tpc::kNumOfPad] / 2.);
      Double_t sTheta = (-1 + st / 180.) * TMath::Pi();
      Double_t dTheta = (360. / tpc::padParameter[l][tpc::kNumOfDivision]) / 180. * TMath::Pi();
      Double_t cRad = tpc::padParameter[l][tpc::kRadius];
      Int_t nPad = static_cast<Int_t>(tpc::padParameter[l][tpc::kNumOfPad]);
      for(Int_t j = 0; j < nPad; ++j) {
        // For Y-Z view: TH2Poly X = actual Z (horizontal), TH2Poly Y = actual Y (vertical)
        // According to tpc::getPosition: Z = radius * cos(theta) + ZTarget, Y = 0
        Double_t theta1 = j * dTheta + sTheta;
        Double_t theta2 = (j + 1) * dTheta + sTheta;
        Double_t r_in = cRad - (pLength / 2.);
        Double_t r_out = cRad + (pLength / 2.);
        
        // Z = radius * cos(theta) + ZTarget
        // TH2Poly X = actual Z (horizontal axis)
        X[1] = r_out * TMath::Cos(theta1) + tpc::ZTarget;
        X[2] = r_out * TMath::Cos(theta2) + tpc::ZTarget;
        X[3] = r_in * TMath::Cos(theta2) + tpc::ZTarget;
        X[4] = r_in * TMath::Cos(theta1) + tpc::ZTarget;
        X[0] = X[4];
        
        // TH2Poly Y = actual Y (vertical axis, use small range for pad thickness)
        Double_t yThickness = 5.0; // Approximate pad thickness in Y direction
        Y[1] = yThickness / 2.0;
        Y[2] = yThickness / 2.0;
        Y[3] = -yThickness / 2.0;
        Y[4] = -yThickness / 2.0;
        Y[0] = Y[4];
        
        h2->AddBin(5, X, Y);
      }
    }
    h2->SetMinimum(0);
    h2->SetMaximum(1000);
    
    // Fill pads with hit de values
    std::map<Int_t, Double_t> pad_de_map; // padid -> max de
    for(Int_t i = 0; i < gEvent.nhTpc; i++) {
      if(i < static_cast<Int_t>(gEvent.raw_padid.size()) && gEvent.raw_padid[i] >= 0) {
        Int_t padid = gEvent.raw_padid[i];
        Double_t de = (i < static_cast<Int_t>(gEvent.raw_de.size())) ? gEvent.raw_de[i] : 0.0;
        if(pad_de_map.find(padid) == pad_de_map.end() || pad_de_map[padid] < de) {
          pad_de_map[padid] = de;
        }
      }
    }
    
    // Set pad values based on de
    for(const auto& pair : pad_de_map) {
      Int_t padid = pair.first;
      Double_t de = pair.second;
      TVector3 padPos = tpc::getPosition(padid);
      // For Y-Z view, Z is horizontal (first arg), Y is vertical (second arg)
      Int_t bin = h2->FindBin(padPos.Z(), 0.0);
      if(bin > 0) {
        h2->SetBinContent(bin, de);
      }
    }
    
    h2->Draw("COLZ");
    
    // Draw clusters on top
    // Note: Z is horizontal (x), Y is vertical (y)
    TGraph* gCl = new TGraph();
    gCl->SetMarkerStyle(21);
    gCl->SetMarkerSize(1.0);
    gCl->SetMarkerColor(kBlue);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      gCl->SetPoint(gCl->GetN(), gEvent.cluster_z[i], gEvent.cluster_y[i]);
    }
    gCl->Draw("P same");
    
    // Draw Hough-selected clusters
    TGraph* gClHough = new TGraph();
    gClHough->SetMarkerStyle(22);
    gClHough->SetMarkerSize(1.2);
    gClHough->SetMarkerColor(kRed);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(gEvent.cluster_houghflag[i] > 0) {
        gClHough->SetPoint(gClHough->GetN(), gEvent.cluster_z[i], gEvent.cluster_y[i]);
      }
    }
    gClHough->Draw("P same");
    
    // Draw tracks
    // Note: Z is horizontal (x), Y is vertical (y)
    for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
      Double_t y0 = gEvent.y0Tpc[itrack];
      Double_t v0 = gEvent.v0Tpc[itrack];
      Double_t y1 = y0 + v0 * zmin;
      Double_t y2 = y0 + v0 * zmax;
      
      TGraph* gTrack = new TGraph();
      gTrack->SetLineWidth(2);
      gTrack->SetLineColor(itrack % 9 + 1);
      gTrack->SetPoint(0, zmin, y1);
      gTrack->SetPoint(1, zmax, y2);
      gTrack->Draw("L same");
    }
  }
  
  gMacroCanvas->Update();
  
  // Print summary
  std::cout << "\n=== Event Summary ===" << std::endl;
  std::cout << "Raw hits: " << gEvent.nhTpc << std::endl;
  std::cout << "Clusters: " << gEvent.nclTpc << std::endl;
  Int_t nHoughClusters = 0;
  for(Int_t i = 0; i < gEvent.nclTpc; i++) {
    if(gEvent.cluster_houghflag[i] > 0) nHoughClusters++;
  }
  std::cout << "Hough-selected clusters: " << nHoughClusters << std::endl;
  std::cout << "Tracks: " << gEvent.ntTpc << std::endl;
  for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
    std::cout << "  Track " << itrack << ": "
         << "x0=" << gEvent.x0Tpc[itrack] << ", y0=" << gEvent.y0Tpc[itrack]
         << ", u0=" << gEvent.u0Tpc[itrack] << ", v0=" << gEvent.v0Tpc[itrack]
         << ", theta=" << gEvent.theta[itrack] << ", nhits=" << gEvent.nhtrack[itrack] << std::endl;
  }
  std::cout << "===================" << std::endl;
}
