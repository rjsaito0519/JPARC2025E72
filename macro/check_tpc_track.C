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
#include <vector>
#include <iostream>

#include "TPCPadHelper.hh"

// Global variables
TFile* gFile = nullptr;
TTreeReader* gReader = nullptr;
TRandom3* gRandom = nullptr;
TCanvas* gCanvas = nullptr;

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
      X[1] = (cRad + (pLength / 2.)) * TMath::Cos(j * dTheta + sTheta);
      X[2] = (cRad + (pLength / 2.)) * TMath::Cos((j + 1) * dTheta + sTheta);
      X[3] = (cRad - (pLength / 2.)) * TMath::Cos((j + 1) * dTheta + sTheta);
      X[4] = (cRad - (pLength / 2.)) * TMath::Cos(j * dTheta + sTheta);
      X[0] = X[4];
      Y[1] = (cRad + (pLength / 2.)) * TMath::Sin(j * dTheta + sTheta);
      Y[2] = (cRad + (pLength / 2.)) * TMath::Sin((j + 1) * dTheta + sTheta);
      Y[3] = (cRad - (pLength / 2.)) * TMath::Sin((j + 1) * dTheta + sTheta);
      Y[4] = (cRad - (pLength / 2.)) * TMath::Sin(j * dTheta + sTheta);
      Y[0] = Y[4];
      for(Int_t k = 0; k < 5; ++k) X[k] += tpc::ZTarget;
      h->AddBin(5, X, Y);
    }
  }
  h->SetMaximum(0x1000);
}

//______________________________________________________________________________
void
set_path(const char* path)
{
  if(gFile) {
    gFile->Close();
    delete gFile;
  }
  if(gReader) {
    delete gReader;
  }
  
  gFile = TFile::Open(path);
  if(!gFile || gFile->IsZombie()) {
    std::cerr << "Error: Cannot open file " << path << std::endl;
    return;
  }
  
  TTree* tree = dynamic_cast<TTree*>(gFile->Get("tpc"));
  if(!tree) {
    std::cerr << "Error: Cannot find tree 'tpc' in file" << std::endl;
    return;
  }
  
  gReader = new TTreeReader(tree);
  
  // Set up branches
  TTreeReaderValue<UInt_t> runnum(*gReader, "run_number");
  TTreeReaderValue<UInt_t> evnum(*gReader, "event_number");
  TTreeReaderValue<Int_t> nhTpc(*gReader, "nhTpc");
  TTreeReaderValue<Int_t> nclTpc(*gReader, "nclTpc");
  TTreeReaderValue<Int_t> ntTpc(*gReader, "ntTpc");
  
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_x(*gReader, "raw_hitpos_x");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_y(*gReader, "raw_hitpos_y");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_z(*gReader, "raw_hitpos_z");
  TTreeReaderValue<std::vector<Int_t>> raw_padid(*gReader, "raw_padid");
  TTreeReaderValue<std::vector<Int_t>> raw_layer(*gReader, "raw_layer");
  TTreeReaderValue<std::vector<Int_t>> raw_row(*gReader, "raw_row");
  
  TTreeReaderValue<std::vector<Double_t>> cluster_x(*gReader, "cluster_x");
  TTreeReaderValue<std::vector<Double_t>> cluster_y(*gReader, "cluster_y");
  TTreeReaderValue<std::vector<Double_t>> cluster_z(*gReader, "cluster_z");
  TTreeReaderValue<std::vector<Double_t>> cluster_de(*gReader, "cluster_de");
  TTreeReaderValue<std::vector<Int_t>> cluster_layer(*gReader, "cluster_layer");
  TTreeReaderValue<std::vector<Int_t>> cluster_row_center(*gReader, "cluster_row_center");
  TTreeReaderValue<std::vector<Int_t>> cluster_houghflag(*gReader, "cluster_houghflag");
  
  TTreeReaderValue<std::vector<Double_t>> x0Tpc(*gReader, "x0Tpc");
  TTreeReaderValue<std::vector<Double_t>> y0Tpc(*gReader, "y0Tpc");
  TTreeReaderValue<std::vector<Double_t>> u0Tpc(*gReader, "u0Tpc");
  TTreeReaderValue<std::vector<Double_t>> v0Tpc(*gReader, "v0Tpc");
  TTreeReaderValue<std::vector<Double_t>> theta(*gReader, "theta");
  TTreeReaderValue<std::vector<Int_t>> nhtrack(*gReader, "nhtrack");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_x(*gReader, "hitpos_x");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_y(*gReader, "hitpos_y");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_z(*gReader, "hitpos_z");
  
  // Read first entry to get tree size
  if(gReader->Next()) {
    gEvent.runnum = *runnum;
    gEvent.evnum = *evnum;
    gEvent.nhTpc = *nhTpc;
    gEvent.nclTpc = *nclTpc;
    gEvent.ntTpc = *ntTpc;
    gEvent.raw_hitpos_x = *raw_hitpos_x;
    gEvent.raw_hitpos_y = *raw_hitpos_y;
    gEvent.raw_hitpos_z = *raw_hitpos_z;
    gEvent.raw_padid = *raw_padid;
    gEvent.raw_layer = *raw_layer;
    gEvent.raw_row = *raw_row;
    gEvent.cluster_x = *cluster_x;
    gEvent.cluster_y = *cluster_y;
    gEvent.cluster_z = *cluster_z;
    gEvent.cluster_de = *cluster_de;
    gEvent.cluster_layer = *cluster_layer;
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
  }
  
  gReader->SetEntry(0);
  
  std::cout << "File opened: " << path << std::endl;
  std::cout << "Total entries: " << tree->GetEntries() << std::endl;
  
  if(!gRandom) gRandom = new TRandom3();
}

//______________________________________________________________________________
void
load_event(Long64_t entry = -1)
{
  if(!gReader) {
    std::cerr << "Error: No file opened. Use set_path() first." << std::endl;
    return;
  }
  
  Long64_t nentries = gReader->GetEntries(false);
  if(nentries == 0) {
    std::cerr << "Error: No entries in tree" << std::endl;
    return;
  }
  
  if(entry < 0) {
    if(!gRandom) gRandom = new TRandom3();
    entry = gRandom->Integer(nentries);
  }
  
  if(entry >= nentries) {
    std::cerr << "Error: Entry " << entry << " is out of range [0, " << nentries << ")" << std::endl;
    return;
  }
  
  gReader->SetEntry(entry);
  
  // Read event data
  TTreeReaderValue<UInt_t> runnum(*gReader, "run_number");
  TTreeReaderValue<UInt_t> evnum(*gReader, "event_number");
  TTreeReaderValue<Int_t> nhTpc(*gReader, "nhTpc");
  TTreeReaderValue<Int_t> nclTpc(*gReader, "nclTpc");
  TTreeReaderValue<Int_t> ntTpc(*gReader, "ntTpc");
  
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_x(*gReader, "raw_hitpos_x");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_y(*gReader, "raw_hitpos_y");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_z(*gReader, "raw_hitpos_z");
  TTreeReaderValue<std::vector<Int_t>> raw_padid(*gReader, "raw_padid");
  TTreeReaderValue<std::vector<Int_t>> raw_layer(*gReader, "raw_layer");
  TTreeReaderValue<std::vector<Int_t>> raw_row(*gReader, "raw_row");
  
  TTreeReaderValue<std::vector<Double_t>> cluster_x(*gReader, "cluster_x");
  TTreeReaderValue<std::vector<Double_t>> cluster_y(*gReader, "cluster_y");
  TTreeReaderValue<std::vector<Double_t>> cluster_z(*gReader, "cluster_z");
  TTreeReaderValue<std::vector<Double_t>> cluster_de(*gReader, "cluster_de");
  TTreeReaderValue<std::vector<Int_t>> cluster_layer(*gReader, "cluster_layer");
  TTreeReaderValue<std::vector<Int_t>> cluster_row_center(*gReader, "cluster_row_center");
  TTreeReaderValue<std::vector<Int_t>> cluster_houghflag(*gReader, "cluster_houghflag");
  
  TTreeReaderValue<std::vector<Double_t>> x0Tpc(*gReader, "x0Tpc");
  TTreeReaderValue<std::vector<Double_t>> y0Tpc(*gReader, "y0Tpc");
  TTreeReaderValue<std::vector<Double_t>> u0Tpc(*gReader, "u0Tpc");
  TTreeReaderValue<std::vector<Double_t>> v0Tpc(*gReader, "v0Tpc");
  TTreeReaderValue<std::vector<Double_t>> theta(*gReader, "theta");
  TTreeReaderValue<std::vector<Int_t>> nhtrack(*gReader, "nhtrack");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_x(*gReader, "hitpos_x");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_y(*gReader, "hitpos_y");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_z(*gReader, "hitpos_z");
  
  gEvent.runnum = *runnum;
  gEvent.evnum = *evnum;
  gEvent.nhTpc = *nhTpc;
  gEvent.nclTpc = *nclTpc;
  gEvent.ntTpc = *ntTpc;
  gEvent.raw_hitpos_x = *raw_hitpos_x;
  gEvent.raw_hitpos_y = *raw_hitpos_y;
  gEvent.raw_hitpos_z = *raw_hitpos_z;
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
  
  if(!gCanvas) {
    gCanvas = new TCanvas("c1", "TPC Track Check", 1400, 1000);
    gCanvas->Divide(2, 2);
  }
  
  gCanvas->cd(1);
  {
    // X-Z view: Raw hits and clusters
    TH2D* h1 = new TH2D("h1", Form("X-Z View (Run %u, Event %u);X [mm];Z [mm]", 
                                    gEvent.runnum, gEvent.evnum),
                         100, -200, 200, 100, -200, 200);
    h1->SetStats(0);
    h1->Draw();
    
    // Draw raw hits (using TPCPadHelper to get pad position if padid is available)
    TGraph* gRaw = new TGraph();
    gRaw->SetMarkerStyle(20);
    gRaw->SetMarkerSize(0.5);
    gRaw->SetMarkerColor(kGray);
    for(Int_t i = 0; i < gEvent.nhTpc; i++) {
      Double_t x = gEvent.raw_hitpos_x[i];
      Double_t z = gEvent.raw_hitpos_z[i];
      // If padid is available, use TPCPadHelper to get more accurate position
      if(i < static_cast<Int_t>(gEvent.raw_padid.size()) && gEvent.raw_padid[i] >= 0) {
        TVector3 padPos = tpc::getPosition(gEvent.raw_padid[i]);
        if(!TMath::IsNaN(padPos.X())) {
          x = padPos.X();
          z = padPos.Z();
        }
      }
      gRaw->SetPoint(gRaw->GetN(), x, z);
    }
    gRaw->Draw("P same");
    
    // Draw clusters
    TGraph* gCl = new TGraph();
    gCl->SetMarkerStyle(21);
    gCl->SetMarkerSize(1.0);
    gCl->SetMarkerColor(kBlue);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      gCl->SetPoint(gCl->GetN(), gEvent.cluster_x[i], gEvent.cluster_z[i]);
    }
    gCl->Draw("P same");
    
    // Draw Hough-selected clusters
    TGraph* gClHough = new TGraph();
    gClHough->SetMarkerStyle(22);
    gClHough->SetMarkerSize(1.2);
    gClHough->SetMarkerColor(kRed);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(gEvent.cluster_houghflag[i] > 0) {
        gClHough->SetPoint(gClHough->GetN(), gEvent.cluster_x[i], gEvent.cluster_z[i]);
      }
    }
    gClHough->Draw("P same");
    
    // Draw tracks
    for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
      Double_t x0 = gEvent.x0Tpc[itrack];
      Double_t u0 = gEvent.u0Tpc[itrack];
      Double_t zmin = -200, zmax = 200;
      Double_t x1 = x0 + u0 * zmin;
      Double_t x2 = x0 + u0 * zmax;
      
      TGraph* gTrack = new TGraph();
      gTrack->SetLineWidth(2);
      gTrack->SetLineColor(itrack % 9 + 1);
      gTrack->SetPoint(0, x1, zmin);
      gTrack->SetPoint(1, x2, zmax);
      gTrack->Draw("L same");
    }
    
    TLegend* leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(gRaw, "Raw hits", "p");
    leg1->AddEntry(gCl, "Clusters", "p");
    leg1->AddEntry(gClHough, "Hough selected", "p");
    leg1->Draw();
  }
  
  gCanvas->cd(2);
  {
    // Y-Z view: Raw hits and clusters
    TH2D* h2 = new TH2D("h2", Form("Y-Z View (Run %u, Event %u);Y [mm];Z [mm]", 
                                    gEvent.runnum, gEvent.evnum),
                         100, -200, 200, 100, -200, 200);
    h2->SetStats(0);
    h2->Draw();
    
    // Draw raw hits (using TPCPadHelper to get pad position if padid is available)
    TGraph* gRaw = new TGraph();
    gRaw->SetMarkerStyle(20);
    gRaw->SetMarkerSize(0.5);
    gRaw->SetMarkerColor(kGray);
    for(Int_t i = 0; i < gEvent.nhTpc; i++) {
      Double_t y = gEvent.raw_hitpos_y[i];
      Double_t z = gEvent.raw_hitpos_z[i];
      // If padid is available, use TPCPadHelper to get more accurate position
      // Note: getPosition returns (x, 0, z), so y is always 0 for pad center
      // We use raw_hitpos_y for y coordinate
      if(i < static_cast<Int_t>(gEvent.raw_padid.size()) && gEvent.raw_padid[i] >= 0) {
        TVector3 padPos = tpc::getPosition(gEvent.raw_padid[i]);
        if(!TMath::IsNaN(padPos.Z())) {
          z = padPos.Z();
        }
      }
      gRaw->SetPoint(gRaw->GetN(), y, z);
    }
    gRaw->Draw("P same");
    
    // Draw clusters
    TGraph* gCl = new TGraph();
    gCl->SetMarkerStyle(21);
    gCl->SetMarkerSize(1.0);
    gCl->SetMarkerColor(kBlue);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      gCl->SetPoint(gCl->GetN(), gEvent.cluster_y[i], gEvent.cluster_z[i]);
    }
    gCl->Draw("P same");
    
    // Draw Hough-selected clusters
    TGraph* gClHough = new TGraph();
    gClHough->SetMarkerStyle(22);
    gClHough->SetMarkerSize(1.2);
    gClHough->SetMarkerColor(kRed);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(gEvent.cluster_houghflag[i] > 0) {
        gClHough->SetPoint(gClHough->GetN(), gEvent.cluster_y[i], gEvent.cluster_z[i]);
      }
    }
    gClHough->Draw("P same");
    
    // Draw tracks
    for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
      Double_t y0 = gEvent.y0Tpc[itrack];
      Double_t v0 = gEvent.v0Tpc[itrack];
      Double_t zmin = -200, zmax = 200;
      Double_t y1 = y0 + v0 * zmin;
      Double_t y2 = y0 + v0 * zmax;
      
      TGraph* gTrack = new TGraph();
      gTrack->SetLineWidth(2);
      gTrack->SetLineColor(itrack % 9 + 1);
      gTrack->SetPoint(0, y1, zmin);
      gTrack->SetPoint(1, y2, zmax);
      gTrack->Draw("L same");
    }
  }
  
  gCanvas->cd(3);
  {
    // X-Y view at target
    TH2D* h3 = new TH2D("h3", Form("X-Y View at Target (Run %u, Event %u);X [mm];Y [mm]", 
                                    gEvent.runnum, gEvent.evnum),
                         100, -100, 100, 100, -100, 100);
    h3->SetStats(0);
    h3->Draw();
    
    // Draw clusters (projected to z=0)
    TGraph* gCl = new TGraph();
    gCl->SetMarkerStyle(21);
    gCl->SetMarkerSize(1.0);
    gCl->SetMarkerColor(kBlue);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      gCl->SetPoint(gCl->GetN(), gEvent.cluster_x[i], gEvent.cluster_y[i]);
    }
    gCl->Draw("P same");
    
    // Draw Hough-selected clusters
    TGraph* gClHough = new TGraph();
    gClHough->SetMarkerStyle(22);
    gClHough->SetMarkerSize(1.2);
    gClHough->SetMarkerColor(kRed);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(gEvent.cluster_houghflag[i] > 0) {
        gClHough->SetPoint(gClHough->GetN(), gEvent.cluster_x[i], gEvent.cluster_y[i]);
      }
    }
    gClHough->Draw("P same");
    
    // Draw track starting points
    TGraph* gTrackStart = new TGraph();
    gTrackStart->SetMarkerStyle(29);
    gTrackStart->SetMarkerSize(2.0);
    gTrackStart->SetMarkerColor(kGreen);
    for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
      gTrackStart->SetPoint(gTrackStart->GetN(), gEvent.x0Tpc[itrack], gEvent.y0Tpc[itrack]);
    }
    gTrackStart->Draw("P same");
  }
  
  gCanvas->cd(4);
  {
    // Layer vs Row: Clustering visualization
    TH2D* h4 = new TH2D("h4", Form("Layer vs Row (Run %u, Event %u);Row;Layer", 
                                    gEvent.runnum, gEvent.evnum),
                         200, 0, 200, 20, 0, 20);
    h4->SetStats(0);
    h4->Draw();
    
    // Draw raw hits
    for(Int_t i = 0; i < gEvent.nhTpc; i++) {
      h4->Fill(gEvent.raw_row[i], gEvent.raw_layer[i]);
    }
    
    // Draw clusters
    TGraph* gCl = new TGraph();
    gCl->SetMarkerStyle(21);
    gCl->SetMarkerSize(1.5);
    gCl->SetMarkerColor(kBlue);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      gCl->SetPoint(gCl->GetN(), gEvent.cluster_row_center[i], gEvent.cluster_layer[i]);
    }
    gCl->Draw("P same");
    
    // Draw Hough-selected clusters
    TGraph* gClHough = new TGraph();
    gClHough->SetMarkerStyle(22);
    gClHough->SetMarkerSize(1.8);
    gClHough->SetMarkerColor(kRed);
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(gEvent.cluster_houghflag[i] > 0) {
        gClHough->SetPoint(gClHough->GetN(), gEvent.cluster_row_center[i], gEvent.cluster_layer[i]);
      }
    }
    gClHough->Draw("P same");
  }
  
  gCanvas->Update();
  
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
