// -*- C++ -*-
// Macro to check TPC tracking results in 3D: interactive 3D visualization
// Usage:
//   root
//   .L check_tpc_track_3d.C
//   set_path("path/to/rootfile.root")
//   event()        // random event
//   event(evnum)   // specific event
//
// Controls:
//   - Left mouse button + drag: Rotate view
//   - Right mouse button + drag: Zoom
//   - Middle mouse button + drag: Pan
//   - Mouse wheel: Zoom in/out

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TMarker3DBox.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TText.h>
#include <vector>
#include <iostream>
#include "../include/TPCPadHelper.hh"

// Global variables
TFile* gMacroFile = nullptr;
TTreeReader* gMacroReader = nullptr;
TRandom3* gMacroRandom = nullptr;
TCanvas* gMacroCanvas = nullptr;
TView* gMacroView = nullptr;

// Event data structure (same as check_tpc_track.C)
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

// Storage for 3D objects
std::vector<TPolyMarker3D*> gRawHits;
std::vector<TPolyMarker3D*> gClusters;
std::vector<TPolyMarker3D*> gHoughClusters;
std::vector<TPolyLine3D*> gTracks;
std::vector<TPolyLine3D*> gTPCFrame;

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
clear_objects()
{
  // Delete existing 3D objects
  for(auto* obj : gRawHits) {
    if(obj) delete obj;
  }
  gRawHits.clear();
  
  for(auto* obj : gClusters) {
    if(obj) delete obj;
  }
  gClusters.clear();
  
  for(auto* obj : gHoughClusters) {
    if(obj) delete obj;
  }
  gHoughClusters.clear();
  
  for(auto* obj : gTracks) {
    if(obj) delete obj;
  }
  gTracks.clear();
  
  for(auto* obj : gTPCFrame) {
    if(obj) delete obj;
  }
  gTPCFrame.clear();
}

//______________________________________________________________________________
void
event(Long64_t evnum = -1)
{
  load_event(evnum);
  
  // Clear previous objects
  clear_objects();
  
  // Create or reuse canvas
  if(!gMacroCanvas) {
    gMacroCanvas = new TCanvas("c3d", "TPC Track 3D View", 1200, 800);
  } else {
    gMacroCanvas->Clear();
  }
  
  // Calculate range for view
  // Note: Coordinate transformation (x,y,z) -> (z,x,y)
  // So we need to find range in original coordinates and transform
  Double_t xmin_orig = -300.0, xmax_orig = 300.0;
  Double_t ymin_orig = -200.0, ymax_orig = 200.0;
  Double_t zmin_orig = -300.0, zmax_orig = 300.0;
  
  // Find actual range from data (in original coordinates)
  if(gEvent.nhTpc > 0) {
    for(Int_t i = 0; i < gEvent.nhTpc; i++) {
      if(i < static_cast<Int_t>(gEvent.raw_hitpos_x.size())) {
        if(gEvent.raw_hitpos_x[i] < xmin_orig) xmin_orig = gEvent.raw_hitpos_x[i];
        if(gEvent.raw_hitpos_x[i] > xmax_orig) xmax_orig = gEvent.raw_hitpos_x[i];
      }
      if(i < static_cast<Int_t>(gEvent.raw_hitpos_y.size())) {
        if(gEvent.raw_hitpos_y[i] < ymin_orig) ymin_orig = gEvent.raw_hitpos_y[i];
        if(gEvent.raw_hitpos_y[i] > ymax_orig) ymax_orig = gEvent.raw_hitpos_y[i];
      }
      if(i < static_cast<Int_t>(gEvent.raw_hitpos_z.size())) {
        if(gEvent.raw_hitpos_z[i] < zmin_orig) zmin_orig = gEvent.raw_hitpos_z[i];
        if(gEvent.raw_hitpos_z[i] > zmax_orig) zmax_orig = gEvent.raw_hitpos_z[i];
      }
    }
  }
  
  // Add some margin
  Double_t xmargin = (xmax_orig - xmin_orig) * 0.1;
  Double_t ymargin = (ymax_orig - ymin_orig) * 0.1;
  Double_t zmargin = (zmax_orig - zmin_orig) * 0.1;
  xmin_orig -= xmargin; xmax_orig += xmargin;
  ymin_orig -= ymargin; ymax_orig += ymargin;
  zmin_orig -= zmargin; zmax_orig += zmargin;
  
  // Transform range: (x,y,z) -> (z,x,y)
  Double_t xmin = zmin_orig, xmax = zmax_orig; // Display X = original Z
  Double_t ymin = xmin_orig, ymax = xmax_orig; // Display Y = original X
  Double_t zmin = ymin_orig, zmax = ymax_orig; // Display Z = original Y
  
  // Create 3D view
  gMacroView = TView::CreateView(1, 0, 0);
  gMacroView->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);
  gMacroView->ShowAxis();
  
  // Draw TPC frame (cylindrical layers)
  // Draw circles for each layer at different Z positions
  const Int_t nCirclePoints = 64;
  for(Int_t layer = 0; layer < tpc::NumOfLayersTPC; layer++) {
    Double_t radius = tpc::padParameter[layer][tpc::kRadius];
    Double_t padLength = tpc::padParameter[layer][tpc::kLength];
    Double_t r_in = radius - padLength / 2.0;
    Double_t r_out = radius + padLength / 2.0;
    
    // Draw inner and outer circles at different Y positions (original X)
    // Coordinate transform: (x,y,z) -> (z,x,y)
    // For frame, we draw circles in X-Y plane at different Z positions
    // But with coordinate transform, we need to draw in transformed space
    
    // Inner circle
    TPolyLine3D* circle_in = new TPolyLine3D(nCirclePoints + 1);
    circle_in->SetLineColor(kGray + 1);
    circle_in->SetLineWidth(1);
    circle_in->SetLineStyle(2);
    for(Int_t i = 0; i <= nCirclePoints; i++) {
      Double_t theta = 2.0 * TMath::Pi() * i / nCirclePoints;
      Double_t x_orig = r_in * TMath::Cos(theta);
      Double_t y_orig = r_in * TMath::Sin(theta);
      Double_t z_orig = 0.0; // Center of TPC in Z
      // Transform: (x,y,z) -> (z,x,y)
      circle_in->SetPoint(i, z_orig, x_orig, y_orig);
    }
    circle_in->Draw();
    gTPCFrame.push_back(circle_in);
    
    // Outer circle
    TPolyLine3D* circle_out = new TPolyLine3D(nCirclePoints + 1);
    circle_out->SetLineColor(kGray + 1);
    circle_out->SetLineWidth(1);
    circle_out->SetLineStyle(1);
    for(Int_t i = 0; i <= nCirclePoints; i++) {
      Double_t theta = 2.0 * TMath::Pi() * i / nCirclePoints;
      Double_t x_orig = r_out * TMath::Cos(theta);
      Double_t y_orig = r_out * TMath::Sin(theta);
      Double_t z_orig = 0.0; // Center of TPC in Z
      // Transform: (x,y,z) -> (z,x,y)
      circle_out->SetPoint(i, z_orig, x_orig, y_orig);
    }
    circle_out->Draw();
    gTPCFrame.push_back(circle_out);
  }
  
  // Draw raw hits
  if(gEvent.nhTpc > 0) {
    TPolyMarker3D* rawHits = new TPolyMarker3D(gEvent.nhTpc);
    rawHits->SetMarkerStyle(20);
    rawHits->SetMarkerSize(0.5);
    rawHits->SetMarkerColor(kGray);
    
    Int_t npoints = 0;
    for(Int_t i = 0; i < gEvent.nhTpc; i++) {
      if(i < static_cast<Int_t>(gEvent.raw_hitpos_x.size()) &&
         i < static_cast<Int_t>(gEvent.raw_hitpos_y.size()) &&
         i < static_cast<Int_t>(gEvent.raw_hitpos_z.size())) {
        // Transform: (x,y,z) -> (z,x,y)
        rawHits->SetPoint(npoints++,
                         gEvent.raw_hitpos_z[i],  // Display X = original Z
                         gEvent.raw_hitpos_x[i],  // Display Y = original X
                         gEvent.raw_hitpos_y[i]); // Display Z = original Y
      }
    }
    if(npoints > 0) {
      rawHits->Draw();
      gRawHits.push_back(rawHits);
    } else {
      delete rawHits;
    }
  }
  
  // Draw clusters
  if(gEvent.nclTpc > 0) {
    TPolyMarker3D* clusters = new TPolyMarker3D(gEvent.nclTpc);
    clusters->SetMarkerStyle(21);
    clusters->SetMarkerSize(1.0);
    clusters->SetMarkerColor(kBlue);
    
    Int_t npoints = 0;
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(i < static_cast<Int_t>(gEvent.cluster_x.size()) &&
         i < static_cast<Int_t>(gEvent.cluster_y.size()) &&
         i < static_cast<Int_t>(gEvent.cluster_z.size())) {
        // Transform: (x,y,z) -> (z,x,y)
        clusters->SetPoint(npoints++,
                          gEvent.cluster_z[i],  // Display X = original Z
                          gEvent.cluster_x[i],  // Display Y = original X
                          gEvent.cluster_y[i]); // Display Z = original Y
      }
    }
    if(npoints > 0) {
      clusters->Draw();
      gClusters.push_back(clusters);
    } else {
      delete clusters;
    }
  }
  
  // Draw Hough-selected clusters
  if(gEvent.nclTpc > 0) {
    TPolyMarker3D* houghClusters = new TPolyMarker3D(gEvent.nclTpc);
    houghClusters->SetMarkerStyle(22);
    houghClusters->SetMarkerSize(1.2);
    houghClusters->SetMarkerColor(kRed);
    
    Int_t npoints = 0;
    for(Int_t i = 0; i < gEvent.nclTpc; i++) {
      if(i < static_cast<Int_t>(gEvent.cluster_houghflag.size()) &&
         gEvent.cluster_houghflag[i] > 0) {
        if(i < static_cast<Int_t>(gEvent.cluster_x.size()) &&
           i < static_cast<Int_t>(gEvent.cluster_y.size()) &&
           i < static_cast<Int_t>(gEvent.cluster_z.size())) {
          // Transform: (x,y,z) -> (z,x,y)
          houghClusters->SetPoint(npoints++,
                                 gEvent.cluster_z[i],  // Display X = original Z
                                 gEvent.cluster_x[i],  // Display Y = original X
                                 gEvent.cluster_y[i]); // Display Z = original Y
        }
      }
    }
    if(npoints > 0) {
      houghClusters->Draw();
      gHoughClusters.push_back(houghClusters);
    } else {
      delete houghClusters;
    }
  }
  
  // Draw tracks
  if(gEvent.ntTpc > 0) {
    for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
      if(itrack < static_cast<Int_t>(gEvent.x0Tpc.size()) &&
         itrack < static_cast<Int_t>(gEvent.y0Tpc.size()) &&
         itrack < static_cast<Int_t>(gEvent.u0Tpc.size()) &&
         itrack < static_cast<Int_t>(gEvent.v0Tpc.size())) {
        
        Double_t x0 = gEvent.x0Tpc[itrack];
        Double_t y0 = gEvent.y0Tpc[itrack];
        Double_t u0 = gEvent.u0Tpc[itrack];
        Double_t v0 = gEvent.v0Tpc[itrack];
        
        // Calculate track points at z boundaries (original coordinates)
        Double_t z1_orig = zmin_orig;
        Double_t z2_orig = zmax_orig;
        Double_t x1_orig = x0 + u0 * z1_orig;
        Double_t y1_orig = y0 + v0 * z1_orig;
        Double_t x2_orig = x0 + u0 * z2_orig;
        Double_t y2_orig = y0 + v0 * z2_orig;
        
        TPolyLine3D* track = new TPolyLine3D(2);
        track->SetLineWidth(2);
        track->SetLineColor(itrack % 9 + 1);
        // Transform: (x,y,z) -> (z,x,y)
        track->SetPoint(0, z1_orig, x1_orig, y1_orig);
        track->SetPoint(1, z2_orig, x2_orig, y2_orig);
        track->Draw();
        gTracks.push_back(track);
      }
    }
  }
  
  // Draw track hits if available
  if(gEvent.ntTpc > 0) {
    for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
      if(itrack < static_cast<Int_t>(gEvent.nhtrack.size()) &&
         gEvent.nhtrack[itrack] > 0) {
        Int_t nhits = gEvent.nhtrack[itrack];
        if(nhits > 0 &&
           itrack < static_cast<Int_t>(gEvent.hitpos_x.size()) &&
           nhits <= static_cast<Int_t>(gEvent.hitpos_x[itrack].size())) {
          TPolyMarker3D* trackHits = new TPolyMarker3D(nhits);
          trackHits->SetMarkerStyle(24);
          trackHits->SetMarkerSize(0.8);
          trackHits->SetMarkerColor(itrack % 9 + 1);
          
          for(Int_t i = 0; i < nhits; i++) {
            // Transform: (x,y,z) -> (z,x,y)
            trackHits->SetPoint(i,
                               gEvent.hitpos_z[itrack][i],  // Display X = original Z
                               gEvent.hitpos_x[itrack][i],  // Display Y = original X
                               gEvent.hitpos_y[itrack][i]); // Display Z = original Y
          }
          trackHits->Draw();
          gRawHits.push_back(trackHits);
        }
      }
    }
  }
  
  gMacroCanvas->Update();
  
  // Print summary
  std::cout << "\n=== Event Summary (3D) ===" << std::endl;
  std::cout << "Raw hits: " << gEvent.nhTpc << std::endl;
  std::cout << "Clusters: " << gEvent.nclTpc << std::endl;
  Int_t nHoughClusters = 0;
  for(Int_t i = 0; i < gEvent.nclTpc; i++) {
    if(i < static_cast<Int_t>(gEvent.cluster_houghflag.size()) &&
       gEvent.cluster_houghflag[i] > 0) nHoughClusters++;
  }
  std::cout << "Hough-selected clusters: " << nHoughClusters << std::endl;
  std::cout << "Tracks: " << gEvent.ntTpc << std::endl;
  for(Int_t itrack = 0; itrack < gEvent.ntTpc; itrack++) {
    if(itrack < static_cast<Int_t>(gEvent.x0Tpc.size())) {
      std::cout << "  Track " << itrack << ": "
           << "x0=" << gEvent.x0Tpc[itrack] << ", y0=" << gEvent.y0Tpc[itrack]
           << ", u0=" << gEvent.u0Tpc[itrack] << ", v0=" << gEvent.v0Tpc[itrack]
           << ", theta=" << (itrack < static_cast<Int_t>(gEvent.theta.size()) ? gEvent.theta[itrack] : 0.0)
           << ", nhits=" << (itrack < static_cast<Int_t>(gEvent.nhtrack.size()) ? gEvent.nhtrack[itrack] : 0) << std::endl;
    }
  }
  std::cout << "=========================" << std::endl;
  std::cout << "\n3D View Controls:" << std::endl;
  std::cout << "  - Left mouse + drag: Rotate" << std::endl;
  std::cout << "  - Right mouse + drag: Zoom" << std::endl;
  std::cout << "  - Middle mouse + drag: Pan" << std::endl;
  std::cout << "  - Mouse wheel: Zoom in/out" << std::endl;
}
