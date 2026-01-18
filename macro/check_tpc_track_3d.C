// -*- C++ -*-
// Macro to check TPC tracking results in 3D: interactive 3D visualization
// Usage:
//   root
//   .L check_tpc_track_3d.C
//   set_path_3d("path/to/rootfile.root")
//   event_3d()        // random event
//   event_3d(evnum)   // specific event
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
TFile* gMacroFile3D = nullptr;
TTreeReader* gMacroReader3D = nullptr;
TRandom3* gMacroRandom3D = nullptr;
TCanvas* gMacroCanvas3D = nullptr;
TView* gMacroView3D = nullptr;

// Event data structure (same as check_tpc_track.C)
struct EventData3D {
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

EventData3D gEvent3D;

// Storage for 3D objects
std::vector<TPolyMarker3D*> gRawHits3D;
std::vector<TPolyMarker3D*> gClusters3D;
std::vector<TPolyMarker3D*> gHoughClusters3D;
std::vector<TPolyLine3D*> gTracks3D;

//______________________________________________________________________________
void
set_path_3d(const char* path)
{
  if(gMacroFile3D) {
    gMacroFile3D->Close();
    delete gMacroFile3D;
  }
  if(gMacroReader3D) {
    delete gMacroReader3D;
  }
  
  gMacroFile3D = TFile::Open(path);
  if(!gMacroFile3D || gMacroFile3D->IsZombie()) {
    std::cerr << "Error: Cannot open file " << path << std::endl;
    return;
  }
  
  TTree* tree = dynamic_cast<TTree*>(gMacroFile3D->Get("tpc"));
  if(!tree) {
    std::cerr << "Error: Cannot find tree 'tpc' in file" << std::endl;
    return;
  }
  
  gMacroReader3D = new TTreeReader(tree);
  
  std::cout << "File opened: " << path << std::endl;
  std::cout << "Total entries: " << tree->GetEntries() << std::endl;
  
  if(!gMacroRandom3D) gMacroRandom3D = new TRandom3();
}

//______________________________________________________________________________
void
load_event_3d(Long64_t entry = -1)
{
  if(!gMacroReader3D) {
    std::cerr << "Error: No file opened. Use set_path_3d() first." << std::endl;
    return;
  }
  
  Long64_t nentries = gMacroReader3D->GetEntries(false);
  if(nentries == 0) {
    std::cerr << "Error: No entries in tree" << std::endl;
    return;
  }
  
  if(entry < 0) {
    if(!gMacroRandom3D) gMacroRandom3D = new TRandom3();
    entry = gMacroRandom3D->Integer(nentries);
  }
  
  if(entry >= nentries) {
    std::cerr << "Error: Entry " << entry << " is out of range [0, " << nentries << ")" << std::endl;
    return;
  }
  
  // Restart the reader to allow creating new TTreeReaderValue objects
  gMacroReader3D->Restart();
  
  // Create TTreeReaderValue objects BEFORE calling SetEntry()
  TTreeReaderValue<UInt_t> runnum(*gMacroReader3D, "run_number");
  TTreeReaderValue<UInt_t> evnum(*gMacroReader3D, "event_number");
  TTreeReaderValue<Int_t> nhTpc(*gMacroReader3D, "nhTpc");
  TTreeReaderValue<Int_t> nclTpc(*gMacroReader3D, "nclTpc");
  TTreeReaderValue<Int_t> ntTpc(*gMacroReader3D, "ntTpc");
  
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_x(*gMacroReader3D, "raw_hitpos_x");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_y(*gMacroReader3D, "raw_hitpos_y");
  TTreeReaderValue<std::vector<Double_t>> raw_hitpos_z(*gMacroReader3D, "raw_hitpos_z");
  TTreeReaderValue<std::vector<Double_t>> raw_de(*gMacroReader3D, "raw_de");
  TTreeReaderValue<std::vector<Int_t>> raw_padid(*gMacroReader3D, "raw_padid");
  TTreeReaderValue<std::vector<Int_t>> raw_layer(*gMacroReader3D, "raw_layer");
  TTreeReaderValue<std::vector<Int_t>> raw_row(*gMacroReader3D, "raw_row");
  
  TTreeReaderValue<std::vector<Double_t>> cluster_x(*gMacroReader3D, "cluster_x");
  TTreeReaderValue<std::vector<Double_t>> cluster_y(*gMacroReader3D, "cluster_y");
  TTreeReaderValue<std::vector<Double_t>> cluster_z(*gMacroReader3D, "cluster_z");
  TTreeReaderValue<std::vector<Double_t>> cluster_de(*gMacroReader3D, "cluster_de");
  TTreeReaderValue<std::vector<Int_t>> cluster_layer(*gMacroReader3D, "cluster_layer");
  TTreeReaderValue<std::vector<Int_t>> cluster_row_center(*gMacroReader3D, "cluster_row_center");
  TTreeReaderValue<std::vector<Int_t>> cluster_houghflag(*gMacroReader3D, "cluster_houghflag");
  
  TTreeReaderValue<std::vector<Double_t>> x0Tpc(*gMacroReader3D, "x0Tpc");
  TTreeReaderValue<std::vector<Double_t>> y0Tpc(*gMacroReader3D, "y0Tpc");
  TTreeReaderValue<std::vector<Double_t>> u0Tpc(*gMacroReader3D, "u0Tpc");
  TTreeReaderValue<std::vector<Double_t>> v0Tpc(*gMacroReader3D, "v0Tpc");
  TTreeReaderValue<std::vector<Double_t>> theta(*gMacroReader3D, "theta");
  TTreeReaderValue<std::vector<Int_t>> nhtrack(*gMacroReader3D, "nhtrack");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_x(*gMacroReader3D, "hitpos_x");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_y(*gMacroReader3D, "hitpos_y");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> hitpos_z(*gMacroReader3D, "hitpos_z");
  
  // Now set the entry after all TTreeReaderValue objects are created
  gMacroReader3D->SetEntry(entry);
  
  gEvent3D.runnum = *runnum;
  gEvent3D.evnum = *evnum;
  gEvent3D.nhTpc = *nhTpc;
  gEvent3D.nclTpc = *nclTpc;
  gEvent3D.ntTpc = *ntTpc;
  gEvent3D.raw_hitpos_x = *raw_hitpos_x;
  gEvent3D.raw_hitpos_y = *raw_hitpos_y;
  gEvent3D.raw_hitpos_z = *raw_hitpos_z;
  gEvent3D.raw_de = *raw_de;
  gEvent3D.raw_padid = *raw_padid;
  gEvent3D.raw_layer = *raw_layer;
  gEvent3D.raw_row = *raw_row;
  gEvent3D.cluster_x = *cluster_x;
  gEvent3D.cluster_y = *cluster_y;
  gEvent3D.cluster_z = *cluster_z;
  gEvent3D.cluster_de = *cluster_de;
  gEvent3D.cluster_layer = *cluster_layer;
  gEvent3D.cluster_row_center = *cluster_row_center;
  gEvent3D.cluster_houghflag = *cluster_houghflag;
  gEvent3D.x0Tpc = *x0Tpc;
  gEvent3D.y0Tpc = *y0Tpc;
  gEvent3D.u0Tpc = *u0Tpc;
  gEvent3D.v0Tpc = *v0Tpc;
  gEvent3D.theta = *theta;
  gEvent3D.nhtrack = *nhtrack;
  gEvent3D.hitpos_x = *hitpos_x;
  gEvent3D.hitpos_y = *hitpos_y;
  gEvent3D.hitpos_z = *hitpos_z;
  
  std::cout << "Event loaded: Run=" << gEvent3D.runnum 
       << ", Event=" << gEvent3D.evnum 
       << ", Entry=" << entry << std::endl;
  std::cout << "  Raw hits: " << gEvent3D.nhTpc << std::endl;
  std::cout << "  Clusters: " << gEvent3D.nclTpc << std::endl;
  std::cout << "  Tracks: " << gEvent3D.ntTpc << std::endl;
}

//______________________________________________________________________________
void
clear_3d_objects()
{
  // Delete existing 3D objects
  for(auto* obj : gRawHits3D) {
    if(obj) delete obj;
  }
  gRawHits3D.clear();
  
  for(auto* obj : gClusters3D) {
    if(obj) delete obj;
  }
  gClusters3D.clear();
  
  for(auto* obj : gHoughClusters3D) {
    if(obj) delete obj;
  }
  gHoughClusters3D.clear();
  
  for(auto* obj : gTracks3D) {
    if(obj) delete obj;
  }
  gTracks3D.clear();
}

//______________________________________________________________________________
void
event_3d(Long64_t evnum = -1)
{
  load_event_3d(evnum);
  
  // Clear previous objects
  clear_3d_objects();
  
  // Create or reuse canvas
  if(!gMacroCanvas3D) {
    gMacroCanvas3D = new TCanvas("c3d", "TPC Track 3D View", 1200, 800);
  } else {
    gMacroCanvas3D->Clear();
  }
  
  // Calculate range for view
  Double_t xmin = -300.0, xmax = 300.0;
  Double_t ymin = -200.0, ymax = 200.0;
  Double_t zmin = -300.0, zmax = 300.0;
  
  // Find actual range from data
  if(gEvent3D.nhTpc > 0) {
    for(Int_t i = 0; i < gEvent3D.nhTpc; i++) {
      if(i < static_cast<Int_t>(gEvent3D.raw_hitpos_x.size())) {
        if(gEvent3D.raw_hitpos_x[i] < xmin) xmin = gEvent3D.raw_hitpos_x[i];
        if(gEvent3D.raw_hitpos_x[i] > xmax) xmax = gEvent3D.raw_hitpos_x[i];
      }
      if(i < static_cast<Int_t>(gEvent3D.raw_hitpos_y.size())) {
        if(gEvent3D.raw_hitpos_y[i] < ymin) ymin = gEvent3D.raw_hitpos_y[i];
        if(gEvent3D.raw_hitpos_y[i] > ymax) ymax = gEvent3D.raw_hitpos_y[i];
      }
      if(i < static_cast<Int_t>(gEvent3D.raw_hitpos_z.size())) {
        if(gEvent3D.raw_hitpos_z[i] < zmin) zmin = gEvent3D.raw_hitpos_z[i];
        if(gEvent3D.raw_hitpos_z[i] > zmax) zmax = gEvent3D.raw_hitpos_z[i];
      }
    }
  }
  
  // Add some margin
  Double_t xmargin = (xmax - xmin) * 0.1;
  Double_t ymargin = (ymax - ymin) * 0.1;
  Double_t zmargin = (zmax - zmin) * 0.1;
  xmin -= xmargin; xmax += xmargin;
  ymin -= ymargin; ymax += ymargin;
  zmin -= zmargin; zmax += zmargin;
  
  // Create 3D view
  gMacroView3D = TView::CreateView(1, 0, 0);
  gMacroView3D->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);
  gMacroView3D->ShowAxis();
  
  // Draw raw hits
  if(gEvent3D.nhTpc > 0) {
    TPolyMarker3D* rawHits = new TPolyMarker3D(gEvent3D.nhTpc);
    rawHits->SetMarkerStyle(20);
    rawHits->SetMarkerSize(0.5);
    rawHits->SetMarkerColor(kGray);
    
    Int_t npoints = 0;
    for(Int_t i = 0; i < gEvent3D.nhTpc; i++) {
      if(i < static_cast<Int_t>(gEvent3D.raw_hitpos_x.size()) &&
         i < static_cast<Int_t>(gEvent3D.raw_hitpos_y.size()) &&
         i < static_cast<Int_t>(gEvent3D.raw_hitpos_z.size())) {
        rawHits->SetPoint(npoints++,
                         gEvent3D.raw_hitpos_x[i],
                         gEvent3D.raw_hitpos_y[i],
                         gEvent3D.raw_hitpos_z[i]);
      }
    }
    if(npoints > 0) {
      rawHits->SetNextPoint(npoints);
      rawHits->Draw();
      gRawHits3D.push_back(rawHits);
    } else {
      delete rawHits;
    }
  }
  
  // Draw clusters
  if(gEvent3D.nclTpc > 0) {
    TPolyMarker3D* clusters = new TPolyMarker3D(gEvent3D.nclTpc);
    clusters->SetMarkerStyle(21);
    clusters->SetMarkerSize(1.0);
    clusters->SetMarkerColor(kBlue);
    
    Int_t npoints = 0;
    for(Int_t i = 0; i < gEvent3D.nclTpc; i++) {
      if(i < static_cast<Int_t>(gEvent3D.cluster_x.size()) &&
         i < static_cast<Int_t>(gEvent3D.cluster_y.size()) &&
         i < static_cast<Int_t>(gEvent3D.cluster_z.size())) {
        clusters->SetPoint(npoints++,
                          gEvent3D.cluster_x[i],
                          gEvent3D.cluster_y[i],
                          gEvent3D.cluster_z[i]);
      }
    }
    if(npoints > 0) {
      clusters->SetNextPoint(npoints);
      clusters->Draw();
      gClusters3D.push_back(clusters);
    } else {
      delete clusters;
    }
  }
  
  // Draw Hough-selected clusters
  if(gEvent3D.nclTpc > 0) {
    TPolyMarker3D* houghClusters = new TPolyMarker3D(gEvent3D.nclTpc);
    houghClusters->SetMarkerStyle(22);
    houghClusters->SetMarkerSize(1.2);
    houghClusters->SetMarkerColor(kRed);
    
    Int_t npoints = 0;
    for(Int_t i = 0; i < gEvent3D.nclTpc; i++) {
      if(i < static_cast<Int_t>(gEvent3D.cluster_houghflag.size()) &&
         gEvent3D.cluster_houghflag[i] > 0) {
        if(i < static_cast<Int_t>(gEvent3D.cluster_x.size()) &&
           i < static_cast<Int_t>(gEvent3D.cluster_y.size()) &&
           i < static_cast<Int_t>(gEvent3D.cluster_z.size())) {
          houghClusters->SetPoint(npoints++,
                                 gEvent3D.cluster_x[i],
                                 gEvent3D.cluster_y[i],
                                 gEvent3D.cluster_z[i]);
        }
      }
    }
    if(npoints > 0) {
      houghClusters->SetNextPoint(npoints);
      houghClusters->Draw();
      gHoughClusters3D.push_back(houghClusters);
    } else {
      delete houghClusters;
    }
  }
  
  // Draw tracks
  if(gEvent3D.ntTpc > 0) {
    for(Int_t itrack = 0; itrack < gEvent3D.ntTpc; itrack++) {
      if(itrack < static_cast<Int_t>(gEvent3D.x0Tpc.size()) &&
         itrack < static_cast<Int_t>(gEvent3D.y0Tpc.size()) &&
         itrack < static_cast<Int_t>(gEvent3D.u0Tpc.size()) &&
         itrack < static_cast<Int_t>(gEvent3D.v0Tpc.size())) {
        
        Double_t x0 = gEvent3D.x0Tpc[itrack];
        Double_t y0 = gEvent3D.y0Tpc[itrack];
        Double_t u0 = gEvent3D.u0Tpc[itrack];
        Double_t v0 = gEvent3D.v0Tpc[itrack];
        
        // Calculate track points at z boundaries
        Double_t z1 = zmin;
        Double_t z2 = zmax;
        Double_t x1 = x0 + u0 * z1;
        Double_t y1 = y0 + v0 * z1;
        Double_t x2 = x0 + u0 * z2;
        Double_t y2 = y0 + v0 * z2;
        
        TPolyLine3D* track = new TPolyLine3D(2);
        track->SetLineWidth(2);
        track->SetLineColor(itrack % 9 + 1);
        track->SetPoint(0, x1, y1, z1);
        track->SetPoint(1, x2, y2, z2);
        track->Draw();
        gTracks3D.push_back(track);
      }
    }
  }
  
  // Draw track hits if available
  if(gEvent3D.ntTpc > 0) {
    for(Int_t itrack = 0; itrack < gEvent3D.ntTpc; itrack++) {
      if(itrack < static_cast<Int_t>(gEvent3D.nhtrack.size()) &&
         gEvent3D.nhtrack[itrack] > 0) {
        Int_t nhits = gEvent3D.nhtrack[itrack];
        if(nhits > 0 &&
           itrack < static_cast<Int_t>(gEvent3D.hitpos_x.size()) &&
           nhits <= static_cast<Int_t>(gEvent3D.hitpos_x[itrack].size())) {
          TPolyMarker3D* trackHits = new TPolyMarker3D(nhits);
          trackHits->SetMarkerStyle(24);
          trackHits->SetMarkerSize(0.8);
          trackHits->SetMarkerColor(itrack % 9 + 1);
          
          for(Int_t i = 0; i < nhits; i++) {
            trackHits->SetPoint(i,
                               gEvent3D.hitpos_x[itrack][i],
                               gEvent3D.hitpos_y[itrack][i],
                               gEvent3D.hitpos_z[itrack][i]);
          }
          trackHits->SetNextPoint(nhits);
          trackHits->Draw();
          gRawHits3D.push_back(trackHits);
        }
      }
    }
  }
  
  gMacroCanvas3D->Update();
  
  // Print summary
  std::cout << "\n=== Event Summary (3D) ===" << std::endl;
  std::cout << "Raw hits: " << gEvent3D.nhTpc << std::endl;
  std::cout << "Clusters: " << gEvent3D.nclTpc << std::endl;
  Int_t nHoughClusters = 0;
  for(Int_t i = 0; i < gEvent3D.nclTpc; i++) {
    if(i < static_cast<Int_t>(gEvent3D.cluster_houghflag.size()) &&
       gEvent3D.cluster_houghflag[i] > 0) nHoughClusters++;
  }
  std::cout << "Hough-selected clusters: " << nHoughClusters << std::endl;
  std::cout << "Tracks: " << gEvent3D.ntTpc << std::endl;
  for(Int_t itrack = 0; itrack < gEvent3D.ntTpc; itrack++) {
    if(itrack < static_cast<Int_t>(gEvent3D.x0Tpc.size())) {
      std::cout << "  Track " << itrack << ": "
           << "x0=" << gEvent3D.x0Tpc[itrack] << ", y0=" << gEvent3D.y0Tpc[itrack]
           << ", u0=" << gEvent3D.u0Tpc[itrack] << ", v0=" << gEvent3D.v0Tpc[itrack]
           << ", theta=" << (itrack < static_cast<Int_t>(gEvent3D.theta.size()) ? gEvent3D.theta[itrack] : 0.0)
           << ", nhits=" << (itrack < static_cast<Int_t>(gEvent3D.nhtrack.size()) ? gEvent3D.nhtrack[itrack] : 0) << std::endl;
    }
  }
  std::cout << "=========================" << std::endl;
  std::cout << "\n3D View Controls:" << std::endl;
  std::cout << "  - Left mouse + drag: Rotate" << std::endl;
  std::cout << "  - Right mouse + drag: Zoom" << std::endl;
  std::cout << "  - Middle mouse + drag: Pan" << std::endl;
  std::cout << "  - Mouse wheel: Zoom in/out" << std::endl;
}
