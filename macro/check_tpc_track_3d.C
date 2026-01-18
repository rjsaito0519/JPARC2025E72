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

// Forward declaration
void DrawAxisLabel(Double_t x, Double_t y, Double_t z, Char_t label, Color_t color, Double_t labelSize);

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
DrawAxisLabel(Double_t x, Double_t y, Double_t z, Char_t label, Color_t color, Double_t labelSize)
{
  // Draw a simple axis label using small lines to form letters
  // X: two crossing lines
  if(label == 'X') {
    TPolyLine3D* x1 = new TPolyLine3D(2);
    x1->SetLineColor(color);
    x1->SetLineWidth(2);
    x1->SetPoint(0, x - labelSize * 0.3, y - labelSize * 0.3, z);
    x1->SetPoint(1, x + labelSize * 0.3, y + labelSize * 0.3, z);
    x1->Draw();
    gTPCFrame.push_back(x1);
    
    TPolyLine3D* x2 = new TPolyLine3D(2);
    x2->SetLineColor(color);
    x2->SetLineWidth(2);
    x2->SetPoint(0, x - labelSize * 0.3, y + labelSize * 0.3, z);
    x2->SetPoint(1, x + labelSize * 0.3, y - labelSize * 0.3, z);
    x2->Draw();
    gTPCFrame.push_back(x2);
  }
  // Y: V shape with vertical line
  else if(label == 'Y') {
    TPolyLine3D* y1 = new TPolyLine3D(2);
    y1->SetLineColor(color);
    y1->SetLineWidth(2);
    y1->SetPoint(0, x, y - labelSize * 0.4, z);
    y1->SetPoint(1, x, y + labelSize * 0.2, z);
    y1->Draw();
    gTPCFrame.push_back(y1);
    
    TPolyLine3D* y2 = new TPolyLine3D(2);
    y2->SetLineColor(color);
    y2->SetLineWidth(2);
    y2->SetPoint(0, x, y + labelSize * 0.2, z);
    y2->SetPoint(1, x - labelSize * 0.3, y + labelSize * 0.4, z);
    y2->Draw();
    gTPCFrame.push_back(y2);
    
    TPolyLine3D* y3 = new TPolyLine3D(2);
    y3->SetLineColor(color);
    y3->SetLineWidth(2);
    y3->SetPoint(0, x, y + labelSize * 0.2, z);
    y3->SetPoint(1, x + labelSize * 0.3, y + labelSize * 0.4, z);
    y3->Draw();
    gTPCFrame.push_back(y3);
  }
  // Z: three lines forming Z shape
  else if(label == 'Z') {
    TPolyLine3D* z1 = new TPolyLine3D(2);
    z1->SetLineColor(color);
    z1->SetLineWidth(2);
    z1->SetPoint(0, x - labelSize * 0.3, y + labelSize * 0.3, z);
    z1->SetPoint(1, x + labelSize * 0.3, y + labelSize * 0.3, z);
    z1->Draw();
    gTPCFrame.push_back(z1);
    
    TPolyLine3D* z2 = new TPolyLine3D(2);
    z2->SetLineColor(color);
    z2->SetLineWidth(2);
    z2->SetPoint(0, x + labelSize * 0.3, y + labelSize * 0.3, z);
    z2->SetPoint(1, x - labelSize * 0.3, y - labelSize * 0.3, z);
    z2->Draw();
    gTPCFrame.push_back(z2);
    
    TPolyLine3D* z3 = new TPolyLine3D(2);
    z3->SetLineColor(color);
    z3->SetLineWidth(2);
    z3->SetPoint(0, x - labelSize * 0.3, y - labelSize * 0.3, z);
    z3->SetPoint(1, x + labelSize * 0.3, y - labelSize * 0.3, z);
    z3->Draw();
    gTPCFrame.push_back(z3);
  }
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
  
  // Fixed range for view (to prevent scale changes between events)
  // Note: Coordinate transformation (x,y,z) -> (z,x,y)
  // Original coordinates: X: -300 to 300, Y: -310 to 310 (TPC frame height), Z: -300 to 300
  const Double_t xmin_orig = -300.0, xmax_orig = 300.0;
  const Double_t ymin_orig = -310.0, ymax_orig = 310.0;  // Fixed to TPC frame height
  const Double_t zmin_orig = -300.0, zmax_orig = 300.0;
  
  // Transform range: (x,y,z) -> (z,x,y)
  const Double_t xmin = zmin_orig, xmax = zmax_orig; // Display X = original Z
  const Double_t ymin = xmin_orig, ymax = xmax_orig; // Display Y = original X
  const Double_t zmin = ymin_orig, zmax = ymax_orig; // Display Z = original Y
  
  // Create 3D view (without axis marks)
  gMacroView = TView::CreateView(1, 0, 0);
  gMacroView->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);
  // Don't call ShowAxis() to avoid axis marks
  
  // Draw small axis arrows in the corner
  // Position: near the minimum corner of the view
  Double_t axisOriginX = xmin + (xmax - xmin) * 0.05;
  Double_t axisOriginY = ymin + (ymax - ymin) * 0.05;
  Double_t axisOriginZ = zmin + (zmax - zmin) * 0.05;
  Double_t axisLength = TMath::Min(TMath::Min(xmax - xmin, ymax - ymin), zmax - zmin) * 0.1;
  
  // X axis (original Z direction) - red
  TPolyLine3D* axisX = new TPolyLine3D(2);
  axisX->SetLineColor(kRed);
  axisX->SetLineWidth(2);
  axisX->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ);
  axisX->SetPoint(1, axisOriginX + axisLength, axisOriginY, axisOriginZ);
  axisX->Draw();
  gTPCFrame.push_back(axisX);
  
  // Arrow head for X axis
  TPolyLine3D* arrowX1 = new TPolyLine3D(3);
  arrowX1->SetLineColor(kRed);
  arrowX1->SetLineWidth(2);
  Double_t arrowSize = axisLength * 0.15;
  arrowX1->SetPoint(0, axisOriginX + axisLength, axisOriginY, axisOriginZ);
  arrowX1->SetPoint(1, axisOriginX + axisLength - arrowSize, axisOriginY - arrowSize * 0.3, axisOriginZ);
  arrowX1->SetPoint(2, axisOriginX + axisLength - arrowSize, axisOriginY + arrowSize * 0.3, axisOriginZ);
  arrowX1->Draw();
  gTPCFrame.push_back(arrowX1);
  
  // Y axis (original X direction) - green
  TPolyLine3D* axisY = new TPolyLine3D(2);
  axisY->SetLineColor(kGreen + 2);
  axisY->SetLineWidth(2);
  axisY->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ);
  axisY->SetPoint(1, axisOriginX, axisOriginY + axisLength, axisOriginZ);
  axisY->Draw();
  gTPCFrame.push_back(axisY);
  
  // Arrow head for Y axis
  TPolyLine3D* arrowY1 = new TPolyLine3D(3);
  arrowY1->SetLineColor(kGreen + 2);
  arrowY1->SetLineWidth(2);
  arrowY1->SetPoint(0, axisOriginX, axisOriginY + axisLength, axisOriginZ);
  arrowY1->SetPoint(1, axisOriginX - arrowSize * 0.3, axisOriginY + axisLength - arrowSize, axisOriginZ);
  arrowY1->SetPoint(2, axisOriginX + arrowSize * 0.3, axisOriginY + axisLength - arrowSize, axisOriginZ);
  arrowY1->Draw();
  gTPCFrame.push_back(arrowY1);
  
  // Z axis (original Y direction) - blue
  TPolyLine3D* axisZ = new TPolyLine3D(2);
  axisZ->SetLineColor(kBlue);
  axisZ->SetLineWidth(2);
  axisZ->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ);
  axisZ->SetPoint(1, axisOriginX, axisOriginY, axisOriginZ + axisLength);
  axisZ->Draw();
  gTPCFrame.push_back(axisZ);
  
  // Arrow head for Z axis
  TPolyLine3D* arrowZ1 = new TPolyLine3D(3);
  arrowZ1->SetLineColor(kBlue);
  arrowZ1->SetLineWidth(2);
  arrowZ1->SetPoint(0, axisOriginX, axisOriginY, axisOriginZ + axisLength);
  arrowZ1->SetPoint(1, axisOriginX - arrowSize * 0.3, axisOriginY, axisOriginZ + axisLength - arrowSize);
  arrowZ1->SetPoint(2, axisOriginX + arrowSize * 0.3, axisOriginY, axisOriginZ + axisLength - arrowSize);
  arrowZ1->Draw();
  gTPCFrame.push_back(arrowZ1);
  
  // Draw axis labels using small markers/lines to form letters
  // X axis label
  DrawAxisLabel(axisOriginX + axisLength * 1.3, axisOriginY, axisOriginZ, 'X', kRed, axisLength);
  // Y axis label
  DrawAxisLabel(axisOriginX, axisOriginY + axisLength * 1.3, axisOriginZ, 'Y', kGreen + 2, axisLength);
  // Z axis label
  DrawAxisLabel(axisOriginX, axisOriginY, axisOriginZ + axisLength * 1.3, 'Z', kBlue, axisLength);
  
  // Draw TPC frame (regular octagonal prism)
  // Frame dimensions:
  //   - Height: 620 in Y direction (center at 0, so -310 to +310)
  //   - Octagon face-to-face distance: 622 (center at 0)
  //   - Coordinate transform: (x,y,z) -> (z,x,y)
  
  const Double_t frameHeight = 620.0;  // Y direction
  const Double_t frameHalfHeight = frameHeight / 2.0;  // 310.0
  const Double_t frameFaceToFace = 622.0;
  const Double_t frameRadius = frameFaceToFace / 2.0;  // 311.0
  const Int_t nOctagonSides = 8;
  
  // Calculate octagon vertices (in X-Z plane, Y is height)
  // For regular octagon, vertices are at angles: 0, 45, 90, 135, 180, 225, 270, 315 degrees
  // Distance from center to vertex = face-to-face / (2 * cos(22.5°))
  const Double_t octagonVertexRadius = frameRadius / TMath::Cos(TMath::Pi() / 8.0);
  
  std::vector<Double_t> octagon_x_orig(nOctagonSides);
  std::vector<Double_t> octagon_z_orig(nOctagonSides);
  for(Int_t i = 0; i < nOctagonSides; i++) {
    Double_t angle = i * TMath::Pi() / 4.0;  // 45 degrees per side
    octagon_x_orig[i] = octagonVertexRadius * TMath::Cos(angle);
    octagon_z_orig[i] = octagonVertexRadius * TMath::Sin(angle);
  }
  
  // Draw top and bottom octagons
  for(Int_t y_sign = -1; y_sign <= 1; y_sign += 2) {
    Double_t y_orig = y_sign * frameHalfHeight;  // -310 or +310
    TPolyLine3D* octagon = new TPolyLine3D(nOctagonSides + 1);
    octagon->SetLineColor(kGray + 1);
    octagon->SetLineWidth(2);
    for(Int_t i = 0; i <= nOctagonSides; i++) {
      Int_t idx = i % nOctagonSides;
      // Transform: (x,y,z) -> (z,x,y)
      octagon->SetPoint(i, octagon_z_orig[idx], octagon_x_orig[idx], y_orig);
    }
    octagon->Draw();
    gTPCFrame.push_back(octagon);
  }
  
  // Draw vertical edges (8 edges connecting top and bottom octagons)
  for(Int_t i = 0; i < nOctagonSides; i++) {
    TPolyLine3D* edge = new TPolyLine3D(2);
    edge->SetLineColor(kGray + 1);
    edge->SetLineWidth(2);
    // Bottom vertex
    edge->SetPoint(0, octagon_z_orig[i], octagon_x_orig[i], -frameHalfHeight);
    // Top vertex
    edge->SetPoint(1, octagon_z_orig[i], octagon_x_orig[i], frameHalfHeight);
    edge->Draw();
    gTPCFrame.push_back(edge);
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
