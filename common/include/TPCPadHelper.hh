// -*- C++ -*-
/*
 * TPCPadHelper for myanalysis
 * - Pad / Layer / Row indices are unified to start from 0
 * - Dependencies are self-contained within myanalysis
 *
 * Based on: include/TPCPadHelper.hh
 * Author: Haein Lee
 * Date  : 2025-12-26
 */

#ifndef TPC_PAD_HELPER_HH
#define TPC_PAD_HELPER_HH

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stdexcept>

#include <TDirectory.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TVector3.h>
#include <TString.h>

namespace tpc
{
// Constants from DetectorID.hh
const Int_t NumOfLayersTPC = 32;
const Int_t NumOfPadTPC = 5768;

const Double_t ZTarget = -143.; // Target from center
const Double_t TargetVtxWindow = 30.;

enum EPadParameter
{
  kLayerID,
  kNumOfPad,
  kRadius,
  kNumOfDivision,
  kDummy,
  kLength,
  NPadParameter
};

//for kinematic fitting
//E42
static const Double_t PositionScale = 1.;

//for clustering
//E42
static const Int_t MaxRowDifTPC = 2;

//for helix tracking
static const Double_t ConstC = 0.299792458; //=c/10^9
  
//_____________________________________________________________________________
//#OfPad #division #radius padLength
static const Double_t padParameter[NumOfLayersTPC][NPadParameter] =
{{0, 48,    14.75, 48, 0,  9.},
 {1, 48,    24.25, 48, 0,  9.},
 {2, 72,    33.75, 72, 0,  9.},
 {3, 96,    43.25, 96, 0,  9.},
 {4, 120,    52.75,120,0,   9.},
 {5, 144,    62.25,144,0,   9.},
 {6, 168,    71.75,168,0,   9.},
 {7, 192,    81.25,192,0,   9.},
 {8, 216,    90.75,216,0,   9.},
 {9, 240,    100.25,240,0,  9.},
 {10,208,    111.5,241, 0,  12.5},
 {11,218,    124.5,271, 0,  12.5},
 {12,230,    137.5,300, 0,  12.5},
 {13,214,    150.5,330, 0,  12.5},
 {14,212,    163.5,360, 0,  12.5},
 {15,214,    176.5,390, 0,  12.5},
 {16,220,    189.5,420, 0,  12.5},
 {17,224,    202.5,449, 0,  12.5},
 {18,232,    215.5,479, 0,  12.5},
 {19,238,    228.5,509, 0,  12.5},
 {20,244,    241.5,539, 0,  12.5},
 {21,232,    254.5,569, 0,  12.5},
 {22,218,    267.5,599, 0,  12.5},
 {23,210,    280.5,628, 0,  12.5},
 {24,206,    293.5,658, 0,  12.5},
 {25,202,    306.5,688, 0,  12.5},
 {26,200,    319.5,718, 0,  12.5},
 {27,196,    332.5,748, 0,  12.5},
 {28,178,    345.5,777, 0,  12.5},
 {29,130,    358.5,807, 0,  12.5},
 {30,108,    371.5,837, 0, 12.5},
 {31,90,     384.5,867, 0, 12.5}};

// Dead channels and frame pads (simplified - only essential arrays)
static const Int_t deadChannel[] =
{
  72,134,135,222,332,408,467,626,809,59,110,185,284,555,726,921,1140,1363,1564,1565,1778,1818,1819,2036,2037,2245,2246,2247,2455,2456,2668,2669,2887,3111,3112,3342,3343,3579,3580,3814,4035,1629,1850,2068,2277,2486,2487,2699,2700,2918,3142,3374,3611,3845,4066,4277,4481,4681,4879,5073
};
static const Int_t NumDeadChannels = sizeof(deadChannel) / sizeof(deadChannel[0]);

//_____________________________________________________________________________
// Simple error handling (replacing Exception.hh)
inline void ErrorMessage(const char* funcName, const char* message)
{
  std::cerr << "[tpc::" << funcName << "] " << message << std::endl;
}

//_____________________________________________________________________________
// Check validity of Layer ID
inline void ValidateLayer(Int_t layer, const char* funcName)
{
  if (layer < 0 || NumOfLayersTPC <= layer) {
    std::string msg = "Invalid layerID " + std::to_string(layer) + 
                      " (Max: " + std::to_string(NumOfLayersTPC - 1) + ")";
    ErrorMessage(funcName, msg.c_str());
    throw std::runtime_error(msg);
  }
}

//_____________________________________________________________________________
// Check validity of Row ID for a specific Layer
inline void ValidateRow(Int_t layer, Int_t row, const char* funcName)
{
  ValidateLayer(layer, funcName);
  Int_t max_row = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  if (row < 0 || max_row <= row) {
    std::string msg = "Invalid rowID " + std::to_string(row) + 
                      " for layer " + std::to_string(layer) + 
                      " (Max: " + std::to_string(max_row - 1) + ")";
    ErrorMessage(funcName, msg.c_str());
    throw std::runtime_error(msg);
  }
}

//_____________________________________________________________________________
// Check validity of Pad ID (Global ID)
inline void ValidatePadID(Int_t padID, const char* funcName)
{
  if (padID < 0 || NumOfPadTPC <= padID) {
    std::string msg = "Invalid padID " + std::to_string(padID) + 
                      " (Max: " + std::to_string(NumOfPadTPC - 1) + ")";
    ErrorMessage(funcName, msg.c_str());
    throw std::runtime_error(msg);
  }
}

//_____________________________________________________________________________
inline Int_t GetPadId(Int_t layerID, Int_t rowID)
{
  if (rowID < 0) return -1;
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layerID, rowID, __func__);
#endif

  Int_t padID = 0;
  for (Int_t layer = 0; layer < layerID; layer++)
    padID += static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  padID += rowID;
  return padID;
}

//_____________________________________________________________________________
inline Int_t getLayerID(Int_t padID)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(padID, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    sum += static_cast<Int_t>(padParameter[layer][kNumOfPad]);
    if (padID < sum) return layer;
  }
  
  std::string msg = "Invalid padID " + std::to_string(padID) + 
                    " (Max: " + std::to_string(NumOfPadTPC - 1) + ")";
  ErrorMessage(__func__, msg.c_str());
  throw std::runtime_error(msg);
}

//_____________________________________________________________________________
inline Int_t getRowID(Int_t padID)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(padID, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t nPad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
    if (padID < sum + nPad) return padID - sum;
    sum += nPad;
  }

  std::string msg = "Logic Error: padID " + std::to_string(padID) + " not found in parameter table.";
  ErrorMessage(__func__, msg.c_str());
  throw std::runtime_error(msg);
}

//_____________________________________________________________________________
inline Double_t getTheta(Int_t padID)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(padID, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t nPad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);

    if (padID < sum + nPad) {
      Int_t row = padID - sum;
      Double_t nDiv = padParameter[layer][kNumOfDivision];
      Double_t sTheta = 180. - (360. / nDiv) * nPad / 2.;
      Double_t theta  = sTheta + (row + 0.5) * 360. / nDiv - 180;
      return theta;
    }
    sum += nPad;
  }

  std::string msg = "Logic Error: padID " + std::to_string(padID) + " not found.";
  ErrorMessage(__func__, msg.c_str());
  throw std::runtime_error(msg);
}

//_____________________________________________________________________________
inline Double_t getTheta(Int_t layer, Double_t m_row)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
  Double_t max_row = padParameter[layer][kNumOfPad];
  Double_t epsilon = -1.e-3;
  if (m_row < epsilon || max_row < m_row) {
    std::string msg = "Invalid m_row " + std::to_string(m_row) + 
                      " for layer " + std::to_string(layer) + 
                      " (Limit: < " + std::to_string(max_row) + ")";
    ErrorMessage(__func__, msg.c_str());
    throw std::runtime_error(msg);
  }
#endif

  Int_t    nPad   = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
  Double_t nDiv   = padParameter[layer][kNumOfDivision];
  Double_t sTheta = 180. - (360. / nDiv) * nPad / 2.;
  Double_t theta  = sTheta + (m_row + 0.5) * 360. / nDiv - 180;

  return theta;
}

//_____________________________________________________________________________
inline Double_t GetRadius(Int_t layer)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif
  return padParameter[layer][kRadius];
}

//_____________________________________________________________________________
inline Double_t getR(Int_t padID)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(padID, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t nPad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
    if (padID < sum + nPad) {
      return padParameter[layer][kRadius];
    }
    sum += nPad;
  }

  std::string msg = "Logic Error: padID " + std::to_string(padID) + " not found.";
  ErrorMessage(__func__, msg.c_str());
  throw std::runtime_error(msg);
}

//_____________________________________________________________________________
inline TVector3 getPosition(Int_t padID)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(padID, __func__);
#endif

  Int_t sum = 0;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Int_t nPad = static_cast<Int_t>(padParameter[layer][kNumOfPad]);

    if (padID < sum + nPad) {
      Int_t row = padID - sum;
      Double_t theta  = getTheta(layer, static_cast<Double_t>(row))*TMath::DegToRad();
      Double_t radius = padParameter[layer][kRadius];
      Double_t x = radius * std::sin(theta);
      Double_t z = radius * std::cos(theta) + ZTarget;
      return TVector3(x, 0., z);
    }
    sum += nPad;
  }

  return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
}

//_____________________________________________________________________________
inline TVector3 getPosition(Int_t layer, Double_t m_row)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif

  Double_t nPad = padParameter[layer][kNumOfPad];
  if (m_row < 0 || nPad < m_row) {
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  }

  Double_t theta  = getTheta(layer, m_row) * TMath::DegToRad();
  Double_t radius = padParameter[layer][kRadius];
  Double_t x = radius * std::sin(theta);
  Double_t z = radius * std::cos(theta) + ZTarget;
  return TVector3(x, 0., z);
}

//_____________________________________________________________________________
inline Int_t getSection(Int_t padID)
{
  TVector3 pad_pos = getPosition(padID);
  Double_t x = pad_pos.x();
  Double_t z = pad_pos.z();

  Int_t section = -1;
  if (std::abs(x) == std::abs(z)) return section;
  if (std::abs(x) <  std::abs(z)) {
    if (z < 0) section = 1;
    else       section = 3;
  } else {
    if (x < 0) section = 2;
    else       section = 4;
  }

  return section;
}

//_____________________________________________________________________________
inline Int_t getSection(Int_t layer, Double_t m_row)
{
  Int_t padid = GetPadId(layer, m_row);
  return getSection(padid);
}

//_____________________________________________________________________________
// Find PadID from global position (z, x)
inline Int_t findPadID(Double_t z, Double_t x)
{
  Double_t radius = std::hypot(x, z-ZTarget);
  Double_t angle  = 180.0 + std::atan2(x, z-ZTarget) * TMath::RadToDeg();
  if (angle >= 360.0) angle -= 360.0;
  if (angle <    0.0) angle += 360.0;

  Int_t hit_layer = -1;
  for (Int_t layer = 0; layer < NumOfLayersTPC; layer++) {
    Double_t r_pad = padParameter[layer][kRadius];
    Double_t l_pad = padParameter[layer][kLength];
    Double_t r_in  = r_pad - l_pad * 0.5;
    Double_t r_out = r_pad + l_pad * 0.5;
    if (r_in <= radius && radius <= r_out) {
      hit_layer = layer;
      break;
    }

    if (layer > 0) {
      Double_t r_prev_out = padParameter[layer-1][kRadius] + padParameter[layer-1][kLength] * 0.5;
      if (r_prev_out < radius && radius < r_in) {
        return -layer; 
      }
    }
  }
  if (hit_layer == -1) return -1000;

  Double_t nPad = padParameter[hit_layer][kNumOfPad];
  Double_t nDiv = padParameter[hit_layer][kNumOfDivision];
  Double_t sTheta = 180. - (360. / nDiv) * nPad / 2.;
  Double_t dTheta = 360. / nDiv;

  Double_t diff = angle - sTheta;
  if (std::isnan(diff) || diff < 0) return -1000;
  Int_t row = static_cast<Int_t>(diff / dTheta);
  if (row < 0 || static_cast<Int_t>(nPad) <= row) return -1000;

  return GetPadId(hit_layer, row);
}

//_____________________________________________________________________________
inline Double_t ArcLength(Int_t layer, Double_t row1, Double_t row2)
{
#ifdef PAD_HELPER_DEBUG
  ValidateLayer(layer, __func__);
#endif

  Double_t radius = padParameter[layer][kRadius];
  Double_t theta1 = getTheta(layer, row1);
  Double_t theta2 = getTheta(layer, row2);
  Double_t diff   = std::abs(theta1 - theta2);
  diff = std::fmod(diff, 360.0);
  if (diff > 180.0) diff = 360.0 - diff;

  return radius * diff * TMath::DegToRad();
}

//_____________________________________________________________________________
inline void 
InitializeHistograms(const char* name)
{
  auto h1 = gDirectory->Get<TH2Poly>(name);
  if (!h1) {
    std::cerr << "[tpc::InitializeHistograms] Error: TH2Poly '" << name << "' not found in gDirectory" << std::endl;
    return;
  }

  Double_t X[5];
  Double_t Y[5];
  for (Int_t layer = 0; layer < NumOfLayersTPC; ++layer) {    
    Int_t nPad      = static_cast<Int_t>(padParameter[layer][kNumOfPad]);
    Double_t nDiv   = padParameter[layer][kNumOfDivision];
    Double_t radius = padParameter[layer][kRadius];
    Double_t length = padParameter[layer][kLength];
    Double_t r_min  = radius - length / 2.0;
    Double_t r_max  = radius + length / 2.0;
    Double_t dTheta = 360.0 / nDiv;
    Double_t sTheta = - dTheta * nPad / 2.0;

    for (Int_t j = 0; j < nPad; ++j) {      
      Double_t theta1 = (sTheta + j*dTheta)     * TMath::DegToRad();
      Double_t theta2 = (sTheta + (j+1)*dTheta) * TMath::DegToRad();

      X[0] = r_max * std::cos(theta1);
      X[1] = r_max * std::cos(theta2);
      X[2] = r_min * std::cos(theta2);
      X[3] = r_min * std::cos(theta1);
      X[4] = X[0];
      for (Int_t k = 0; k < 5; ++k) X[k] += ZTarget; 
            
      Y[0] = r_max * std::sin(theta1);
      Y[1] = r_max * std::sin(theta2);
      Y[2] = r_min * std::sin(theta2);
      Y[3] = r_min * std::sin(theta1);
      Y[4] = Y[0];

      h1->AddBin(5, X, Y);
    }
  }
}

//_____________________________________________________________________________
inline Bool_t IsDead(Int_t padID)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(padID, __func__);
#endif

  Bool_t is_dead_channel = std::find(deadChannel, deadChannel + NumDeadChannels, padID) 
                           != (deadChannel + NumDeadChannels);
  return is_dead_channel;
}

//_____________________________________________________________________________
inline Bool_t IsDead(Int_t layer, Int_t row){
  Int_t padID = GetPadId(layer, row);
  return IsDead(padID);
}

} // namespace tpc

#endif
