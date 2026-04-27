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
#include <cstddef>
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
static const Int_t padOnCenterFrame[] =
{

//Pads on the frame
964,965,966,967,968,969,970,971,972,973,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1394,1395,1396,1397,1398,1399,1400,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1435,1436,1437,1438,1439,1440,1441,1454,1455,1456,1457,1458,1459,1460,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1495,1496,1497,1498,1499,1500,1501,1594,1595,1596,1597,1603,1604,1605,1606,1607,1646,1647,1648,1649,1650,1651,1657,1658,1659,1662,1663,1664,1670,1671,1672,1673,1674,1675,1714,1715,1716,1717,1718,1724,1725,1726,1727,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1878,1879,1880,1881,1888,1889,1890,1891,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,2017,2018,2019,2024,2025,2026,2027,2099,2100,2101,2102,2103,2110,2111,2112,2113,2114,2186,2187,2188,2189,2194,2195,2196,2218,2219,2220,2221,2222,2223,2224,2225,2226,2227,2308,2309,2310,2316,2317,2322,2323,2329,2330,2331,2412,2413,2414,2415,2416,2417,2418,2419,2420,2421,2426,2427,2428,2429,2517,2518,2519,2520,2521,2522,2523,2524,2525,2526,2539,2540,2541,2542,2543,2544,2545,2546,2547,2548,2636,2637,2638,2639,2731,2732,2738,2739,2760,2761,2767,2768,2949,2950,2951,2952,2953,2954,2955,2956,2957,2986,2987,2988,2989,2990,2991,2992,2993,2994,3173,3174,3175,3176,3177,3178,3179,3180,3181,3218,3219,3220,3221,3222,3223,3224,3225,3226,3404,3405,3406,3407,3408,3409,3410,3411,3412,3457,3458,3459,3460,3461,3462,3463,3464,3465,3641,3642,3643,3644,3645,3646,3647,3648,3649,3702,3703,3704,3705,3706,3707,3708,3709,3710,3876,3877,3878,3879,3880,3881,3882,3883,3944,3945,3946,3947,3948,3949,3950,3951,4097,4098,4104,4173,4179,4180,4307,4308,4309,4310,4311,4312,4313,4314,4315,4390,4391,4392,4393,4394,4395,4396,4397,4398,4511,4512,4519,4602,4609,4610,4712,4713,4714,4715,4716,4717,4718,4719,4810,4811,4812,4813,4814,4815,4816,4817,4909,4910,4911,4912,4913,4914,4915,4916,5015,5016,5017,5018,5019,5020,5021,5022,5103,5104,5107,5110,5111,5216,5217,5220,5223,5224,5286,5287,5288,5289,5290,5291,5292,5293,5294,5407,5408,5409,5410,5411,5412,5413,5414,5415,5440,5441,5442,5443,5444,5565,5566,5567,5568,5569,

//Empty Pads
1016,1017,1393,1401,1402,1403,1432,1433,1434,1461,1462,1463,1492,1493,1494,1502,1598,1599,1600,1601,1602,1652,1653,1654,1655,1656,1665,1666,1667,1668,1669,1719,1720,1721,1722,1723,1876,1877,1882,1883,1884,1885,1886,1887,1892,1893,2020,2021,2022,2023,2104,2105,2106,2107,2108,2109,2190,2191,2192,2193,2311,2312,2313,2314,2315,2324,2325,2326,2327,2328,2733,2734,2735,2736,2737,2762,2763,2764,2765,2766,4099,4100,4101,4102,4103,4174,4175,4176,4177,4178,4513,4514,4515,4516,4517,4518,4603,4604,4605,4606,4607,4608,5105,5106,5108,5109,5218,5219,5221,5222

};

static const std::size_t NumPadOnCenterFrame = sizeof(padOnCenterFrame) / sizeof(padOnCenterFrame[0]);

//_____________________________________________________________________________
static const Int_t padOnSectionFrame[] =
{
//Gem section 1 frame 0
1503 ,1504, 1505 ,1506 ,1507 ,1508 ,1509 ,1510 ,1511 ,1512 ,1513 ,1514 ,1515 ,1728 ,1729 ,1730,

//Gem section 1 frame 1
72 ,134 ,135 ,222 ,332 , 408,  467 ,626 ,809 ,1016, 59 ,110 ,185 ,284,555 ,726, 921 ,1140 ,1363, 1564 ,1565 ,1778,

//Gem section 2 frame 0, frame 1 + adjacent pads
1628,1629,1630,1817,1818,1819,1820,1849,1850,1851,2035,2036,2037,2038,2067,2068,2069,2244,2245,2246,2247,2248,2276,2277,2278,2454,2455,2456,2457,2485,2486,2487,2488,2667,2668,2669,2670,2698,2699,2700,2701,2886,2887,2888,2917,2918,2919,3110,3111,3112,3113,3141,3142,3143,3341,3342,3343,3344,3373,3374,3375,3578,3579,3580,3581,3610,3611,3612,3813,3814,3815,3844,3845,3846,4034,4035,4036,4065,4066,4067,4276,4277,4278,4480,4481,4482,4680,4681,4682,4878,4879,4880,5072,5073,5074,

//Gem section 3 frame 0
4720 ,4721 ,4722 ,4723 ,4724 ,4725 ,4726 ,4727 ,4728 ,4929 ,4930 ,4931 ,4932 ,4933 ,4934 ,4935 ,4936 ,5134 ,5135 ,5136 ,5137 ,5138 ,5139 ,5140 ,5328 ,5329 ,5330 ,5331 ,5332 ,5333 ,5487 ,5488 ,5489 ,5490 ,5491 ,5492 ,5612 ,5613 ,5614 ,5615 ,5616 ,5716 ,5717 ,5718 ,5719 ,5720, 5532 ,5533, 5612 ,5613 ,5614 ,5615 ,5616 ,5654 ,5655 ,5716 ,5717 ,5718 ,5719 ,5720 ,5757 ,5758,

//Gem section 3 frame 1
3413 ,3414 ,3415 ,3416 ,3417 ,3418 ,3419 ,3660 ,3661 ,3662 ,3663 ,3664 ,3665 ,3904 ,3905 ,3906 ,3907 ,3908 ,4134 ,4135 ,4136 ,4137 ,4353 ,4354 ,4355 ,4356 ,4565 ,4566 ,4567 ,4568 ,4774 ,4775 ,4776 ,4981 ,4982 ,5183 ,5184 ,5372 ,5373 ,5374 ,5530 ,5531 ,5532 ,5533 ,5654 ,5655 ,5757 ,5758,

//Gem section 4 frame 0
4409 ,4410 ,4411 ,4412 ,4413 ,4414 ,4415 ,4416 ,4417 ,4418 ,4419 ,4420 ,4421 ,4422 ,4423 ,4424 ,4425 ,4426 ,4427 ,4428 ,4429 ,4430 ,4431 ,4432 ,4433 ,4434 ,4435 ,4436 ,4437 ,4438 ,4439 ,4440 ,4441 ,4442 ,4443 ,4444 ,4445 ,4446 ,4447 ,4448 ,4449 ,4450 ,4451 ,4452 ,4611 ,4612 ,4613 ,4614 ,4615 ,4616 ,4617 ,4618 ,4619,

//Gem section 4 frame 1
2794 ,2795 ,2796 ,2797 ,2798 ,2799 ,2800 ,2801 ,2802 ,2803 ,2804 ,2805 ,2806 ,2807 ,2808 ,2809 ,3001 ,3002 ,3003 ,3004 ,3005 ,3006 ,3007 ,3008 ,3009 ,3010 ,3011 ,3012 ,3013 ,3014 ,3015 ,3016 ,3017, 3018, 3037, 3038, 3039 ,3040 ,3041 ,3042 ,3043 ,3044 ,3045 ,3046 ,3047 ,3048 ,3049 ,3050 ,3051 ,3052 ,3053 ,3054 ,3227 ,3228 ,3229 ,3230 ,3231 ,3289 ,3290 ,3291 ,3292 ,3293 ,3294 ,3295 ,3296 ,3297 ,3540 ,3541 ,3542 ,3543 ,3544 ,3545 ,3546
};

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
inline Bool_t IsPadOnCenterFrame(Int_t padID)
{
#ifdef PAD_HELPER_DEBUG
  ValidatePadID(padID, __func__);
#endif
  return std::find(padOnCenterFrame, padOnCenterFrame + NumPadOnCenterFrame, padID)
         != (padOnCenterFrame + NumPadOnCenterFrame);
}

inline Bool_t IsPadOnCenterFrame(Int_t layer, Int_t row)
{
  return IsPadOnCenterFrame(GetPadId(layer, row));
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

// ASAD index 0..30 (same convention as main analyzer TPCPadHelper / HistTools)
static const Int_t NumOfAsadTPC = 31;

//_____________________________________________________________________________
inline Int_t GetASADId(Int_t layer, Int_t row) // 0~30
{
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer, row, __func__);
#endif
  Int_t flag = layer / 4;
  Int_t section;
  if (flag == 0)
    section = 0; // layer 0~3
  else if (flag == 1)
    section = 1; // layer 4~7
  else if (layer == 30 || layer == 31)
    section = 3; // layer 30~31
  else
    section = 2; // layer 8~29

  Int_t half = static_cast<Int_t>(padParameter[layer][kNumOfPad] / 2);
  Int_t division1 = static_cast<Int_t>(padParameter[layer][kNumOfPad] / 6);
  Int_t division2 = static_cast<Int_t>(padParameter[layer][kNumOfPad] * 5 / 6);

  switch (section) {
  case 0:
    if (row < half)
      return 0;
    else
      return 1;
  case 1: {
    const Int_t dummy = (layer % 4 < 2) ? 2 : 5;
    if (row < division1 || division2 <= row)
      return dummy;
    else if (division1 <= row && row < half)
      return dummy + 1;
    else
      return dummy + 2;
  }
  case 3:
    return 30;
  default:
    if (row < half)
      return layer - layer % 2;
    else
      return layer - layer % 2 + 1;
  }
}

//_____________________________________________________________________________
inline Int_t GetCoBoId(Int_t layer, Int_t row)
{
#ifdef PAD_HELPER_DEBUG
  ValidateRow(layer, row, __func__);
#endif
  switch (layer) {
  case 4:
    if (60 <= row && row <= 99)
      return 1;
    else
      return 0;
  case 5:
    if (72 <= row && row <= 119)
      return 1;
    else
      return 0;
  default:
    return layer / 4;
  }
}

//_____________________________________________________________________________
inline void 
InitializeHistograms(TH2Poly* h1)
{
  if (!h1) {
    std::cerr << "[tpc::InitializeHistograms] Error: TH2Poly pointer is null" << std::endl;
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
inline void 
InitializeHistograms(const char* name)
{
  auto h1 = gDirectory->Get<TH2Poly>(name);
  if (!h1) {
    std::cerr << "[tpc::InitializeHistograms] Error: TH2Poly '" << name << "' not found in gDirectory" << std::endl;
    return;
  }
  InitializeHistograms(h1);
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
