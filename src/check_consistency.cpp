#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>

// ROOT Headers
#include "TVector3.h"
#include "TRandom3.h"
#include "Math/Rotation3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Vector3D.h"

// ---------------------------------------------------------
// 1. オリジナルの実装 (Manual Calculation - Inverse)
// ---------------------------------------------------------
TVector3 Local2GlobalPos_Original(TVector3 LocalPos, TVector3 DetectorGlobalPos,
                                  Double_t tilt, Double_t RA1, Double_t RA2) {
    Double_t ct0 = TMath::Cos(tilt*TMath::DegToRad());
    Double_t st0 = TMath::Sin(tilt*TMath::DegToRad());
    Double_t ct1 = TMath::Cos(RA1*TMath::DegToRad());
    Double_t st1 = TMath::Sin(RA1*TMath::DegToRad());
    Double_t ct2 = TMath::Cos(RA2*TMath::DegToRad());
    Double_t st2 = TMath::Sin(RA2*TMath::DegToRad());

    Double_t dxds =  ct0*ct2+st0*st1*st2;
    Double_t dxdt = -st0*ct2+ct0*st1*st2;
    Double_t dxdu =  ct1*st2;

    Double_t dyds =  st0*ct1;
    Double_t dydt =  ct0*ct1;
    Double_t dydu = -st1;

    Double_t dzds = -ct0*st2+st0*st1*ct2;
    Double_t dzdt =  st0*st2+ct0*st1*ct2;
    Double_t dzdu =  ct1*ct2;

    // 行列の転置（M^T）と平行移動の逆操作
    Double_t x = dxds*LocalPos.x() + dxdt*LocalPos.y() + dxdu*LocalPos.z() + DetectorGlobalPos.x();
    Double_t y = dyds*LocalPos.x() + dydt*LocalPos.y() + dydu*LocalPos.z() + DetectorGlobalPos.y();
    Double_t z = dzds*LocalPos.x() + dzdt*LocalPos.y() + dzdu*LocalPos.z() + DetectorGlobalPos.z();

    return TVector3(x, y, z);
}

// ---------------------------------------------------------
// 2. 新しい実装 (ROOT::Math::GenVector - Inverse)
// ---------------------------------------------------------
TVector3 Local2GlobalPos_New(TVector3 LocalPos, TVector3 DetectorGlobalPos,
                             Double_t tilt, Double_t RA1, Double_t RA2) {

    constexpr double DegToRad = TMath::DegToRad();

    // 1. Input -> XYZVector
    ROOT::Math::XYZVector localVec(LocalPos.X(), LocalPos.Y(), LocalPos.Z());

    // 2. Rotation Matrices (Positive angles for inverse)
    ROOT::Math::RotationY rY(RA2 * DegToRad);
    ROOT::Math::RotationX rX(RA1 * DegToRad);
    ROOT::Math::RotationZ rZ(tilt * DegToRad);

    // 3. Compose Matrix (Inverse Order: Ry * Rx * Rz)
    // キャストを入れてリンクエラーを回避
    ROOT::Math::Rotation3D rotMatrix = ROOT::Math::Rotation3D(rY) 
                                     * ROOT::Math::Rotation3D(rX) 
                                     * ROOT::Math::Rotation3D(rZ);
    
    ROOT::Math::XYZVector rotatedVec = rotMatrix * localVec;

    // 4. Translation
    ROOT::Math::XYZVector globalVec(
        rotatedVec.X() + DetectorGlobalPos.X(),
        rotatedVec.Y() + DetectorGlobalPos.Y(),
        rotatedVec.Z() + DetectorGlobalPos.Z()
    );

    return TVector3(globalVec.X(), globalVec.Y(), globalVec.Z());
}

// ---------------------------------------------------------
// Main Function
// ---------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <x_local> <y_local> <z_local>" << std::endl;
        return 1;
    }

    // コマンドライン引数 (Local座標)
    Double_t inX = std::atof(argv[1]);
    Double_t inY = std::atof(argv[2]);
    Double_t inZ = std::atof(argv[3]);
    TVector3 localPos(inX, inY, inZ);

    TRandom3 randGen(0); 

    // ランダムなパラメータ設定
    TVector3 detectorPos(
        randGen.Uniform(-100.0, 100.0), 
        randGen.Uniform(-100.0, 100.0), 
        randGen.Uniform(0.0, 2000.0)
    );
    Double_t tilt = randGen.Uniform(-30.0, 30.0);
    Double_t RA1  = randGen.Uniform(-30.0, 30.0);
    Double_t RA2  = randGen.Uniform(-30.0, 30.0);

    // 実行
    TVector3 resOrig = Local2GlobalPos_Original(localPos, detectorPos, tilt, RA1, RA2);
    TVector3 resNew  = Local2GlobalPos_New(localPos, detectorPos, tilt, RA1, RA2);

    // 結果比較
    std::cout << "========================================" << std::endl;
    std::cout << "Consistency Check: Local -> Global" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::cout << "[Inputs (Local)]" << std::endl;
    std::cout << " Local Pos    : (" << inX << ", " << inY << ", " << inZ << ")" << std::endl;
    std::cout << " Detector Pos : (" << detectorPos.X() << ", " << detectorPos.Y() << ", " << detectorPos.Z() << ")" << std::endl;
    std::cout << " Angles [deg] : Tilt=" << tilt << ", RA1=" << RA1 << ", RA2=" << RA2 << std::endl;
    
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::fixed << std::setprecision(10);

    std::cout << "[Result: Original]" << std::endl;
    std::cout << " X: " << resOrig.X() << std::endl;
    std::cout << " Y: " << resOrig.Y() << std::endl;
    std::cout << " Z: " << resOrig.Z() << std::endl;

    std::cout << "[Result: New]" << std::endl;
    std::cout << " X: " << resNew.X() << std::endl;
    std::cout << " Y: " << resNew.Y() << std::endl;
    std::cout << " Z: " << resNew.Z() << std::endl;

    TVector3 diff = resOrig - resNew;
    Double_t mag = diff.Mag();

    std::cout << "----------------------------------------" << std::endl;
    std::cout << " Difference Magnitude: " << mag << std::endl;

    if (mag < 1e-9) {
        std::cout << "\n>>> SUCCESS: Results are consistent! <<<" << std::endl;
    } else {
        std::cout << "\n>>> WARNING: Results differ! <<<" << std::endl;
    }

    return 0;
}