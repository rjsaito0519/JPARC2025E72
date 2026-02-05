//_____________________________________________________________________________
/**
 * DstTPCBcOutTracking の出力 ROOT から、
 * TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock (補正前) と
 * TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d (補正後) を描画し、PDF に出力する。
 *
 * オプションで TpcPhase.root も読み込み、RawClock ヒストグラムに
 * フィット曲線を重ね描きする。
 *
 * 2D ヒストが重いため、各ページを PNG 等のラスター形式で出力し、
 * ImageMagick の convert で PDF に埋め込む（軽量 PDF）。
 * convert がなければ PNG のみ出力し、手動で convert *.png out.pdf 可能。
 *
 * Usage:
 *   tpc_phase_plot <bcout.root> [output.pdf] [--phase <TpcPhase.root>] [--vdrift V]
 *
 * オプション:
 *   --phase <file>  : TpcPhase.root を指定してフィット曲線を重ね描き
 *   --vdrift V      : drift velocity [mm/ns] (デフォルト 0.055)
 *
 * 出力 PDF 構成:
 *   Page 1: CoBo RawClock (補正前、2x4) + フィット曲線（--phase 指定時）
 *   Page 2: CoBo 補正後 (2x4) - 補正効果の確認用
 *   Page 3: Asad RawClock (補正前、4x8)
 *   Page 4: Asad 補正後 (4x8) - 補正効果の確認用
 */
//_____________________________________________________________________________

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

//_____________________________________________________________________________
static const Int_t NumOfSegCOBO = 8;
static const Int_t NumOfAsadTPC = 31;

// PNG 解像度（ピクセル）。ラスター PDF 用なので控えめに。
static const Int_t kCoboW = 1200, kCoboH = 600;
static const Int_t kAsadW = 1600, kAsadH = 800;

//_____________________________________________________________________________
static TH2D* GetHist(TFile* f, const TString& hname)
{
  const TString paths[] = { hname, TString("hist/") + hname };
  for (const TString& path : paths) {
    TObject* obj = f->Get(path.Data());
    if (obj) {
      TH2D* h = dynamic_cast<TH2D*>(obj);
      if (h) return h;
    }
  }
  return nullptr;
}

//_____________________________________________________________________________
static void SetupStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleSize(0.06, "t");
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(50);
}

//_____________________________________________________________________________
static void DrawCoboPage(TFile* f, TCanvas* c, Bool_t corrected, 
                         TFile* phase_file = nullptr, Double_t vdrift = 0.055)
{
  c->Clear();
  c->Divide(4, 2, 0.001, 0.001);

  const char* suffix = corrected ? "_Corrected" : "";
  const char* label = corrected ? "Corrected" : "Raw";

  for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
    // _RawClock ヒストグラムを使用（補正前）、または補正後はsuffix無し
    TString hname;
    if (!corrected) {
      hname = Form("TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock", cobo);
    } else {
      hname = Form("TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d", cobo);
    }
    TH2D* h = GetHist(f, hname);
    c->cd(cobo + 1);
    gPad->SetRightMargin(0.12);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    if (h) {
      h->SetTitle(Form("CoBo %d (%s);Clock Time [ns];Residual Y [mm]", cobo, label));
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetLabelSize(0.045);
      h->GetYaxis()->SetLabelSize(0.045);
      h->Draw("colz");
      
      // フィット曲線を重ね描き（RawClockかつphase_file指定時のみ）
      if (!corrected && phase_file) {
        TString phase_name = Form("TpcPhase_Cobo%d", cobo);
        TGraph* g_phase = (TGraph*)phase_file->Get(phase_name);
        if (g_phase) {
          // TGraph (Δclock vs clock) を ResY 空間に変換
          Int_t npoints = g_phase->GetN();
          TGraph* g_resy = new TGraph(npoints);
          for (Int_t i = 0; i < npoints; ++i) {
            Double_t clk, delta_clk;
            g_phase->GetPoint(i, clk, delta_clk);
            Double_t resy = -delta_clk * vdrift;  // ResY = -Δclock * vdrift
            g_resy->SetPoint(i, clk, resy);
          }
          g_resy->SetLineColor(kRed);
          g_resy->SetLineWidth(3);
          g_resy->Draw("L same");
        }
      }
    } else {
      TH2D* dummy = new TH2D(Form("dummy_cobo%d_%s", cobo, label),
                             Form("CoBo %d (%s) - not found", cobo, label),
                             1, -60, 50, 1, -20, 20);
      dummy->Draw();
    }
  }
  c->Update();
}

//_____________________________________________________________________________
static void DrawAsadPage(TFile* f, TCanvas* c, Bool_t corrected)
{
  c->Clear();
  c->Divide(8, 4, 0.001, 0.001);

  const char* suffix = corrected ? "_Corrected" : "";
  const char* label = corrected ? "Corrected" : "Raw";

  for (Int_t asad = 0; asad < NumOfAsadTPC; ++asad) {
    // _RawClock ヒストグラムを使用（補正前）、または補正後はsuffix無し
    TString hname;
    if (!corrected) {
      hname = Form("TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d_RawClock", asad);
    } else {
      hname = Form("TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d", asad);
    }
    TH2D* h = GetHist(f, hname);
    c->cd(asad + 1);
    gPad->SetRightMargin(0.12);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    if (h) {
      h->SetTitle(Form("Asad %d (%s)", asad, label));
      h->GetXaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelSize(0.05);
      h->Draw("colz");
    } else {
      TH2D* dummy = new TH2D(Form("dummy_asad%d_%s", asad, label),
                             Form("Asad %d - not found", asad),
                             1, -60, 50, 1, -20, 20);
      dummy->Draw();
    }
  }
  c->Update();
}

//_____________________________________________________________________________
static int run_convert(const std::vector<std::string>& pngs, const std::string& pdf)
{
  std::ostringstream cmd;
  cmd << "convert";
  for (const auto& p : pngs) cmd << " \"" << p << "\"";
  cmd << " \"" << pdf << "\"";
  return std::system(cmd.str().c_str());
}

//_____________________________________________________________________________
int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <bcout.root> [output.pdf] [--phase <TpcPhase.root>] [--vdrift V]\n"
              << "  Draws 2D hists as PNG, embeds into PDF via 'convert' (ImageMagick).\n"
              << "  Options:\n"
              << "    --phase <file> : overlay fit curves from TpcPhase.root\n"
              << "    --vdrift V     : drift velocity [mm/ns] (default: 0.055)\n"
              << "  If convert is missing, PNGs are kept; merge with:\n"
              << "    convert <base>_page*.png <output.pdf>\n";
    return 1;
  }

  std::string inpath = argv[1];
  std::string outpath;
  std::string phase_path;
  Double_t vdrift = 0.055;
  
  // パース引数
  Int_t argidx = 2;
  while (argidx < argc) {
    std::string arg = argv[argidx];
    if (arg == "--phase" && argidx + 1 < argc) {
      phase_path = argv[++argidx];
      ++argidx;
    } else if (arg == "--vdrift" && argidx + 1 < argc) {
      vdrift = std::atof(argv[++argidx]);
      if (vdrift <= 0.0) vdrift = 0.055;
      ++argidx;
    } else if (arg.find("--") == 0) {
      ++argidx;
    } else {
      outpath = arg;
      ++argidx;
      break;
    }
  }

  if (outpath.empty()) {
    fs::path p(inpath);
    std::string base = p.stem().string();
    outpath = (p.parent_path() / (base + "_TPC_Phase.pdf")).string();
  }

  fs::path out(outpath);
  fs::path outdir = out.parent_path();
  std::string outbase = outdir.empty()
    ? out.stem().string()
    : (outdir / out.stem()).string();

  TFile* f = TFile::Open(inpath.c_str());
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << inpath << std::endl;
    return 1;
  }
  
  // Phase ファイルを開く（オプション）
  TFile* f_phase = nullptr;
  if (!phase_path.empty()) {
    f_phase = TFile::Open(phase_path.c_str());
    if (!f_phase || f_phase->IsZombie()) {
      std::cerr << "Warning: cannot open phase file " << phase_path << ", skipping fit overlay\n";
      f_phase = nullptr;
    } else {
      std::cout << "Phase file loaded: " << phase_path << " (vdrift=" << vdrift << " mm/ns)\n";
    }
  }

  SetupStyle();
  gROOT->SetBatch(kTRUE);

  fs::path tmpdir = fs::temp_directory_path() / ("tpc_phase_plot_" + std::to_string(getpid()));
  std::error_code ec;
  if (!fs::create_directories(tmpdir, ec)) {
    std::cerr << "Error: cannot create temp dir " << tmpdir << std::endl;
    f->Close();
    delete f;
    return 1;
  }

  TCanvas* c_cobo = new TCanvas("c_cobo", "CoBo", kCoboW, kCoboH);
  TCanvas* c_asad = new TCanvas("c_asad", "Asad", kAsadW, kAsadH);

  std::vector<std::string> pngs;
  auto save = [&](const char* tag) {
    std::string png = (tmpdir / (std::string("page_") + tag + ".png")).string();
    pngs.push_back(png);
  };

  // Page 1: CoBo RawClock (補正前) + フィット曲線
  DrawCoboPage(f, c_cobo, false, f_phase, vdrift);
  save("cobo_raw");
  c_cobo->Print(pngs.back().c_str());

  // Page 2: CoBo 補正後（補正効果の確認用）
  DrawCoboPage(f, c_cobo, true, nullptr, vdrift);
  save("cobo_corrected");
  c_cobo->Print(pngs.back().c_str());

  // Page 3: Asad RawClock (補正前)
  DrawAsadPage(f, c_asad, false);
  save("asad_raw");
  c_asad->Print(pngs.back().c_str());

  // Page 4: Asad 補正後（補正効果の確認用）
  DrawAsadPage(f, c_asad, true);
  save("asad_corrected");
  c_asad->Print(pngs.back().c_str());

  delete c_cobo;
  delete c_asad;
  f->Close();
  delete f;
  if (f_phase) {
    f_phase->Close();
    delete f_phase;
  }

  int ret = run_convert(pngs, outpath);
  for (const auto& p : pngs)
    fs::remove(p, ec);

  if (ret != 0) {
    std::cerr << "Warning: 'convert' failed (ImageMagick?). Saving PNGs next to output.\n";
    if (!outdir.empty()) fs::create_directories(outdir, ec);
    TCanvas* c_cobo2 = new TCanvas("c_cobo", "CoBo", kCoboW, kCoboH);
    TCanvas* c_asad2 = new TCanvas("c_asad", "Asad", kAsadW, kAsadH);
    TFile* f2 = TFile::Open(inpath.c_str());
    TFile* f_phase2 = nullptr;
    if (!phase_path.empty()) {
      f_phase2 = TFile::Open(phase_path.c_str());
    }
    if (f2 && !f2->IsZombie()) {
      DrawCoboPage(f2, c_cobo2, false, f_phase2, vdrift);
      c_cobo2->Print((outbase + "_page1_cobo_raw.png").c_str());
      DrawAsadPage(f2, c_asad2, false);
      c_asad2->Print((outbase + "_page2_asad_raw.png").c_str());
      f2->Close();
      delete f2;
    }
    if (f_phase2) {
      f_phase2->Close();
      delete f_phase2;
    }
    delete c_cobo2;
    delete c_asad2;
    std::cerr << "PNGs: " << outbase << "_page{1,2}_*.png\n"
              << "Merge: convert " << outbase << "_page*.png " << outpath << "\n";
    fs::remove_all(tmpdir, ec);
    return 0;
  }

  fs::remove(tmpdir, ec);
  std::cout << "Saved: " << outpath << " (raster PDF via convert)\n";
  return 0;
}
