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
 *   --mode hit|trk  : ヒスト名の優先順（trk: TPCTrk_ResY_vs_ClockTime_... を先に、デフォルト hit）
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
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TPRegexp.h>

#include "paths.h"
#include "ana_helper.h"

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
static const Int_t kCoboW = 1400, kCoboH = 760;
static const Int_t kAsadW = 1900, kAsadH = 980;

//_____________________________________________________________________________
// 検索するヒストグラム名のフォーマット候補 (優先順)
// フォーマット文字列には channel (%d または %02d) と suffix ("_RawClock" か empty) が入る
static const std::vector<std::string> HIST_FMTS_COBO_HIT = {
  "TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d%s",
  "TPCHit_ResY_vs_ClockTime_CoBo%d%s",
  "TPCCl_ResY_vs_ClockTime_CoBo%d%s",
  "TPCTrk_ResY_vs_ClockTime_CoBo%d%s",
};
static const std::vector<std::string> HIST_FMTS_COBO_TRK = {
  "TPCTrk_ResY_vs_ClockTime_CoBo%d%s",
  "TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d%s",
  "TPCHit_ResY_vs_ClockTime_CoBo%d%s",
  "TPCCl_ResY_vs_ClockTime_CoBo%d%s",
};

static const std::vector<std::string> HIST_FMTS_ASAD_HIT = {
  "TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d%s",
  "TPCHit_ResY_vs_ClockTime_Asad%02d%s",
  "TPCCl_ResY_vs_ClockTime_Asad%02d%s",
  "TPCTrk_ResY_vs_ClockTime_Asad%02d%s",
};
static const std::vector<std::string> HIST_FMTS_ASAD_TRK = {
  "TPCTrk_ResY_vs_ClockTime_Asad%02d%s",
  "TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d%s",
  "TPCHit_ResY_vs_ClockTime_Asad%02d%s",
  "TPCCl_ResY_vs_ClockTime_Asad%02d%s",
};

static const std::vector<std::string>& HistFmtsCobo(const TString& mode)
{
  if (mode == "trk" || mode == "track" || mode == "TPCTrk")
    return HIST_FMTS_COBO_TRK;
  return HIST_FMTS_COBO_HIT;
}
static const std::vector<std::string>& HistFmtsAsad(const TString& mode)
{
  if (mode == "trk" || mode == "track" || mode == "TPCTrk")
    return HIST_FMTS_ASAD_TRK;
  return HIST_FMTS_ASAD_HIT;
}
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
// TGraphErrors (clk, dclk[ns]) -> TGraph (clk, ResY[mm]); max_points<=0 は全点
static TGraph* MakeResyGraphFromDclock(TGraphErrors* gsrc, Double_t vdrift, Int_t max_points)
{
  if (!gsrc) return nullptr;
  const Int_t n = gsrc->GetN();
  if (n <= 0) return nullptr;
  Int_t step = 1;
  if (max_points > 0 && n > max_points)
    step = (n + max_points - 1) / max_points;
  Int_t nout = 0;
  for (Int_t i = 0; i < n; i += step) ++nout;
  TGraph* g = new TGraph(nout);
  Int_t j = 0;
  for (Int_t i = 0; i < n; i += step) {
    Double_t clk = 0.0, dclk = 0.0;
    gsrc->GetPoint(i, clk, dclk);
    g->SetPoint(j++, clk, -dclk * vdrift);
  }
  return g;
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
                         const std::vector<std::string>& hist_fmts_cobo,
                         TFile* phase_file = nullptr, Double_t vdrift = 0.055)
{
  c->Clear();
  c->Divide(4, 2, 0.001, 0.001);

  const char* suffix = corrected ? "_Corrected" : "";
  const char* label = corrected ? "Corrected" : "Raw";

  for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
    // 候補フォーマットを順に試す
    TH2D* h = nullptr;
    TString found_name = "";
    
    // suffix: RawClock or empty (corrected)
    const char* hist_suffix = !corrected ? "_RawClock" : "";

    for (const auto& fmt : hist_fmts_cobo) {
      TString hname = Form(fmt.c_str(), cobo, hist_suffix);
      h = GetHist(f, hname);
      if (h) {
        found_name = hname;
        break;
      }
    }

    c->cd(cobo + 1);
    gPad->SetRightMargin(0.18);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);
    if (h) {
      h->SetTitle(Form("CoBo %d (%s);Clock Time [ns];Residual Y [mm]", cobo, label));
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetLabelSize(0.045);
      h->GetYaxis()->SetLabelSize(0.045);
      h->Draw("colz");
      
      // フィット曲線と補助線を重ね描き（RawClockかつphase_file指定時のみ）
      if (!corrected && phase_file) {
        TString phase_name = Form("TpcPhase_Cobo%d", cobo);
        TGraph* g_phase = (TGraph*)phase_file->Get(phase_name);
        TF1* f_phase = (TF1*)phase_file->Get(Form("TpcPhase_Fit_Cobo%d", cobo));
        TF1* f_init  = (TF1*)phase_file->Get(Form("TpcPhase_Init_Cobo%d", cobo));

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
          g_resy->DrawClone("L same");
          delete g_resy;
        }

        // --- Profile overlay (fit デバッグ): 全 profile は細い折れ線のみ（点マーカーなし）
        TGraphErrors* g_profile = (TGraphErrors*)phase_file->Get(Form("TpcPhase_Profile_Cobo%d", cobo));
        if (g_profile) {
          TGraph* g_resy = MakeResyGraphFromDclock(g_profile, vdrift, 0);
          if (g_resy) {
            g_resy->SetLineColor(kAzure + 2);
            g_resy->SetLineWidth(1);
            g_resy->SetLineStyle(1);
            g_resy->DrawClone("L same");
            delete g_resy;
          }
        }

        // Fit-used core: プロファイル（青線）と同じく全点を描く（間引きしない）
        TGraphErrors* g_fitused = (TGraphErrors*)phase_file->Get(Form("TpcPhase_ProfileFitUsed_Cobo%d", cobo));
        if (g_fitused) {
          TGraph* g_resy = MakeResyGraphFromDclock(g_fitused, vdrift, 0);
          if (g_resy) {
            g_resy->SetMarkerStyle(24);
            g_resy->SetMarkerSize(0.75);
            g_resy->SetMarkerColor(kOrange + 7);
            g_resy->SetLineColor(kOrange + 7);
            g_resy->DrawClone("P same");
            delete g_resy;
          }
        }

        if (f_phase) {
          // フィット範囲と t0(ステップ位置) の縦線を描画
          const Double_t fxmin = f_phase->GetXmin();
          const Double_t fxmax = f_phase->GetXmax();
          const Int_t npar = f_phase->GetNpar();
          const Double_t t0_fit = (npar >= 3) ? f_phase->GetParameter(2) : 0.0;

          const Double_t y_min = h->GetYaxis()->GetXmin();
          const Double_t y_max = h->GetYaxis()->GetXmax();

          // フィットに使った X 範囲
          TLine* l_xmin = new TLine(fxmin, y_min, fxmin, y_max);
          TLine* l_xmax = new TLine(fxmax, y_min, fxmax, y_max);
          l_xmin->SetLineColor(kBlue);
          l_xmax->SetLineColor(kBlue);
          l_xmin->SetLineStyle(3);
          l_xmax->SetLineStyle(3);
          l_xmin->Draw("same");
          l_xmax->Draw("same");

          // フィット後のステップ位置 t0
          TLine* l_t0 = new TLine(t0_fit, y_min, t0_fit, y_max);
          l_t0->SetLineColor(kGreen+2);
          l_t0->SetLineStyle(2);
          l_t0->Draw("same");
        }

        if (f_init) {
          // 初期ステップ位置 t0_initial と初期ステップ関数のイメージを描画
          const Int_t npar_i = f_init->GetNpar();
          const Double_t t0_init = (npar_i >= 3) ? f_init->GetParameter(2) : 0.0;
          const Double_t y_min = h->GetYaxis()->GetXmin();
          const Double_t y_max = h->GetYaxis()->GetXmax();

          TLine* l_t0_init = new TLine(t0_init, y_min, t0_init, y_max);
          l_t0_init->SetLineColor(kMagenta+1);
          l_t0_init->SetLineStyle(7);
          l_t0_init->Draw("same");

          // 初期ステップ関数を Δclock→ResY に変換して細い線で描画
          const Double_t x1 = f_init->GetXmin();
          const Double_t x2 = f_init->GetXmax();
          const Int_t nseg = 200;
          TGraph* g_init_resy = new TGraph(nseg);
          for (Int_t i = 0; i < nseg; ++i) {
            const Double_t x = x1 + (x2 - x1) * (Double_t)i / (Double_t)(nseg - 1);
            const Double_t dclk = f_init->Eval(x);
            const Double_t resy = -dclk * vdrift;
            g_init_resy->SetPoint(i, x, resy);
          }
          g_init_resy->SetLineColor(kMagenta+1);
          g_init_resy->SetLineStyle(7);
          g_init_resy->SetLineWidth(2);
          g_init_resy->Draw("L same");
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
static void DrawAsadPage(TFile* f, TCanvas* c, Bool_t corrected,
                         const std::vector<std::string>& hist_fmts_asad)
{
  c->Clear();
  c->Divide(8, 4, 0.001, 0.001);

  const char* suffix = corrected ? "_Corrected" : "";
  const char* label = corrected ? "Corrected" : "Raw";

  for (Int_t asad = 0; asad < NumOfAsadTPC; ++asad) {
    // 候補フォーマットを順に試す
    TH2D* h = nullptr;
    TString found_name = "";

    // suffix: RawClock or empty (corrected)
    const char* hist_suffix = !corrected ? "_RawClock" : "";

    for (const auto& fmt : hist_fmts_asad) {
      TString hname = Form(fmt.c_str(), asad, hist_suffix);
      h = GetHist(f, hname);
      if (h) {
        found_name = hname;
        break;
      }
    }
    
    c->cd(asad + 1);
    gPad->SetRightMargin(0.18);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);
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
    std::cerr << "Usage: " << argv[0]
              << " <bcout.root> [output.pdf] [--phase <TpcPhase.root>] [--vdrift V] [--mode hit|trk]\n"
              << "  Draws 2D hists as PNG, embeds into PDF via 'convert' (ImageMagick).\n"
              << "  Options:\n"
              << "    --phase <file> : overlay fit curves from TpcPhase.root\n"
              << "    --vdrift V     : drift velocity [mm/ns] (default: 0.055)\n"
              << "    --mode hit|trk : histogram name priority (default: hit)\n"
              << "  If convert is missing, PNGs are kept; merge with:\n"
              << "    convert <base>_page*.png <output.pdf>\n";
    return 1;
  }

  std::string inpath = argv[1];
  std::string outpath;
  std::string phase_path;
  Double_t vdrift = 0.055;
  TString mode = "hit";

  // パース引数
  Int_t argidx = 2;
  while (argidx < argc) {
    std::string arg = argv[argidx];
    if (arg == "--phase" && argidx + 1 < argc) {
      phase_path = argv[++argidx];
      ++argidx;
    } else if (arg == "--mode" && argidx + 1 < argc) {
      mode = argv[++argidx];
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
    Int_t run_number = -1;
    TFile* f_tmp = TFile::Open(inpath.c_str());
    if (f_tmp && !f_tmp->IsZombie()) {
      TTree *tree = (TTree*)f_tmp->Get("tpc");
      if (tree) {
        try {
          TTreeReader reader(tree);
          TTreeReaderValue<UInt_t> rv(reader, "run_number");
          if (reader.Next()) {
            run_number = *rv;
            std::cout << "Auto-detected run number from TTree (TTreeReader): " << run_number << std::endl;
          }
        } catch (...) {
          // Ignore errors and try fallback
        }
      }
      f_tmp->Close();
      delete f_tmp;
    }

    if (run_number < 0) {
      TString filename = inpath.c_str();
      TPRegexp run_regex("run([0-9]+)");
      TString run_str = filename(run_regex);
      if (!run_str.IsNull()) {
        run_str.ReplaceAll("run", "");
        run_number = run_str.Atoi();
        std::cout << "Auto-detected run number from filename: " << run_number << std::endl;
      }
    }

    fs::path p(inpath);
    std::string base = p.stem().string();
    if (run_number >= 0) {
      TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
      outpath = std::string(img_dir.Data()) + "/" + base + "_TPC_Phase.pdf";
    } else {
      outpath = (p.parent_path() / (base + "_TPC_Phase.pdf")).string();
    }
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

  const std::vector<std::string>& hist_fmts_cobo = HistFmtsCobo(mode);
  const std::vector<std::string>& hist_fmts_asad = HistFmtsAsad(mode);
  std::cout << "Histogram mode: " << mode << std::endl;

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
  DrawCoboPage(f, c_cobo, false, hist_fmts_cobo, f_phase, vdrift);
  save("cobo_raw");
  c_cobo->Print(pngs.back().c_str());

  // Page 2: CoBo 補正後（補正効果の確認用）
  DrawCoboPage(f, c_cobo, true, hist_fmts_cobo, nullptr, vdrift);
  save("cobo_corrected");
  c_cobo->Print(pngs.back().c_str());

  // Page 3: Asad RawClock (補正前)
  DrawAsadPage(f, c_asad, false, hist_fmts_asad);
  save("asad_raw");
  c_asad->Print(pngs.back().c_str());

  // Page 4: Asad 補正後（補正効果の確認用）
  DrawAsadPage(f, c_asad, true, hist_fmts_asad);
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
      DrawCoboPage(f2, c_cobo2, false, hist_fmts_cobo, f_phase2, vdrift);
      c_cobo2->Print((outbase + "_page1_cobo_raw.png").c_str());
      DrawAsadPage(f2, c_asad2, false, hist_fmts_asad);
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
