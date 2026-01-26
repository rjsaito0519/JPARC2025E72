//_____________________________________________________________________________
/**
 * DstTPCTracking の TTree 'tpc' を読んで、CoBo ごとに
 * Residual Y vs Clock Time の TProfile を描画し PDF に出力する。
 *
 * イベントを先頭から --chunk 個ずつ区切り、各区間のプロファイルを
 * Draw("HIST SAME") で重ね書きし、区間ごとに色を変える。凡例で色と
 * チャンク（ev 範囲）の対応を表示。巻末は同じ CoBo でまとめ、各 CoBo 1 ページに
 * その CoBo の全チャンクの 2D ヒスト（Clock Time vs Res Y、COLZ）を並べる。
 *
 * Usage:
 *   tpc_resy_vs_cobotime_from_tree <rootfile> [output.pdf] [--chunk N] [--no-appendix]
 *   --no-appendix : 巻末の 2D ヒストを出さない（PDF 軽量化）
 *
 * Ref: TPCPadHelper.hh GetCoBoId, DstTPCTracking tree
 */
//_____________________________________________________________________________

#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TProfile.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TROOT.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "progress_bar.h"

//_____________________________________________________________________________
// CoBo mapping (TPCPadHelper.hh GetCoBoId)
static const int NumOfSegCOBO = 8;

// プログレスバー更新間隔（定数 N エントリごとに表示）
static const Long64_t kProgressStep = 10000;

// イベント区切り数（デフォルト。--chunk で上書き）
static const Long64_t kChunkSizeDefault = 1000;

// 区間ごとの描画色（重ね書きで区別）
static const int kChunkColors[] = {
  kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kCyan+1, kOrange+7, kYellow+2,
  kPink+1, kSpring+2, kTeal+1, kAzure+2, kViolet+1, kGray+2
};
static const int kChunkColorsN = sizeof(kChunkColors) / sizeof(kChunkColors[0]);

static Int_t GetCoBoId(Int_t layer, Int_t row)
{
  switch (layer) {
  case 4:
    return (60 <= row && row <= 99) ? 1 : 0;
  case 5:
    return (72 <= row && row <= 119) ? 1 : 0;
  default:
    return layer / 4;
  }
}

// DstTPCTracking / HistTools の TPC_ResidualY_vs_ClockTime_CoBo と同一 (X はそのまま)
// (HistTools trk_bins_2d_clocktime_resy: 200,-60,50, 200,-50,50)
static const Double_t TMIN = -60.0,  TMAX = 50.0;   // Clock Time [ns]
static const Double_t YMIN = -10.0,  YMAX = 10.0;   // Residual Y [mm] 描画範囲
static const Int_t   TBINS = 200;                   // TProfile の X ビン数
static const Int_t   TBINS2D = 100,  YBINS2D = 100; // 巻末 2D のビン数（PDF 軽量化のため粗め）

//_____________________________________________________________________________
int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <rootfile> [output.pdf] [--chunk N] [--no-appendix]\n";
    return 1;
  }

  std::string rootpath = argv[1];
  std::string outpath;
  Long64_t chunk_size = kChunkSizeDefault;
  Bool_t no_appendix = kFALSE;

  for (int i = 2; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--chunk" && i + 1 < argc) {
      chunk_size = std::atol(argv[++i]);
      if (chunk_size < 1) chunk_size = 1;
    } else if (a == "--no-appendix") {
      no_appendix = kTRUE;
    } else if (a.find("--") != 0) {
      outpath = a;
    }
  }

  if (outpath.empty()) {
    std::string base = rootpath;
    size_t s = base.find_last_of("/\\");
    if (s != std::string::npos) base = base.substr(s + 1);
    s = base.find_last_of('.');
    if (s != std::string::npos) base = base.substr(0, s);
    outpath = base + "_TPC_ResY_vs_CoBoTime_from_tree.pdf";
  }

  TFile* f = TFile::Open(rootpath.c_str());
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << rootpath << std::endl;
    return 1;
  }

  TTree* tree = (TTree*)f->Get("tpc");
  if (!tree) {
    std::cerr << "Error: tree 'tpc' not found in " << rootpath << std::endl;
    f->Close();
    return 1;
  }

  // TTreeReader (DstTPCTracking tree structure)
  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> rv_ntTpc(reader, "ntTpc");
  TTreeReaderValue<std::vector<Int_t>> rv_nhtrack(reader, "nhtrack");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_hitlayer(reader, "hitlayer");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_track_cluster_row_center(reader, "track_cluster_row_center");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_residual_y(reader, "residual_y");
  TTreeReaderValue<std::vector<Double_t>> rv_clkTpc(reader, "clkTpc");

  Long64_t nent = tree->GetEntries();
  Long64_t n_chunks = (nent + chunk_size - 1) / chunk_size;
  if (n_chunks < 1) n_chunks = 1;

  // TProfile per (chunk, CoBo)。区間ごとに Fill して後で Draw SAME で重ねる
  std::vector<std::vector<TProfile*>> prof(n_chunks, std::vector<TProfile*>(NumOfSegCOBO, nullptr));
  // TH2D per (chunk, CoBo)。巻末で同じ CoBo のチャンクを COLZ で並べる（--no-appendix で省略可）
  std::vector<std::vector<TH2D*>> hist2d;
  if (!no_appendix) hist2d.resize(n_chunks, std::vector<TH2D*>(NumOfSegCOBO, nullptr));
  for (Long64_t ch = 0; ch < n_chunks; ++ch) {
    Long64_t ev_lo = ch * chunk_size + 1;
    Long64_t ev_hi = (ch + 1) * chunk_size;
    if (ev_hi > nent) ev_hi = nent;
    for (int c = 0; c < NumOfSegCOBO; ++c) {
      prof[ch][c] = new TProfile(Form("prof_CoBo%d_chunk%lld", c, (long long)ch),
                                 Form("CoBo %d;Clock Time [ns];Residual Y [mm]", c),
                                 TBINS, TMIN, TMAX, YMIN, YMAX, "s");
      if (!no_appendix)
        hist2d[ch][c] = new TH2D(Form("h2_CoBo%d_chunk%lld", c, (long long)ch),
                                 Form("Chunk %lld (ev %lld-%lld);Clock Time [ns];Residual Y [mm]", (long long)ch, (long long)ev_lo, (long long)ev_hi),
                                 TBINS2D, TMIN, TMAX, YBINS2D, YMIN, YMAX);
    }
  }

  Long64_t nfilled = 0, nskip = 0;
  Long64_t ie = 0;

  while (reader.Next()) {
    ++ie;
    if (ie % kProgressStep == 0 || ie == nent)
      displayProgressBar(static_cast<Int_t>(ie), static_cast<Int_t>(nent));

    Int_t ntTpc = *rv_ntTpc;
    const auto& nhtrack = *rv_nhtrack;
    const auto& hitlayer = *rv_hitlayer;
    const auto& track_cluster_row_center = *rv_track_cluster_row_center;
    const auto& residual_y = *rv_residual_y;
    const auto& clkTpc = *rv_clkTpc;

    if ((Int_t)clkTpc.size() < NumOfSegCOBO) { ++nskip; continue; }

    for (Int_t it = 0; it < ntTpc; ++it) {
      Int_t nht = (it < (Int_t)nhtrack.size()) ? nhtrack[it] : 0;
      if (it >= (Int_t)hitlayer.size() || it >= (Int_t)track_cluster_row_center.size() || it >= (Int_t)residual_y.size())
        continue;

      for (Int_t ih = 0; ih < nht; ++ih) {
        if (ih >= (Int_t)hitlayer[it].size() || ih >= (Int_t)track_cluster_row_center[it].size() || ih >= (Int_t)residual_y[it].size())
          continue;

        Double_t ly = hitlayer[it][ih];
        Double_t rc = track_cluster_row_center[it][ih];
        Double_t ry = residual_y[it][ih];
        Int_t layer = (Int_t)std::round(ly);
        Int_t row   = (Int_t)std::round(rc);

        Int_t cobo = GetCoBoId(layer, row);
        if (cobo < 0 || cobo >= NumOfSegCOBO) { ++nskip; continue; }
        if (!TMath::Finite(clkTpc[cobo])) { ++nskip; continue; }

        Double_t t = clkTpc[cobo];
        if (!TMath::Finite(t)) { ++nskip; continue; }
        if (!TMath::Finite(ry)) { ++nskip; continue; }

        Long64_t chunk_id = (ie - 1) / chunk_size;
        if (chunk_id >= n_chunks) chunk_id = n_chunks - 1;
        prof[chunk_id][cobo]->Fill(t, ry);
        if (!no_appendix) hist2d[chunk_id][cobo]->Fill(t, ry);
        ++nfilled;
      }
    }
  }

  std::cout << "Filled " << nfilled << " hits, skipped " << nskip << ", entries " << nent
            << ", chunks " << n_chunks << " (size " << chunk_size << ")" << std::endl;

  // Draw: 各 CoBo ごとに 1 ページ。同じ CoBo ならチャンクを Draw SAME で重ね、チャンクごとに色分け
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);
  TCanvas* c1 = new TCanvas("c1", "TPC ResY vs CoBo Clock Time from tree", 800, 600);

  c1->Print((outpath + "(").c_str());   // PDF を開く

  for (int c = 0; c < NumOfSegCOBO; ++c) {
    c1->Clear();
    for (Long64_t ch = 0; ch < n_chunks; ++ch) {
      prof[ch][c]->SetLineColor(kChunkColors[ch % kChunkColorsN]);
      prof[ch][c]->SetLineWidth(1);
      prof[ch][c]->GetYaxis()->SetRangeUser(YMIN, YMAX);
      if (ch == 0)
        prof[ch][c]->Draw("HIST");
      else
        prof[ch][c]->Draw("HIST SAME");
    }
    // 凡例: 色とチャンクの対応（小さめ・右上に寄せる）
    TLegend* leg = new TLegend(0.78, 0.80, 0.97, 0.96);
    leg->SetTextSize(0.022);
    for (Long64_t ch = 0; ch < n_chunks; ++ch) {
      Long64_t ev_lo = ch * chunk_size + 1;
      Long64_t ev_hi = (ch + 1) * chunk_size;
      if (ev_hi > nent) ev_hi = nent;
      leg->AddEntry(prof[ch][c], Form("Chunk %lld (ev %lld-%lld)", (long long)ch, (long long)ev_lo, (long long)ev_hi), "L");
    }
    leg->Draw();
    c1->Print(outpath.c_str());         // 1 ページ = 1 CoBo（チャンクは重ね描き＋凡例）
  }

  // 巻末: 同じ CoBo でまとめ、1 CoBo 1 ページに全チャンクの 2D ヒスト（COLZ）を並べる
  if (!no_appendix) {
    Int_t ncol = (Int_t)std::ceil(std::sqrt((Double_t)n_chunks));
    if (ncol < 1) ncol = 1;
    Int_t nrow = (Int_t)std::ceil((Double_t)n_chunks / (Double_t)ncol);
    if (nrow < 1) nrow = 1;
    const Int_t kPadPx = 260;  // 1 パッドのピクセル（1:1 になる）
    Int_t cw = ncol * kPadPx, ch_h = nrow * kPadPx;
    for (int c = 0; c < NumOfSegCOBO; ++c) {
      c1->Clear();
      c1->SetCanvasSize(cw, ch_h);
      c1->Divide(ncol, nrow);
      for (Long64_t ch = 0; ch < n_chunks; ++ch) {
        c1->cd((Int_t)ch + 1);
        hist2d[ch][c]->Draw("COLZ");
      }
      c1->Print(outpath.c_str());
    }
  }

  c1->Print((outpath + ")").c_str());   // PDF を閉じる

  std::cout << "Saved: " << outpath << std::endl;

  // Cleanup
  for (Long64_t ch = 0; ch < n_chunks; ++ch) {
    for (int c = 0; c < NumOfSegCOBO; ++c) {
      delete prof[ch][c];
      if (!no_appendix) delete hist2d[ch][c];
    }
  }
  delete c1;
  f->Close();
  delete f;

  return 0;
}
