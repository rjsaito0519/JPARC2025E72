// c++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <map>
#include <limits>
#include <initializer_list>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TLine.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"

Config& conf = Config::getInstance();

// HTOF calibration in one pass (TENTATIVE): TDC, ADC (pedestal/MIP), PHC and the
// absolute TOF offset.
//
// Inputs: one or more (Hodo, DST) pairs; several pairs merge runs (the run numbers of each pair
// must agree). Output is labeled with the first pair's run number (representative run).
//   Hodo = run{run}_Hodo.root (UserHodoscope), DST = run{run}_TPCHelixHTOF.root (DstTPCHelixHTOF)
//
// - TDC : all-event HTOF_TDC_seg{i}{U|D|S} histograms of Hodo (summed over runs), ana_helper::tdc_fit.
// - ADC : pedestal from the all-event HTOF_ADC_seg{i}{U|D|S} histograms of Hodo; MIP from the ADC of
//         HTOF-matched TPC tracks with the dE/dx pion bit set and the proton bit not set (kaon-
//         ambiguous kept; beam/accidental tracks are not excluded), ana_helper::htof_adc_fit.
//         Sample: vertex tracks + no-vertex tracks with extrap_seg == match_cl_seg (dseg = 0).
//         Beam-window short segments, weak side (seg1 D, seg2 U, seg3 D, seg4 U): refitted on tracks
//         whose strong side is near its MIP with a spike + hump model, ana_helper::htof_adc_fit_weak.
// - PHC : U/D pairs from the Hodo "hodo" tree; dE recomputed from the raw ADC with the ADC fit above
//         (dE = (adc - ped)/(mip - ped)), so no re-decoding is needed in between; times are the
//         decode-time pre-PHC htof_time_{u|d}. U is fitted on Y_U = (time0 - t_U) + (time0 - ct_D)
//         vs dE_U (mean time: the hit-position correlation cancels), D likewise, alternating; p2 of
//         each side re-centred per channel (the mean time only fixes the U+D sum). Two stages:
//         all hits, then the pion-tagged selection with p1 fixed to stage 1. ana_helper::htof_phc_fit.
//         Function and parameters as the analyzer's HodoPHCParam Type1: ctime = time - p0/sqrt|dE-p1| + p2.
// - TOF offset: common shift of p2 (U and D) per segment so that
//         dt_pi = (ct_mean - time0) - (t_beam + tof_calc_pi) has median 0 for non-beam pion-tagged
//         tracks with a valid TOF and p_vtx >= tof_offset_min_p in the DST (expectation independent
//         of HTOF parameters). Low-momentum pion-tagged tracks contain electrons, which arrive early
//         and would pull the offset; above the cut the e/pi TOF difference is below the resolution.
//
// Outputs (read by update_param.py htofprm / htofphc):
//   run{run}_HTOF_HDPRM_All.root : tdc_p0_val/err
//   run{run}_HTOF_ADC_Pi.root    : adc_p0/p1_val/err, adc_flag/chi2/ndf, raw/selected histograms
//   run{run}_HTOF_PHC_Pi.root    : p0/p1/p2_val/err (final), p*_all_val (stage 1), fit_lo/hi,
//                                  tof_offset(_n), Y histograms, dt_pi after correction
//   run{run}_HTOF_Calib_Pi.pdf   : one page per segment, rows U/D/S, columns
//     ADC raw | ADC selected | ADC close-up | TDC | PHC stage 1 | PHC final | dt_pi after correction
//     (S row, PHC columns: value summary)
//
// Exit status: 0 on success, 1 on any input/processing error (run_hodo.py relies on it so that
// update_param.py never runs on stale output files).

static TFile* open_output(Int_t run_num, const char* name) {
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString path = Form("%s/run%05d_HTOF_%s.root", root_base_dir.Data(), run_num, name);
    return new TFile(path.Data(), "RECREATE");
}

// True if every reader value of r is set up (branch exists with the expected type).
// A missing/renamed branch would otherwise make the loop read nothing and the tool
// "succeed" with an empty sample. Rewinds the reader so the caller's Next() loop
// starts from the first entry.
static Bool_t all_branches_found(TTreeReader& r,
                                 std::initializer_list<ROOT::Internal::TTreeReaderValueBase*> vals,
                                 const TString& file) {
    r.SetEntry(0);   // triggers the branch setup
    Bool_t ok = kTRUE;
    for (auto* v : vals) {
        if (v->GetSetupStatus() < 0) {
            std::cerr << "Error: branch '" << v->GetBranchName() << "' missing or of unexpected type in "
                      << file << " (tree " << r.GetTree()->GetName() << ")" << std::endl;
            ok = kFALSE;
        }
    }
    r.Restart();
    return ok;
}

Bool_t analyze(std::vector<TString> hodo_paths, std::vector<TString> dst_paths, TString particle) {
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0/256, 160.0/256, 44.0/256);

    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x");
    gStyle->SetTitleSize(0.06, "y");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // +------------+
    // | load files |
    // +------------+
    if (hodo_paths.empty() || hodo_paths.size() != dst_paths.size()) {
        std::cerr << "Error: hodo_paths/dst_paths must be non-empty and of equal length" << std::endl;
        return kFALSE;
    }
    const std::size_t n_runs = hodo_paths.size();

    std::vector<TFile*> f_hodo_list(n_runs), f_dst_list(n_runs);
    for (std::size_t k = 0; k < n_runs; k++) {
        f_hodo_list[k] = new TFile(hodo_paths[k].Data());
        if (!f_hodo_list[k] || f_hodo_list[k]->IsZombie()) {
            std::cerr << "Error: Could not open file : " << hodo_paths[k] << std::endl;
            return kFALSE;
        }
        f_dst_list[k] = new TFile(dst_paths[k].Data());
        if (!f_dst_list[k] || f_dst_list[k]->IsZombie()) {
            std::cerr << "Error: Could not open file : " << dst_paths[k] << std::endl;
            return kFALSE;
        }
    }

    // Representative run number for output naming/dirs = first pair. The Hodo and DST of each
    // pair must be the same run (they are joined by event_number below).
    Int_t run_num = -1;
    for (std::size_t k = 0; k < n_runs; k++) {
        TTreeReader r_dst("htof", f_dst_list[k]);
        TTreeReaderValue<unsigned int> rn_dst(r_dst, "run_number");
        TTreeReader r_hodo("hodo", f_hodo_list[k]);
        TTreeReaderValue<unsigned int> rn_hodo(r_hodo, "run_number");
        if (r_dst.SetEntry(0) != TTreeReader::kEntryValid || r_hodo.SetEntry(0) != TTreeReader::kEntryValid) {
            std::cerr << "Error: cannot read run_number from " << hodo_paths[k] << " / " << dst_paths[k] << std::endl;
            return kFALSE;
        }
        std::cout << "#D HTOF_Calib: input " << k << ": run" << *rn_dst
                   << " (" << hodo_paths[k] << ", " << dst_paths[k] << ")" << std::endl;
        if (*rn_dst != *rn_hodo) {
            std::cerr << "Error: run number mismatch in pair " << k << ": Hodo run" << *rn_hodo
                      << " vs DST run" << *rn_dst << std::endl;
            return kFALSE;
        }
        if (k == 0) run_num = *rn_dst;
    }
    conf.run_num = run_num;
    if (n_runs > 1) {
        std::cout << "#D HTOF_Calib: merged " << n_runs
                   << " runs; output labeled with representative run" << run_num << std::endl;
    }

    const Int_t n_seg = conf.num_of_ch.at("htof");
    const char* side_name[3] = {"U", "D", "S"};

    // Sum an all-event per-seg/side histogram over all hodo files.
    auto merge_hist = [&](const TString& name, const TString& new_name) -> TH1* {
        TH1* h_sum = nullptr;
        for (std::size_t k = 0; k < n_runs; k++) {
            auto* h = (TH1*)f_hodo_list[k]->Get(name);
            if (!h) {
                std::cerr << "Error: Missing " << name << " in " << hodo_paths[k] << std::endl;
                return nullptr;
            }
            if (k == 0) {
                h_sum = (TH1*)h->Clone(new_name);
                h_sum->SetDirectory(nullptr);
            } else {
                h_sum->Add(h);
            }
        }
        return h_sum;
    };

    // +-----+
    // | TDC |
    // +-----+
    std::vector<TH1D*> h_tdc[3];
    for (Int_t s = 0; s < 3; s++) h_tdc[s].resize(n_seg, nullptr);
    for (Int_t s = 0; s < 3; s++) {
        for (Int_t i = 0; i < n_seg; i++) {
            h_tdc[s][i] = (TH1D*)merge_hist(Form("HTOF_TDC_seg%d%s", i, side_name[s]),
                                            Form("HTOF_TDC_seg%d%s_merged", i, side_name[s]));
            if (!h_tdc[s][i]) return kFALSE;
        }
    }
    // search ranges: individual (U+D) and sum (S) channels
    TH1D *h_sum_tdc[2];
    h_sum_tdc[0] = (TH1D*)h_tdc[0][0]->Clone("h_sum_tdc_indiv");
    h_sum_tdc[0]->Reset();
    h_sum_tdc[1] = (TH1D*)h_tdc[2][0]->Clone("h_sum_tdc_sum");
    h_sum_tdc[1]->Reset();
    for (Int_t i = 0; i < n_seg; i++) {
        h_sum_tdc[0]->Add(h_tdc[0][i]);
        h_sum_tdc[0]->Add(h_tdc[1][i]);
        h_sum_tdc[1]->Add(h_tdc[2][i]);
    }

    // +-----+
    // | ADC |
    // +-----+
    // "Raw" = all-event histograms (pedestal); "selected" = TPC dE/dx pion-tagged (MIP)
    std::vector<TH1D*> h_adc_raw[3], h_adc_selected[3];
    for (Int_t s = 0; s < 3; s++) { h_adc_raw[s].resize(n_seg, nullptr); h_adc_selected[s].resize(n_seg, nullptr); }
    for (Int_t s = 0; s < 3; s++) {
        for (Int_t i = 0; i < n_seg; i++) {
            h_adc_raw[s][i] = (TH1D*)merge_hist(Form("HTOF_ADC_seg%d%s", i, side_name[s]),
                                                Form("HTOF_ADC_seg%d%s_merged", i, side_name[s]));
            if (!h_adc_raw[s][i]) return kFALSE;
            h_adc_selected[s][i] = new TH1D(
                Form("HTOF_ADC_seg%d%s_%s", i, side_name[s], particle.Data()),
                Form("HTOF ADC seg%d %s (%s);ADC;Counts", i, side_name[s], particle.Data()),
                conf.adc_bin_num, conf.adc_min, conf.adc_max);
            h_adc_selected[s][i]->SetDirectory(nullptr);
        }
    }

    // per run: event_number -> bitmask of HTOF segments matched to a pion-tagged,
    // non-proton TPC track with a vertex (used to select the PHC hits below)
    std::vector<std::unordered_map<UInt_t, ULong64_t>> sel_seg_mask(n_runs);
    // (ADC U, ADC D) of the ADC-sample tracks on the beam-window short segments 1-4 (weak-side refit)
    std::vector<std::pair<Double_t, Double_t>> win_adc[5];
    // per run: event_number -> (segment -> E), E = t_beam + tof_calc_pi of the first non-beam,
    // non-accidental, pion-tagged (proton bit not set) track with a valid TOF and
    // p_vtx >= tof_offset_min_p matched to that segment: the expected BH2 -> vertex -> HTOF time
    // for a pion (independent of HTOF params).
    // Used for the absolute TOF offset (see the PHC section).
    std::vector<std::unordered_map<UInt_t, std::map<Int_t, Float_t>>> pi_expect(n_runs);
    const Double_t tof_offset_min_p = 0.20; // [GeV/c]
    Long64_t n_entries = 0, n_matched = 0, n_pi_htof = 0, n_pi_htof_novtx = 0, n_pi_htof_novtx_rej = 0, n_pi_tof = 0;
    for (std::size_t k = 0; k < n_runs; k++) {
        TTreeReader reader("htof", f_dst_list[k]);
        TTreeReaderValue<UInt_t> event_number(reader, "event_number");
        TTreeReaderValue<std::vector<Int_t>> pid(reader, "pid");
        TTreeReaderValue<std::vector<Int_t>> seg_match(reader, "seg_match");
        TTreeReaderValue<std::vector<Int_t>> vertex_source(reader, "vertex_source");
        TTreeReaderValue<std::vector<Double_t>> extrap_seg(reader, "extrap_seg");
        TTreeReaderValue<std::vector<Double_t>> match_cl_seg(reader, "match_cl_seg");
        TTreeReaderValue<std::vector<Double_t>> htof_adc_u(reader, "htof_adc_u");
        TTreeReaderValue<std::vector<Double_t>> htof_adc_d(reader, "htof_adc_d");
        TTreeReaderValue<std::vector<Double_t>> htof_adc_s(reader, "htof_adc_s");
        TTreeReaderValue<std::vector<Int_t>> is_beam(reader, "is_beam");
        TTreeReaderValue<std::vector<Int_t>> is_accidental(reader, "is_accidental");
        TTreeReaderValue<std::vector<Double_t>> t_beam(reader, "t_beam");
        TTreeReaderValue<std::vector<Double_t>> tof_calc_pi(reader, "tof_calc_pi");
        TTreeReaderValue<std::vector<Double_t>> dt_pi(reader, "dt_pi");
        TTreeReaderValue<std::vector<Double_t>> p_vtx(reader, "p_vtx");
        if (!all_branches_found(reader, {&event_number, &pid, &seg_match, &vertex_source, &extrap_seg, &match_cl_seg, &htof_adc_u, &htof_adc_d,
                                &htof_adc_s, &is_beam, &is_accidental, &t_beam, &tof_calc_pi, &dt_pi, &p_vtx},
                                dst_paths[k])) return kFALSE;
        while (reader.Next()) {
            n_entries++;
            const auto& v_pid   = *pid;
            const auto& v_seg_match = *seg_match;
            const auto& v_seg   = *extrap_seg;
            const auto& v_adc_u = *htof_adc_u;
            const auto& v_adc_d = *htof_adc_d;
            const auto& v_adc_s = *htof_adc_s;
            for (std::size_t it = 0, n_tr = v_pid.size(); it < n_tr; it++) {
                if (it >= v_seg_match.size() || v_seg_match[it] != 1) continue;       // valid HTOF extrapolation+cluster match
                if (it >= v_seg.size() || TMath::IsNaN(v_seg[it])) continue;
                const Int_t seg = static_cast<Int_t>(std::lround(v_seg[it]));
                if (seg < 0 || seg >= n_seg) continue;
                n_matched++;

                // HypTPCdEdxPID bitmask: bit0=pi, bit1=K, bit2=proton.
                // Require pion bit set and proton bit NOT set; kaon-ambiguous (pid==3) is kept.
                // NOTE: beam/accidental tracks are not excluded here (they dominate the
                // downstream segments).
                if (!(v_pid[it] & 1) || (v_pid[it] & 4)) continue;
                n_pi_htof++;
                // ADC MIP sample: vertex tracks, plus (DstTPCHelixHTOF ExtrapAll=1 only) no-vertex
                // tracks whose extrapolated segment equals the matched cluster segment (dseg = 0).
                // No-vertex |dseg| = 1 matches give scattered ADC peaks, so they
                // are dropped. Timing samples (PHC proton-excluded mask, absolute offset): vertex
                // tracks only. With ExtrapAll=0 every matched track has a vertex, so both samples are the same.
                const Bool_t has_vtx = (it < vertex_source->size() && (*vertex_source)[it] != 0);
                if (!has_vtx) {
                    const Bool_t dseg0 = (it < match_cl_seg->size() && !TMath::IsNaN((*match_cl_seg)[it])
                                          && TMath::Abs(v_seg[it] - (*match_cl_seg)[it]) < 0.5);
                    if (!dseg0) { n_pi_htof_novtx_rej++; continue; }
                    n_pi_htof_novtx++;
                }
                if (has_vtx) sel_seg_mask[k][*event_number] |= (1ULL << seg);
                const Bool_t beamlike = (it < is_beam->size() && (*is_beam)[it]) || (it < is_accidental->size() && (*is_accidental)[it]);
                if (has_vtx && !beamlike && it < dt_pi->size() && !TMath::IsNaN((*dt_pi)[it])
                    && it < p_vtx->size() && (*p_vtx)[it] >= tof_offset_min_p) {
                    auto& m = pi_expect[k][*event_number];
                    if (!m.count(seg)) { m[seg] = static_cast<Float_t>((*t_beam)[it] + (*tof_calc_pi)[it]); n_pi_tof++; }
                }
                if (it < v_adc_u.size() && !TMath::IsNaN(v_adc_u[it])) h_adc_selected[0][seg]->Fill(v_adc_u[it]);
                if (it < v_adc_d.size() && !TMath::IsNaN(v_adc_d[it])) h_adc_selected[1][seg]->Fill(v_adc_d[it]);
                if (it < v_adc_s.size() && !TMath::IsNaN(v_adc_s[it])) h_adc_selected[2][seg]->Fill(v_adc_s[it]);
                if (seg >= 1 && seg <= 4 && it < v_adc_u.size() && it < v_adc_d.size())
                    win_adc[seg].emplace_back(v_adc_u[it], v_adc_d[it]);
            }
        }
    }
    std::cout << "#D HTOF_Calib: " << n_pi_tof << " non-beam pion-tagged tracks with a valid TOF and p_vtx >= " << tof_offset_min_p << " GeV/c (absolute offset sample)" << std::endl;
    std::cout << "#D HTOF_Calib: " << n_entries << " DST events, " << n_matched << " HTOF-matched hits, "
              << n_pi_htof << " pion-tagged (proton excluded) among them; ADC MIP sample "
              << (n_pi_htof - n_pi_htof_novtx_rej) << " = " << (n_pi_htof - n_pi_htof_novtx - n_pi_htof_novtx_rej)
              << " with vertex (also used for the PHC selection) + " << n_pi_htof_novtx
              << " without vertex, dseg = 0 (ADC only); " << n_pi_htof_novtx_rej
              << " without vertex, dseg != 0 (not used)" << std::endl;

    // +---------------------------------+
    // | fit TDC and ADC, one canvas/seg |
    // +---------------------------------+
    // pads (7 columns): row r (U/D/S) -> 7r+1 ADC raw, 7r+2 ADC selected, 7r+3 ADC close-up,
    //       7r+4 TDC, 7r+5 PHC stage 1 (all hits), 7r+6 PHC final (pion-tagged),
    //       7r+7 dt_pi after correction (U/D rows: vs dE; S row: 1D before/after)
    const Int_t rows = 3, cols = 7;
    std::vector<TCanvas*> c_seg(n_seg);
    std::vector<FitResult> tdc_result[3], adc_result[3];
    for (Int_t i = 0; i < n_seg; i++) {
        c_seg[i] = new TCanvas(Form("htof_calib_seg%d", i), "", 4200, 1500);
        c_seg[i]->Divide(cols, rows);
        for (Int_t s = 0; s < 3; s++) {
            const Int_t pad0 = cols*s + 1;
            // rebin is chosen automatically by htof_adc_fit (n_rebin=0)
            adc_result[s].push_back(ana_helper::htof_adc_fit(h_adc_raw[s][i], h_adc_selected[s][i], c_seg[i], pad0, 0,
                                                             Form("htof-%d-%c", i, "uds"[s])));
            ana_helper::set_tdc_search_range(s < 2 ? h_sum_tdc[0] : h_sum_tdc[1]);
            tdc_result[s].push_back(ana_helper::tdc_fit(h_tdc[s][i], c_seg[i], pad0+3));
        }
    }

    // Beam-window short segments (plane 0, seg 1-4): hits concentrate next to the window and the far
    // readout end sees little light, so the weak-side ADC of MIP-like tracks is a spike just above
    // the pedestal plus a MIP hump. Refit the weak side on tracks whose strong-side
    // ADC is within win_cond_k (MIP - ped) of its MIP, spike + hump model (MIP = hump maximum).
    // Replaces the MIP of these channels (used for the PHC dE below and written to the ADC output).
    const std::pair<Int_t, Int_t> win_weak[4] = { {1, 1}, {2, 0}, {3, 1}, {4, 0} }; // (seg, side 0 U / 1 D)
    const Double_t win_cond_k        = 0.3;
    const Double_t win_min_hump_frac = 0.5; // below: flag 64 (hump not separated, MIP uncertain)
    std::vector<TH1D*> h_adc_cond;
    for (const auto& [seg, w] : win_weak) {
        const Int_t st = 1 - w;
        const FitResult& rs = adc_result[st][seg];
        const Double_t ped_s = rs.par[1], mip_s = rs.par[4], u = mip_s - ped_s;
        if (u <= 0.0 || (static_cast<Int_t>(rs.additional[0]) & 8)) {
            std::cerr << "Warning: HTOF_Calib: seg" << seg << side_name[st] << " (strong side) MIP unusable; "
                      << "seg" << seg << side_name[w] << " keeps the standard ADC fit" << std::endl;
            continue;
        }
        TH1D* hc = new TH1D(Form("HTOF_ADC_seg%d%s_%s_cond", seg, side_name[w], particle.Data()),
                            Form("HTOF ADC seg%d %s (%s, %s near MIP);ADC;Counts", seg, side_name[w], particle.Data(), side_name[st]),
                            conf.adc_bin_num, conf.adc_min, conf.adc_max);
        hc->SetDirectory(nullptr);
        for (const auto& [a_u, a_d] : win_adc[seg]) {
            const Double_t a[2] = {a_u, a_d};
            if (TMath::IsNaN(a[w]) || TMath::IsNaN(a[st])) continue;
            if (std::fabs(a[st] - mip_s) < win_cond_k*u) hc->Fill(a[w]);
        }
        const Double_t mip_std = adc_result[w][seg].par[4];
        adc_result[w][seg] = ana_helper::htof_adc_fit_weak(hc, adc_result[w][seg], u, c_seg[seg], cols*w + 2, win_min_hump_frac,
                                                           Form("htof-%d-%c", seg, "uds"[w]));
        h_adc_cond.push_back(hc);
        std::cout << "#D HTOF_Calib: seg" << seg << side_name[w] << " weak-side refit: MIP " << adc_result[w][seg].par[4]
                  << " (standard fit " << mip_std << "), flag " << static_cast<Int_t>(adc_result[w][seg].additional[0])
                  << ", " << hc->GetEntries() << " tracks with " << side_name[st] << " near its MIP" << std::endl;
    }

    // +-----------------------------------------------------+
    // | PHC: mean-time based, with DeltaE from the ADC fit above |
    // +-----------------------------------------------------+
    // A per-channel fit of (time0 - t_U) vs dE_U also absorbs the hit-position correlation
    // along the bar (light propagation time vs attenuation), which does not cancel in the
    // analyzer's mean time and over-corrects it. So U is fitted with
    //   Y_U = (time0 - t_U) + (time0 - ct_D)   ( = 2 x (time0 - mean time) )
    // vs dE_U, where ct_D is D corrected with the current D parameters (propagation cancels),
    // then D likewise, alternating phc_n_alt times. Y_U vs dE_U then follows the U walk
    // (same Type1 function and parameter meaning as a per-channel fit).
    // yu = time0 - t_U, yd = time0 - t_D; E = expected pion BH2->vtx->HTOF time (NaN if none)
    struct PhcPair { Float_t yu, yd, du, dd, E; Bool_t sel; };
    std::vector<std::vector<PhcPair>> phc_pairs(n_seg);
    Long64_t n_hodo_events = 0, n_pair_all = 0, n_pair_sel = 0;
    for (std::size_t k = 0; k < n_runs; k++) {
        TTreeReader reader("hodo", f_hodo_list[k]);
        TTreeReaderValue<UInt_t> event_number(reader, "event_number");
        TTreeReaderValue<Double_t> time0(reader, "time0");
        TTreeReaderValue<std::vector<Double_t>> raw_seg(reader, "htof_raw_seg");
        TTreeReaderValue<std::vector<Double_t>> adc_u(reader, "htof_adc_u");
        TTreeReaderValue<std::vector<Double_t>> adc_d(reader, "htof_adc_d");
        TTreeReaderValue<std::vector<Double_t>> hit_seg(reader, "htof_hit_seg");
        TTreeReaderValue<std::vector<std::vector<Double_t>>> time_u(reader, "htof_time_u");
        TTreeReaderValue<std::vector<std::vector<Double_t>>> time_d(reader, "htof_time_d");
        if (!all_branches_found(reader, {&event_number, &time0, &raw_seg, &adc_u, &adc_d, &hit_seg, &time_u, &time_d},
                                hodo_paths[k])) return kFALSE;
        while (reader.Next()) {
            n_hodo_events++;
            const auto it_mask = sel_seg_mask[k].find(*event_number);
            // 0 if the event has no pion-tagged non-proton track
            const ULong64_t mask = (it_mask == sel_seg_mask[k].end()) ? 0ULL : it_mask->second;
            const auto it_exp = pi_expect[k].find(*event_number);
            const Double_t t0 = *time0;
            if (TMath::IsNaN(t0)) continue;
            const auto& v_raw_seg = *raw_seg;
            const auto& v_hit_seg = *hit_seg;
            for (std::size_t ih = 0; ih < v_hit_seg.size(); ih++) {
                const Int_t seg = static_cast<Int_t>(std::lround(v_hit_seg[ih]));
                if (seg < 0 || seg >= n_seg) continue;
                const Bool_t is_sel = (mask & (1ULL << seg)); // segment matched to such a track
                // raw ADC of the same segment (adc_u/d are indexed like htof_raw_seg)
                Int_t ir = -1;
                for (std::size_t j = 0; j < v_raw_seg.size(); j++) {
                    if (static_cast<Int_t>(std::lround(v_raw_seg[j])) == seg) { ir = static_cast<Int_t>(j); break; }
                }
                if (ir < 0) continue;
                Double_t de[2];
                Bool_t de_ok = kTRUE;
                const Double_t adc[2] = { (*adc_u)[ir], (*adc_d)[ir] };
                for (Int_t s = 0; s < 2; s++) {
                    const Double_t ped = adc_result[s][seg].par[1];
                    const Double_t mip = adc_result[s][seg].par[4];
                    if (TMath::IsNaN(adc[s]) || mip - ped <= 0.0) { de_ok = kFALSE; break; }
                    de[s] = (adc[s] - ped) / (mip - ped);
                }
                if (!de_ok) continue;
                Float_t e_pi = std::numeric_limits<Float_t>::quiet_NaN();
                if (it_exp != pi_expect[k].end()) {
                    const auto je = it_exp->second.find(seg);
                    if (je != it_exp->second.end()) e_pi = je->second;
                }
                // U/D leading times of a hit are already paired (same index) by the analyzer
                const auto& vu = (*time_u)[ih];
                const auto& vd = (*time_d)[ih];
                for (std::size_t j = 0, nj = std::min(vu.size(), vd.size()); j < nj; j++) {
                    if (TMath::IsNaN(vu[j]) || TMath::IsNaN(vd[j])) continue;
                    phc_pairs[seg].push_back({ static_cast<Float_t>(t0 - vu[j]), static_cast<Float_t>(t0 - vd[j]),
                                               static_cast<Float_t>(de[0]), static_cast<Float_t>(de[1]), e_pi, is_sel });
                    n_pair_all++;
                    if (is_sel) n_pair_sel++;
                }
            }
        }
    }
    std::cout << "#D HTOF_Calib: " << n_hodo_events << " hodo events, "
              << n_pair_all << " U/D-paired (hit x multi-hit) PHC entries for all hits, "
              << n_pair_sel << " for the proton-excluded selection" << std::endl;

    // Binning taken from the analyzer's all-event HTOF_seg0U_TOF_vs_DeltaE histogram.
    auto* h_tmpl = (TH2D*)f_hodo_list[0]->Get("HTOF_seg0U_TOF_vs_DeltaE");
    if (!h_tmpl) {
        std::cerr << "Error: Missing HTOF_seg0U_TOF_vs_DeltaE (binning template) in " << hodo_paths[0] << std::endl;
        return kFALSE;
    }
    // Type1 walk term as fitted: (time0 - time) = -p0/sqrt|dE-p1| + p2
    auto walk = [](const std::vector<Double_t>& p, Double_t de) {
        return -p[0]/std::sqrt(std::fabs(de - p[1])) + p[2];
    };
    // Y histogram for side s (0: U, 1: D) of segment i, the other side corrected with p_other.
    // x binning from the analyzer template; y in 0.1 ns bins (the band core is ~0.4 ns wide).
    auto fill_y = [&](Int_t s, Int_t i, Bool_t sel_only, const std::vector<Double_t>& p_other, const TString& name) {
        TH2D* h = new TH2D(name, "", h_tmpl->GetNbinsX(), h_tmpl->GetXaxis()->GetXmin(), h_tmpl->GetXaxis()->GetXmax(),
                           200, h_tmpl->GetYaxis()->GetXmin(), h_tmpl->GetYaxis()->GetXmax());
        h->SetDirectory(nullptr);
        h->SetTitle(Form("HTOF seg%d %s (time0-t_{%s})+(time0-ct_{%s}) %s;#DeltaE_{%s} [MIP];[ns]", i, side_name[s],
                         side_name[s], side_name[1-s], sel_only ? "(proton excluded)" : "(all hits)", side_name[s]));
        for (const auto& p : phc_pairs[i]) {
            if (sel_only && !p.sel) continue;
            if (s == 0) h->Fill(p.du, p.yu + (p.yd - walk(p_other, p.dd)));
            else        h->Fill(p.dd, p.yd + (p.yu - walk(p_other, p.du)));
        }
        return h;
    };

    // dE window [MIP] for the per-channel p2 re-centring (median residual = 0)
    const Double_t phc_recenter_de_lo = 0.5, phc_recenter_de_hi = 2.5;
    // The mean time only constrains the U+D sum of the offsets p2 (p2_U + c, p2_D - c leaves
    // it unchanged, but shifts U-D time differences), so after every alternation each
    // side's p2 is shifted to make its median per-channel residual (time0 - t) - walk(dE)
    // zero (hits with 0.5 < dE < 2.5 MIP), i.e. each channel aligned to time0 as in a
    // per-channel fit, without changing the mean-time shape.
    auto recenter = [&](Int_t i, Bool_t sel_only, std::vector<Double_t> p_cur[2]) {
        std::vector<Double_t> res[2];
        for (const auto& p : phc_pairs[i]) {
            if (sel_only && !p.sel) continue;
            if (p.du > phc_recenter_de_lo && p.du < phc_recenter_de_hi) res[0].push_back(p.yu - walk(p_cur[0], p.du));
            if (p.dd > phc_recenter_de_lo && p.dd < phc_recenter_de_hi) res[1].push_back(p.yd - walk(p_cur[1], p.dd));
        }
        for (Int_t s = 0; s < 2; s++) {
            if (res[s].empty()) continue;
            std::nth_element(res[s].begin(), res[s].begin() + res[s].size()/2, res[s].end());
            p_cur[s][2] += res[s][res[s].size()/2];
        }
    };

    // Two stages: (1) all hits, (2) proton-excluded hits with the parameters in
    // phc_fix_mask fixed to the stage-1 values (the proton-excluded sample alone is
    // concentrated near 1 MIP, which leaves the curvature (p0 vs p1) unconstrained).
    // bit i -> p_i; TENTATIVE choice: fix p1 (pole position) only.
    const Int_t phc_fix_mask = (1 << 1);
    // fit-range threshold (TENTATIVE 5% / 4%): low enough that the fit reaches into the
    // sparse low-dE region.
    const Double_t phc_dense_ratio_all = 0.05, phc_dense_ratio_sel = 0.04;
    // U/D alternations per stage (the mean-time fit of one side depends on the other side)
    const Int_t phc_n_alt = 4;
    // last-iteration Y histograms (drawn/saved): stage 1 (all hits) and final
    std::vector<TH2D*> h_phc_all[2], h_phc[2];
    for (Int_t s = 0; s < 2; s++) { h_phc_all[s].resize(n_seg, nullptr); h_phc[s].resize(n_seg, nullptr); }
    // absolute TOF offset per segment (shift added to p2 of U and D) and its sample size
    std::vector<Double_t> phc_tof_offset(n_seg, 0.0);
    std::vector<Int_t> phc_tof_offset_n(n_seg, 0);
    const std::size_t phc_tof_offset_min_n = 100;
    std::vector<FitResult> phc_result_all[2], phc_result[2];
    TCanvas* c_phc_tmp = new TCanvas("c_phc_tmp", "", 400, 300);
    for (Int_t i = 0; i < n_seg; i++) {
        std::vector<Double_t> p_cur[2] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };  // start: no correction
        FitResult r_stage[2][2];  // [stage][side]
        TH2D* r_hist[2] = { nullptr, nullptr };
        for (Int_t stage = 0; stage < 2; stage++) {
            const Bool_t sel_only = (stage == 1);
            for (Int_t it = 0; it < phc_n_alt; it++) {
                const Bool_t last = (it == phc_n_alt - 1);
                for (Int_t s = 0; s < 2; s++) {
                    TH2D* h = fill_y(s, i, sel_only, p_cur[1-s],
                                     Form("HTOF_seg%d%s_PHC_Y_%s%s", i, side_name[s], sel_only ? "sel" : "all", last ? "" : "_tmp"));
                    TCanvas* c = last ? c_seg[i] : c_phc_tmp;
                    // pad 0 = the (undivided) scratch canvas itself; cd(1) on an undivided
                    // canvas does not change gPad and would draw into the page's last pad
                    const Int_t pad = last ? cols*s + 5 + stage : 0;
                    FitResult r = (stage == 0)
                        ? ana_helper::htof_phc_fit(h, c, pad, nullptr, 0, phc_dense_ratio_all)
                        : ana_helper::htof_phc_fit(h, c, pad, &r_stage[0][s], phc_fix_mask, phc_dense_ratio_sel);
                    p_cur[s] = { r.par[0], r.par[1], r.par[2] };
                    r_stage[stage][s] = r;
                    r_hist[s] = h;
                    if (last) {
                        (stage == 0 ? h_phc_all : h_phc)[s][i] = h;
                    }
                }
                recenter(i, sel_only, p_cur);
                for (Int_t s = 0; s < 2; s++) r_stage[stage][s].par[2] = p_cur[s][2];
                if (!last) {
                    c_phc_tmp->Clear();
                    for (Int_t s = 0; s < 2; s++) delete r_hist[s];
                }
            }
        }
        // absolute TOF offset: common shift of p2 (U and D, so the U-D difference is kept)
        // making the median of (ct_mean - time0) - E zero for the pion sample, i.e. dt_pi = 0.
        // The walk shape (p0, p1) and U/D balance are those of the stages above.
        {
            std::vector<Double_t> r;
            for (const auto& p : phc_pairs[i]) {
                if (std::isnan(p.E)) continue;
                r.push_back(0.5*((p.yu + p.E - walk(p_cur[0], p.du)) + (p.yd + p.E - walk(p_cur[1], p.dd))));
            }
            phc_tof_offset_n[i] = static_cast<Int_t>(r.size());
            if (r.size() >= phc_tof_offset_min_n) {
                std::nth_element(r.begin(), r.begin() + r.size()/2, r.end());
                phc_tof_offset[i] = r[r.size()/2];
                for (Int_t s = 0; s < 2; s++) { p_cur[s][2] += phc_tof_offset[i]; r_stage[1][s].par[2] = p_cur[s][2]; }
            } else {
                std::cerr << "Warning: HTOF seg" << i << ": only " << r.size()
                          << " pion TOF entries, absolute TOF offset not applied" << std::endl;
            }
        }
        for (Int_t s = 0; s < 2; s++) {
            phc_result_all[s].push_back(r_stage[0][s]);
            phc_result[s].push_back(r_stage[1][s]);
        }
        // S row, PHC columns: value summary
        c_seg[i]->cd(cols*2 + 5);
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.06);
        tx.DrawLatex(0.08, 0.90, Form("HTOF seg%d (run%05d%s)", i, run_num, n_runs > 1 ? Form(", %zu runs", n_runs) : ""));
        for (Int_t s = 0; s < 3; s++) {
            const FitResult& a = adc_result[s][i];
            tx.DrawLatex(0.08, 0.78 - 0.11*s, Form("%s: ped %.1f  MIP %.1f  flag %d", side_name[s],
                                                  a.par[1], a.par[4], static_cast<Int_t>(a.additional[0])));
        }
        for (Int_t s = 0; s < 3; s++) {
            tx.DrawLatex(0.08, 0.42 - 0.11*s, Form("%s: TDC %.1f", side_name[s], tdc_result[s][i].par[1]));
        }
        c_seg[i]->cd(cols*2 + 6);
        TLatex tp;
        tp.SetNDC();
        tp.SetTextSize(0.055);
        tp.DrawLatex(0.05, 0.90, "PHC  -p0/#sqrt{|dE-p1|} + p2, fitted on Y = (t0-t)+(t0-ct_{other})");
        for (Int_t s = 0; s < 2; s++) {
            const FitResult& ra = phc_result_all[s][i];
            const FitResult& r  = phc_result[s][i];
            tp.DrawLatex(0.05, 0.76 - 0.30*s, Form("%s all : %.3f %.3f %.3f", side_name[s], ra.par[0], ra.par[1], ra.par[2]));
            tp.DrawLatex(0.05, 0.66 - 0.30*s, Form("%s final: %.3f %.3f %.3f", side_name[s], r.par[0], r.par[1], r.par[2]));
        }
        tp.DrawLatex(0.05, 0.22, Form("TOF offset (pion dt=0): %+.3f ns from %d entries", phc_tof_offset[i], phc_tof_offset_n[i]));
        tp.DrawLatex(0.05, 0.12, Form("final: fix mask %d; range thr. %.1f%% / %.1f%%; %d U/D alternations",
                                      phc_fix_mask, 100*phc_dense_ratio_all, 100*phc_dense_ratio_sel, phc_n_alt));
    }
    delete c_phc_tmp;

    // -- after correction: pion TOF sample (non-beam, pion-tagged, valid TOF), final parameters
    //    applied as the analyzer does (Type1 per channel, then mean time):
    //    dt_pi = (ct_mean - time0) - E, should be 0 and flat in dE -----
    std::vector<TH2D*> h_dtpi2d[2];
    for (Int_t s = 0; s < 2; s++) h_dtpi2d[s].resize(n_seg, nullptr);
    std::vector<TH1D*> h_dtpi_before(n_seg, nullptr), h_dtpi_after(n_seg, nullptr);
    for (Int_t i = 0; i < n_seg; i++) {
        const std::vector<Double_t> pu = { phc_result[0][i].par[0], phc_result[0][i].par[1], phc_result[0][i].par[2] };
        const std::vector<Double_t> pd = { phc_result[1][i].par[0], phc_result[1][i].par[1], phc_result[1][i].par[2] };
        for (Int_t s = 0; s < 2; s++) {
            h_dtpi2d[s][i] = new TH2D(Form("HTOF_seg%d_dTPi_vs_DeltaE%s_calib", i, side_name[s]),
                                     Form("HTOF seg%d dt_{#pi} after correction vs #DeltaE_{%s};#DeltaE_{%s} [MIP];(ct_{U}+ct_{D})/2 - time0 - t_{exp,#pi} [ns]",
                                          i, side_name[s], side_name[s]),
                                     h_tmpl->GetNbinsX(), h_tmpl->GetXaxis()->GetXmin(), h_tmpl->GetXaxis()->GetXmax(), 200, -10.0, 10.0);
            h_dtpi2d[s][i]->SetDirectory(nullptr);
        }
        h_dtpi_before[i] = new TH1D(Form("HTOF_seg%d_dTPi_before_calib", i),
                                  Form("HTOF seg%d dt_{#pi} (gray: no PHC/offset, blue: final);mean time - time0 - t_{exp,#pi} [ns];Counts", i), 200, -10.0, 10.0);
        h_dtpi_after[i]  = new TH1D(Form("HTOF_seg%d_dTPi_after_calib", i), "", 200, -10.0, 10.0);
        h_dtpi_before[i]->SetDirectory(nullptr);
        h_dtpi_after[i]->SetDirectory(nullptr);
        for (const auto& p : phc_pairs[i]) {
            if (std::isnan(p.E)) continue;
            const Double_t y_before = -0.5*(p.yu + p.yd) - p.E;
            const Double_t y_after  = -0.5*((p.yu - walk(pu, p.du)) + (p.yd - walk(pd, p.dd))) - p.E;
            h_dtpi2d[0][i]->Fill(p.du, y_after);
            h_dtpi2d[1][i]->Fill(p.dd, y_after);
            h_dtpi_before[i]->Fill(y_before);
            h_dtpi_after[i]->Fill(y_after);
        }
        // residual dE dependence: TOF peak (mean within +-0.8 ns of the local maximum) in
        // dE slices of 0.4-2.1 MIP; "spread" = max - min of the peaks (0 if fully corrected)
        const Double_t edges[] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.1};
        for (Int_t s = 0; s < 2; s++) {
            TGraph* g = new TGraph();
            Double_t lo = 1e9, hi = -1e9;
            for (Int_t k = 0; k < 7; k++) {
                TH1D* pj = h_dtpi2d[s][i]->ProjectionY(Form("%s_pj%d", h_dtpi2d[s][i]->GetName(), k),
                                                      h_dtpi2d[s][i]->GetXaxis()->FindBin(edges[k]),
                                                      h_dtpi2d[s][i]->GetXaxis()->FindBin(edges[k+1]) - 1);
                if (pj->Integral() >= 30) {
                    const Double_t c = pj->GetBinCenter(pj->GetMaximumBin());
                    pj->GetXaxis()->SetRangeUser(c - 0.8, c + 0.8);
                    const Double_t y = pj->GetMean();
                    g->SetPoint(g->GetN(), 0.5*(edges[k] + edges[k+1]), y);
                    lo = std::min(lo, y); hi = std::max(hi, y);
                }
                delete pj;
            }
            c_seg[i]->cd(cols*s + 7);
            gPad->SetLogz(0);
            h_dtpi2d[s][i]->GetYaxis()->SetRangeUser(-3.0, 3.0);
            h_dtpi2d[s][i]->Draw("colz");
            TLine* l0 = new TLine(h_dtpi2d[s][i]->GetXaxis()->GetXmin(), 0.0, h_dtpi2d[s][i]->GetXaxis()->GetXmax(), 0.0);
            l0->SetLineStyle(2);
            l0->SetLineColor(kGray+2);
            l0->Draw("same");
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.8);
            g->SetMarkerColor(kRed);
            if (g->GetN() > 0) g->Draw("P same");
            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.05);
            if (g->GetN() >= 3) tl.DrawLatex(0.20, 0.85, Form("peak spread (0.4-2.1 MIP) %.2f ns", hi - lo));
        }
        // S row: 1D mean time before / after
        c_seg[i]->cd(cols*2 + 7);
        h_dtpi_before[i]->SetLineColor(kGray+2);
        h_dtpi_after[i]->SetLineColor(kBlue);
        h_dtpi_before[i]->SetMaximum(1.15*std::max(h_dtpi_before[i]->GetMaximum(), h_dtpi_after[i]->GetMaximum()));
        h_dtpi_before[i]->Draw("HIST");
        h_dtpi_after[i]->Draw("HIST same");
        TLatex tm;
        tm.SetNDC();
        tm.SetTextSize(0.05);
        tm.SetTextColor(kGray+2);
        tm.DrawLatex(0.18, 0.85, Form("before: RMS %.2f ns", h_dtpi_before[i]->GetRMS()));
        tm.SetTextColor(kBlue);
        tm.DrawLatex(0.18, 0.78, Form("after : RMS %.2f ns", h_dtpi_after[i]->GetRMS()));
    }

    // +-----+
    // | PDF |
    // +-----+
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_HTOF_Calib_%s.pdf", img_base_dir.Data(), run_num, particle.Data());
    c_seg[0]->Print(pdf_path + "[");
    for (Int_t i = 0; i < n_seg; i++) c_seg[i]->Print(pdf_path);
    c_seg[0]->Print(pdf_path + "]");

    // +-------+
    // | Write |
    // +-------+
    // -- TDC -----
    {
        TFile* fout = open_output(run_num, "HDPRM_All");
        TTree* tree = new TTree("tree", "");
        Int_t ch;
        std::vector<Double_t> tdc_p0_val, tdc_p0_err;
        tree->Branch("ch", &ch, "ch/I");
        tree->Branch("tdc_p0_val", &tdc_p0_val);
        tree->Branch("tdc_p0_err", &tdc_p0_err);
        for (Int_t i = 0; i < n_seg; i++) {
            ch = i;
            tdc_p0_val.clear(); tdc_p0_err.clear();
            for (Int_t s = 0; s < 3; s++) {
                tdc_p0_val.push_back(tdc_result[s][i].par[1]);
                tdc_p0_err.push_back(tdc_result[s][i].err[1]);
            }
            tree->Fill();
        }
        tree->Write();
        fout->Close();
    }
    // -- ADC -----
    {
        TFile* fout = open_output(run_num, Form("ADC_%s", particle.Data()));
        TTree* tree = new TTree("tree", "");
        Int_t ch;
        std::vector<Double_t> adc_p0_val, adc_p1_val, adc_p0_err, adc_p1_err, adc_chi2;
        // fit quality (U, D, S): flag bits documented at ana_helper::htof_adc_fit
        std::vector<Int_t> adc_flag, adc_ndf;
        tree->Branch("ch", &ch, "ch/I");
        tree->Branch("adc_p0_val", &adc_p0_val);
        tree->Branch("adc_p1_val", &adc_p1_val);
        tree->Branch("adc_p0_err", &adc_p0_err);
        tree->Branch("adc_p1_err", &adc_p1_err);
        tree->Branch("adc_flag", &adc_flag);
        tree->Branch("adc_chi2", &adc_chi2);
        tree->Branch("adc_ndf", &adc_ndf);
        for (Int_t i = 0; i < n_seg; i++) {
            ch = i;
            adc_p0_val.clear(); adc_p1_val.clear(); adc_p0_err.clear(); adc_p1_err.clear();
            adc_flag.clear(); adc_chi2.clear(); adc_ndf.clear();
            for (Int_t s = 0; s < 3; s++) {
                const FitResult& r = adc_result[s][i];
                adc_p0_val.push_back(r.par[1]); adc_p0_err.push_back(r.err[1]);   // pedestal
                adc_p1_val.push_back(r.par[4]); adc_p1_err.push_back(r.err[4]);   // mip
                adc_flag.push_back(static_cast<Int_t>(r.additional[0]));
                adc_chi2.push_back(r.chi_square);
                adc_ndf.push_back(r.ndf);
            }
            tree->Fill();
        }
        // Keep the underlying histograms too, for offline inspection without rerunning.
        for (Int_t s = 0; s < 3; s++) {
            for (Int_t i = 0; i < n_seg; i++) {
                h_adc_raw[s][i]->Write();
                h_adc_selected[s][i]->Write();
            }
        }
        for (TH1D* hc : h_adc_cond) hc->Write();
        tree->Write();
        fout->Close();
    }
    // -- PHC -----
    {
        TFile* fout = open_output(run_num, Form("PHC_%s", particle.Data()));
        TTree* tree = new TTree("tree", "");
        Int_t ch;
        std::vector<Double_t> p0_val, p1_val, p2_val, p0_err, p1_err, p2_err;
        // stage-1 (all hits) values, for reference only (update_param.py reads p*_val)
        std::vector<Double_t> p0_all_val, p1_all_val, p2_all_val;
        // final-fit dE range [MIP]
        std::vector<Double_t> fit_lo, fit_hi;
        // absolute TOF offset included in p2 (same for U and D) and its pion sample size
        Double_t tof_offset; Int_t tof_offset_n;
        tree->Branch("ch", &ch, "ch/I");
        tree->Branch("p0_val", &p0_val);
        tree->Branch("p1_val", &p1_val);
        tree->Branch("p2_val", &p2_val);
        tree->Branch("p0_err", &p0_err);
        tree->Branch("p1_err", &p1_err);
        tree->Branch("p2_err", &p2_err);
        tree->Branch("tof_offset", &tof_offset, "tof_offset/D");
        tree->Branch("tof_offset_n", &tof_offset_n, "tof_offset_n/I");
        tree->Branch("fit_lo", &fit_lo);
        tree->Branch("fit_hi", &fit_hi);
        tree->Branch("p0_all_val", &p0_all_val);
        tree->Branch("p1_all_val", &p1_all_val);
        tree->Branch("p2_all_val", &p2_all_val);
        for (Int_t i = 0; i < n_seg; i++) {
            ch = i;
            p0_val.clear(); p1_val.clear(); p2_val.clear();
            p0_err.clear(); p1_err.clear(); p2_err.clear();
            p0_all_val.clear(); p1_all_val.clear(); p2_all_val.clear();
            tof_offset = phc_tof_offset[i]; tof_offset_n = phc_tof_offset_n[i];
            fit_lo.clear(); fit_hi.clear();
            for (Int_t s = 0; s < 2; s++) {
                const FitResult& ra = phc_result_all[s][i];
                p0_all_val.push_back(ra.par[0]); p1_all_val.push_back(ra.par[1]); p2_all_val.push_back(ra.par[2]);
                const FitResult& r = phc_result[s][i];
                fit_lo.push_back(r.additional[0]); fit_hi.push_back(r.additional[1]);
                p0_val.push_back(r.par[0]); p0_err.push_back(r.err[0]);
                p1_val.push_back(r.par[1]); p1_err.push_back(r.err[1]);
                p2_val.push_back(r.par[2]); p2_err.push_back(r.err[2]);
            }
            tree->Fill();
        }
        for (Int_t s = 0; s < 2; s++) {
            for (Int_t i = 0; i < n_seg; i++) { h_phc_all[s][i]->Write(); h_phc[s][i]->Write(); h_dtpi2d[s][i]->Write(); }
        }
        for (Int_t i = 0; i < n_seg; i++) { h_dtpi_before[i]->Write(); h_dtpi_after[i]->Write(); }
        tree->Write();
        fout->Close();
    }
    return kTRUE;
}

Int_t main(int argc, char** argv) {

    // -- check arguments -----
    // <Hodo1.root> <DST1.root> [<Hodo2.root> <DST2.root> ...] <particle>
    // (a single pair is the normal single-run case; more pairs merge runs)
    if (argc < 4 || (argc - 2) % 2 != 0) {
        std::cerr << "Usage: " << argv[0]
                   << " <Hodo1.root> <DST1.root> [<Hodo2.root> <DST2.root> ...] <particle>" << std::endl;
        std::cerr << "  Hodo.root: all-event HTOF TDC/ADC histograms and the hodo tree (PHC)" << std::endl;
        std::cerr << "  DST.root : DstTPCHelixHTOF output (TPC dE/dx pion tag for the ADC MIP)" << std::endl;
        std::cerr << "  particle : Pi (only supported value; used for output naming)" << std::endl;
        return 1;
    }
    TString particle = argv[argc-1];
    if (particle != "Pi") {
        std::cerr << "Error: Unsupported particle '" << particle
                  << "'. Only 'Pi' is implemented (dE/dx pion-bit selection, proton excluded)." << std::endl;
        return 1;
    }
    std::vector<TString> hodo_paths, dst_paths;
    for (Int_t k = 0; k < (argc - 2) / 2; k++) {
        hodo_paths.push_back(argv[1 + 2*k]);
        dst_paths.push_back(argv[2 + 2*k]);
    }

    conf.detector = "htof";
    conf.phc_de_range_min = 0.2; // lower dE bound [MIP] of the PHC fit range (used by htof_phc_fit)
    const Bool_t ok = analyze(hodo_paths, dst_paths, particle);
    if (!ok) std::cerr << "Error: HTOF_Calib failed; outputs are not valid" << std::endl;
    // ROOT 6.40.04's own static destructor (RConcurrentHashColl) crashes on
    // normal exit after unloading libRIO.so; bypass it. Not this file's bug.
    // Same workaround as the analyzer's example/DstTPCHelixHTOF.cc. Output is already written above.
    std::cout.flush();
    std::cerr.flush();
    std::_Exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
