#ifndef ANALYSIS_LINKING_DATA_H
#define ANALYSIS_LINKING_DATA_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include <algorithm>
#include <functional>

#include "../includes/MAMBA_presetting.h"
#include "../includes/STRAW_presetting.h"
#include "../includes/histos.h"
#include "../includes/profiling.h"

class linking_data {
 public:
    // TDC conversion for hits: 1024 is 25ns
    Double_t conversion;

    Int_t mam_nEntries;
    Int_t N;

    Int_t count_processed;
    Int_t count_good;
    Int_t count_very_good;  // have a good track

    Long_t run;

    Double_t Max_Short, Max_Long;
    Double_t a_L, b_L, c_L;
    Double_t a_S, b_S, c_S;
    Double_t a_L_err, b_L_err, c_L_err;
    Double_t a_S_err, b_S_err, c_S_err;

    TF1 *fit_prof_L;
    TF1 *fit_prof_S;

    Double_t vertex_err_L, vertex_err_S;
    Double_t Discrim2_L, Discrim2_S;
    vector<Double_t> timing_L, timing_S, timing_L_after_cut, timing_S_after_cut;
    vector<Double_t> StrawPoint_L, StrawPoint_S, StrawPoint_L_after_cut, StrawPoint_S_after_cut;

    Double_t range_S[3] = {};
    Double_t range_L[3] = {};
    Double_t left_limit_S, left_limit_L, right_limit_S, right_limit_L;

    linking_data(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree, Long_t run_num);

    virtual ~linking_data();
    virtual void init();
    virtual void merging(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree);
    virtual void track_ana(MAMBA_presetting *mamba, STRAW_presetting *straw);
    virtual void makingProfile(histos *hist);
    virtual void StrawResolution_L(histos *hist);
    virtual void StrawResolution_S(Double_t a, Double_t b, Double_t c, TH1F *histo);
    virtual void filling_hists(histos *hist, STRAW_presetting *straw, MAMBA_presetting *mamba);
    virtual void limits_S(histos *hist, Double_t &limits);
    virtual void limits_L(histos *hist, Double_t &limits);
    virtual void draft(histos *hist, profiling *profile);
};

#endif //ANALYSIS_LINKING_DATA_H
