#ifndef ANALYSIS_HISTOS_H
#define ANALYSIS_HISTOS_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include <algorithm>
#include <functional>
#include "TGraphErrors.h"
#include "TProfile.h"

using namespace std;

class histos {
 public:

// V-shapes (the r(t) relation)
// short straw = S
    TH2F *vshapeUS;
    TH2F *vshapeVS;
// long straw = L
    TH2F *vshapeUL;
    TH2F *vshapeVL;

// V-shapes (the r(t) relation) CLEAR
// short straw = S
    TH2F *vshapeUS_clear;
    TH2F *vshapeUS_cleaned;
    TH2F *vshapeVS_clear;
// long straw = L
    TH2F *vshapeUL_clear;
    TH2F *vshapeUL_cleaned;
    TH2F *vshapeVL_clear;

    TProfile *straw_L;
    TProfile *straw_S;

    TH1F *Resolution_L;
    TH1F *Resolution_S;

    TH1F *n_hits_long_straw;
    TH1F *n_hits_short_straw;

    TH1D *vshapeUS_proj;
    TH1D *vshapeUL_proj;

    TH1D *check_prof;

    TH1D *vshapeUS_proj_Y;
    TH1D *vshapeUL_proj_Y;

    histos();
    virtual ~histos();
    virtual void Init_histos();
    virtual void Drawing_histos();
};

#endif //ANALYSIS_HISTOS_H
