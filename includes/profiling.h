#ifndef NEW_DATA_TYPE_ANA_PROFILING_H
#define NEW_DATA_TYPE_ANA_PROFILING_H

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
#include <TStyle.h>

#include "../includes/histos.h"
class profiling {
 public:
    Double_t X_Axis_Point;
    Double_t Y_Axis_Point;
    Double_t X_Axis_Point_Err;
    Double_t Y_Axis_Point_Err;

    Double_t second;
    Double_t first;

    TH1D *bin_proj;
    Double_t fit_range_L;
    Double_t fit_range_R;

    profiling();

    virtual ~ profiling();
    virtual void init();
    virtual void bins_filling(TH2F *histo, Int_t i);
    virtual void finding_fit_range();
    virtual void fitting_bin_hist(histos *hist, Int_t i);
    virtual void main_algorithm(Int_t start, Int_t end, TH2F *histo, histos *hist);
};

#endif //NEW_DATA_TYPE_ANA_PROFILING_H
