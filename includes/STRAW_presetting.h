#ifndef ANALYSIS_STRAW_PRESETTING_H
#define ANALYSIS_STRAW_PRESETTING_H

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

using namespace std;

class STRAW_presetting {
 public:
    Int_t PCclock;
    Int_t pattern;
    Int_t TDCtrig;
    Double_t t;
    Int_t nt00;
    Int_t tt00;
    Int_t nt01;
    Int_t tt01;
    Int_t nt02;
    Int_t tt02;
    Int_t nt03;
    Int_t tt03;
    Int_t nt04;
    Int_t tt04;
    Int_t nt05;
    Int_t tt05;
    Int_t nt06;
    Int_t tt06;
    Int_t nt08;
    Int_t tt08;

    STRAW_presetting(TTree *STRAW_EVENT_tree);
    virtual ~STRAW_presetting();
    virtual void Init_STRAW_presetting();
    virtual void SetUpSTRAWTree(TTree *STRAW_EVENT_tree);
};

#endif //ANALYSIS_STRAW_PRESETTING_H
