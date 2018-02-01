#ifndef ANALYSIS_MAMBA_PRESETTING_H
#define ANALYSIS_MAMBA_PRESETTING_H

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

class MAMBA_presetting {
 private:
    static const Int_t length_ntrk = 10000;
    static const Int_t length_inter = 70000;
    static const Int_t length_sp = 140;
 public:
    Int_t nevent;
    ULong64_t timestamp;
    Int_t nDet;
///--------
    Int_t ntrk;
//    Double_t chi2[length_ntrk];
//    Double_t ndof[length_ntrk];
//    Int_t fitStatus[length_ntrk];
//    Double_t ItpX[length_ntrk];
//    Double_t ItpY[length_ntrk];
//    Double_t SlpX[length_ntrk];
//    Double_t SlpY[length_ntrk];
//    Double_t XCov00[length_ntrk];
//    Double_t XCov01[length_ntrk];
//    Double_t XCov10[length_ntrk];
//    Double_t XCov11[length_ntrk];
//    Double_t YCov00[length_ntrk];
//    Double_t YCov01[length_ntrk];
//    Double_t YCov10[length_ntrk];
//    Double_t YCov11[length_ntrk];
//    Double_t XrecoT0[length_ntrk];
//    Double_t YrecoT0[length_ntrk];
//    Double_t XrecoT1[length_ntrk];
//    Double_t YrecoT1[length_ntrk];
//    Double_t XrecoT2[length_ntrk];
//    Double_t YrecoT2[length_ntrk];
//    Double_t XrecoT3[length_ntrk];
//    Double_t YrecoT3[length_ntrk];
//    Double_t XfitT0[length_ntrk];
//    Double_t YfitT0[length_ntrk];
//    Double_t XfitT1[length_ntrk];
//    Double_t YfitT1[length_ntrk];
//    Double_t XfitT2[length_ntrk];
//    Double_t YfitT2[length_ntrk];
//    Double_t XfitT3[length_ntrk];
//    Double_t YfitT3[length_ntrk];
/////--------
    Int_t nIntersect;
//    Int_t intTrkID[length_inter];
//    Int_t intLayer[length_inter];
//    Int_t intInside[length_inter];
//    Double_t intXPos[length_inter];
//    Double_t intYPos[length_inter];
//    Double_t intZPos[length_inter];
//    Double_t intUPos[length_inter];
//    Double_t intVPos[length_inter];
///--------
    Int_t nSP;
    Int_t spLayer[length_sp];
    Double_t spXPos[length_sp];
    Double_t spYPos[length_sp];
    Double_t spZPos[length_sp];
    Double_t spXErr[length_sp];
    Double_t spYErr[length_sp];
    Double_t spZErr[length_sp];
    Double_t spUPos[length_sp];
    Double_t spVPos[length_sp];
    Bool_t spIsOnTrk[length_sp];

    MAMBA_presetting(TTree *MAMBA_EVENT_tree);
    virtual ~MAMBA_presetting();
    virtual void Init_MAMBA_presetting();
    virtual void SetUpMAMBATree(TTree *MAMBA_EVENT_tree);
};

#endif //ANALYSIS_MAMBA_PRESETTING_H
