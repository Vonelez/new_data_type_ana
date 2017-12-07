#include "../includes/MAMBA_presetting.h"

MAMBA_presetting::MAMBA_presetting(TTree *MAMBA_EVENT_tree) {
    Init_MAMBA_presetting();
    SetUpMAMBATree(MAMBA_EVENT_tree);
}

MAMBA_presetting::~MAMBA_presetting() = default;

void MAMBA_presetting::Init_MAMBA_presetting() {
    nevent = 0;
    timestamp = 0;
    nDet = 0;
    ntrk = 0;
    nIntersect = 0;
    nSP = 0;
}

void MAMBA_presetting::SetUpMAMBATree(TTree *MAMBA_EVENT_tree) {
    MAMBA_EVENT_tree->SetBranchAddress("nevent", &nevent);
    MAMBA_EVENT_tree->SetBranchAddress("timestamp", &timestamp);
    MAMBA_EVENT_tree->SetBranchAddress("nDet", &nDet);
///--------
    MAMBA_EVENT_tree->SetBranchAddress("ntrk", &ntrk);
//    MAMBA_EVENT_tree->SetBranchAddress("chi2", &chi2);
//    MAMBA_EVENT_tree->SetBranchAddress("ndof", &ndof);
//    MAMBA_EVENT_tree->SetBranchAddress("fitStatus", &fitStatus);
//    MAMBA_EVENT_tree->SetBranchAddress("ItpX", &ItpX);
//    MAMBA_EVENT_tree->SetBranchAddress("ItpY", &ItpY);
//    MAMBA_EVENT_tree->SetBranchAddress("SlpX", &SlpX);
//    MAMBA_EVENT_tree->SetBranchAddress("SlpY", &SlpY);
//    MAMBA_EVENT_tree->SetBranchAddress("XCov00", &XCov00);
//    MAMBA_EVENT_tree->SetBranchAddress("XCov01", &XCov01);
//    MAMBA_EVENT_tree->SetBranchAddress("XCov10", &XCov10);
//    MAMBA_EVENT_tree->SetBranchAddress("XCov11", &XCov11);
//    MAMBA_EVENT_tree->SetBranchAddress("YCov00", &YCov00);
//    MAMBA_EVENT_tree->SetBranchAddress("YCov01", &YCov01);
//    MAMBA_EVENT_tree->SetBranchAddress("YCov10", &YCov10);
//    MAMBA_EVENT_tree->SetBranchAddress("YCov11", &YCov11);
//    MAMBA_EVENT_tree->SetBranchAddress("XrecoT0", &XrecoT0);
//    MAMBA_EVENT_tree->SetBranchAddress("YrecoT0", &YrecoT0);
//    MAMBA_EVENT_tree->SetBranchAddress("XrecoT1", &XrecoT1);
//    MAMBA_EVENT_tree->SetBranchAddress("YrecoT1", &YrecoT1);
//    MAMBA_EVENT_tree->SetBranchAddress("XrecoT2", &XrecoT2);
//    MAMBA_EVENT_tree->SetBranchAddress("YrecoT2", &YrecoT2);
//    MAMBA_EVENT_tree->SetBranchAddress("XrecoT3", &XrecoT3);
//    MAMBA_EVENT_tree->SetBranchAddress("YrecoT3", &YrecoT3);
//    MAMBA_EVENT_tree->SetBranchAddress("XfitT0", &XfitT0);
//    MAMBA_EVENT_tree->SetBranchAddress("YfitT0", &YfitT0);
//    MAMBA_EVENT_tree->SetBranchAddress("XfitT1", &XfitT1);
//    MAMBA_EVENT_tree->SetBranchAddress("YfitT1", &YfitT1);
//    MAMBA_EVENT_tree->SetBranchAddress("XfitT2", &XfitT2);
//    MAMBA_EVENT_tree->SetBranchAddress("YfitT2", &YfitT2);
//    MAMBA_EVENT_tree->SetBranchAddress("XfitT3", &XfitT3);
//    MAMBA_EVENT_tree->SetBranchAddress("YfitT3", &YfitT3);
///--------
//    MAMBA_EVENT_tree->SetBranchAddress("nIntersect", &nIntersect);
//    MAMBA_EVENT_tree->SetBranchAddress("intTrkID", &intTrkID);
//    MAMBA_EVENT_tree->SetBranchAddress("intLayer", &intLayer);
//    MAMBA_EVENT_tree->SetBranchAddress("intInside", &intInside);
//    MAMBA_EVENT_tree->SetBranchAddress("intXPos", &intXPos);
//    MAMBA_EVENT_tree->SetBranchAddress("intYPos", &intYPos);
//    MAMBA_EVENT_tree->SetBranchAddress("intZPos", &intZPos);
//    MAMBA_EVENT_tree->SetBranchAddress("intUPos", &intUPos);
//    MAMBA_EVENT_tree->SetBranchAddress("intVPos", &intVPos);
///--------
    MAMBA_EVENT_tree->SetBranchAddress("nSP", &nSP);
    MAMBA_EVENT_tree->SetBranchAddress("spLayer", &spLayer);
    MAMBA_EVENT_tree->SetBranchAddress("spXPos", &spXPos);
    MAMBA_EVENT_tree->SetBranchAddress("spYPos", &spYPos);
    MAMBA_EVENT_tree->SetBranchAddress("spZPos", &spZPos);
    MAMBA_EVENT_tree->SetBranchAddress("spXErr", &spXErr);
    MAMBA_EVENT_tree->SetBranchAddress("spYErr", &spYErr);
    MAMBA_EVENT_tree->SetBranchAddress("spZErr", &spZErr);
    MAMBA_EVENT_tree->SetBranchAddress("spUPos", &spUPos);
    MAMBA_EVENT_tree->SetBranchAddress("spVPos", &spVPos);
    MAMBA_EVENT_tree->SetBranchAddress("spIsOnTrk", &spIsOnTrk);
}
