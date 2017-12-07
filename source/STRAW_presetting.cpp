#include "../includes/STRAW_presetting.h"

STRAW_presetting::STRAW_presetting(TTree *STRAW_EVENT_tree) {
    Init_STRAW_presetting();
    SetUpSTRAWTree(STRAW_EVENT_tree);
}

STRAW_presetting::~STRAW_presetting() = default;

void STRAW_presetting::Init_STRAW_presetting() {

    PCclock = 0;
    pattern = 0;
    TDCtrig = 0;
    t = 0;
    nt00 = 0;
    tt00 = 0;
    nt01 = 0;
    tt01 = 0;
    nt02 = 0;
    tt02 = 0;
    nt03 = 0;
    tt03 = 0;
    nt04 = 0;
    tt04 = 0;
    nt05 = 0;
    tt05 = 0;
    nt06 = 0;
    tt06 = 0;
    nt08 = 0;
    tt08 = 0;
}

void STRAW_presetting::SetUpSTRAWTree(TTree *STRAW_EVENT_tree) {
    STRAW_EVENT_tree->SetBranchAddress("PCclock", &PCclock);
    STRAW_EVENT_tree->SetBranchAddress("pattern", &pattern);
    STRAW_EVENT_tree->SetBranchAddress("TDCtrig", &TDCtrig);
    STRAW_EVENT_tree->SetBranchAddress("t", &t);
    STRAW_EVENT_tree->SetBranchAddress("nt00", &nt00);
    STRAW_EVENT_tree->SetBranchAddress("tt00", &tt00);
    STRAW_EVENT_tree->SetBranchAddress("nt01", &nt01);
    STRAW_EVENT_tree->SetBranchAddress("tt01", &tt01);
    STRAW_EVENT_tree->SetBranchAddress("nt02", &nt02);
    STRAW_EVENT_tree->SetBranchAddress("tt02", &tt02);
    STRAW_EVENT_tree->SetBranchAddress("nt03", &nt03);
    STRAW_EVENT_tree->SetBranchAddress("tt03", &tt03);
    STRAW_EVENT_tree->SetBranchAddress("nt04", &nt04);
    STRAW_EVENT_tree->SetBranchAddress("tt04", &tt04);
    STRAW_EVENT_tree->SetBranchAddress("nt05", &nt05);
    STRAW_EVENT_tree->SetBranchAddress("tt05", &tt05);
    STRAW_EVENT_tree->SetBranchAddress("nt06", &nt06);
    STRAW_EVENT_tree->SetBranchAddress("tt06", &tt06);
    STRAW_EVENT_tree->SetBranchAddress("nt08", &nt08);
    STRAW_EVENT_tree->SetBranchAddress("tt08", &tt08);
}