#include <iostream>
#include "includes/linking_data.h"

int main() {
    TString dir ("/Users/andrew_zelenov/Documents/SHiP/DATA/SOFT/");
    TString Name ("Synchronized_data_run_222.root");

    TString RootName (dir + Name);

    TFile *AnaFile = TFile::Open(RootName,"read");
    TTree *mt;
    TTree *st;
    AnaFile->GetObject("tree", mt);
    AnaFile->GetObject("ADC1", st);

    linking_data(st, mt);

    AnaFile->TFile::Close();
    return 0;
}