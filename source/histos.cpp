#include <TStyle.h>
#include "../includes/histos.h"

histos::histos() {
    Init_histos();
}

histos::~histos() = default;

void histos::Init_histos() {
// V-shapes (the r(t) relation)
// short straw = S
    vshapeUS = new TH2F("vshapeU_Short", "vshapeU_Short", 150, 0.0, 0.0, 150, 0.0, 0.0);
    vshapeUS->GetXaxis()->SetTitle("U (mm)");
    vshapeUS->GetYaxis()->SetTitle("T (ns)");
    vshapeVS = new TH2F("vshapeV_Short", "vshapeV_Short", 150, 0.0, 0.0, 150, 0.0, 0.0);
// long straw = L
    vshapeUL = new TH2F("vshapeU_Long", "vshapeU_Long", 150, 0.0, 0.0, 150, 0.0, 0.0);
    vshapeUL->GetXaxis()->SetTitle("U (mm)");
    vshapeUL->GetYaxis()->SetTitle("T (ns)");
    vshapeVL = new TH2F("vshapeV_Long", "vshapeV_Long", 150, 0.0, 0.0, 150, 0.0, 0.0);

// V-shapes (the r(t) relation) CLEAR
// short straw = S
    vshapeUS_clear = new TH2F("vshapeU_Short_clear", "vshapeU_Short_clear", 150, 0.0, 0.0, 150, 0.0, 0.0);
    vshapeUS_clear->GetXaxis()->SetTitle("U (mm)");
    vshapeUS_clear->GetYaxis()->SetTitle("T (ns)");
    vshapeVS_clear = new TH2F("vshapeV_Short_clear", "vshapeV_Short_clear", 150, 0.0, 0.0, 150, 0.0, 0.0);
// long straw = L
    vshapeUL_clear = new TH2F("vshapeU_Long_clear", "vshapeU_Long_clear", 150, 0.0, 0.0, 150, 0.0, 0.0);
    vshapeUL_clear->GetXaxis()->SetTitle("U (mm)");
    vshapeUL_clear->GetYaxis()->SetTitle("T (ns)");
    vshapeVL_clear = new TH2F("vshapeV_Long_clear", "vshapeV_Long_clear", 150, 0.0, 0.0, 150, 0.0, 0.0);

    Resolution_L = new TH1F("Straw Resolution L", "Straw Resolution L", 200, -4, 4);
    Resolution_L->GetXaxis()->SetTitle("Residual (mm)");
    Resolution_S = new TH1F("Straw Resolution S", "Straw Resolution S", 200, -4, 4);
    Resolution_S->GetXaxis()->SetTitle("Residual (mm)");

    n_hits_long_straw = new TH1F("Hits in L straw", "Hits in L straw", 10, 0.0, 0.0);
    n_hits_short_straw = new TH1F("Hits in S straw", "Hits in S straw", 10, 0.0, 0.0);
}

void histos::Drawing_histos() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);

    TFile myfile("MY_OUT.root", "RECREATE");
    vshapeUL->Write("Before L");
    vshapeUS->Write("Before S");
    vshapeUL_clear->Write("After L");
    vshapeUS_clear->Write("After S");
    straw_L->Write("Profile L");
    straw_S->Write("Profile S");

    n_hits_short_straw->Write("hits in S straw");
    n_hits_long_straw->Write("hits in L straw");

    vshapeUS_proj->Write("Proj_x S");
    vshapeUL_proj->Write("Proj_x L");
    vshapeUS_proj_Y->Write("Proj_y S");
    vshapeUL_proj_Y->Write("Proj_y L");

    Resolution_L->Write("Resolution L");
    Resolution_S->Write("Resolution S");

    TCanvas *shapes = new TCanvas("Shapes", "Shapes", 1440, 1000);
    shapes->Divide(2, 2);
    shapes->cd(1);
    vshapeUS->Draw("COLZ");
    shapes->cd(2);
    vshapeUL->Draw("COLZ");
    shapes->cd(3);
    vshapeUS_clear->Draw("COLZ");
    shapes->cd(4);
    vshapeUL_clear->Draw("COLZ");

    TCanvas *profiles = new TCanvas("Profiles", "Profiles", 1440, 500);
    profiles->Divide(2, 1);
    profiles->cd(1);
    straw_S->Draw();
    profiles->cd(2);
    straw_L->Draw();

    TCanvas *projections = new TCanvas("Projections", "Projections", 1440, 1000);
    projections->Divide(2, 2);
    projections->cd(1);
    vshapeUS_cleaned->Draw("COLZ");
    projections->cd(2);
    vshapeUS_proj_Y->Draw();
    projections->cd(3);
    vshapeUL_cleaned->Draw("COLZ");
    projections->cd(4);
    vshapeUL_proj_Y->Draw();

    TCanvas *resol = new TCanvas("Resolution", "Resolution", 1440, 500);
    resol->Divide(2, 1);
    resol->cd(1);
    Resolution_S->Draw();
    resol->cd(2);
    Resolution_L->Draw();

    TString way("/Users/andrew_zelenov/Documents/SHiP/new_data_type_ana/img/");
    TString shape_name("Shapes.pdf");
    TString prof_name("Profiles.pdf");
    TString proj_name("Projections.pdf");
    TString resol_name("Resolution.pdf");

    shapes->SaveAs(way + shape_name, "Q");
    profiles->SaveAs(way + prof_name, "Q");
    projections->SaveAs(way + proj_name, "Q");
    resol->SaveAs(way + resol_name, "Q");
}
