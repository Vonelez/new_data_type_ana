#include <TStyle.h>
#include "../includes/histos.h"

histos::histos(Long_t run_num) {
    run = run_num;
    Init_histos();
}

histos::~histos() = default;

void histos::Init_histos() {
// V-shapes (the r(t) relation)
// short straw = S
    vshapeUS = new TH2F("vshapeU_Short", "vshapeU_Short", 600, 0.0, 0.0, 900, 0.0, 0.0);
    vshapeUS->GetXaxis()->SetTitle("U (mm)");
    vshapeUS->GetYaxis()->SetTitle("T (ns)");
    vshapeVS = new TH2F("vshapeV_Short", "vshapeV_Short", 150, 0.0, 0.0, 300, 0.0, 0.0);
// long straw = L
    vshapeUL = new TH2F("vshapeU_Long", "vshapeU_Long", 150, 0.0, 0.0, 300, 0.0, 0.0);
    vshapeUL->GetXaxis()->SetTitle("U (mm)");
    vshapeUL->GetYaxis()->SetTitle("T (ns)");
    vshapeVL = new TH2F("vshapeV_Long", "vshapeV_Long", 150, 0.0, 0.0, 300, 0.0, 0.0);

// V-shapes (the r(t) relation) CLEAR
// short straw = S
    vshapeUS_clear = new TH2F("vshapeU_Short_clear", "vshapeU_Short_clear", 150, 0.0, 0.0, 300, 0.0, 0.0);
    vshapeUS_clear->GetXaxis()->SetTitle("U (mm)");
    vshapeUS_clear->GetYaxis()->SetTitle("T (ns)");
    vshapeVS_clear = new TH2F("vshapeV_Short_clear", "vshapeV_Short_clear", 150, 0.0, 0.0, 300, 0.0, 0.0);
// long straw = L
    vshapeUL_clear = new TH2F("vshapeU_Long_clear", "vshapeU_Long_clear", 150, 0.0, 0.0, 300, 0.0, 0.0);
    vshapeUL_clear->GetXaxis()->SetTitle("U (mm)");
    vshapeUL_clear->GetYaxis()->SetTitle("T (ns)");
    vshapeVL_clear = new TH2F("vshapeV_Long_clear", "vshapeV_Long_clear", 150, 0.0, 0.0, 300, 0.0, 0.0);

    Resolution_L = new TH1F("Straw Resolution L", "Straw Resolution L", 1000, -4, 4);
    Resolution_L->GetXaxis()->SetTitle("Residual (mm)");
    Resolution_S = new TH1F("Straw Resolution S", "Straw Resolution S", 1000, -4, 4);
    Resolution_S->GetXaxis()->SetTitle("Residual (mm)");
    Resolution_S_bins = new TH1F("Straw Resolution S (bins)", "Straw Resolution S (bins)", 1000, -4, 4);
    Resolution_S_bins->GetXaxis()->SetTitle("Residual (mm)");

    n_hits_long_straw = new TH1F("Hits in L straw", "Hits in L straw", 10, 0.0, 0.0);
    n_hits_short_straw = new TH1F("Hits in S straw", "Hits in S straw", 10, 0.0, 0.0);

    g = new TGraphErrors();
    g->GetXaxis()->SetTitle("U (mm)");
    g->GetYaxis()->SetTitle("T (ns)");
    g->SetTitle("T = f(U)");
    sigma = new TGraphErrors();
    sigma->GetXaxis()->SetTitle("U (mm)");
    sigma->GetYaxis()->SetTitle("#sigma_{T} (ns)");
    sigma->SetTitle("#sigma_{T} = f(U)");
    U_res = new TGraphErrors();
    U_res->GetXaxis()->SetTitle("U (mm)");
    U_res->GetYaxis()->SetTitle("#sigma_{U} (mm)");
    U_res->SetTitle("#sigma_{U} = f(#sigma_{T})");
    resolgeom = new TGraphErrors();
    resolgeom->GetXaxis()->SetTitle("U (mm)");
    resolgeom->GetYaxis()->SetTitle("#sigma_{U} (mm)");
    resolgeom->SetTitle("#sigma_{U} = f(#sigma_{T})");
    chi_ndf = new TGraph();
    chi_ndf->GetXaxis()->SetTitle("U (mm)");
    chi_ndf->GetYaxis()->SetTitle("#chi^{2} / NDF");
    ndf = new TGraph();
    ndf->GetXaxis()->SetTitle("U (mm)");
    ndf->GetYaxis()->SetTitle("NDF");
}

void histos::Drawing_histos() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);

    TString Name("MY_OUT_");
    TString suf(".root");
    TFile myfile(Name + run + suf, "RECREATE");
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
    Resolution_S_bins->Write("Resolution S (bins)");

    g->Write("Parabola");
    gnew->Write("Parabola2");
    sigma->Write("Sigma(x)");
    U_res->Write("RESOL");
    resolgeom->Write("RESOLgeom");
    chi_ndf->Write("chi/ndf");
    ndf->Write("ndf");

    TCanvas *shapes = new TCanvas("Shapes", "Shapes", 1400, 1000);
    shapes->cd();
    vshapeUS->Draw("COLZ");

    resolgeom->GetYaxis()->SetRangeUser(0.0 , 0.5);
    resolgeom->SetMarkerSize(1.5);
    resolgeom->SetMarkerColor(kGreen+2);
    resolgeom->SetMarkerStyle(47);
    TCanvas *resol = new TCanvas("Resolution", "Resolution", 1400, 1000);
    resol->cd(0);
    resolgeom->Draw("AP");


    gnew->SetLineColor(kGreen+2);
    //making legends for hists
    auto leg1 = new TLegend(0.55, 0.8, 0.9, 0.9);
    leg1->AddEntry(gnew, "bin-by-bin Gaussian fit mean ", "l");
    leg1->AddEntry(straw_S, "VShape TProfile", "l");

    TCanvas *prof_graph = new TCanvas("TProfile vs bin by bin fit", "TProfile vs bin by bin fit", 720, 500);
    prof_graph->cd();
    gnew->Draw("AP");
    straw_S->Draw("SAME");
    leg1->Draw();
    prof_graph->Write("2plots");

    chi_ndf->SetMarkerColor(kGreen+2);
    chi_ndf->SetMarkerStyle(47);
    chi_ndf->SetMarkerSize(1.5);
    ndf->SetMarkerColor(kRed+2);
    ndf->SetMarkerStyle(22);
    ndf->SetMarkerSize(1.5);

    TCanvas *chi_ndf_canv = new TCanvas("#chi^{2} / NDF", "#chi^{2} / NDF", 1200, 1400);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    Float_t small = 1e-5;
    chi_ndf_canv->Divide(1,2,small,small);
    chi_ndf_canv->cd(1);
    gPad->SetBottomMargin(small);
    gPad->SetLogy();
    chi_ndf->Draw("AP");
    chi_ndf_canv->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    ndf->Draw("AP");

    TString path("/Users/andrew_zelenov/Documents/SHiP/new_data_type_ana/img");
    TString way(path + run + "/");
    TString shape_name("Shapes.pdf");
    TString resol_name("Resolution.pdf");
    TString chi("chi_ndf.pdf");
    TString resolution("Resol_RUN_");

    shapes->SaveAs(way + shape_name, "Q");
    resol->SaveAs(way + resol_name, "Q");
    chi_ndf_canv->SaveAs(way + chi, "Q");
    resolgeom->SaveAs(way + resolution + run + ".root");
}
