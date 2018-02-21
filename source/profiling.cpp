
#include "../includes/profiling.h"

profiling::profiling() {
    init();
}

profiling::~profiling() = default;

void profiling::init() {
    X_Axis_Point = 0;
    X_Axis_Point_Err = 0;
    Y_Axis_Point = 0;
    Y_Axis_Point_Err = 0;
    first = 0;
    second = 0;
    fit_range_L = 0;
    fit_range_R = 0;
}

void profiling::bins_filling(TH2F *histo, Int_t i) {
    X_Axis_Point = histo->GetXaxis()->GetBinCenter(i);
    TString name = {"bin_XAxis_"};
    TString num;
    num.Form("%f", X_Axis_Point);
    TString full = name + num;
    bin_proj = new TH1D();
    bin_proj = histo->ProjectionY(full, i, i+1);
}

void profiling::finding_fit_range() {
    Int_t max_bin = bin_proj->GetMaximumBin();
    left_bin = 0;
    right_bin = 290;
    Int_t i = 1;
    Int_t j = 1;
    do {
        right_bin = max_bin + i;
        i++;
    } while (bin_proj->GetBinContent(max_bin) * 0.2 < bin_proj->GetBinContent(max_bin + i));

    do {
        left_bin = max_bin - j;
        j++;
    } while (bin_proj->GetBinContent(max_bin) * 0.2 < bin_proj->GetBinContent(max_bin - j));

    fit_range_L = bin_proj->GetBinCenter(left_bin - 1);
    fit_range_R = bin_proj->GetBinCenter(right_bin + 1);
}

void profiling::fitting_bin_hist() {
    bin_proj->Fit("gaus", "QR", "SAME", fit_range_L, fit_range_R);
    TF1 *fit_func = bin_proj->GetFunction("gaus");
    Y_Axis_Point = fit_func->GetParameter(1);
    sigma = fit_func->GetParameter(2);
    Y_Axis_Point_Err = sigma / sqrt(bin_proj->GetEntries());
    sigma_err = fit_func->GetParError(2);
    fit = gMinuit->fCstatu;
}

void profiling::main_algorithm(Int_t start, Int_t end, TH2F *histo, histos *hist) {
    gStyle->SetOptFit(1111);
    first = histo->GetXaxis()->GetBinCenter(1);
    second = histo->GetXaxis()->GetBinCenter(2);
    X_Axis_Point_Err = (second - first) / 2.0;
    Int_t j = 0;
    for (int i = 0; i < end - start; ++i) {
        bins_filling(histo, start + i);
        finding_fit_range();
        fitting_bin_hist();
        if (strncmp(fit.Data(), res.Data(), 5) != 0) continue;
        (hist->g)->SetPoint(j, X_Axis_Point, Y_Axis_Point);
        (hist->g)->SetPointError(j, X_Axis_Point_Err, Y_Axis_Point_Err);
        (hist->sigma)->SetPoint(j, X_Axis_Point, sigma);
        (hist->sigma)->SetPointError(j, X_Axis_Point_Err, sigma_err);
        if (j % 10 == 0) {
            TCanvas *c1 = new TCanvas();
            TString name = {"img/test/h_x_"};
            Long_t num = j;
            TString end = {".pdf"};
            TString full = name + num + end;
            c1->cd();
            bin_proj->Draw();
            c1->SaveAs(full);
        }
        j++;
    }
}



