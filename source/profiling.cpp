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
    bin_proj = NULL;
}

void profiling::bins_filling(TH2F *histo, Int_t i) {
    X_Axis_Point = histo->GetXaxis()->GetBinCenter(i);
    bin_proj = histo->ProjectionY("Bin proj", i, i+1);
}

void profiling::finding_fit_range() {
    Int_t max_bin = bin_proj->GetMaximumBin();
    Int_t left_bin = 0;
    Int_t right_bin = 0;
    Int_t i = 1;
    Int_t j = 1;
    while (bin_proj->GetBinContent(max_bin) * 0.2 < bin_proj->GetBinContent(max_bin + i)) {
        right_bin = max_bin + i;
        i++;
    }

    while (bin_proj->GetBinContent(max_bin) * 0.2 < bin_proj->GetBinContent(max_bin - j)) {
        left_bin = max_bin - j;
        j++;
    }

    fit_range_L = bin_proj->GetBinCenter(left_bin);
    fit_range_R = bin_proj->GetBinCenter(right_bin);
}

void profiling::fitting_bin_hist(histos *hist, Int_t i) {
    bin_proj->Fit("gaus", "QR", "SAME", fit_range_L, fit_range_R);
    TF1 *fit_func = bin_proj->GetFunction("gaus");
    Y_Axis_Point = fit_func->GetParameter(1);
    Y_Axis_Point_Err = fit_func->GetParameter(2);

    (hist->g)->SetPoint(i, X_Axis_Point, Y_Axis_Point);
    (hist->g)->SetPointError(i, X_Axis_Point_Err, Y_Axis_Point_Err);
}

void profiling::main_algorithm(Int_t start, Int_t end, TH2F *histo, histos *hist) {
    gStyle->SetOptFit(1111);
    first = histo->GetXaxis()->GetBinCenter(1);
    second = histo->GetXaxis()->GetBinCenter(2);
    X_Axis_Point_Err = (second - first) / 2.0;

    for (int i = 0; i < end - start; ++i) {
        bins_filling(histo, start + i);
        finding_fit_range();
        fitting_bin_hist(hist, i);
        if (i % 10 == 0) {
            TCanvas *c1 = new TCanvas();
            char *histname = new char[10];
            sprintf(histname, "img/h_x_%d.pdf", i);
            c1->cd();
            bin_proj->Draw();
            c1->SaveAs(histname);
        }
    }
}



