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
    right_bin = 150;
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
    chi = fit_func->GetChisquare();
    ndf = fit_func->GetNDF();
    fit = gMinuit->fCstatu;
}

void profiling::finding_range(histos *hist) {
    Int_t N = (hist->g)->GetN();
    Double_t *x, *y;
    x = (hist->g)->GetX();
    y = (hist->g)->GetY();
    Int_t max_bin = TMath::LocMax(N, y);
    Max1 = x[max_bin];
    Max2 = 0;

    if (Max1 < 0.) {
        for (int i = max_bin; i < N; i++) {

            if (abs(y[i] / y[i+1]) > 5.) {
                Max2 = x[i];
//                cout << x[i] << "  " << abs(y[i] / y[i+1]) << endl;
            }
            if (abs(x[i] - x[max_bin]) > 21.) break;
        }
    } else {
        for (int i = max_bin; i < N; --i) {
            if (abs(y[i] / y[i-1]) > 5.) {
                Max2 = x[i];
                break;
            }
            if (abs(x[i] - x[max_bin]) > 21.) break;
        }
    }

    if (Max1 > Max2) {
        Double_t c = Max1;
        Max1 = Max2;
        Max2 = c;
    }
}

void profiling::main_algorithm(Int_t start, Int_t end, TH2F *histo, histos *hist) {
    gStyle->SetOptFit(1111);
    first = histo->GetXaxis()->GetBinCenter(1);
    second = histo->GetXaxis()->GetBinCenter(2);
    X_Axis_Point_Err = (second - first) / 2.0;
    Int_t j = 0;
    for (int i = 0; i < end - start; ++i) {
        bins_filling(histo, start + i);
        if (bin_proj->GetEntries() == 0) continue;
        finding_fit_range();
        fitting_bin_hist();
        if (strncmp(fit.Data(), res.Data(), 5) != 0) continue;
//        if (Y_Axis_Point > 800.) continue; ///I'm KOSTYL' for run 222!!!!
        (hist->g)->SetPoint(j, X_Axis_Point, Y_Axis_Point);
        (hist->g)->SetPointError(j, X_Axis_Point_Err, Y_Axis_Point_Err);
        (hist->sigma)->SetPoint(j, X_Axis_Point, sigma);
        (hist->sigma)->SetPointError(j, X_Axis_Point_Err, sigma_err);
        if (j % 10 == 0) {
            TCanvas *c1 = new TCanvas();
            TString folder("img");
            TString path("/test/h_x_");
            TString name = folder + hist->run + path;
            Long_t num = j;
            TString end = {".pdf"};
            TString full = name + num + end;
            c1->cd();
            bin_proj->Draw();
            c1->SaveAs(full);
        }
        if (ndf == 0) {
            (hist->chi_ndf)->SetPoint(j, X_Axis_Point, 50);
            (hist->ndf)->SetPoint(j, X_Axis_Point, ndf);
        } else {
            Double_t chindf = chi / ndf;
            (hist->chi_ndf)->SetPoint(j, X_Axis_Point, chindf);
            (hist->ndf)->SetPoint(j, X_Axis_Point, ndf);
        }
        j++;
    }
}

void profiling::maean_and_wm(histos *hist) {
    Int_t N = (hist->gnew)->GetN();
    Double_t *x, *y;
    x = (hist->gnew)->GetX();
    y = (hist->gnew)->GetY();
    Double_t min = 1000.;
    Double_t coord = 0;

    for (int i = 0; i < N; ++i) {
        if (x[i] >= Max1 && x[i] <= Max2) {
            if (y[i] < min) {
                min = y[i];
                coord = x[i];
            }
        }
    }

    cout << Max1 << " --- " << Max2 << " --- " << coord << endl;

    ofstream myfile;
    TString folder("img");
    TString path("/run_");
    TString suf(".txt");
    TString name = folder + hist->run + path + hist->run + suf;
    myfile.open(name);
    myfile << "R-L \t L \t R \t C \n";
    myfile << Max2 - Max1 << "\t" << Max1 << "\t" << Max2 << "\t" << coord << "\t" << endl;
    myfile.close();

}

void profiling::geometric_resol(histos *hist) {
    Double_t *x = (hist->sigma)->GetX();
    Double_t *y = (hist->sigma)->GetY();
    Double_t *y_2 = (hist->gnew)->GetY();

    Int_t Num = (hist->sigma)->GetN();
    int j = 0;

    for (int i = 1; i < Num - 1; ++i) {
        Double_t derivative = (1./2.) * ((y_2[i+1] - y_2[i]) / (x[i+1] - x[i]) + (y_2[i] - y_2[i-1]) / (x[i] - x[i-1]));
        Double_t resol = y[i] / abs(derivative);
        if (resol > 100) continue;

        (hist->resolgeom)->SetPoint(j, x[i], resol);
        j++;
    }
}



