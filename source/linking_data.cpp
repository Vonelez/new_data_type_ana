#include "../includes/linking_data.h"

linking_data::linking_data(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree) {
    init();
    merging(STRAW_EVENT_tree, MAMBA_EVENT_tree);
}

linking_data::~linking_data() = default;

void linking_data::init() {
    conversion = 25.0 / 1024;
    mam_nEntries = 0;
    N = 0;
    count_processed = 0;
    count_good = 0;
    count_very_good = 0;
    Max_Short = -1;
    Max_Long = -1;
    a_L = 0;
    a_L_err = 0;
    b_L = 0;
    b_L_err = 0;
    c_L = 0;
    c_L_err = 0;
    a_S = 0;
    a_S_err = 0;
    b_S = 0;
    b_S_err = 0;
    c_S = 0;
    c_S_err = 0;

    fit_prof_L = new TF1("fit_prof_L", "[0]+[1]*x+[2]*x*x");
    fit_prof_S = new TF1("fit_prof_S", "[0]+[1]*x+[2]*x*x");

    vertex_err_L = 0;
    vertex_err_S = 0;

    Discrim2_L = 0;
    Discrim2_S = 0;
    timing_L.clear();
    timing_S.clear();
    StrawPoint_L.clear();
    StrawPoint_S.clear();
}

void linking_data::merging(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree) {
    auto *straw = new STRAW_presetting(STRAW_EVENT_tree);
    auto *mamba = new MAMBA_presetting(MAMBA_EVENT_tree);
    auto *hist = new histos();
    auto *profile = new profiling();

    mam_nEntries = (Int_t) MAMBA_EVENT_tree->GetEntries();
    N = mam_nEntries;
    cout << "requesting " << N << " events..." << endl;

    cout << "Begin analysis..." << endl;
    for (Int_t i = 0; i < N; i++) {
/// progress bar
        int barWidth = 70;
        cout << "\033[0;32m[\033[0m";
        auto pos = (int) (barWidth * ((float) i / (float) N));
        for (int bar = 0; bar < barWidth; ++bar) {
            if (bar < pos) cout << "\033[0;32m=\033[0m";
            else if (bar == pos) cout << "\033[0;32m>\033[0m";
            else cout << " ";
        }
        cout << "\033[0;32m]\033[0m " << "\033[0;31m" << int(((float) i / (float) N) * 100.0 + 1) << "%\r \033[0m";
        cout.flush();
/// progress bar

        STRAW_EVENT_tree->GetEntry(i);
        MAMBA_EVENT_tree->GetEntry(i);

        filling_hists(hist, straw, mamba);
        track_ana(mamba, straw);

        count_processed++;
    }
    cout << endl;

    (hist->vshapeUS_proj) = (hist->vshapeUS)->ProjectionX();
    // need to make copies of this hist because range will be changed in Limits_S/L
    (hist->vshapeUL_proj) = (hist->vshapeUL)->ProjectionX();

    limits_S(hist, *range_S);
    limits_L(hist, *range_L);
    left_limit_S = range_S[0] - 1;
    left_limit_L = range_L[0];
    right_limit_S = range_S[1] + 3;
    right_limit_L = range_L[1] - 3;

    cout << "------- " << StrawPoint_S.size() << endl;

    for (int i = 0; i < StrawPoint_S.size(); i++) {
        if (StrawPoint_S.at(i) >= left_limit_S && StrawPoint_S.at(i) <= right_limit_S) {
            (hist->vshapeUS_clear)->Fill(StrawPoint_S.at(i), timing_S.at(i));
            if (StrawPoint_S.at(i) >= left_limit_S + 1.5 && StrawPoint_S.at(i) <= right_limit_S - 1.5) {
                StrawPoint_S_after_cut.push_back(StrawPoint_S.at(i));
                timing_S_after_cut.push_back(timing_S.at(i));
            }
            if (StrawPoint_S.at(i) >= -5.0 && StrawPoint_S.at(i) <= -4.8) {
                (hist->check_prof)->Fill(timing_S.at(i));
            }
        }
    }


    for (int i = 0; i < StrawPoint_L.size(); i++) {
        if (StrawPoint_L.at(i) >= left_limit_L && StrawPoint_L.at(i) <= right_limit_L) {
            (hist->vshapeUL_clear)->Fill(StrawPoint_L.at(i), timing_L.at(i));
            if (StrawPoint_L.at(i) >= left_limit_L + 0.5 && StrawPoint_L.at(i) <= right_limit_L - 0.5) {
                StrawPoint_L_after_cut.push_back(StrawPoint_L.at(i));
                timing_L_after_cut.push_back(timing_L.at(i));
            }
        }
    }

    makingProfile(hist);
    draft(hist, profile);
    StrawResolution_S(hist);
    StrawResolution_L(hist);

    cout << "Linkdata: Events processed: " << count_processed << endl;
    cout << "Linkdata: Events with a track: " << count_good << " of which " << count_very_good << " with 8 clusters"
         << endl;

    hist->Drawing_histos();
}

void linking_data::track_ana(MAMBA_presetting *mamba, STRAW_presetting *straw) {
    if (mamba->ntrk >= 1) {
        count_good++;
        if (mamba->ntrk == 1) {
            count_very_good++;

            if ((straw->nt00 == 1) && (straw->nt08 == 1)) {
                timing_S.push_back(conversion * (straw->tt00 - straw->tt08));
                StrawPoint_S.push_back(mamba->spUPos[0] * 10);
            }

            if ((straw->nt01 == 1) && (straw->nt08 == 1)) {
                timing_L.push_back(conversion * (straw->tt01 - straw->tt08));
                StrawPoint_L.push_back(mamba->spUPos[0] * 10);
            }
        }
    } else {
        return;
    }
}

void linking_data::makingProfile(histos *hist) {
    (hist->vshapeUS_cleaned) = (TH2F*)(hist->vshapeUS_clear)->Clone("cleaned_S");
    (hist->vshapeUL_cleaned) = (TH2F*)(hist->vshapeUL_clear)->Clone("cleaned_L");

    for (int i = 0; i < (hist->vshapeUL_cleaned)->GetSize(); i++) {
        if ((hist->vshapeUL_cleaned)->GetBinContent(i) > Max_Long) {
            Max_Long = (hist->vshapeUL_cleaned)->GetBinContent(i);
        }
    }
    for (int i = 0; i < (hist->vshapeUS_cleaned)->GetSize(); i++) {
        if ((hist->vshapeUS_cleaned)->GetBinContent(i) > Max_Short) {
            Max_Short = (hist->vshapeUS_cleaned)->GetBinContent(i);
        }
    }

    Max_Long *= 0.05;
    Max_Short *= 0.01;

    for (int i = 0; i < (hist->vshapeUL_cleaned)->GetSize(); i++) {
        if ((hist->vshapeUL_cleaned)->GetBinContent(i) < Max_Long && (hist->vshapeUL_cleaned)->GetBinContent(i) != 0) {
            (hist->vshapeUL_cleaned)->SetBinContent(i, 0);
        }
    }
    for (int i = 0; i < (hist->vshapeUS_cleaned)->GetSize(); i++) {
        if ((hist->vshapeUS_cleaned)->GetBinContent(i) < Max_Short && (hist->vshapeUS_cleaned)->GetBinContent(i) != 0) {
            (hist->vshapeUS_cleaned)->SetBinContent(i, 0);
        }
    }

    (hist->vshapeUS_proj_Y) = (hist->vshapeUS_cleaned)->ProjectionY();
    (hist->vshapeUL_proj_Y) = (hist->vshapeUL_cleaned)->ProjectionY();

    (hist->straw_L) = (hist->vshapeUL_clear)->ProfileX("L Straw Profile", -1, -1, "o");
    (hist->straw_S) = (hist->vshapeUS_clear)->ProfileX("S Straw Profile", -1, -1, "o");

    (hist->straw_L)->Fit(fit_prof_L, "QR", "SAME", left_limit_L + 0.5, right_limit_L - 0.5);
    (hist->straw_S)->Fit(fit_prof_S, "QR", "SAME", left_limit_S + 1.5, right_limit_S - 1.5);

    a_L = fit_prof_L->GetParameter(2);
    b_L = fit_prof_L->GetParameter(1);
    c_L = fit_prof_L->GetParameter(0);
    a_L_err = fit_prof_L->GetParError(2);
    b_L_err = fit_prof_L->GetParError(1);
    c_L_err = fit_prof_L->GetParError(0);

    a_S = fit_prof_S->GetParameter(2);
    b_S = fit_prof_S->GetParameter(1);
    c_S = fit_prof_S->GetParameter(0);
    a_S_err = fit_prof_S->GetParError(2);
    b_S_err = fit_prof_S->GetParError(1);
    c_S_err = fit_prof_S->GetParError(0);

    vertex_err_L = (1.0 / 2.0) * (abs((1 / a_L) * b_L_err) + abs((b_L / (a_L * a_L)) * a_L_err));
    vertex_err_S = (1.0 / 2.0) * (abs((1 / a_S) * b_S_err) + abs((b_S / (a_S * a_S)) * a_S_err));

    cout << "Parabolic parameters (For Long Straw): " << a_L << " " << b_L << " " << c_L << endl;
    cout << "Parabolic parameters (For Short Straw): " << a_S << " " << b_S << " " << c_S << endl;
    cout << "Parabolic parameters errors (For Long Straw): " << a_L_err << " " << b_L_err << " " << c_L_err << endl;
    cout << "Parabolic parameters errors (For Short Straw): " << a_S_err << " " << b_S_err << " " << c_S_err << endl;

    cout << "=====> Vertex (Long Straw)= " << -b_L / (2.0 * a_L) << " " << vertex_err_L << endl;
    cout << "=====> Vertex (Short Straw)= " << -b_S / (2.0 * a_S) << " " << vertex_err_S << endl;

}
void linking_data::StrawResolution_L(histos *hist) {
    Double_t mindev = 99999, maxdev = -9999999;

    for (int i = 0; i < timing_L_after_cut.size(); ++i) {
        Discrim2_L = (b_L * b_L) - 4.0 * a_L * (c_L - timing_L_after_cut.at(i));
        if (Discrim2_L < 0) continue;
        Double_t discr = sqrt(Discrim2_L);
        Double_t tekdev_plus = (-b_L + discr) / (2.0 * a_L);
        Double_t tekdev_minus = (-b_L - discr) / (2.0 * a_L);
        Double_t tekdev;

        if (fabs(tekdev_minus - StrawPoint_L_after_cut.at(i)) > fabs(tekdev_plus - StrawPoint_L_after_cut.at(i))) {
            tekdev = tekdev_plus;
        } else {
            tekdev = tekdev_minus;
        }
        if (tekdev < mindev) {
            mindev = tekdev;
        }
        if (tekdev > maxdev) {
            maxdev = tekdev;
        }
        (hist->Resolution_L)->Fill(tekdev - StrawPoint_L_after_cut.at(i));
    }

    (hist->Resolution_L)->Fit("gaus", "QR", "SAME", -0.5, 0.5);
}
void linking_data::StrawResolution_S(histos *hist) {
    Double_t mindev = 99999, maxdev = -9999999;

    Double_t b_S_check = -19.88;
    Double_t c_S_check = 9.149;
    Double_t a_S_check = 6.771;

    for (int i = 0; i < timing_S_after_cut.size(); ++i) {
//        Discrim2_S = (b_S * b_S) - 4.0 * a_S * (c_S - timing_S_after_cut.at(i));
//        if (Discrim2_S < 0) continue;
//        Double_t discr = sqrt(Discrim2_S);
//        Double_t tekdev_plus = (-b_S + discr) / (2.0 * a_S);
//        Double_t tekdev_minus = (-b_S - discr) / (2.0 * a_S);

        Discrim2_S = (b_S_check * b_S_check) - 4.0 * a_S_check * (c_S_check - timing_S_after_cut.at(i));
        if (Discrim2_S < 0) continue;
        Double_t discr = sqrt(Discrim2_S);
        Double_t tekdev_plus = (-b_S_check + discr) / (2.0 * a_S_check);
        Double_t tekdev_minus = (-b_S_check - discr) / (2.0 * a_S_check);

        Double_t tekdev;

        if (fabs(tekdev_minus - StrawPoint_S_after_cut.at(i)) > fabs(tekdev_plus - StrawPoint_S_after_cut.at(i))) {
            tekdev = tekdev_plus;
        } else {
            tekdev = tekdev_minus;
        }
        if (tekdev < mindev) {
            mindev = tekdev;
        }
        if (tekdev > maxdev) {
            maxdev = tekdev;
        }
        (hist->Resolution_S)->Fill(tekdev - StrawPoint_S_after_cut.at(i));
    }

    (hist->Resolution_S)->Fit("gaus", "QR", "SAME", -0.5, 0.5);
}
void linking_data::filling_hists(histos *hist, STRAW_presetting *straw, MAMBA_presetting *mamba) {
    if (mamba->ntrk == 1) {
        if ((straw->nt00 == 1) && (straw->nt08 == 1)) {
            (hist->vshapeUS)->Fill(mamba->spUPos[0] * 10,
                                   conversion * (straw->tt00 - straw->tt08)); //conversion = 25ms / 1024bits
            (hist->vshapeVS)->Fill(mamba->spVPos[0] * 10, conversion * (straw->tt00 - straw->tt08));
        }

        if ((straw->nt01 == 1) && (straw->nt08 == 1)) {
            (hist->vshapeUL)->Fill(mamba->spUPos[0] * 10, conversion * (straw->tt01 - straw->tt08));
            (hist->vshapeVL)->Fill(mamba->spVPos[0] * 10, conversion * (straw->tt01 - straw->tt08));
        }
        if (straw->nt00 == 1) (hist->n_hits_short_straw)->Fill(straw->nt08);
        if (straw->nt01 == 1) (hist->n_hits_long_straw)->Fill(straw->nt08);
    }
}

void linking_data::limits_S(histos *hist, Double_t &limits) {
    Double_t values[3] = {};
    (&limits)[0] = (hist->vshapeUS_proj)->GetXaxis()->GetBinCenter((hist->vshapeUS_proj)->GetMaximumBin());
    values[0] = (hist->vshapeUS_proj)->GetBinContent((hist->vshapeUS_proj)->GetMaximumBin());
    (hist->vshapeUS_proj)->SetAxisRange((&limits)[0] - 25, (&limits)[0] - 19);
    (&limits)[1] = (hist->vshapeUS_proj)->GetXaxis()->GetBinCenter((hist->vshapeUS_proj)->GetMaximumBin());
    values[1] = (hist->vshapeUS_proj)->GetBinContent((hist->vshapeUS_proj)->GetMaximumBin());
    (hist->vshapeUS_proj)->SetAxisRange((&limits)[0] + 19, (&limits)[0] + 25);
    (&limits)[2] = (hist->vshapeUS_proj)->GetXaxis()->GetBinCenter((hist->vshapeUS_proj)->GetMaximumBin());
    values[2] = (hist->vshapeUS_proj)->GetBinContent((hist->vshapeUS_proj)->GetMaximumBin());

    if (values[1] < values[2]) {
        (&limits)[1] = (&limits)[2];
    }
    if ((&limits)[1] < (&limits)[0]) {
        Double_t temp = 0;
        temp = (&limits)[0];
        (&limits)[0] = (&limits)[1];
        (&limits)[1] = temp;
    }
}

void linking_data::limits_L(histos *hist, Double_t &limits) {
    Double_t values[3] = {};
    (&limits)[0] = (hist->vshapeUL_proj)->GetXaxis()->GetBinCenter((hist->vshapeUL_proj)->GetMaximumBin());
    values[0] = (hist->vshapeUL_proj)->GetBinContent((hist->vshapeUL_proj)->GetMaximumBin());
    (hist->vshapeUL_proj)->SetAxisRange((&limits)[0] - 15, (&limits)[0] - 9);
    (&limits)[1] = (hist->vshapeUL_proj)->GetXaxis()->GetBinCenter((hist->vshapeUL_proj)->GetMaximumBin());
    values[1] = (hist->vshapeUL_proj)->GetBinContent((hist->vshapeUL_proj)->GetMaximumBin());
    (hist->vshapeUL_proj)->SetAxisRange((&limits)[0] + 9, (&limits)[0] + 15);
    (&limits)[2] = (hist->vshapeUL_proj)->GetXaxis()->GetBinCenter((hist->vshapeUL_proj)->GetMaximumBin());
    values[2] = (hist->vshapeUL_proj)->GetBinContent((hist->vshapeUL_proj)->GetMaximumBin());

    if (values[1] < values[2]) {
        (&limits)[1] = (&limits)[2];
    }
    if ((&limits)[1] < (&limits)[0]) {
        Double_t temp = 0;
        temp = (&limits)[0];
        (&limits)[0] = (&limits)[1];
        (&limits)[1] = temp;
    }
}

void linking_data::draft(histos *hist, profiling *profile) {
    profile->main_algorithm((hist->vshapeUS_clear)->FindFirstBinAbove(), (hist->vshapeUS_clear)->FindLastBinAbove(), hist->vshapeUS_clear, hist);
}