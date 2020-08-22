double abcModel(double* x, double* par) {
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double E = x[0];
    return TMath::Sqrt(a/E+b+c/E/E);
}

void compare_centeruni_gamma()
{
    double scale = 1363;

    const int N = 4;
    double center_pe[N] = {8.98451e+02, 1.42441e+03, 2.08021e+03, 9.19555e+03};
    double center_peerr[N] = {2.25358e-01, 1.44615e+00, 3.77970e-01, 9.13540e-01};
    double center_sigma[N] = {3.22243e+01, 4.08414e+01, 5.45291e+01, 1.28307e+02};
    double center_sigmaerr[N] = {1.57736e-01, 1.31273e+00, 2.64676e-01, 6.09010e-01};

    std::string par[N] = {"Cs137", "Ge68", "K40", "O16"};
    TH1D* pe_hist[20]; // maximum 20 energies
    Float_t m_evis; TBranch* b_evis;
    Float_t m_edepZ; TBranch* b_edepZ;
    Float_t m_edepR; TBranch* b_edepR;
    for(int i=0; i<N; i++) {
        string hist_name = "evis"+(par[i])+"MeV";
        pe_hist[i] = new TH1D(hist_name.c_str(), "", 2400, 0, 12);
        pe_hist[i]->SetName(hist_name.c_str());
        string filename = "/junofs/users/yumiao/simulation/energy_model/uniformity/gamma/"+par[i]+".root";
        TFile* file = new TFile(filename.c_str(), "read");
        if(!file) continue;
        cout << "Processing " << filename << endl;
        TTree* tree = (TTree*)file->Get("petree");
        tree->SetBranchAddress("evis", &m_evis, &b_evis);
        tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
        tree->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
        for(int j=0; j<tree->GetEntries(); j++) {
            tree->GetEntry(j);
            if(m_edepR > 17200) continue;
            double delta_R = 18*18*18./20.;
            int ir = int( TMath::Power(m_edepR/1000, 3) / delta_R);
            double cos_theta = m_edepZ/m_edepR;
            double delta_cos_theta = 2./10.;
            int itheta = int((cos_theta+1)/delta_cos_theta);
            int ibin = ir*10+itheta;
            //m_evis /= corr[ibin];  // secondary correction
            pe_hist[i]->Fill(m_evis);
            //if(m_edepR<15000) {hTotPE_FV->Fill(m_evis);}
        }
    }

    TGraphErrors* ge_uni = new TGraphErrors();
    TGraphErrors* ge_center = new TGraphErrors();
    TGraphErrors* ge_diff = new TGraphErrors();
    for(int i=0; i<N; i++) {
        pe_hist[i]->Fit("gaus", "Q0");
        TF1* f1 = (TF1*)pe_hist[i]->GetFunction("gaus");
        double mu1 = f1->GetParameter(1);
        double mu1_err = f1->GetParError(1);
        double sigma1 = f1->GetParameter(2);
        double sigma1_err = f1->GetParError(2);
        double rsl1 = sigma1/mu1;
        double rsl1_err = TMath::Sqrt(sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);

        double mu2 = center_pe[i]/scale;
        double mu2_err = center_peerr[i]/scale;
        double sigma2 = center_sigma[i]/scale;
        double sigma2_err = center_sigmaerr[i]/scale;
        double rsl2 = sigma2/mu2;
        double rsl2_err = TMath::Sqrt(sigma2_err*sigma2_err/mu2/mu2 + mu2_err*mu2_err*sigma2*sigma2/mu2/mu2/mu2/mu2);

        ge_uni->SetPoint(i, mu1, rsl1);
        ge_uni->SetPointError(i, mu1_err, rsl1_err);
        ge_center->SetPoint(i, mu2, rsl2);
        ge_center->SetPointError(i, mu2_err, rsl2_err);

        cout << i << " " << center_pe[i] << " " << rsl1 << " " << rsl2 << endl;

        //ge_diff->SetPoint(i, mu1, rsl1*rsl1-rsl2*rsl2);
        //ge_diff->SetPointError(i, 0, TMath::Sqrt( TMath::Power(2*rsl1*rsl1_err,2) + TMath::Power(2*rsl2*rsl2_err, 2) ));

        ge_diff->SetPoint(i, mu1, rsl1-rsl2);
        ge_diff->SetPointError(i, 0, TMath::Sqrt(rsl1_err*rsl1_err + rsl2_err*rsl2_err));
    }


    //TF1* fAbcModel_uni = new TF1("fAbcModel_uni", abcModel, 0, 12000, 3);
    //fAbcModel_uni->SetParameters(8.3e-4, 3.5e-5, -8.6e-5);
    //ge_uni->Fit(fAbcModel_uni);
    //TF1* fAbcModel_center = new TF1("fAbcModel_center", abcModel, 0, 12000, 3);
    //fAbcModel_center->SetParameters(8.3e-4, 3.5e-5, -8.6e-5);
    //ge_center->Fit(fAbcModel_center);
    //for(int i=0; i<700; i++) {
    //    double npe2 = 1.0+0.01*i;
    //    //gExtra->SetPoint(i, npe2, fAbcModel_uni->Eval(npe2)- fAbcModel_center->Eval(npe2));
    //    double fAbcModel_center_err = TMath::Sqrt(  TMath::Power(fAbcModel_center->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_center->GetParError(1),2) +  TMath::Power(fAbcModel_center->GetParError(2)/npe2/npe2, 2)  );
    //    double fAbcModel_uni_err = TMath::Sqrt(  TMath::Power(fAbcModel_uni->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_uni->GetParError(1),2) +  TMath::Power(fAbcModel_uni->GetParError(2)/npe2/npe2, 2)  );
    //    //gExtra->SetPoint(i, npe2, fAbcModel_uni->Eval(npe2)-fAbcModel_center->Eval(npe2));
    //    ge_diff->SetPoint(i, npe2, (fAbcModel_uni->Eval(npe2)*fAbcModel_uni->Eval(npe2) - fAbcModel_center->Eval(npe2)*fAbcModel_center->Eval(npe2)));
    //    ge_diff->SetPointError(i, 0, TMath::Sqrt(fAbcModel_center_err*fAbcModel_center_err + fAbcModel_uni_err*fAbcModel_uni_err));
    //    //double fAbcModel_center_err_abs = TMath::Sqrt(  TMath::Power(fAbcModel_center->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_center->GetParError(1),2) +  TMath::Power(fAbcModel_center->GetParError(2)/npe2/npe2, 2)  ) / TMath::Sqrt(fAbcModel_center->Eval(npe2));
    //    //double fAbcModel_uni_err_abs = TMath::Sqrt(  TMath::Power(fAbcModel_uni->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_uni->GetParError(1),2) +  TMath::Power(fAbcModel_uni->GetParError(2)/npe2/npe2, 2)  ) / TMath::Sqrt(fAbcModel_uni->Eval(npe2));
    //    //ge_diff->SetPoint(i, npe2, fAbcModel_uni->Eval(npe2) - fAbcModel_center->Eval(npe2) );
    //    //ge_diff->SetPointError(i, 0, TMath::Sqrt(fAbcModel_center_err_abs*fAbcModel_center_err_abs + fAbcModel_uni_err_abs*fAbcModel_uni_err_abs));
    //}


    TCanvas* cc1 = new TCanvas("resol", "", 800, 600);
    auto mg = new TMultiGraph(); 
    ge_uni->SetLineColor(38);
    ge_uni->SetLineWidth(2);
    ge_uni->SetMarkerColor(38);
    ge_uni->SetMarkerStyle(24);
    mg->Add(ge_uni);
    ge_center->SetLineColor(48);
    ge_center->SetLineWidth(2);
    ge_center->SetMarkerColor(48);
    ge_center->SetMarkerStyle(25);
    mg->Add(ge_center);
    mg->Draw("AP");
    TLegend* ll1 = new TLegend();
    ll1->AddEntry(ge_uni, "uniform");
    ll1->AddEntry(ge_center, "center");
    ll1->Draw("SAME");

    //TCanvas* cc2 = new TCanvas("diff", "", 800, 600);
    //ge_diff->SetLineColor(48);
    //ge_diff->SetLineWidth(2);
    ////ge_diff->SetMarkerColor(48);
    ////ge_diff->SetMarkerStyle(24);
    //ge_diff->SetFillStyle(3003);
    //ge_diff->SetFillColor(44);
    ////ge_diff->SetTitle("; Evis/MeV;  extra term");
    //ge_diff->SetTitle("; Evis/MeV; absolute resolution increase");
    //ge_diff->Draw("A3 l");

}

