void compare_resol()
{
    double scale = 3350/2.22;

    std::vector<string> filename;
    filename.push_back("1000keV/sample_detsim_user.root");
    filename.push_back("2000keV/sample_detsim_user.root");
    filename.push_back("3000keV/sample_detsim_user.root");
    filename.push_back("5000keV/sample_detsim_user.root");
    filename.push_back("8000keV/sample_detsim_user.root");

    string path = "/junofs/users/yumiao/simulation/energy_model/production/electron/cerenkov/";

    double nbins[5] = {200, 200, 300, 300, 400};
    double start[5] = {1200,2600,4000,7000,11400};
    double stop[5]  = {1800,3300,5000,8000,12800};
    TGraphErrors *gResol = new TGraphErrors();
    TGraphErrors *gResol_noCer = new TGraphErrors();


    Int_t m_totpe; TBranch* b_totpe;
    Float_t m_ratio; TBranch* b_ratio;

    for(int i=0; i<filename.size(); i++) {
        cout << "Processing " << filename[i] << endl;
        TFile* file = TFile::Open((path+filename[i]).c_str());
        TTree* evt = (TTree*)file->Get("evt");
        evt->SetBranchAddress("totalPE", &m_totpe, &b_totpe);
        evt->SetBranchAddress("CerenkovRatio", &m_ratio, &b_ratio);

        TH1D* h1 = new TH1D("h1", "", nbins[i], start[i], stop[i]);
        TH1D* h2 = new TH1D("h2", "", nbins[i], start[i], stop[i]);

        for(int j=0; j<evt->GetEntries(); j++){
            evt->GetEntry(j);
            h1->Fill(m_totpe);
            h2->Fill(m_totpe*(1-m_ratio));
        }

        // npe fitting: 
        h1->Fit("gaus", "0");
        TF1* f1 = (TF1*)h1->GetFunction("gaus");
        double mu1 = f1->GetParameter(1);
        double mu1_err = f1->GetParError(1);
        double sigma1 = f1->GetParameter(2);
        double sigma1_err = f1->GetParError(2);
        double rsl1 = sigma1/mu1;
        double rsl1_err = TMath::Sqrt(sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);

        h2->Fit("gaus", "0");
        TF1* f2 = (TF1*)h2->GetFunction("gaus");
        double mu2 = f2->GetParameter(1);
        double mu2_err = f2->GetParError(1);
        double sigma2 = f2->GetParameter(2);
        double sigma2_err = f2->GetParError(2);
        if(i==4) sigma2_err = 1.6;
        double rsl2 = sigma2/mu2;
        double rsl2_err = TMath::Sqrt(sigma2_err*sigma2_err/mu2/mu2 + mu2_err*mu2_err*sigma2*sigma2/mu2/mu2/mu2/mu2);

        gResol->SetPoint(i, mu1/scale, rsl1);
        gResol->SetPointError(i, mu1_err/scale, rsl1_err);
        gResol_noCer->SetPoint(i, mu2/scale, rsl2);
        gResol_noCer->SetPointError(i, mu2_err/scale, rsl2_err);

        delete h1; 
        delete h2; 
        delete evt;
        delete file;
    }


    TCanvas* cc = new TCanvas();  cc->SetGrid();
    gResol->SetMarkerStyle(24);
    gResol->SetMarkerColor(36);
    gResol->SetLineColor(36);
    gResol->SetLineWidth(2);
    gResol_noCer->SetMarkerStyle(25);
    gResol_noCer->SetMarkerColor(42);
    gResol_noCer->SetLineColor(42);
    gResol_noCer->SetLineWidth(2);
    gResol->SetTitle("Resolution; Evis/MeV; resolution");
    gResol->GetYaxis()->SetRangeUser(0,0.03);
    gResol->Draw("APL");
    gResol_noCer->Draw("PL SAME");

    TLegend* ll = new TLegend();
    ll->AddEntry(gResol, "w/ Cerenkov", "L");
    ll->AddEntry(gResol_noCer, "w/o Cerenkov", "L");
    ll->Draw("SAME");
}
