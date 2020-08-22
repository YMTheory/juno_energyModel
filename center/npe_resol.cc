void npe_resol()
{
    std::vector<string> filename;
    string path = "/junofs/users/yumiao/simulation/energy_model/production/electron/cerenkov/";
    for(int i=0; i<7; i++) {
        string No = to_string(i*1000+1000);
        string name = path+No+"keV/sample_detsim_user.root";
        filename.push_back(name);
    }

    TGraphErrors* gResol = new TGraphErrors();

    Int_t m_totpe; TBranch* b_totpe;

    for(int i=0; i<filename.size(); i++) {
        cout << "Processing " << filename[i] << endl;
        TFile* file = TFile::Open((filename[i]).c_str());
        TTree* evt = (TTree*)file->Get("evt");
        evt->SetBranchAddress("totalPE", &m_totpe, &b_totpe);

        TH1D* h1 = new TH1D("h1", "", 2400, 0, 12000);

        for(int j=0; j<evt->GetEntries(); j++){
            evt->GetEntry(j);
            h1->Fill(m_totpe);
        }

        // npe fitting: 
        h1->Fit("gaus", "Q0");
        TF1* f1 = (TF1*)h1->GetFunction("gaus");
        double mu1 = f1->GetParameter(1);
        double mu1_err = f1->GetParError(1);
        double sigma1 = f1->GetParameter(2);
        double sigma1_err = f1->GetParError(2);
        double rsl1 = sigma1/mu1;
        double rsl1_err = TMath::Sqrt(sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);


        gResol->SetPoint(i, mu1, rsl1);
        gResol->SetPointError(i, mu1_err, rsl1_err);
        //gResol->SetPoint(i, mu1, sigma1);
        //gResol->SetPointError(i, mu1_err, sigma1_err);

        delete h1; 
        delete evt;
        delete file;
    }


    gResol->SetMarkerColor(32);
    gResol->SetMarkerStyle(26);
    gResol->SetLineColor(32);
    gResol->SetLineWidth(2);

    TF1 *func = new TF1("func", "TMath::Sqrt([0]/x+[1]+[2]/x/x)", 0, 12000);
    //TF1 *func = new TF1("func", "[0]+[1]*x+[2]*x*x", 0, 12000);
    gResol->Fit(func);
    gResol->Draw("APL");
}
