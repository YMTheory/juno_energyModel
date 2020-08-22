void create_average_map() {

    gStyle->SetOptStat(0);

    const int N = 5;
    Int_t mom_array[N] = {1, 2, 4, 6, 8};
    double rsd_ratio[N][200];

    for(int index=0; index<1; index++) {

        Float_t m_evis; TBranch* b_evis;
        Float_t m_edepX; TBranch* b_edepX;
        Float_t m_edepY; TBranch* b_edepY;
        Float_t m_edepZ; TBranch* b_edepZ;
        Float_t m_edepR; TBranch* b_edepR;

        string path = "../../uniformity/electron/pe"; string suffix = "MeV.root";
        string name = path+to_string(mom_array[index])+suffix;
        
        name = "../../uniformity/gamma/K40.root";

        TFile* ff = TFile::Open(name.c_str());
        TTree* tree = (TTree*)ff->Get("petree");
        tree->SetBranchAddress("evis", &m_evis, &b_evis);
        tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
        tree->SetBranchAddress("edepX", &m_edepX, &b_edepX);
        tree->SetBranchAddress("edepY", &m_edepY, &b_edepY);
        tree->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
        TProfile2D* prof = new TProfile2D("prof", "", 20, 0, 18*18*18, 10, -1, 1, 0., 12);
        TH1D* hTotPE[200];
        for(int i=0; i<200; i++) {
            TString name = Form("pe%d",i);
            hTotPE[i] = new TH1D(name, "", 2000, 0, 12 );
        }
        double delta_R = 18*18*18./20.;
        for(int i=0; i<tree->GetEntries(); i++) {
            tree->GetEntry(i);

            // vertex smearing:
            //double sigma_vertex = 230./TMath::Sqrt(m_evis);  // mm unit               
            //m_edepX = gRandom->Gaus(m_edepX, sigma_vertex);
            //m_edepY = gRandom->Gaus(m_edepY, sigma_vertex);
            //m_edepZ = gRandom->Gaus(m_edepZ, sigma_vertex);

            m_edepR = TMath::Sqrt(m_edepX*m_edepX+m_edepY*m_edepY+m_edepZ*m_edepZ);
            prof->Fill(m_edepR/1000*m_edepR/1000*m_edepR/1000, m_edepZ/m_edepR, m_evis);
            int rbin = TMath::Power(m_edepR/1000,3)/(18*18*18/20);
            int thetabin = (m_edepZ/m_edepR+1) / (2./10);
            int nbin = rbin*10+thetabin;
            hTotPE[nbin]->Fill(m_evis);
        }

        int idx = 0;
        for(int i=0; i<prof->GetNbinsX(); i++) {
            for(int j=0; j<prof->GetNbinsY(); j++) {
                rsd_ratio[index][idx] = prof->GetBinContent(i+1, j+1);
                idx++;
            }
        }


        TFile* out = new TFile("K40.root", "recreate");
        for(int i=0; i<200; i++) {
            TString n1 = Form("sub%d", i);
            //TCanvas* c1 = new TCanvas();
            //c1->SetName(n1);
            hTotPE[i]->SetName(n1);
            //hTotPE[i]->Draw();
            //c1->Print("smear_sub_100mm.pdf(");
            hTotPE[i]->Write();
        }
        out->Close();
        //TCanvas* cc = new TCanvas();
        //cc->Print("smear_sub_100mm.pdf)");
    }

    double factor[N] = {rsd_ratio[0][0], rsd_ratio[1][0], rsd_ratio[2][0], rsd_ratio[3][0], rsd_ratio[4][0]};
    double average[200] = {0.};   
    for(int idx=0; idx<200; idx++) {
        for(int j=0; j<1; j++) {
            average[idx] += rsd_ratio[j][idx]/factor[j];
        }
        average[idx] /= 1;
    }

    TGraph* gg = new TGraph(); gg->SetName("average");
    for(int i=0; i<200; i++)  {
        gg->SetPoint(i, i, average[i]);
    }

    gg->SetMarkerColor(28);
    gg->SetMarkerStyle(24);
    gg->SetMarkerSize(0.5);
    gg->SetLineColor(28);
    gg->SetLineWidth(2);
    TCanvas* c1= new TCanvas();
    gg->Draw("APL");

    //TFile* ff = new TFile("subDet-corr-smear500mm.root", "recreate");
    //gg->Write();
    //ff->Close();


}
