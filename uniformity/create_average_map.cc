void create_average_map() {

    gStyle->SetOptStat(0);

    const int N = 5;
    Int_t mom_array[N] = {1, 2, 4, 6, 8};
    double rsd_ratio[N][200];

    for(int index=0; index<N; index++) {

        Float_t m_evis; TBranch* b_evis;
        Float_t m_edepZ; TBranch* b_edepZ;
        Float_t m_edepR; TBranch* b_edepR;

        string path = "../../uniformity/electron/pe"; string suffix = "MeV.root";
        string name = path+to_string(mom_array[index])+suffix;

        TFile* ff = TFile::Open(name.c_str());
        TTree* tree = (TTree*)ff->Get("petree");
        tree->SetBranchAddress("evis", &m_evis, &b_evis);
        tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
        tree->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
        TProfile2D* prof = new TProfile2D("prof", "", 20, 0, 18*18*18, 10, -1, 1, 0., 12);
        TH1D* hTotPE_FV = new TH1D("hTotPE_FV", "", 2400, 0, 12);
        double delta_R = 18*18*18./20.;
        for(int i=0; i<tree->GetEntries(); i++) {
            tree->GetEntry(i);
            prof->Fill(m_edepR/1000*m_edepR/1000*m_edepR/1000, m_edepZ/m_edepR, m_evis);
            if(m_edepR<15000) {hTotPE_FV->Fill(m_evis);}
        }

        int idx = 0;
        for(int i=0; i<prof->GetNbinsX(); i++) {
            for(int j=0; j<prof->GetNbinsY(); j++) {
                rsd_ratio[index][idx] = prof->GetBinContent(i+1, j+1);
                idx++;
            }
        }

    }

    double factor[N] = {rsd_ratio[0][0], rsd_ratio[1][0], rsd_ratio[2][0], rsd_ratio[3][0], rsd_ratio[4][0]};
    double average[200] = {0.};   
    for(int idx=0; idx<200; idx++) {
        for(int j=0; j<N; j++) {
            average[idx] += rsd_ratio[j][idx]/factor[j];
        }
        average[idx] /= N;
    }

    TGraph* gg = new TGraph();
    for(int i=0; i<200; i++)  {
        gg->SetPoint(i, i, average[i]);
    }

    TCanvas* c1= new TCanvas();
    gg->SetMarkerColor(28);
    gg->SetMarkerStyle(24);
    gg->SetMarkerSize(0.5);
    gg->SetLineColor(28);
    gg->SetLineWidth(2);
    gg->Draw("APL");


}
