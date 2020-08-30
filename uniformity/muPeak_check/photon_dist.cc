void photon_dist()
{
    gStyle->SetOptStat(0);

    TFile* file = TFile::Open("./0_0_12m/normal1.root"); 
    TTree* evt = (TTree*)file->Get("evt");
    Int_t m_Npmt; TBranch* b_Npmt;
    Int_t m_pmtID[200000]; TBranch* b_pmtID;
    Int_t m_nPE[200000]; TBranch* b_nPE;
    evt->SetBranchAddress("nPMTs", &m_Npmt, &b_Npmt);
    evt->SetBranchAddress("PMTID_byPMT", m_pmtID, &b_pmtID);
    evt->SetBranchAddress("nPE_byPMT", m_nPE, &b_nPE);

    TH2D* hmap = new TH2D("pemap","", 720, -180, 180, 360, -90, 90);

    ifstream in; in.open("/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc830/Pre-Release/J20v1r0-Pre2/offline/Simulation/DetSimV2/DetSimOptions/data/PMTPos_Acrylic_with_chimney.csv");
    const int num = 17612; int index = 0;
    string line ; double pmt_id[num], pmt_x[num], pmt_y[num], pmt_z[num], pmt_theta[num], pmt_phi[num];
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> pmt_id[index] >> pmt_x[index] >> pmt_y[index] >> pmt_z[index] >> pmt_theta[index] >> pmt_phi[index]; 
        pmt_theta[index] = 90-pmt_theta[index];
        //cout << pmt_theta[index] << " " << pmt_phi[index] <<endl;
        index++;
    }

    int entries = evt->GetEntries();
    for(int i=0; i<30000; i++) {
    //for(int i=0; i<entries; i++) {
        //double tmp_pe[100] = {0.};
        evt->GetEntry(i); cout << "Processing " << i << "th Event" << endl;
        for(int ipmt=0; ipmt<m_Npmt; ipmt++) {
            int idx = m_pmtID[ipmt];
            if(idx>=17612) continue;
            for(int num=0; num<m_nPE[ipmt]; num++) {
                double plotx = (pmt_phi[idx]-180)*TMath::Cos(pmt_theta[idx]/180*TMath::Pi());
                double ploty = pmt_theta[idx];
                hmap->Fill(plotx, ploty);
            }
        } 
    }

    hmap->Scale(1./30000.);
    hmap->Draw("COLZ");
    //hmap->SaveAs("00_9m.pdf");
}
