extern double rsd_ratio[500]={0.}; 
double residual(int ibin);
void read_corr();

void toyMC_study()
{
    //TFile* f1 = TFile::Open("../../uniformity/electron/analysis8MeV_test.root");
    //TGraphErrors* gTotPE = (TGraphErrors*)f1->Get("pe_graph");
    //TGraphErrors* gResol = (TGraphErrors*)f1->Get("rsl_graph");
    //double* pos = gTotPE->GetX();
    //double* pe  = gTotPE->GetY();
    //double* peerr = gTotPE->GetEY();
    //double* rsl = gResol->GetY();
    //double* rslerr = gResol->GetEY();
    //TGraphErrors* gTotPE_15m = new TGraphErrors();
    //TGraphErrors* gResol_15m = new TGraphErrors();
    //for(int i=0; i<gTotPE->GetN(); i++) {
    //    if(pos[i]<15*15*15) {
    //        gTotPE_15m->SetPoint(i, pos[i], pe[i]);
    //        gTotPE_15m->SetPointError(i, 0., peerr[i]);
    //        gResol_15m->SetPoint(i, pos[i], rsl[i]);
    //        gResol_15m->SetPointError(i, 0, rslerr[i]);
    //    }
    //}

    //TF1* func = new TF1("func", "[0]+[1]*x", 0, 15*15*15);
    //gResol_15m->Fit(func, "RE");
    //gResol_15m->SetMarkerStyle(25);
    //gResol_15m->SetMarkerColor(46);
    //gResol_15m->SetLineColor(46);
    //gResol_15m->SetLineWidth(2);
    //gResol_15m->Draw("APL");

    // totpe info w/o correction : 

    gStyle->SetOptStat(0);
    gRandom->SetSeed(43);

    TFile* f1 = TFile::Open("../../uniformity/electron/totalpe2r2theta_40rbins20thetabins_1MeV.root");
    TProfile2D* prof2d = (TProfile2D*)f1->Get("pe2r2theta");
    const int rbins = prof2d->GetNbinsX();
    const int thetabins = prof2d->GetNbinsY();
    const int totalbins = rbins*thetabins;
    double pe_array[totalbins];
    int idx = 0;
    for(int irbin=0; irbin<rbins; irbin++) {
        for(int ithetabin=0; ithetabin<thetabins; ithetabin++) {
            pe_array[idx] = prof2d->GetBinContent(irbin+1, ithetabin+1);
            //cout << irbin << " " << ithetabin << " " << pe_array[idx] <<endl;
            idx++;
        }
    }

    // MC sampling
    // read totpe mean value from histogram, and use abc model fitting result to estimate sigma at each positions
    double p0 = -3.69505e+02;
    double p1 = 1.30168e+00;
    double p2 = 1.90109e-05;
    double sigma_array[totalbins];
    for(int idx=0; idx<totalbins; idx++) {
        sigma_array[idx] = TMath::Sqrt(pe_array[idx]*pe_array[idx]*p2 + pe_array[idx]*p1 + p0);
    }

    read_corr();

    TH1D* hTotPE = new TH1D("hTotPE", "", 200, 1000,2000);
    TH1D* hTotPE_center = new TH1D("hTotPE_center" ,"" ,200, 1000, 2000);
    idx = 0;
    double scale = pe_array[0]; double scale_sigma = sigma_array[0];
    while(idx<50000) {
        int ibin = int(gRandom->Uniform(0,totalbins));
        if(ibin >=180 ) continue;
        int totpe = int(gRandom->Gaus(pe_array[ibin], sigma_array[ibin]));
        //int totpe = int(gRandom->Gaus(scale*rsd_ratio[ibin]/rsd_ratio[0], sigma_array[ibin]*scale/pe_array[ibin]*rsd_ratio[ibin]/rsd_ratio[0]));
        hTotPE->Fill(totpe);
        int totpe_center = int(gRandom->Gaus(scale, scale_sigma));
        hTotPE_center->Fill(totpe_center);
        idx++;
    }

    ////prof2d->Draw("COLZ");
    cout << " >>>>> Fitting MC sampling for whole CD" << endl;
    hTotPE->Fit("gaus", "0");
    cout << " >>>>> Fitting MC sampling for CD center" << endl;
    hTotPE_center->Fit("gaus","0");

    hTotPE->SetLineColor(42);
    hTotPE->SetLineWidth(2);
    hTotPE_center->SetLineColor(36);
    hTotPE_center->SetLineWidth(2);
    //hTotPE_center->Draw("");
    hTotPE->Draw("");

    TFile* f2 = TFile::Open("../../uniformity/electron/analysis1MeV_test.root");
    TH1D* hTotPEData = (TH1D*)f2->Get("hTotPEorigin_15m");
    double s1 = hTotPE->GetEntries();
    double s2 = hTotPEData->GetEntries();
    hTotPEData->Scale(s1/s2);
    hTotPEData->SetLineColor(37);
    hTotPEData->SetLineWidth(2);
    hTotPEData->Draw("SAME");

    TLegend* ll = new TLegend();
    ll->AddEntry(hTotPE, "toyMC", "L");
    ll->AddEntry(hTotPEData, "simul", "L");
    ll->Draw("SAME");
}

// residual non-uniformity effect
double residual(int ibin) {
    return 1 - 0.003/12*ibin;
}

void read_corr()
{
    Float_t m_evis; TBranch* b_evis;
    Float_t m_edepZ; TBranch* b_edepZ;
    Float_t m_edepR; TBranch* b_edepR;

    TFile* ff = TFile::Open("../../uniformity/electron/pe1MeV.root");
    TTree* tree = (TTree*)ff->Get("petree");
    tree->SetBranchAddress("evis", &m_evis, &b_evis);
    tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
    tree->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
    TProfile2D* prof = new TProfile2D("prof", "", 20, 0, 18*18*18, 10, -1, 1, 0.8, 1.4);
    double delta_R = 18*18*18./20.;
    for(int i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        prof->Fill(m_edepR/1000*m_edepR/1000*m_edepR/1000, m_edepZ/m_edepR, m_evis);
    }

    int idx = 0;
    for(int i=0; i<prof->GetNbinsX(); i++) {
        for(int j=0; j<prof->GetNbinsY(); j++) {
            rsd_ratio[idx] = prof->GetBinContent(i+1, j+1);
            idx++;
        }
    }
    //TCanvas* cc = new TCanvas();
    //prof->Draw("COLZ");
    //cc->SaveAs("corr.pdf");

}

