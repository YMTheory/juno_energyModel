extern double rsd_ratio[500]={0.}; 
double residual(int ibin);
void read_corr(int index);
void customed_corr();
TGraphErrors* gResol_data;
const double energy_scale = 3350/2.22;
const int mom_num = 5;
int mom_array[mom_num] = {1,2,4,6,8};

void toyMC_study()
{
    gStyle->SetOptStat(0);
    gRandom->SetSeed(52);

    //TFile* f3 = TFile::Open("../../uniformity/electron/analysis1MeV_test.root");
    //TH1D* hR0 = (TH1D*)f3->Get("radius0");
    //TH1D* hCorr15m = (TH1D*)f3->Get("hTotPE15m");
    //double r1 = hR0->GetEntries();
    //double r2 = hCorr15m->GetEntries();
    //hCorr15m->Scale(r1/r2);
    //hR0->SetLineColor(36);
    //hCorr15m->SetLineColor(42);
    //TCanvas* cg1 = new TCanvas();
    //hR0->Draw();
    //hCorr15m->Draw("SAME");
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

    int nbins[5] = {200, 400, 400, 1000, 1000};
    double start[5] = {1000, 2000, 5000, 8000, 10000};
    double stop[5] = {2000, 4000, 7000, 13000, 15000};
    TGraphErrors* gResol = new TGraphErrors();
    TGraphErrors* gResol_ideal = new TGraphErrors();
    TGraphErrors* gResol_diff = new TGraphErrors();
    gResol_data = new TGraphErrors();

    for(int i=0; i<mom_num; i++){
        string path = "../../uniformity/electron/totalpe2r2theta_40rbins20thetabins_"; string suffix = "MeV.root";
        string name = path+to_string(mom_array[i])+suffix;
        TFile* ft = TFile::Open(name.c_str()); if(!ft) {cout << "No such file " << name <<endl; break; }
        cout << "Processing file " << name << endl;
        TProfile2D* prof2d = (TProfile2D*)ft->Get("pe2r2theta");
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

        read_corr(i);
        //customed_corr();

        TH1D* hTotPE = new TH1D("hTotPE", "", nbins[i], start[i], stop[i]);
        TH1D* hTotPE_ideal = new TH1D("hTotPE_ideal", "", nbins[i], start[i], stop[i]);
        TH1D* hTotPE_center = new TH1D("hTotPE_center" ,"" ,nbins[i], start[i], stop[i]);
        idx = 0;
        double scale = pe_array[0]; double scale_sigma = sigma_array[0];
        while(idx<100000) {
            int ibin = int(gRandom->Uniform(0,totalbins));
            if(ibin >=120 ) continue;
            int totpe1 = int(gRandom->Gaus(scale, sigma_array[ibin]*scale/pe_array[ibin]));
            hTotPE_ideal->Fill(totpe1);
            //int totpe = int(gRandom->Gaus(pe_array[ibin], sigma_array[ibin]));
            int totpe = int(gRandom->Gaus(scale*rsd_ratio[ibin]/rsd_ratio[0], sigma_array[ibin]*scale/pe_array[ibin]*rsd_ratio[ibin]/rsd_ratio[0]));
            hTotPE->Fill(totpe);
            int totpe_center = int(gRandom->Gaus(scale, scale_sigma));
            hTotPE_center->Fill(totpe_center);
            idx++;
        }


        hTotPE->Fit("gaus", "0Q");
        TF1* f1 = (TF1*)hTotPE->GetFunction("gaus");
        double mu1 = f1->GetParameter(1);
        double mu1_err = f1->GetParError(1);
        double sigma1 = f1->GetParameter(2);
        double sigma1_err = f1->GetParError(2);
        double rsl1 = sigma1/mu1;
        double rsl1_err = TMath::Sqrt(sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);

        hTotPE_ideal->Fit("gaus", "Q0");
        TF1* f2 = (TF1*)hTotPE_ideal->GetFunction("gaus");
        double mu2 = f2->GetParameter(1);
        double mu2_err = f2->GetParError(1);
        double sigma2 = f2->GetParameter(2);
        double sigma2_err = f2->GetParError(2);
        if(i==4) sigma2_err = 1.6;
        double rsl2 = sigma2/mu2;
        double rsl2_err = TMath::Sqrt(sigma2_err*sigma2_err/mu2/mu2 + mu2_err*mu2_err*sigma2*sigma2/mu2/mu2/mu2/mu2);

        gResol->SetPoint(i, mu1/energy_scale, rsl1);
        gResol->SetPointError(i, mu1_err/energy_scale, rsl1_err);
        gResol_ideal->SetPoint(i, mu2/energy_scale, rsl2);
        gResol_ideal->SetPointError(i, mu2_err/energy_scale, rsl2_err);

        //gResol_diff->SetPoint(i, mu2/energy_scale, (rsl1-rsl2));
        

        //hTotPE_center->Fit("gaus","0Q");

    }
    //TCanvas* cg2 = new TCanvas(); cg2->cd();
    //hTotPE_center->GetXaxis()->SetTitle("NPE");
    //hTotPE->SetLineColor(42);
    //hTotPE->SetLineWidth(2);
    //hTotPE_center->SetLineColor(36);
    //hTotPE_center->SetLineWidth(2);
    //hTotPE_center->Draw("");
    //hTotPE->Draw("SAME");

    //TFile* f2 = TFile::Open("../../uniformity/electron/analysis1MeV_test.root");
    //TH1D* hTotPEData = (TH1D*)f2->Get("hTotPEorigin_15m");
    //double s1 = hTotPE->GetEntries();
    //double s2 = hTotPEData->GetEntries();
    //hTotPEData->Scale(s1/s2);
    //hTotPEData->SetLineColor(37);
    //hTotPEData->SetLineWidth(2);
    //hTotPEData->Draw("SAME");

    //TLegend* ll = new TLegend();
    //ll->AddEntry(hTotPE, "toyMC uniform", "L");
    //ll->AddEntry(hTotPE_center, "toyMC center", "L");
    //ll->Draw("SAME");

    TCanvas* cg3 = new TCanvas(); cg3->SetGrid();
    gResol->SetMarkerStyle(24);
    gResol->SetMarkerColor(36);
    gResol->SetLineColor(36);
    gResol->SetLineWidth(2);
    gResol_ideal->SetMarkerStyle(25);
    gResol_ideal->SetMarkerColor(42);
    gResol_ideal->SetLineColor(42);
    gResol_ideal->SetLineWidth(2);
    gResol->SetTitle("Resolution; Evis/MeV; resolution");
    gResol->GetYaxis()->SetRangeUser(0,0.03);
    gResol->Draw("APL");
    //gResol_ideal->Draw("PL SAME");


    /// scale evis as the same: 
    double *x1 = gResol->GetX();
    for(int ii=0; ii<mom_num; ii++) {
        gResol_data->SetPointX(ii, x1[ii]);
    }

    /////////
    gResol_data->SetMarkerStyle(20);
    gResol_data->SetMarkerColor(42);
    gResol_data->SetLineColor(42);
    gResol_data->SetLineWidth(2);
    gResol_data->Draw("PL SAME");

    TLegend* ll = new TLegend();
    ll->AddEntry(gResol, "simul residual non-uniformity toyMC", "L");
    //ll->AddEntry(gResol, "2\% linear residual non-uniformity", "L");
    ll->AddEntry(gResol_ideal, "full simulation", "L");
    ll->Draw("SAME");
    
    double *y1 = gResol->GetY();
    double *y2 = gResol_data->GetY();
    for(int i1=0; i1<mom_num; i1++) {
        gResol_diff->SetPoint(i1, x1[i1], (y2[i1]-y1[i1]));
    }

    TCanvas* cg4 = new TCanvas(); cg4->SetGrid();
    gResol_diff->SetMarkerStyle(25);
    gResol_diff->SetMarkerColor(42);
    gResol_diff->SetLineColor(42);
    gResol_diff->SetLineWidth(2);
    gResol_diff->SetTitle("resolution difference; Evis/MeV; absolute increase");
    gResol_diff->Draw("APL");
}

// residual non-uniformity effect
double residual(int ibin) {
    return 1 - 0.003/12*ibin;
}

void read_corr(int index)
{
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
            rsd_ratio[idx] = prof->GetBinContent(i+1, j+1);
            idx++;
        }
    }

    hTotPE_FV->Fit("gaus", "Q");
    TF1* f1 = (TF1*)hTotPE_FV->GetFunction("gaus");
    double mu1 = f1->GetParameter(1);
    double mu1_err = f1->GetParError(1);
    double sigma1 = f1->GetParameter(2);
    double sigma1_err = f1->GetParError(2);
    double rsl1 = sigma1/mu1;
    double rsl1_err = TMath::Sqrt(sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);

    cout << index << " " << mu1 << " " << rsl1 << endl;

    gResol_data->SetPoint(index, mu1, rsl1);
    gResol_data->SetPointError(index, mu1_err, rsl1_err);

    //TCanvas* cc = new TCanvas();
    //prof->Draw("COLZ");
    //cc->SaveAs("corr.pdf");
    
    //ff->Close();
    //delete hTotPE_FV;
    //delete tree;
    //delete ff;

}

void customed_corr()
{
    int idx = 0;
    for(int i=0; i<20; i++) {
        for(int j=0; j<10; j++) {
            rsd_ratio[idx] = 1 - 0.02/20.*i; idx++;
        }
    }
}

