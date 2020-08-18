extern double rsd_ratio[500]={0.}; 
TGraphErrors* gResol_data;
const double energy_scale = 3350/2.22;

void draw_NU(double *f);

double NPE_smear(double* x, double* par) {
    double val = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];

     sigma2 fitting 
    if(val <= 0)
        return 0;
    else {
        double fval = TMath::Sqrt(val);    
        return fval;
    }

    //if(val < 0) return 0;
    //return val;
}

double abcModel(double* x, double* par) {
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double E = x[0];
    return TMath::Sqrt(a/E+b+c/E/E);
}


void read_corr()
{
    Float_t m_evis; TBranch* b_evis;
    Float_t m_edepZ; TBranch* b_edepZ;
    Float_t m_edepR; TBranch* b_edepR;

    TString name = "../../uniformity/electron/pe1MeV.root";

    TFile* ff = TFile::Open(name, "read");
    TTree* tree = (TTree*)ff->Get("petree");
    tree->SetBranchAddress("evis", &m_evis, &b_evis);
    tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
    tree->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
    TProfile2D* prof = new TProfile2D("prof", "", 20, 0, 18*18*18, 10, -1, 1, 0., 12);
    //TH1D* hTotPE_FV = new TH1D("hTotPE_FV", "", 2400, 0, 12);
    double delta_R = 18*18*18./20.;
    for(int i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        prof->Fill(m_edepR/1000*m_edepR/1000*m_edepR/1000, m_edepZ/m_edepR, m_evis);
        //if(m_edepR<15000) {hTotPE_FV->Fill(m_evis);}
    }

    int idx = 0;
    for(int i=0; i<prof->GetNbinsX(); i++) {
        for(int j=0; j<prof->GetNbinsY(); j++) {
            rsd_ratio[idx] = prof->GetBinContent(i+1, j+1);
            idx++;
        }
    }


}


void toyMC_fullRange()
{
    TF1* fNPEsmear = new TF1("fNPEsmear", NPE_smear, 0, 12000, 3);
    fNPEsmear->SetParameters(-3.69505e2, 1.30168, 1.90109e-05);   // sigma2 fitting relation
    //fNPEsmear->SetParameters(2.06408e+01, 1.40885e-02, -4.04637e-07);  // sigma fitting relation

    //draw NPE smear :
    //TGraph* gNPEsmear = new TGraph();
    //for(int i=1; i<1000; i++) {
    //    int j = i*10;
    //    gNPEsmear->SetPoint(i-1, j, fNPEsmear->Eval(j));
    //}

    //gNPEsmear->SetLineWidth(2);
    //gNPEsmear->SetLineColor(42);
    //gNPEsmear->SetTitle("NPE smear; NPE; NPE sigma");
    //gNPEsmear->Draw("AL");


    TString name = "../../uniformity/electron/totalpe2r2theta_40rbins20thetabins_1MeV.root";
    TFile* ft = new TFile(name, "read");
    TProfile2D* prof2d = (TProfile2D*)ft->Get("pe2r2theta");
    const int rbins = prof2d->GetNbinsX();
    const int thetabins = prof2d->GetNbinsY();
    const int totalbins = rbins*thetabins;
    cout << totalbins << endl;
    double pe_array[totalbins];
    int idx = 0;
    for(int irbin=0; irbin<rbins; irbin++) {
        for(int ithetabin=0; ithetabin<thetabins; ithetabin++) {
            pe_array[idx] = prof2d->GetBinContent(irbin+1, ithetabin+1);
            //cout << irbin << " " << ithetabin << " " << pe_array[idx] <<endl;
            idx++;
        }
    }
    ft->Close();


    read_corr();

    auto mg1 = new TMultiGraph();
    auto led1 = new TLegend();

    const int N = 5;
    Color_t color1[N] = {36, 38, 40, 42, 44};
    Color_t color2[1] = {48};
    TString label1[5] = {"scale factor:1", "scale factor: 1.1", "scale factor: 1.2", "scale factor: 1.3", "scale factor: 1.4"};
    double factor[N] = {1., 1.10, 1.20, 1.30, 1.4};
    //draw_NU(factor);
    TGraphErrors* gResol_real[N];
    TGraphErrors* gResol_ideal;
    TGraphErrors* gResol_center;
    TH1D* hTotPE_real = new TH1D("hTotPE_real", "", 2400, 0, 12000);
    TH1D* hTotPE_ideal = new TH1D("hTotPE_ideal", "", 2400, 0, 12000);
    TH1D* hTotPE_center = new TH1D("hTotPE_center", "", 2400, 0, 12000);
    for(int iloop=0; iloop<N; iloop++) {
        cout <<  "LOOPI ================> " << iloop << endl;
        gResol_real[iloop] = new TGraphErrors();
        if( iloop==0 ) { gResol_ideal = new TGraphErrors(); gResol_center = new TGraphErrors(); }
        for(int jloop=0; jloop<20; jloop++) {
            hTotPE_real->Reset();
            if(iloop==0) { hTotPE_ideal->Reset(); hTotPE_center->Reset(); }
            double npe = (jloop+2)*500;
            //cout << " LOOPJ ===> " <<  jloop << "   NPE ===> " << npe << endl;
            idx = 0;
            while(idx<100000) {  // sampling uniformly
                int ibin = int(gRandom->Uniform(0,totalbins));
                if(ibin <=100 and ibin>=170 ) continue;  // FV cut
                int npe1 = npe*((rsd_ratio[ibin]-rsd_ratio[0])/rsd_ratio[0]*factor[iloop]+1);  //npe in current bin
                int totpe_real = int(gRandom->Gaus(npe1, fNPEsmear->Eval(npe1)));
                hTotPE_real->Fill(totpe_real);
                if( iloop==0 ) {
                    int totpe_ideal = int(gRandom->Gaus(npe, fNPEsmear->Eval(npe/pe_array[ibin]*pe_array[0])));
                    hTotPE_ideal->Fill(totpe_ideal);

                    int totpe_center = int(gRandom->Gaus(npe, fNPEsmear->Eval(npe)));
                    hTotPE_center->Fill(totpe_center);
                }
                //cout << ibin << " " << npe1 << " " << fNPEsmear->Eval(npe1) << " " << npe << " " << fNPEsmear->Eval(npe/pe_array[ibin]*pe_array[0]) << endl; 
                idx++;
            }

            hTotPE_real->Fit("gaus", "Q0");
            TF1* f1 = (TF1*)hTotPE_real->GetFunction("gaus");
            double mu1 = f1->GetParameter(1);
            double mu1_err = f1->GetParError(1);
            double sigma1 = f1->GetParameter(2);
            double sigma1_err = f1->GetParError(2);
            double rsl1 = sigma1/mu1;
            double rsl1_err = TMath::Sqrt(sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);
            if(iloop == 0) {
                hTotPE_ideal->Fit("gaus", "Q0");
                TF1* f2 = (TF1*)hTotPE_ideal->GetFunction("gaus");
                double mu2 = f2->GetParameter(1);
                double mu2_err = f2->GetParError(1);
                double sigma2 = f2->GetParameter(2);
                double sigma2_err = f2->GetParError(2);
                double rsl2 = sigma2/mu2;
                double rsl2_err = TMath::Sqrt(sigma2_err*sigma2_err/mu2/mu2 + mu2_err*mu2_err*sigma2*sigma2/mu2/mu2/mu2/mu2);
                gResol_ideal->SetPoint(jloop, mu2, rsl2);
                gResol_ideal->SetPointError(jloop, 0, rsl2_err);

                hTotPE_center->Fit("gaus", "Q0");
                TF1* f3 = (TF1*)hTotPE_center->GetFunction("gaus");
                double mu3 = f3->GetParameter(1);
                double mu3_err = f3->GetParError(1);
                double sigma3 = f3->GetParameter(2);
                double sigma3_err = f3->GetParError(2);
                double rsl3 = sigma3/mu3;
                double rsl3_err = TMath::Sqrt(sigma3_err*sigma3_err/mu3/mu3 + mu3_err*mu3_err*sigma3*sigma3/mu3/mu3/mu3/mu3);
                gResol_center->SetPoint(jloop, mu3, rsl3);
                gResol_center->SetPointError(jloop, 0, rsl3_err);
            }

            //cout << jloop << hTotPE_ideal->GetMean() << " " << hTotPE_ideal->GetStdDev() << endl;

            gResol_real[iloop]->SetPoint(jloop, mu1, rsl1);
            gResol_real[iloop]->SetPointError(jloop, 0, rsl1_err);

        }
        if(iloop == 0) {
            gResol_ideal->SetMarkerColor(color2[iloop]);
            gResol_ideal->SetMarkerStyle(25);
            gResol_ideal->SetLineColor(color2[iloop]);
            gResol_ideal->SetLineWidth(2);
            mg1->Add(gResol_ideal);
            led1->AddEntry(gResol_ideal,"ideal", "PL");

o           gResol_center->SetMarkerColor(25);

            gResol_center->SetMarkerStyle(26);
            gResol_center->SetLineColor(25);
            gResol_center->SetLineWidth(2);
            mg1->Add(gResol_center);
            led1->AddEntry(gResol_center,"center", "PL");
        }
        gResol_real[iloop]->SetMarkerColor(color1[iloop]);
        gResol_real[iloop]->SetMarkerStyle(24);
        gResol_real[iloop]->SetLineColor(color1[iloop]);
        gResol_real[iloop]->SetLineWidth(2);
        mg1->Add(gResol_real[iloop]);
        led1->AddEntry(gResol_real[iloop], label1[iloop], "LP");
    }

    delete hTotPE_real;
    delete hTotPE_ideal;

    TCanvas* cc = new TCanvas(); 
    mg1->SetTitle("ToyMC Resolution; NPE; resolution");
    mg1->Draw("APL");
    led1->Draw("SAME");
    
    // fitting with abc model:
    //TF1* fAbcModel1 = new TF1("fAbcModel1", abcModel, 0, 12000, 3);
    //fAbcModel1->SetParameters(0.98*0.98, 6.62*6.62*1e-6, 0);
    //gResol_ideal->Fit(fAbcModel1);
    //double *par1 = fAbcModel1->GetParameters();
    //TF1* fAbcModel2 = new TF1("fAbcModel2", abcModel, 0, 12000, 3);
    //fAbcModel2->SetParameters(0.98*0.98, 6.62*6.62*1e-6, 0);
    //gResol_real[0]->Fit(fAbcModel2);
    //double *par2 = fAbcModel2->GetParameters();
    //TF1* fAbcModel3 = new TF1("fAbcModel3", abcModel, 0, 12000, 3);
    //fAbcModel3->SetParameters(0.98*0.98, 6.62*6.62*1e-6, 0);
    //gResol_center->Fit(fAbcModel3);
    //TGraph* gExtra = new TGraph();
    //for(int i=0; i<1000; i++) {
    //    double npe2 = 12*(i+100);
    //    gExtra->SetPoint(i, npe2, fAbcModel2->Eval(npe2)- fAbcModel1->Eval(npe2));
    //    //gExtra->SetPoint(i, npe2, TMath::Sqrt(fAbcModel2->Eval(npe2)*fAbcModel2->Eval(npe2) - fAbcModel1->Eval(npe2)*fAbcModel1->Eval(npe2)));
    //}
    //TCanvas* c1 = new TCanvas(); c1->cd();
    //gExtra->SetLineWidth(2);
    //gExtra->SetLineColor(45);
    //gExtra->SetTitle("extra resolution term; NPE; extra term");
    //gExtra->Draw("AL");
    //mg1->Draw("APL");
    //led1->Draw("SAME");
}



void draw_NU(double *f) {
    auto mg = new TMultiGraph();
    auto led = new TLegend();
    TGraph* gNU[5];
    Color_t color1[5] = {36, 38, 40, 42, 44};
    TString label2[5] = {"scale factor:1", "scale factor: 1.1", "scale factor: 1.2", "scale factor: 1.3", "scale factor: 1.4"};
    for(int i=0; i<5; i++) {
        gNU[i] = new TGraph();
        for(int j=0; j<200; j++) {
            gNU[i]->SetPoint(j, j, (rsd_ratio[j]-rsd_ratio[0])/rsd_ratio[0]*f[i]+1);
        }
        gNU[i]->SetMarkerColor(color1[i]);
        gNU[i]->SetMarkerStyle(24);
        gNU[i]->SetLineColor(color1[i]);
        gNU[i]->SetLineWidth(2);
        mg->Add(gNU[i]);
        led->AddEntry(gNU[i], label2[i], "LP");
    }

    TCanvas* cc = new TCanvas(); cc->cd();
    mg->SetTitle("Residual Non-Uniformity; sub-detecotr id; ratio");
    mg->Draw("APL");
    led->Draw("SAME");
    cc->SaveAs("non-uniform-scale.pdf");
}
