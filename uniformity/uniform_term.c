/*************************************************************************
	> File Name: uniform_term.c
	> Author: MiaoYu 
	> Mail: miaoyu@ihep.ac.cn
	> Created Time: Mon Jul 27 15:01:49 2020
 ************************************************************************/

#include<stdio.h>


Double_t pesmear(Double_t totpe) {
    Double_t p0 = -369.5;
    Double_t p1 = 1.320;
    Double_t p2 = 1.901e-5;
    return p0+p1*totpe+p2*totpe*totpe;
   
}

void uniform_term()
{
    gStyle->SetOptFit(1111);
    /*
    TFile* f = TFile::Open("./out6MeV.root");
    TGraph* gtot = (TGraph*)f->Get("totpe");
    Double_t *pos = gtot->GetX();
    Double_t *evis = gtot->GetY();
    Double_t totpe[40]; Double_t totpeSigma[40];
    Double_t scale = 1473*6;


    // parameterisatrion:
    Double_t k0 = 0.931; Double_t k1=44.15; 

    for(int i=0; i<gtot->GetN(); i++)  {
        totpe[i] = evis[i] / evis[0] * scale;
        //Double_t sctpe = totpe[i] * (k0+k1/totpe[i]);
        //Double_t cerpe = totpe[i] - sctpe;
        totpeSigma[i] = TMath::Sqrt( pesmear(totpe[i]) );
        cout << i << " " << totpe[i] << " " << totpeSigma[i] << endl;
    }

    // MC sampling
    TH1D* hTotPE = new TH1D("totPE", "", 3000, 0, 12000);
    for(int i=0; i<50000; i++) {
        int idx = int(gRandom->Uniform(0,40));
        double sample_pe = gRandom->Gaus(totpe[idx], totpeSigma[idx]/totpe[idx]*scale);
        //double sample_pe = gRandom->Gaus(scale, totpeSigma[idx]/totpe[idx]*scale);
        hTotPE->Fill(sample_pe);
    }
    hTotPE->Fit("gaus");
    hTotPE->Draw();
    */

    const Int_t rbins = 40; Double_t r3max = 17.7*17.7*17.7;
    Double_t pos[rbins]; Double_t pos_err[rbins];
    for(int ibin=0; ibin<rbins; ibin++) {
        pos[ibin] = r3max/rbins*ibin;
        pos_err[ibin] = 0.;
    }

    Double_t mean[40]; Double_t sigma[40]; Double_t rsl[rbins];
    Double_t mean_err[rbins]; Double_t rsl_err[rbins];
    TFile* f = TFile::Open("analysis8MeV.root");
    for(int i=0; i<40; i++) {
        TString Name = Form("radius%d", i);
        TH1D* h1 = (TH1D*)f->Get(Name);
        //mean[i]  = h1->GetMean(); 
        //sigma[i] = h1->GetStdDev();
        h1->Fit("gaus", "Q0");
        TF1* func = (TF1*)h1->GetFunction("gaus");
        double totpe = func->GetParameter(1);
        double totpe_err = func->GetParError(1);
        double sigma1 = func->GetParameter(2);
        double sigma_err = func->GetParError(2);
        mean_err[i] = totpe_err;
        rsl_err[i] = TMath::Sqrt(sigma_err*sigma_err/totpe/totpe + totpe_err*totpe_err*sigma1*sigma1/totpe/totpe/totpe/totpe);

        mean[i] = totpe; sigma[i] = sigma1;
        //mean[i] = totpe*1.022; sigma[i] = sigma1*1.022;
        rsl[i] = sigma[i]/mean[i];
        cout << i << " " << mean[i] <<" " << sigma[i] << endl;
        TCanvas* cc = new TCanvas(); cc->SetName(Name);
        h1->GetXaxis()->SetRangeUser(0.0, 12.0);
        h1->Draw();
        cc->Print("fitting6MeVe-_17m_allpmt.pdf(");
        delete cc;
        delete func;
        delete h1;
    
    }

    // MC sampling
    TH1D* hTotPE = new TH1D("totPE", "", 2400, 0, 12);
    //TH1D* hTotPE = new TH1D("totPE", "", 3000, 0, 12000);
    for(int i=0; i<50000; i++) {
        int idx = int(gRandom->Uniform(0,29));
        double sample_pe = gRandom->Gaus(mean[idx], sigma[idx]);
        //double sample_pe = gRandom->Gaus(scale, totpeSigma[idx]/totpe[idx]*scale);
        hTotPE->Fill(sample_pe);
    }
    
    TCanvas* cMC = new TCanvas(); cMC->SetName("MC");
    hTotPE->Fit("gaus");
    hTotPE->SetTitle("MC pe");
    hTotPE->Draw();

    TCanvas* c1 = new TCanvas(); c1->SetName("canvas_pe");
    double rr = 0.01;
    TGraphErrors* gPE = new TGraphErrors(40, pos, mean, pos_err, mean_err); gPE->SetName("pegraph");
    TGraph* gZone = new TGraph();
    gZone->SetPoint(0, 0, gPE->GetMean(2)*(1-rr));
    gZone->SetPoint(1, 0, gPE->GetMean(2)*(1+rr));
    gZone->SetPoint(2, pos[39], gPE->GetMean(2)*(1+rr));
    gZone->SetPoint(3, pos[39], gPE->GetMean(2)*(1-rr));
    gPE->SetMarkerStyle(20);
    gPE->SetMarkerColor(kViolet+1);
    gPE->SetLineColor(kViolet+1);
    gPE->SetLineWidth(2);
    gPE->GetYaxis()->SetRangeUser(0.0, 12.0);
    gPE->SetTitle(";r^3/m^3; Evis/MeV");
    gPE->Draw("APL");
    gZone->SetFillColor(kCyan+1);
    gZone->SetFillStyle(3003);
    gZone->Draw("F SAME");
    c1->Print("fitting6MeVe-_17m_allpmt.pdf(");


    TCanvas* c2 = new TCanvas(); c2->SetName("canvas_rsl");
    TGraphErrors* gSigma = new TGraphErrors(40, pos, rsl, pos_err, rsl_err); gSigma->SetName("rslgraph");
    TGraph* gZone1 = new TGraph();
    double ratio = 0.035;
    gZone1->SetPoint(0, 0, gSigma->GetMean(2)*(1-ratio));
    gZone1->SetPoint(1, 0, gSigma->GetMean(2)*(1+ratio));
    gZone1->SetPoint(2, pos[39], gSigma->GetMean(2)*(1+ratio));
    gZone1->SetPoint(3, pos[39], gSigma->GetMean(2)*(1-ratio));
    gSigma->SetMarkerStyle(20);
    gSigma->SetMarkerColor(kViolet+1);
    gSigma->SetLineColor(kViolet+1);
    gSigma->SetLineWidth(2);
    //gSigma->GetYaxis()->SetRangeUser(6.0, 7.0);
    gSigma->SetTitle(";r^3/m^3; Evis/MeV");
    gSigma->Draw("APL");
    gZone1->SetFillColor(kCyan+1);
    gZone1->SetFillStyle(3003);
    //gZone1->Draw("F SAME");
    c2->Print("fitting6MeVe-_17m_allpmt.pdf(");

    TCanvas* cc = new TCanvas();
    cc->Print("fitting6MeVe-_17m_allpmt.pdf)");

    TFile* fout = new TFile("fitting6MeVe_17m_allpmt.root", "recreate");
    gPE->Write();
    gSigma->Write();
    hTotPE->Write();
    hTotPE->Write();
    fout->Close();
}
