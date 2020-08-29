void residual_uniform_rsl_cd()
{
    TFile *ff = TFile::Open("results1MeV.root");
    double rmax = 17.7*17.7*17.7;
    TGraphErrors* gTotPE = new TGraphErrors();
    TGraphErrors* gRes = new TGraphErrors();

    double pemax = -1; double pemin= 100;
    int nPoints = 40;
    for(int i=0; i<29; i++) {
        TString Name = Form("totpe_r%d", i);
        TH1D* h1 = (TH1D*)ff->Get(Name);
        h1->Fit("gaus", "Q0");
        TF1 *gaus = (TF1*)h1->GetFunction("gaus");
        //cout << gaus->GetParameter(1) << " " << gaus->GetParError(1) << " " << gaus->GetParameter(2) << " " << gaus->GetParError(2) << endl;
        double totpe = gaus->GetParameter(1);
        double totpe_err = gaus->GetParError(1);
        double sigma = gaus->GetParameter(2);
        double sigma_err = gaus->GetParError(2);
        double resol = sigma/totpe;
        double resol_err = TMath::Sqrt(totpe_err*totpe_err*sigma*sigma/totpe/totpe/totpe/totpe+sigma_err*sigma_err/totpe/totpe);
        if(pemax<totpe) pemax = totpe;
        if(pemin>totpe) pemin = totpe;
        gTotPE->SetPoint(i, rmax/nPoints*i, totpe);
        gTotPE->SetPointError(i, 0, totpe_err);
        gRes->SetPoint(i, rmax/nPoints*i, resol);
        gRes->SetPointError(i, 0, resol_err);
        delete gaus;
        delete h1;
    }

    cout << "Max: " << pemax << " Min: " <<pemin <<endl;
    ff->Close();
    delete ff;

    TFile *ff1 = TFile::Open("results1MeV_nocorr.root");
    TGraphErrors* gTotPE1 = new TGraphErrors();
    for(int i=0; i<nPoints; i++) {
        TString Name = Form("totpe_r%d", i);
        TH1D* h1 = (TH1D*)ff1->Get(Name);
        double totpe = h1->GetMean();
        gTotPE1->SetPoint(i, rmax/nPoints*i, totpe);
        delete h1;
    }
    ff1->Close();
    delete ff1;


    gTotPE->SetName("totpe_corr");
    gRes->SetName("rsl");
    gTotPE1->SetName("totpe");
    TFile* out = new TFile("out1MeV.root", "recreate");
    gTotPE->Write();
    gTotPE1->Write();
    gRes->Write();
    out->Close();


    TGraph* gFV = new TGraph();
    gFV->SetPoint(0, 17.2*17.2*17.2, 0);
    gFV->SetPoint(1, 17.2*17.2*17.2, 10);
    gFV->SetPoint(2, 17.7*17.7*17.7, 10);
    gFV->SetPoint(3, 17.7*17.7*17.7, 0);
    gFV->SetFillStyle(3005);
    gFV->SetFillColor(kCyan+1);

    TGraph* gMaxMin = new TGraph();
    pemax = gTotPE->GetMean(2)*1.01;
    pemin = gTotPE->GetMean(2)*0.99;
    gMaxMin->SetPoint(0, 0, pemax);
    gMaxMin->SetPoint(1, rmax, pemax);
    gMaxMin->SetPoint(2, rmax, pemin);
    gMaxMin->SetPoint(3, 0, pemin);
    gMaxMin->SetFillStyle(3004);
    gMaxMin->SetFillColor(kYellow);

    TGraph* gTRZ = new TGraph();
    gTRZ->SetPoint(0, 15.8*15.8*15.8, 0);
    gTRZ->SetPoint(1, 15.8*15.8*15.8, 10);
    gTRZ->SetLineColor(kOrange+1);
    gTRZ->SetLineWidth(2);

    TCanvas* c1 = new TCanvas();
    gTotPE->SetLineColor(kBlue+1);
    gTotPE->SetMarkerStyle(20);
    gTotPE->SetMarkerColor(kBlue+1);
    gTotPE->SetLineWidth(2);
    gTotPE->SetTitle("edep v.s. edpe radius; radius^3/m^3; evis/MeV");
    //gTotPE->GetYaxis()->SetRangeUser(1,1.2);
    gTotPE->Draw("APL");
    gFV->Draw("F SAME");
    gMaxMin->Draw("F SAME");
    gTRZ->Draw("L SAME");

    TCanvas* c2 = new TCanvas();
    gRes->SetLineColor(kBlue+1);
    gRes->SetMarkerStyle(20);
    gRes->SetMarkerColor(kBlue+1);
    gRes->SetLineWidth(2);
    gRes->SetTitle("resolution v.s. edpe radius; radius^3/m^3; resolution");
    gRes->Draw("APL");
    gFV->Draw("F SAME");
    gTRZ->Draw("L SAME");

    TCanvas* c3 = new TCanvas();
    gTotPE1->SetLineColor(kPink+1);
    gTotPE1->SetMarkerStyle(20);
    gTotPE1->SetMarkerColor(kPink+1);
    gTotPE1->SetLineWidth(2);
    gTotPE1->Draw("APL");
}
