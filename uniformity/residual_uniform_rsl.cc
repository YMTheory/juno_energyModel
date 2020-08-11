void  residual_uniform_rsl()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

    std::vector<string> filename;
    filename.push_back("analysis1MeV.root");
    filename.push_back("analysis2MeV.root");
    filename.push_back("analysis4MeV.root");
    filename.push_back("analysis6MeV.root");
    filename.push_back("analysis8MeV.root");
    
    TFile* file;
    TH1D* hEvis;
    double evis_data[5];
    TGraphErrors* hResCorr = new TGraphErrors();
    for(int i=0; i<filename.size(); i++) {
        file = TFile::Open(filename[i].c_str());
        hEvis = (TH1D*)file->Get("hTotPE15m");
        cout << hEvis->GetMean() << endl;
        hEvis->Fit("gaus", "0");
        TF1 *gaus = (TF1*)hEvis->GetFunction("gaus");
        cout << gaus->GetParameter(1) << " " << gaus->GetParError(1) << " " << gaus->GetParameter(2) << " " << gaus->GetParError(2) << endl;
        double totpe = gaus->GetParameter(1); evis_data[i]= totpe;
        double totpe_err = gaus->GetParError(1);
        double sigma = gaus->GetParameter(2);
        double sigma_err = gaus->GetParError(2);
        double resol = sigma/totpe;
        double resol_err = TMath::Sqrt(totpe_err*totpe_err*sigma*sigma/totpe/totpe/totpe/totpe+sigma_err*sigma_err/totpe/totpe);
        hResCorr->SetPoint(i, totpe, resol); 
        hResCorr->SetPointError(i, 0, resol_err);
    }
    TCanvas* c2 = new TCanvas();
    hResCorr->SetLineColor(36);
    hResCorr->SetMarkerStyle(20);
    hResCorr->SetMarkerColor(36);
    hResCorr->SetLineWidth(2);
    hResCorr->SetTitle(";Evis/MeV;resolution");
    hResCorr->Draw("APL");

    // center samples ~
    double scale = 1363;
    double datax[4] = {1473/scale, 3000/scale, 6061./scale, 9115/scale};
    double datay[4] = {40.38/1473, 60.91/3000, 92.38/6061, 118.2/9115}; /// all pmt pe
    TGraph* hRes = new TGraph(4, datax, datay);
    hRes->SetLineColor(kOrange+2);
    hRes->SetMarkerStyle(20);
    hRes->SetMarkerColor(kOrange+2);
    hRes->SetLineWidth(2);
    //hRes->Draw("PL SAME");

    // calculation results ~
    //double calcx[5] = {1472.65/scale, 2945/scale, 5891.21/scale, 8837.02/scale};
    //double calcy[4] = {0.02978/1.108, 0.04482/2.256, 0.06813/4.558, 0.08694/6.8581};  // only lpmt p.e.
    double calcx[5] = {1.114, 2.268, 4.581, 6.892, 9.203};
    double calcy[5] = {0.02948/1.114, 0.04463/2.268, 0.06788/4.581, 0.08733/6.892, 0.1036/9.203};
    TGraph* hResCalc = new TGraph(5, calcx, calcy);
    hResCalc->SetLineColor(45);
    hResCalc->SetMarkerStyle(21);
    hResCalc->SetMarkerColor(45);
    hResCalc->SetLineWidth(2);
    //hResCalc->SetLineStyle(kDashed);
    hResCalc->Draw("PL SAME");

    TF1* f1 = new TF1("f1", "TMath::Sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/x/x)", 0, 10);
    //hResCorr->Fit(f1, "R0");
    //TGraph* gFit1 = new TGraph();
    //for(int i=0; i<1000; i++) {
    //    gFit1->SetPoint(i, 10./1000*i, f1->Eval(10./1000*i));
    //}
    //gFit1->SetLineColor(kGreen+2);
    //gFit1->SetLineWidth(2);
    //gFit1->Draw("L SAME");


    TF1* f2 = new TF1("f2", "TMath::Sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/x/x)", 0, 10);
    //hResCalc->Fit(f2, "R0");
    //TGraph* gFit2 = new TGraph();
    //for(int i=0; i<1000; i++) {
    //    gFit2->SetPoint(i, 10./1000*i, f1->Eval(10./1000*i));
    //}
    //gFit2->SetLineColor(kViolet+1);
    //gFit2->SetLineWidth(2);
    //gFit2->Draw("L SAME");


    TLegend* ll = new TLegend();
    //ll->AddEntry(hRes, "electron @center", "l");
    ll->AddEntry(hResCorr, "electron uniformly within 15.5m", "l");
    ll->AddEntry(hResCalc, "electron uniformly prediction within 15.5m", "l");
    ll->Draw("SAME");

}
