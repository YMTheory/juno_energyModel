void average_corr_gamma()
{
    auto mg = new TMultiGraph();

    const int N = 4;
    string source[N] = {"Cs137", "K40", "O16", "Ge68"};
    Color_t color[N] = {36, 40, 44, 46};
    TString label[N] = {"Cs137", "K40", "O16", "Ge68"};

    TGraphErrors* ge[N]; TLegend* ll = new TLegend();
    string path = "/junofs/users/yumiao/simulation/energy_model/uniformity/gamma/analysis_";
    string suffix = ".root";
    for(int i=0; i<N; i++) {
        string filename = path + source[i] + suffix;
        TFile *ff = new TFile(filename.c_str(), "read");
        ff->GetObject("pe_graph", ge[i]);
        double *y = ge[i]->GetY();
        double* yerr = ge[i]->GetEY();
        double scale = 1/y[0];
        //double scale = 1/ge[i]->GetMean(2);
        for(int j=0; j<ge[i]->GetN(); j++) {
            ge[i]->SetPointY(j, y[j]*scale);
            ge[i]->SetPointError(j, 0, yerr[i]*scale);
        }
        ge[i]->SetLineWidth(2);
        ge[i]->SetLineColor(color[i]);
        ge[i]->SetMarkerStyle(24);
        ge[i]->SetMarkerColor(color[i]);
        //ge[i]->SetName(graph_name[i].c_str());
        ge[i]->SetLineStyle(kDashed);

        mg->Add(ge[i]);
        ll->AddEntry(ge[i], label[i], "L");
    }

    // Remove points outside FV:
    //for(int i=0; i<ge[0]->GetN(); i++) {
    //    for(int j=0; j<N; j++) {
    //        ge[j]->RemovePoint(i);
    //    }
    //} 

    // average curve :
    TGraphErrors* gTot = new TGraphErrors();
    double *xx = ge[0]->GetX();
    for(int i=0; i<ge[0]->GetN(); i++) {
        double tmp_y = 0;
        double tmp_yerr = 0;
        for(int j=0; j<N; j++) {
            //cout << i << " " << j << " " <<  (ge[i]->GetY())[i] << endl;
            tmp_y += (ge[j]->GetY())[i];
            tmp_yerr += (ge[j]->GetEY())[i];
        }
        tmp_y /= N; tmp_yerr /= N;
        gTot->SetPoint(i, xx[i], tmp_y );
        gTot->SetPointError(i, 18*18*18/40., tmp_yerr);
    }
    gTot->SetLineWidth(2);
    gTot->SetLineColor(46);
    gTot->SetMarkerColor(46);
    gTot->SetMarkerStyle(20);
    //mg->Add(gTot);
    //ll->AddEntry(gTot, "average", "L");

    TCanvas* cc = new TCanvas();
    mg->Draw("APL");
    ll->Draw("SAME");
}
