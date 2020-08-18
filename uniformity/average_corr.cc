void average_corr()
{
    auto mg = new TMultiGraph();
    
    const int N = 5;
    int mom[5] = {1, 2, 4, 6, 8};
    string graph_name[N] = {"KE1", "KE2", "KE4", "KE6", "KE8"};
    Color_t color[N] = {36, 38, 40, 42, 44};
    TString label[N] = {"KE=1MeV", "KE=2MeV", "KE=4MeV", "KE=6MeV", "KE=8MeV"};

    TGraphErrors* ge[N]; TLegend* ll = new TLegend();
    string path = "/junofs/users/yumiao/simulation/energy_model/uniformity/electron/analysis";
    string suffix = "MeV_tmp.root";
    for(int i=0; i<N; i++) {
        string filename = path+to_string(mom[i])+suffix;
        cout << "Processing " << filename << endl;
        TFile *ff = new TFile(filename.c_str(), "read");
        ff->GetObject("pe_graph", ge[i]);
        double *y = ge[i]->GetY();
        double* yerr = ge[i]->GetEY();
        double scale = 1/y[0];
        for(int j=0; j<ge[i]->GetN(); j++) {
            ge[i]->SetPointY(j, y[j]*scale);
            ge[i]->SetPointError(j, 0, yerr[i]*scale);
        }
        ge[i]->SetLineWidth(2);
        ge[i]->SetLineColor(color[i]);
        ge[i]->SetMarkerStyle(24);
        ge[i]->SetMarkerColor(color[i]);
        ge[i]->SetName(graph_name[i].c_str());
        ge[i]->SetLineStyle(kDashed);

        mg->Add(ge[i]);
        ll->AddEntry(ge[i], label[i], "LP");
    }

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

    // gamma sources 
    const int N1 = 4;
    string source[N1] = {"Cs137", "K40", "O16", "Ge68"};
    Color_t color1[N1] = {28, 30, 32, 34};
    TString label1[N1] = {"Cs137", "K40", "O16", "Ge68"};

    TGraphErrors* ge1[N1];
    string path1 = "/junofs/users/yumiao/simulation/energy_model/uniformity/gamma/analysis_";
    string suffix1 = ".root";
    for(int i=0; i<N1; i++) {
        string filename = path1 + source[i] + suffix1;
        cout << "Processing " << filename << endl;
        TFile *ff = new TFile(filename.c_str(), "read");
        ff->GetObject("pe_graph", ge1[i]);
        double *y = ge1[i]->GetY();
        double* yerr = ge1[i]->GetEY();
        double scale = 1/y[0];
        //double scale = 1/ge1[i]->GetMean(2);
        for(int j=0; j<ge1[i]->GetN(); j++) {
            //cout << j << " " <<  y[j]*scale << endl;
            ge1[i]->SetPointY(j, y[j]*scale);
            ge1[i]->SetPointError(j, 0, yerr[i]*scale);
        }
        ge1[i]->SetLineWidth(2);
        ge1[i]->SetLineColor(color[i]);
        ge1[i]->SetMarkerStyle(25);
        ge1[i]->SetMarkerColor(color[i]);
        ge1[i]->SetLineStyle(kDotted);

        mg->Add(ge1[i]);
        ll->AddEntry(ge1[i], label1[i], "LP");
    }

    // average curve :
    //TGraphErrors* gTot1 = new TGraphErrors();
    //double *xx1 = ge1[0]->GetX();
    //for(int i=0; i<ge1[0]->GetN(); i++) {
    //    double tmp_y = 0;
    //    double tmp_yerr = 0;
    //    for(int j=0; j<N; j++) {
    //        //cout << i << " " << j << " " <<  (ge[i]->GetY())[i] << endl;
    //        tmp_y += (ge1[j]->GetY())[i];
    //        tmp_yerr += (ge1[j]->GetEY())[i];
    //    }
    //    tmp_y /= N; tmp_yerr /= N;
    //    gTot1->SetPoint(i, xx1[i], tmp_y );
    //    gTot1->SetPointError(i, 18*18*18/40., tmp_yerr);
    //}
    //gTot1->SetLineWidth(2);
    //gTot1->SetLineColor(46);
    //gTot1->SetMarkerColor(46);
    //gTot1->SetMarkerStyle(20);
    //mg->Add(gTot);
    //ll->AddEntry(gTot, "average", "L");




    TCanvas* cc = new TCanvas();
    mg->Draw("APL");
    ll->Draw("SAME");
}
