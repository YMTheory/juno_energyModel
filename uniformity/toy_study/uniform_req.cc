void uniform_req()
{
    
    const int N = 8;
    double rsd_uniform[N] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10};
    double rsl[N] = {2.622, 2.634, 2.672, 2.744, 2.828, 2.938, 3.216, 3.726};

    TCanvas* cc = new TCanvas(); cc->cd();

    TGraph* g1 = new TGraph(N, rsd_uniform, rsl);
    g1->SetMarkerStyle(24);
    g1->SetMarkerColor(35);
    g1->SetLineColor(35);
    g1->SetLineWidth(2);
    g1->SetTitle(";residual non-uniformity; resolution(\%)");
    g1->Draw("APL");

    TGraph *g2 = new TGraph();
    g2->SetPoint(0, 0, 2.682);
    g2->SetPoint(1,0.1,2.682);
    g2->SetLineColor(43);
    g2->SetLineWidth(2);
    g2->Draw("PL");

    TText* t1 = new TText(0.06, 2.7, "center samples"); 
    t1->SetTextColor(43);
    t1->Draw("L SAME");

    TGraph* gZone = new TGraph();
    gZone->SetPoint(0, 0, 0);
    gZone->SetPoint(1, 0.02, 0);
    gZone->SetPoint(2, 0.02, 4);
    gZone->SetPoint(3, 0, 4);
    gZone->SetFillColor(kCyan);
    gZone->SetFillStyle(3004);
    gZone->Draw("F SAME");

//    cc->SaveAs("./plots/uniform_on_rsl.pdf");

/*
    TCanvas* cg2 = new TCanvas(); cg2->Draw(); cg2->SetGrid();
    double scale = 3350/2.22;
    const int N1 = 5;
    //double evis[N1] = {};
    double etrue[N1] = {1,2,4,6,8};
    double rsl_center[N1] = {2.682/100, 2.011/100, 1.479/100, 1.240/100, 1.098/100};
    double rsl_uni[N1] = {2.634/100, 1.971/100, 1.462/100, 1.234/100, 1.102/100};
    double rsl_ideal[N1] = {2.622/100, 1.958/100, 1.441/100, 1.206/100, 1.069/100};
    double rsl_more[N1];
    for(int i=0; i<N1; i++) {
        rsl_more[i] = TMath::Sqrt(rsl_uni[i]*rsl_uni[i] - rsl_ideal[i]*rsl_ideal[i]);
        //rsl_more[i] = TMath::Sqrt(rsl_uni[i]*rsl_uni[i] - rsl_ideal[i]*rsl_ideal[i]);
    }
    TGraph* g3 = new TGraph(N1, etrue, rsl_center);
    TGraph* g4 = new TGraph(N1, etrue, rsl_uni);
    TGraph* g5 = new TGraph(N1, etrue, rsl_ideal);
    TGraph* g6 = new TGraph(N1, etrue, rsl_more);
    g3->SetMarkerStyle(24);
    g3->SetMarkerColor(36);
    g3->SetLineColor(36);
    g3->SetLineWidth(2);
    g4->SetMarkerStyle(25);
    g4->SetMarkerColor(42);
    g4->SetLineColor(42);
    g4->SetLineWidth(2);
    g5->SetMarkerStyle(26);
    g5->SetMarkerColor(30);
    g5->SetLineColor(30);
    g5->SetLineWidth(2);
    g3->SetTitle(";Etrue/MeV;resolution");
    g3->Draw("APL");
    g4->Draw("PL SAME");
    g5->Draw("PL SAME");
    
    TLegend* ll = new TLegend();
    ll->AddEntry(g3, "center", "L");
    ll->AddEntry(g4, "uniform w/ 1\% non-uniformity", "L");
    ll->AddEntry(g5, "uniform w/o non-uniformity", "L");
    ll->Draw("SAME");


    TCanvas* cg3 = new TCanvas(); cg3->SetGrid();
    g6->SetTitle("sqrt(uniform_rsl**2-uniform_idea**2); Etrue/MeV;");
    g6->SetMarkerStyle(24);
    g6->SetMarkerColor(36);
    g6->SetLineColor(36);
    g6->SetLineWidth(2);
    g6->Draw("APL");
    */
}
