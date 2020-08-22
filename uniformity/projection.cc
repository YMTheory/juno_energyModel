void projection()
{
    TFile* ff = new TFile("./subDet-corr.root");
    TGraph* g1 = (TGraph*)ff->Get("average");
    TH1D* h1= new TH1D("h1", "", 100, 0.95, 1.05);
    for(int i=0; i<180; i++) {
        h1->Fill(g1->GetPointY(i));
    }

    h1->Draw();
}
