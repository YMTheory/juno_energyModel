Double_t fit_func1(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t val = 0;
    Double_t mu = par[0];
    Double_t A = par[1];
    Double_t sigma = par[2];
    val = A*TMath::Exp(-(xx-mu)*(xx-mu)/2/sigma/sigma);
    return val;
}   

Double_t fit_func2(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t val = 0;
    Double_t mu = par[0];
    Double_t A = par[1];
    Double_t sigma1 = par[2];
    Double_t sigma2 = par[3];
    if( xx<mu ) {
        val = A*TMath::Exp(-(xx-mu)*(xx-mu)/2/sigma1/sigma1);
    } else if( xx>=mu) {
        val = A*TMath::Exp(-(xx-mu)*(xx-mu)/2/sigma2/sigma2);
    }
    return val;
}   

void evis_tf1fit(int opt)
{
    gStyle->SetOptFit(1111);

    TFile* f1 = new TFile("../../uniformity/electron/pe1MeV.root"); 
    TTree *t1; f1->GetObject("petree", t1);
    Float_t m_evis, m_edepR;
    t1->SetBranchAddress("evis", &m_evis);
    t1->SetBranchAddress("edepR", &m_edepR);

    TH1D* h_evis = new TH1D("h_evis", "", 200, 0.8, 1.4);
    for(int i=0; i<t1->GetEntries(); i++) {
        t1->GetEntry(i);
        if( m_edepR<17200 ) {
            h_evis->Fill(m_evis);
        }
    }
    
    TF1* func;
    if(opt == 0) { // single gaussian fitting
        func = new TF1("fit_func", fit_func1, 0.8, 1.4, 3);
    } else if(opt == 1) {  // double gaussin fitting
        func = new TF1("fit_func", fit_func2, 0.8, 1.4, 4);
    } else if (opt == 2) {
        func = new TF1("fit_func", fit_func1, 1.114, 1.4, 3);
    } else if (opt ==3 ) {
        func = new TF1("fit_func", fit_func1, 0.8, 1.114, 3);
    }

    func->SetParameters(1.1, 100, 0.028, 0.026);
    h_evis->Fit(func, "RE");
    h_evis->Draw();

}
