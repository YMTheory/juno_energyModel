extern const int Nbins = 200;
extern double rsd_ratio[500]={0.}; 
extern double rsd_sigma[500]={0.};

void read_corr()
{
    TH1D* hTotPE[200];
    TFile* in  = TFile::Open("./smear500mm.root", "read");
    for(int i=0; i<200; i++) {
            TString name = Form("sub%d",i);
            hTotPE[i] = (TH1D*)in->Get(name);
            rsd_ratio[i] = hTotPE[i]->GetMean();
            rsd_sigma[i] = hTotPE[i]->GetStdDev();

    }
}

void vertex_study()
{
    read_corr();
    for(int i=0; i<Nbins; i++) {
        rsd
    }
}
