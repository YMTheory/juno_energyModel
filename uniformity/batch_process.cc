#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

int main(int argc, char* argv[])
//void batch_process()
{
    if(argc!=3) {
        throw "Wrong Input Parameters Number!";
        return 0;
    }
    int doFit = 1;
    cout << "do fitting in current process? " << doFit << endl;

    double R_cut = 15000;
    const int nbins = 40;
    double rmax = 17700; // mm
    double delta_R = TMath::Power(rmax/1000., 3)/nbins ;
    TH1D* hTotPE[nbins];
    TH1D* hTotPE15m = new TH1D("hTotPE15m", "", 2400, 0, 12);
    for(int i=0; i<nbins; i++) {
        hTotPE[i] = new TH1D("", "", 2400, 0, 12);
        TString Name = Form("radius%d", i);
        hTotPE[i]->SetName(Name);
    }


    t1->B
    Float_t m_evis; TBranch* b_evis;
    Float_t m_edepR; TBranch* b_edepR;

    TFile *ff = TFile::Open(argv[1]);
    if(!ff) { cout << "No such file " << argv[1] << endl; }
    TTree* tree = (TTree*)ff->Get("petree");
    tree->SetBranchAddress("evis", &m_evis, &b_evis);
    tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
    for(int i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        int ibin = int(TMath::Power(m_edepR/1000, 3)/delta_R); 
        //cout << m_edepR << " " << ibin  << " "<< m_evis << endl;
        hTotPE[ibin]->Fill(m_evis);
        if(m_edepR<R_cut)hTotPE15m->Fill(m_evis);
    }

    TGraphErrors* gTotPE = new TGraphErrors(); gTotPE->SetName("pe_graph");
    TGraphErrors* gResol = new TGraphErrors(); gResol->SetName("rsl_graph");
    if(doFit) {
        for(int i=0; i<nbins; i++) {
            hTotPE[i]->Fit("gaus" );
            TF1* func = (TF1*)hTotPE[i]->GetFunction("gaus");
            double mean = func->GetParameter(1);
            double mean_err = func->GetParError(1);
            double sigma = func->GetParameter(2);
            double sigma_err = func->GetParError(2);
            double resol = sigma/mean;
            double resol_err = TMath::Sqrt(mean_err*mean_err*sigma*sigma/mean/mean/mean/mean+sigma_err*sigma_err/mean/mean);
            gTotPE->SetPoint(i, i*delta_R, mean);
            gTotPE->SetPointError(i, 0, mean_err);
            gResol->SetPoint(i, i*delta_R, resol);
            gResol->SetPointError(i, 0, resol_err);
        }

        gTotPE->SetLineColor(kViolet+1);
        gTotPE->SetLineWidth(2);
        gTotPE->SetMarkerColor(kViolet+1);
        gTotPE->SetMarkerStyle(25);
        gTotPE->SetTitle("; Evis/MeV; totpe");
        gResol->SetLineColor(kViolet+1);
        gResol->SetLineWidth(2);
        gResol->SetMarkerColor(kViolet+1);
        gResol->SetMarkerStyle(25);
        gResol->SetTitle("; Evis/MeV; resolution");

    }

    TFile* out = new TFile(argv[2], "recreate");
    for(int i=0; i<nbins; i++) {
        hTotPE[i]->Write();
    }
    gTotPE->Write();
    gResol->Write();
    hTotPE15m->Write();
    out->Close();


    ff->Close();
    //
    return 1;

}
