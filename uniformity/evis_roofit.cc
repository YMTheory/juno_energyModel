// CB fitting evis distribution

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooPlot.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include <iostream>
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
//#include <fstream.h>

//using namespace std;
using namespace RooFit;

void evis_roofit()
{
    TFile* f1 = TFile::Open("../../uniformity/electron/analysis8MeV_test.root");
    //TH1D* hTotPE15m = (TH1D*)f1->Get("radius39");
    TH1D* hTotPE15m = (TH1D*)f1->Get("hTotPEall");

    double par[4];
    par[0] = 0.5;
    par[1] = 2.5;
    par[2] = 9;
    par[3] = 0.5;

    RooRealVar evis("evis","",9, 8,10);
    std::cout<<"LINE = "<<__LINE__<<std::endl;
    RooRealVar a1("CB_a1", "", par[0], -2 ,4);
    RooRealVar n1("CB_n1", "", par[1], 0, 10);
    RooRealVar mean("CB_mean", "", par[2], 8, 10);
    RooRealVar sigma("CB_sigma","", par[3], 0, 1);
    RooCBShape CB_pdf("cb1", "", evis, mean, sigma, a1, n1);
    RooGaussian gauss_pdf("gauss_pdf", "gaussian PDF", evis, mean, sigma);

    RooPlot* xframe = evis.frame(Title("CB p.d.f"));
    RooDataHist dh("dh","dh", evis, Import(*hTotPE15m));
    dh.plotOn(xframe);

    //gauss_pdf.fitTo(dh);
    //gauss_pdf.plotOn(xframe);
    //gauss_pdf.paramOn(xframe, Layout(0.45));
    CB_pdf.fitTo(dh);
    CB_pdf.plotOn(xframe);
    CB_pdf.paramOn(xframe, Layout(0.45));

    double chi2 = xframe->chiSquare();
    double dof = 2400./12.*2-4;
    cout << "chis/ndf : " << chi2/dof << endl;

    TCanvas* c = new TCanvas("roofit", "roofit", 800,600);
    xframe->Draw();
}
