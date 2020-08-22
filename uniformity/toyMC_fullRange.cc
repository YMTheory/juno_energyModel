#include "data_valid.h"

extern double rsd_ratio[500]={0.}; 
extern double rsd_sigma[500]={0.};
extern double rsd_ratio_sec[500] = {0.};
//TGraphErrors* gResol_data;
const double energy_scale = 1300; //3350/2.22;

bool doVrtSmear = true;

void draw_NU(double *f);

TGraph* gSecCorr = new TGraph();
void read_SecCorr();

double NPE_smear(double* x, double* par) {
    double val = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];

     //sigma2 fitting 
    if(val <= 0)
        return 0;
    else {
        double fval = TMath::Sqrt(val);    
        return fval;
    }

    //if(val < 0) return 0;
    //return val;
}

double abcModel(double* x, double* par) {
    double a = par[0];
    double b = par[1];
    double c = par[2];
    double E = x[0];
    return TMath::Sqrt(a/E+b+c/E/E);
}


void read_corr()
{
    
    //Float_t m_evis; TBranch* b_evis;
    //Float_t m_edepZ; TBranch* b_edepZ;
    //Float_t m_edepR; TBranch* b_edepR;

    //TString name = "../../uniformity/electron/pe1MeV.root";  //1MeV e- correction results
    //TString name = "../../uniformity/electron/smear1MeV-0mm.root";  //1MeV e- correction results, smear 0mm
    //TString name = "../../uniformity/gamma/Ge68.root";  // gamma source correction results, no smear

    //TFile* ff = TFile::Open(name, "read");
    //TTree* tree = (TTree*)ff->Get("petree");
    //tree->SetBranchAddress("evis", &m_evis, &b_evis);
    //tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
    //tree->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
    //TProfile2D* prof = new TProfile2D("prof", "", 20, 0, 18*18*18, 10, -1, 1, 0., 12);
    ////TH1D* hTotPE_FV = new TH1D("hTotPE_FV", "", 2400, 0, 12);
    //double delta_R = 18*18*18./20.;
    //for(int i=0; i<tree->GetEntries(); i++) {
    //    tree->GetEntry(i);
    //    prof->Fill(m_edepR/1000*m_edepR/1000*m_edepR/1000, m_edepZ/m_edepR, m_evis);
    //    //if(m_edepR<15000) {hTotPE_FV->Fill(m_evis);}
    //}

    //int idx = 0;
    //for(int i=0; i<prof->GetNbinsX(); i++) {
    //    for(int j=0; j<prof->GetNbinsY(); j++) {
    //        rsd_ratio[idx] = prof->GetBinContent(i+1, j+1);
    //        idx++;
    //    }
    //}


    

    TH1D* hTotPE[200];
    TFile* in  = TFile::Open("./smear0mm.root", "read");
    for(int i=0; i<200; i++) {
            TString name = Form("sub%d",i);
            hTotPE[i] = (TH1D*)in->Get(name);
            rsd_ratio[i] = hTotPE[i]->GetMean();
            rsd_sigma[i] = hTotPE[i]->GetStdDev();

    }
}


void toyMC_fullRange()
{
    gRandom->SetSeed(39);

    TF1* fNPEsmear = new TF1("fNPEsmear", NPE_smear, 0, 12000, 3);
    fNPEsmear->SetParameters(-3.69505e2, 1.30168, 1.90109e-05);   // sigma2 fitting relation
    //fNPEsmear->SetParameters(2.06408e+01, 1.40885e-02, -4.04637e-07);  // sigma fitting relation

    /*********************************************************/
    //draw NPE smear :
    //TGraph* gNPEsmear = new TGraph();
    //for(int i=1; i<1000; i++) {
    //    int j = i*10;
    //    gNPEsmear->SetPoint(i-1, j, fNPEsmear->Eval(j));
    //}

    //gNPEsmear->SetLineWidth(2);
    //gNPEsmear->SetLineColor(42);
    //gNPEsmear->SetTitle("NPE smear; NPE; NPE sigma");
    //gNPEsmear->Draw("AL");
    /*********************************************************/


    auto mg1 = new TMultiGraph();
    auto led1 = new TLegend();


    /***************** center sample data *******************/
    const int N = 5;
    double center_pe[N] = {1.47416e+03, 3.00112e+03, 6.06205e+03, 9.11628e+03, 1.21718e+04};
    double center_peerr[N] = {5.71360e-01, 8.55662e-01, 1.29060e+00, 1.67041e+00, 1.98174e+00};
    double center_sigma[N] = {3.98369e+01, 5.94883e+01, 8.91297e+01, 1.15738e+02, 1.36322e+02};
    double center_sigmaerr[N] = {3.97852e-01, 6.18851e-01, 9.36577e-01, 1.26519e+00, 1.46872e+00};
    TGraphErrors* gResol_center_data = new TGraphErrors();
    for(int i=0; i<N; i++ ) {
        double mu2 = center_pe[i];
        double mu2_err = center_peerr[i];
        double sigma2 = center_sigma[i];
        double sigma2_err = center_sigmaerr[i];
        double rsl2 = sigma2/mu2;
        double rsl2_err = TMath::Sqrt(sigma2_err*sigma2_err/mu2/mu2 + mu2_err*mu2_err*sigma2*sigma2/mu2/mu2/mu2/mu2);

        gResol_center_data->SetPoint(i, mu2, rsl2);
        gResol_center_data->SetPointError(i, 0, rsl2_err);
    }
    gResol_center_data->SetLineColor(9);
    gResol_center_data->SetLineWidth(2);
    gResol_center_data->SetMarkerStyle(21);
    gResol_center_data->SetMarkerColor(9);
    mg1->Add(gResol_center_data);
    led1->AddEntry(gResol_center_data, "center simData", "PL");
    /*********************************************************/





    TString name = "../../uniformity/electron/totalpe2r2theta_40rbins20thetabins_1MeV.root";
    TFile* ft = new TFile(name, "read");
    TProfile2D* prof2d = (TProfile2D*)ft->Get("pe2r2theta");
    const int rbins = prof2d->GetNbinsX();
    const int thetabins = prof2d->GetNbinsY();
    const int totalbins = rbins*thetabins;
    cout << totalbins << endl;
    double pe_array[totalbins];
    int idx = 0;
    for(int irbin=0; irbin<rbins; irbin++) {
        for(int ithetabin=0; ithetabin<thetabins; ithetabin++) {
            pe_array[idx] = prof2d->GetBinContent(irbin+1, ithetabin+1);
            //cout << irbin << " " << ithetabin << " " << pe_array[idx] <<endl;
            idx++;
        }
    }
    ft->Close();


    read_corr();
    read_SecCorr();
    
    double subID[200] ; for(int i=0; i<200; i++) {subID[i]=i;}
    TGraph* gRNU = new TGraph(180, subID,rsd_ratio);
    TCanvas* tempc = new TCanvas();
    gRNU->Draw("AL");

    //data_valid valid;
    //int const Nmom = 5;
    //int mom[Nmom] = {1,2,4,6,8};
    //valid.read_data(Nmom, mom, rsd_ratio_sec);
    //TGraphErrors* gResol_data = valid.fit_data();
    //for(int i=0; i<gResol_data->GetN(); i++) {  // do energy scale 
    //    double origin_x = gResol_data->GetPointX(i);
    //    gResol_data->SetPointX(i, origin_x*energy_scale);
    //}


    //const int N = 5;
    Color_t color1[N] = {36, 38, 40, 42, 44};
    Color_t color2[1] = {48};
    TString label1[5] = {"scale factor: 1.0", "scale factor: 2.0", "scale factor: 3.0", "scale factor: 4.0", "scale factor: 5.0"};
    double factor[N] = {1, 2, 3, 4, 5};
    //draw_NU(factor);
    TGraphErrors* gResol_real[N];
    TGraphErrors* gResol_smear[N];
    TGraphErrors* gResol_ideal;
    TGraphErrors* gResol_center;
    TH1D* hTotPE_real = new TH1D("hTotPE_real", "", 2400, 0, 12000);
    TH1D* hTotPE_ideal = new TH1D("hTotPE_ideal", "", 2400, 0, 12000);
    TH1D* hTotPE_center = new TH1D("hTotPE_center", "", 2400, 0, 12000);
    for(int iloop=0; iloop<1; iloop++) {
        cout <<  "LOOPI ================> " << iloop << endl;
        gResol_real[iloop] = new TGraphErrors();
        gResol_smear[iloop] = new TGraphErrors();
        if( iloop==0 ) { gResol_ideal = new TGraphErrors(); gResol_center = new TGraphErrors(); }
        for(int jloop=0; jloop<20; jloop++) {
            hTotPE_real->Reset();
            if(iloop==0) { hTotPE_ideal->Reset(); hTotPE_center->Reset(); }
            double npe = (jloop+2)*500;
            //cout << " LOOPJ ===> " <<  jloop << "   NPE ===> " << npe << endl;
            idx = 0;
            while(idx<100000) {  // sampling uniformly

                int ibin = int(gRandom->Uniform(0,totalbins));
                if(ibin>=180 ) continue;  // FV cut
                //if(ibin<=100 or ibin>=180 ) continue;  // FV cut
                int npe1 = npe*((rsd_ratio[ibin]-rsd_ratio[0])/rsd_ratio[0]*factor[iloop]+1);  //npe in current bin
                //int npe1 = npe*((rsd_ratio[ibin]-rsd_ratio[0])/rsd_ratio[0]*factor[iloop]+1)/gSecCorr->GetPointY(ibin);  // do secondary correction
                int totpe_real = int(gRandom->Gaus(npe1, fNPEsmear->Eval(npe1)));
                //totpe_real = int(totpe_real/gSecCorr->GetPointY(ibin));  // secondary correction
                //cout << ibin << " " << gSecCorr->GetPointY(ibin) << endl;
                hTotPE_real->Fill(totpe_real);
                if( iloop==0 ) {
                    //int totpe_ideal = int(gRandom->Gaus(npe/pe_array[0]*pe_array[ibin], fNPEsmear->Eval(npe/pe_array[0]*pe_array[ibin])));
                    //totpe_ideal = int(totpe_ideal / pe_array[ibin] * pe_array[0]);
                    int totpe_ideal = int(gRandom->Gaus(npe, fNPEsmear->Eval(npe/pe_array[0]*pe_array[ibin])/pe_array[ibin]*pe_array[0] ));
                    hTotPE_ideal->Fill(totpe_ideal);

                    int totpe_center = int(gRandom->Gaus(npe, fNPEsmear->Eval(npe)));
                    hTotPE_center->Fill(totpe_center);
                }
                //cout << ibin << " " << npe1 << " " << fNPEsmear->Eval(npe1) << " " << npe << " " << fNPEsmear->Eval(npe/pe_array[ibin]*pe_array[0]) << endl; 
                idx++;
            }

            hTotPE_real->Fit("gaus", "Q0");
            TF1* f1 = (TF1*)hTotPE_real->GetFunction("gaus");
            double mu1 = f1->GetParameter(1);
            double mu1_err = f1->GetParError(1);
            double sigma1 = f1->GetParameter(2);
            double sigma1_err = f1->GetParError(2);
            double rsl1 = sigma1/mu1;
            double rsl1_err = (sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);
            if(iloop == 0) {
                hTotPE_ideal->Fit("gaus", "Q0");
                TF1* f2 = (TF1*)hTotPE_ideal->GetFunction("gaus");
                double mu2 = f2->GetParameter(1);
                double mu2_err = f2->GetParError(1);
                double sigma2 = f2->GetParameter(2);
                double sigma2_err = f2->GetParError(2);
                double rsl2 = sigma2/mu2;
                double rsl2_err = TMath::Sqrt(sigma2_err*sigma2_err/mu2/mu2 + mu2_err*mu2_err*sigma2*sigma2/mu2/mu2/mu2/mu2);
                gResol_ideal->SetPoint(jloop, mu2, rsl2);
                gResol_ideal->SetPointError(jloop, 0, rsl2_err);

                hTotPE_center->Fit("gaus", "Q0");
                TF1* f3 = (TF1*)hTotPE_center->GetFunction("gaus");
                double mu3 = f3->GetParameter(1);
                double mu3_err = f3->GetParError(1);
                double sigma3 = f3->GetParameter(2);
                double sigma3_err = f3->GetParError(2);
                double rsl3 = sigma3/mu3;
                double rsl3_err = TMath::Sqrt(sigma3_err*sigma3_err/mu3/mu3 + mu3_err*mu3_err*sigma3*sigma3/mu3/mu3/mu3/mu3);
                gResol_center->SetPoint(jloop, mu3, rsl3);
                gResol_center->SetPointError(jloop, 0, rsl3_err);
            }

            //cout << jloop << hTotPE_ideal->GetMean() << " " << hTotPE_ideal->GetStdDev() << endl;

            gResol_real[iloop]->SetPoint(jloop, mu1, rsl1);
            gResol_real[iloop]->SetPointError(jloop, 0, rsl1_err);

        }
        if(iloop == 0) {
            gResol_ideal->SetMarkerColor(color2[iloop]);
            gResol_ideal->SetMarkerStyle(25);
            gResol_ideal->SetLineColor(color2[iloop]);
            gResol_ideal->SetLineWidth(2);
            mg1->Add(gResol_ideal);
            led1->AddEntry(gResol_ideal,"ideal", "PL");

            gResol_center->SetMarkerColor(25);
            gResol_center->SetMarkerStyle(26);
            gResol_center->SetLineColor(25);
            gResol_center->SetLineWidth(2);
            mg1->Add(gResol_center);
            led1->AddEntry(gResol_center,"center", "PL");
        }
        gResol_real[iloop]->SetMarkerColor(color1[iloop]);
        gResol_real[iloop]->SetMarkerStyle(24);
        gResol_real[iloop]->SetLineColor(color1[iloop]);
        gResol_real[iloop]->SetLineWidth(2);
        mg1->Add(gResol_real[iloop]);
        led1->AddEntry(gResol_real[iloop], label1[iloop], "LP");
    }

    delete hTotPE_real;
    delete hTotPE_ideal;

    TCanvas* cc = new TCanvas(); 
    mg1->SetTitle("ToyMC Resolution; NPE; resolution");
    mg1->Draw("APL");
    led1->Draw("SAME");

    //TCanvas* cg = new TCanvas();
    //gResol_data->SetMarkerColor(kViolet+1);
    //gResol_data->SetMarkerStyle(20);
    //gResol_data->SetLineColor(kViolet+1);
    //gResol_data->SetLineWidth(2);
    ////gResol_data->Draw("APL");
    //gResol_data->Draw("PL SAME");

    // fitting with abc model:
    gStyle->SetOptFit(1111);
    TF1* fAbcModel_ideal = new TF1("fAbcModel_ideal", abcModel, 0, 12000, 3);
    fAbcModel_ideal->SetParameters(0.98*0.98, 6.62*6.62*1e-6, 0);
    gResol_ideal->Fit(fAbcModel_ideal);
    TF1* fAbcModel_real = new TF1("fAbcModel_real", abcModel, 0, 12000, 3);
    fAbcModel_real->SetParameters(1.3, 6.5e-5, -300);
    gResol_real[0]->Fit(fAbcModel_real);
    TF1* fAbcModel_center = new TF1("fAbcModel_center", abcModel, 0, 12000, 3);
    fAbcModel_center->SetParameters(0.98*0.98, 6.62*6.62*1e-6, 0);
    gResol_center->Fit(fAbcModel_center);
    TGraphErrors* gExtra = new TGraphErrors();
    //double *err_real = gResol_real[0]->GetEY();
    //double *err_center = gResol_center->GetEY();
    for(int i=0; i<1000; i++) {
        double npe2 = 1500+10*i;
        //gExtra->SetPoint(i, npe2, fAbcModel_real->Eval(npe2)- fAbcModel_center->Eval(npe2));
        double fAbcModel_center_err = TMath::Sqrt(  TMath::Power(fAbcModel_center->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_center->GetParError(1),2) +  TMath::Power(fAbcModel_center->GetParError(2)/npe2/npe2, 2)  );
        double fAbcModel_real_err = TMath::Sqrt(  TMath::Power(fAbcModel_real->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_real->GetParError(1),2) +  TMath::Power(fAbcModel_real->GetParError(2)/npe2/npe2, 2)  );
        gExtra->SetPoint(i, npe2, (fAbcModel_real->Eval(npe2)*fAbcModel_real->Eval(npe2) - fAbcModel_center->Eval(npe2)*fAbcModel_center->Eval(npe2)));
        gExtra->SetPointError(i, 0, TMath::Sqrt(TMath::Power(2*fAbcModel_center_err, 2) + TMath::Power(fAbcModel_real_err,2)) );
        //double fAbcModel_center_err_abs = TMath::Sqrt(  TMath::Power(fAbcModel_center->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_center->GetParError(1),2) +  TMath::Power(fAbcModel_center->GetParError(2)/npe2/npe2, 2)  ) / TMath::Sqrt(fAbcModel_center->Eval(npe2));
        //double fAbcModel_real_err_abs = TMath::Sqrt(  TMath::Power(fAbcModel_real->GetParError(0)/npe2, 2) + TMath::Power(fAbcModel_real->GetParError(1),2) +  TMath::Power(fAbcModel_real->GetParError(2)/npe2/npe2, 2)  ) / TMath::Sqrt(fAbcModel_real->Eval(npe2));
        //gExtra->SetPoint(i, npe2, fAbcModel_real->Eval(npe2)-fAbcModel_center->Eval(npe2));
        //gExtra->SetPointError(i, 0, TMath::Sqrt(fAbcModel_center_err_abs*fAbcModel_center_err_abs + fAbcModel_real_err_abs*fAbcModel_real_err_abs));
    }

    TFile* of = new TFile("e1MeV_pred_bterm_smear0mm.root", "recreate");
    gExtra->SetName("pred");
    gExtra->Write();
    of->Close();

    //TCanvas* c1 = new TCanvas(); c1->cd();
    //gExtra->SetLineWidth(2);
    //gExtra->SetLineColor(45);
    ////gExtra->SetTitle("extra resolution term; NPE; absolute increase value");
    //gExtra->SetTitle("extra resolution term; NPE;  extra term");
    ////gExtra->GetYaxis()->SetRangeUser(1e-5, 1e-4);
    //gExtra->SetFillColor(48);
    //gExtra->SetFillStyle(3002);
    //gExtra->Draw("a3 L");
    //gExtra->Draw("AL");
    //mg1->Draw("APL");
    //led1->Draw("SAME");
}



void draw_NU(double *f) {
    auto mg = new TMultiGraph();
    auto led = new TLegend();
    TGraph* gNU[5];
    Color_t color1[5] = {36, 38, 40, 42, 44};
    TString label2[5] = {"scale factor: 1.0", "scale factor: 2.0", "scale factor: 3.0", "scale factor: 4.0", "scale factor: 5.0"};
    for(int i=0; i<5; i++) {
        gNU[i] = new TGraph();
        for(int j=0; j<200; j++) {
            gNU[i]->SetPoint(j, j, (rsd_ratio[j]-rsd_ratio[0])/rsd_ratio[0]*f[i]+1);
        }
        gNU[i]->SetMarkerColor(color1[i]);
        gNU[i]->SetMarkerStyle(24);
        gNU[i]->SetMarkerSize(0.5);
        gNU[i]->SetLineColor(color1[i]);
        gNU[i]->SetLineWidth(2);
        mg->Add(gNU[i]);
        led->AddEntry(gNU[i], label2[i], "LP");
    }

    TCanvas* cc = new TCanvas(); cc->cd();
    mg->SetTitle("Residual Non-Uniformity; sub-detecotr id; ratio");
    mg->Draw("APL");
    led->Draw("SAME");
    cc->SaveAs("non-uniform-scale.pdf");
}


void read_SecCorr()
{
    TFile* file = new TFile("./subDet-corr.root", "read");
    file->GetObject("average", gSecCorr);
    for(int i=0; i<gSecCorr->GetN(); i++) {
        rsd_ratio_sec[i] = gSecCorr->GetPointY(i);
    }
}


