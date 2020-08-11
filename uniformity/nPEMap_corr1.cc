#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TStopwatch.h"

using namespace std;

bool onlyLpmt ;
const int lpmt_num = 17612;
const int spmt_num = 25600;
const int MaxBin = 1440;
const double MapThetaStep = TMath::Pi()/1440.;
const double spmt_r = 19434-50;
double lpmt_x[lpmt_num]; double lpmt_y[lpmt_num]; double lpmt_z[lpmt_num]; int lpmt_id[lpmt_num];
double spmt_x[spmt_num]; double spmt_y[spmt_num]; double spmt_z[spmt_num]; int spmt_id[spmt_num];
void load_pmt();
TH2D* LMu2D[1440];
TH2D* SMu2D[1440];
void load_nPEMap();
double pred_nPE_Lpmt(double R_source, double theta_source, double theta_pmt);
double pred_nPE_Spmt(double R_source, double theta_source, double theta_pmt);


int main(int argc, char* argv[])
{
    // if only use lpmt?
    onlyLpmt = false;

    if(argc!=3) {cout << "Wrong Parameters Number! Expected 2!" << endl; return 0;} 

    double r3_max = 17.7*17.7*17.7; //15.5*15.5*15.5; //18*18*18;
    const int nbin_r = 40;
    TH1D* hTotPE[nbin_r] ;
    for(int i=0; i<nbin_r; i++) {
        TString Name = Form("totpe_r%d",i);  
        hTotPE[i] = new TH1D(Name, Name, 3000, 0, 12);
    }

    load_pmt();
    load_nPEMap();

    Int_t m_totpe;   TBranch* b_totpe;
    Float_t m_edep;  TBranch* b_edep;
    Float_t m_edepX; TBranch* b_edepX;
    Float_t m_edepY; TBranch* b_edepY;
    Float_t m_edepZ; TBranch* b_edepZ;
    Float_t m_evis; 
    Float_t m_evis_lpmt; 
    Float_t m_edepR;
    Float_t m_theta;
    Int_t m_totpe_lpmt;
    Float_t  m_time;
    Int_t m_npmt;    TBranch* b_npmt;
    Int_t m_pmtid_bypmt[200000]; TBranch* b_pmtid_bypmt;
    Int_t m_npe_bypmt[200000]; TBranch* b_npe_bypmt;
    
    TFile* f1 = new TFile(argv[2], "recreate");
    TTree* petree = new TTree("petree", "my tree");
    petree->Branch("edepX", &m_edepX, "edepX/F");
    petree->Branch("edepY", &m_edepY, "edepY/F");
    petree->Branch("edepZ", &m_edepZ, "edepZ/F");
    petree->Branch("edepR", &m_edepR, "edepR/F");
    petree->Branch("edep", &m_edep, "edep/F");
    petree->Branch("theta", &m_theta, "theta/F");
    petree->Branch("totpe", &m_totpe, "totpe/I");
    petree->Branch("evis", &m_evis, "evis/F");
    petree->Branch("totpe_lpmt", &m_totpe_lpmt, "totpe_lpmt/I");
    petree->Branch("evis_lpmt", &m_evis_lpmt, "evis_lpmt/F");
    petree->Branch("m_time", &m_time, "time/F");


    int startNo = atoi(argv[1]);
    int nFiles = 1;//10;
    for(int i=startNo; i<startNo+nFiles; i++) {
        //string no = to_string(i);
        stringstream ss;
        ss << i;
        string no = ss.str();
        string filename = "/junofs/users/yumiao/simulation/energy_model/production/electron/cerenkov/2000keV/sample_detsim_user.root";
        //string filename = "./uniform/8MeV/user-detsim-"+no+".root";
        cout << "Processing " << filename << endl;
        TFile* ff; TTree* evt;
        ff = TFile::Open(filename.c_str());
        evt = (TTree*)ff->Get("evt");
        evt->SetBranchAddress("totalPE", &m_totpe, &b_totpe);
        evt->SetBranchAddress("edep", &m_edep, &b_edep);
        evt->SetBranchAddress("edepX", &m_edepX, &b_edepX);
        evt->SetBranchAddress("edepY", &m_edepY, &b_edepY);
        evt->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
        evt->SetBranchAddress("nPMTs", &m_npmt, &b_npmt);
        evt->SetBranchAddress("PMTID_byPMT", m_pmtid_bypmt, &b_pmtid_bypmt);
        evt->SetBranchAddress("nPE_byPMT", m_npe_bypmt, &b_npe_bypmt);
        //cout << evt->GetEntries() <<endl;
        for(int iEntry=0; iEntry<evt->GetEntries(); iEntry++) {
            cout << "Processing " << iEntry << endl;
            evt->GetEntry(iEntry);
            TStopwatch timer;
            timer.Start();
            double R_source = TMath::Sqrt(m_edepX*m_edepX+m_edepY*m_edepY+m_edepZ*m_edepZ);
            if(R_source>17700) R_source = 17700-0.001;
            //if(R_source>15500) continue;   // 15.5m fidicial volume cut
            //int rbinid = int(R_source*R_source*R_source/1000000000/(r3_max/nbin_r)-0.00001);
            double theta_source = TMath::ACos(m_edepZ/R_source);
            m_edepR = R_source;
            m_theta = theta_source;

            // select lpmt only
            int totpe_lpmt = 0;
            for(int iPmt=0; iPmt<m_npmt; iPmt++) {
                if(m_pmtid_bypmt[iPmt] > lpmt_num) continue;
                else {
                    totpe_lpmt += m_npe_bypmt[iPmt];
                }
            } 

            // uniformity correction
            double totScale = 0;  double totScale_lpmt = 0;
            double totpe_corr = 0; double totpe_lpmt_corr = 0;
            for(int idx=0; idx<lpmt_num; idx++) {
                double R_pmt = TMath::Sqrt(lpmt_x[idx]*lpmt_x[idx]+lpmt_y[idx]*lpmt_y[idx]+lpmt_z[idx]*lpmt_z[idx]);
                double theta_pmt = TMath::ACos((m_edepX*lpmt_x[idx]+m_edepY*lpmt_y[idx]+m_edepZ*lpmt_z[idx])/R_source/R_pmt);
                //cout << R_source << " " << theta_source << " " << R_pmt << " " << theta_pmt << endl;
                double scale = pred_nPE_Lpmt(R_source, theta_source, theta_pmt);           
                totScale_lpmt += scale;
                totScale += scale;
            }
            for(int idx=0; idx<spmt_num; idx++) {
                double theta_pmt = TMath::ACos((m_edepX*spmt_x[idx]+m_edepY*spmt_y[idx]+m_edepZ*spmt_z[idx])/R_source/spmt_r);
                double scale = pred_nPE_Spmt(R_source, theta_source, theta_pmt);
                totScale += scale;
            }
            //totScale = 1.;   // no uniform correction
            //totScale = 1362;   // no uniform correction
    
            totpe_corr      = m_totpe/totScale*1.022;
            totpe_lpmt_corr = totpe_lpmt / totScale_lpmt*1.022;

            m_evis = totpe_corr;
            m_evis_lpmt = totpe_lpmt_corr;

            timer.Stop();
            m_time = timer.RealTime();

            petree->Fill();

            //if(totScale!=0 and onlyLpmt) {
            //    totpe_corr = totpe_lpmt / totScale * 1.022;
            //    h_totpe->Fill(totpe_lpmt);
            //    h_totpe_corr->Fill(totpe_corr);
            //    if(R_source<15500) h_totpe_corr_15m->Fill(totpe_corr);  // 15.5m FV
            //    hTotPE[rbinid]->Fill(totpe_corr);
            //}
            //if(totScale!=0 and !onlyLpmt) {
            //    totpe_corr = m_totpe / totScale * 1.022;
            //    h_totpe->Fill(m_totpe);
            //    h_totpe_corr->Fill(totpe_corr);
            //    if(R_source<15500) h_totpe_corr_15m->Fill(totpe_corr);  // 15.5m FV
            //    hTotPE[rbinid]->Fill(totpe_corr);
            //}

        }
        //ff->Close();
        evt->Delete();
        ff->Delete();
    }
    
    f1->cd();
    petree->Write();
    f1->Close();

    //TFile* f1 = new TFile(argv[2], "recreate");
    //h_totpe->Write();
    //h_totpe_corr->Write();
    //h_totpe_corr_15m->Write();
    //for(int i=0; i<nbin_r; i++) {
    //    hTotPE[i]->Write();
    //}
    //f1->Close();
    return 1;
}


void load_pmt()
{
    ifstream in; in.open("/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc830/Pre-Release/J20v1r0-Pre2/offline/Simulation/DetSimV2/DetSimOptions/data/PMTPos_Acrylic_with_chimney.csv");
    const int num = 17612; int index = 0;
    string line;
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> lpmt_id[index] >> lpmt_x[index] >> lpmt_y[index] >> lpmt_z[index] ;
        //pmt_theta[index] = 90-pmt_theta[index];
        //cout << pmt_theta[index] << " " << pmt_phi[index] <<endl;
        index++;
    }
    in.close();

    double theta; double phi; index = 0;
    in.open("/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc830/Pre-Release/J20v1r0-Pre2/offline/Simulation/DetSimV2/DetSimOptions/data/3inch_pos.csv");
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> spmt_id[index] >> theta >> phi;
        spmt_x[index] = spmt_r*TMath::Sin(theta/180.*TMath::Pi())*TMath::Cos(phi/180.*TMath::Pi());
        spmt_y[index] = spmt_r*TMath::Sin(theta/180.*TMath::Pi())*TMath::Sin(phi/180.*TMath::Pi());
        spmt_z[index] = spmt_r*TMath::Cos(theta/180.*TMath::Pi());
        //pmt_theta[index] = 90-pmt_theta[index];
        //cout << pmt_theta[index] << " " << pmt_phi[index] <<endl;
        index++;
    }
    in.close();
}

void load_nPEMap()
{
    TFile* lfile = TFile::Open("/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc830/Pre-Release/J20v1r0-Pre2/data/Reconstruction/OMILREC/RecMap/nPEMap/LnPEMapFile_Truth.root");
    for(int i=0; i<1440; i++) {
        LMu2D[i] = (TH2D*)lfile->Get(Form("hLMu2D_%d", i));
    }
    //lfile->Close();

    TFile* sfile = TFile::Open("/cvmfs/juno.ihep.ac.cn/sl6_amd64_gcc830/Pre-Release/J20v1r0-Pre2/data/Reconstruction/OMILREC/RecMap/nPEMap/SnPEMapFile_Truth.root");
    for(int i=0; i<1440; i++) {
        SMu2D[i] = (TH2D*)sfile->Get(Form("hSMu2D_%d", i));
    }
    //sfile->Close();
}


double pred_nPE_Lpmt(double R_source, double theta_source, double theta_pmt) {
    int Rbin = int(R_source/10);
    int theta_source_bin = int(theta_source/(TMath::Pi()/18.));

    int theta_pmt_bin1 = theta_pmt / MapThetaStep;
    int theta_pmt_bin2 = theta_pmt_bin1 + 1;

    if(theta_pmt_bin1==MaxBin) { theta_pmt_bin1 = MaxBin-1; theta_pmt_bin2 = theta_pmt_bin1;}
    else if (theta_pmt_bin2==MaxBin) {theta_pmt_bin2 = theta_pmt_bin1;}

    //cout << "LPMT: "<< theta_pmt_bin1 << " " << theta_pmt_bin2 << endl;
        
    double scale1 = LMu2D[theta_pmt_bin1]->Interpolate(R_source, theta_source);
    double scale2 = LMu2D[theta_pmt_bin2]->Interpolate(R_source, theta_source);
    double ThetaAFrac = theta_pmt - double(theta_pmt_bin1)*MapThetaStep;
    double Weight = ThetaAFrac / MapThetaStep;
    double average_scale = (1.-Weight)*scale1 + Weight*scale2;

    //double scale = LMu2D[theta_pmt_bin1]->GetBinContent(Rbin+1, theta_source_bin+1);
    //cout << theta_pmt_bin1 << " " << scale << " " << average_scale << endl;

    return average_scale;
}


double pred_nPE_Spmt(double R_source, double theta_source, double theta_pmt) {
    //int Rbin = int(R_source/10);
    //int theta_source_bin = int(theta_source/(TMath::Pi()/18.));
    ////cout << "theta_pmt_bin : " << theta_pmt_bin << endl;
    //int theta_pmt_bin = int(theta_pmt/(TMath::Pi()/1440.)-0.0001);
    //double scale = SMu2D[theta_pmt_bin]->GetBinContent(Rbin+1, theta_source_bin+1);

    int theta_pmt_bin1 = theta_pmt / MapThetaStep;
    int theta_pmt_bin2 = theta_pmt_bin1 + 1;

    if(theta_pmt_bin1==MaxBin) { theta_pmt_bin1 = MaxBin-1; theta_pmt_bin2 = theta_pmt_bin1;}
    else if (theta_pmt_bin2==MaxBin) {theta_pmt_bin2 = theta_pmt_bin1;}
    //cout << "SPMT: "<< theta_pmt_bin1 << " " << theta_pmt_bin2 << endl;

    double scale1 = SMu2D[theta_pmt_bin1]->Interpolate(R_source, theta_source);
    double scale2 = SMu2D[theta_pmt_bin2]->Interpolate(R_source, theta_source);
    double ThetaAFrac = theta_pmt - double(theta_pmt_bin1)*MapThetaStep;
    double Weight = ThetaAFrac / MapThetaStep;
    double average_scale = (1.-Weight)*scale1 + Weight*scale2;

    //double scale = LMu2D[theta_pmt_bin]->GetBinContent(Rbin+1, theta_source_bin+1);
    //cout << theta_pmt_bin << " " << scale << endl;

    return average_scale;
}
