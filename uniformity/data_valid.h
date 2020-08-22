#ifndef data_valid_h
#define data_valid_h

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH1D.h"

class data_valid{
    public:
        data_valid() {}
        ~data_valid() {}

    public:

        void read_data(int N, int* par) {
            m_histN = N;
            Float_t m_evis; TBranch* b_evis;
            Float_t m_edepZ; TBranch* b_edepZ;
            Float_t m_edepR; TBranch* b_edepR;
            for(int i=0; i<m_histN; i++) {
                string hist_name = "evis"+to_string(par[i])+"MeV";
                pe_hist[i] = new TH1D(hist_name.c_str(), "", 2400, 0, 12);
                pe_hist[i]->SetName(hist_name.c_str());
                string filename = "/junofs/users/yumiao/simulation/energy_model/uniformity/electron/pe"+to_string(par[i])+"MeV.root";
                TFile* file = new TFile(filename.c_str(), "read");
                if(!file) continue;
                cout << "Processing " << filename << endl;
                TTree* tree = (TTree*)file->Get("petree");
                tree->SetBranchAddress("evis", &m_evis, &b_evis);
                tree->SetBranchAddress("edepR", &m_edepR, &b_edepR);
                tree->SetBranchAddress("edepZ", &m_edepZ, &b_edepZ);
                for(int j=0; j<tree->GetEntries(); j++) {
                    tree->GetEntry(j);
                    if(m_edepR > 17200) continue;
                    double delta_R = 18*18*18./20.;
                    int ir = int( TMath::Power(m_edepR/1000, 3) / delta_R);
                    double cos_theta = m_edepZ/m_edepR;
                    double delta_cos_theta = 2./10.;
                    int itheta = int((cos_theta+1)/delta_cos_theta);
                    int ibin = ir*10+itheta;
                    //m_evis /= corr[ibin];  // secondary correction
                    pe_hist[i]->Fill(m_evis);
                    //if(m_edepR<15000) {hTotPE_FV->Fill(m_evis);}
                }
            }
        }


        TGraphErrors* fit_data() {
            TGraphErrors* ge = new TGraphErrors();
            for(int i=0; i<m_histN; i++) {
                pe_hist[i]->Fit("gaus", "Q0");
                TF1* f1 = (TF1*)pe_hist[i]->GetFunction("gaus");
                double mu1 = f1->GetParameter(1);
                double mu1_err = f1->GetParError(1);
                double sigma1 = f1->GetParameter(2);
                double sigma1_err = f1->GetParError(2);
                double rsl1 = sigma1/mu1;
                double rsl1_err = TMath::Sqrt(sigma1_err*sigma1_err/mu1/mu1 + mu1_err*mu1_err*sigma1*sigma1/mu1/mu1/mu1/mu1);
                
                ge->SetPoint(i, mu1, rsl1);
                ge->SetPointError(i, mu1_err, rsl1_err);
            }
            return ge;
        }

    private:
        int m_histN;
        TH1D* pe_hist[20]; // maximum 20 energies

};

#endif
