void muPeak_analysis()
{

    TH1D* hist_arr[227];
    TFile* ff = TFile::Open("/junofs/users/huanggh/EnergyRec/GenCalibData/ACU_CLS_MAP/nPEMap_Truth/J20v1r1-Pre0/e+2MeVScaleCompFile.root");
    for(int i=0; i<227; i++) {
        TString Name=Form("hComp_pos%d",i);
        hist_arr[i] = (TH1D*)ff->Get(Name);
    }


    int idx = 0;
    TGraph* graph = new TGraph();
    for(int ipos=0; ipos<227; ipos++) {
        int peak_pos = -1;
        int peak_value = -1;
        for(int ibin=0; ibin<hist_arr[ipos]->GetNbinsX(); ibin++) {
            if(hist_arr[ipos]->GetBinContent(ibin+1)>peak_value and hist_arr[ipos]->GetBinContent(ibin+1)>1.2) {
                peak_value = hist_arr[ipos]->GetBinContent(ibin+1);
                peak_pos = hist_arr[ipos]->GetBinCenter(ibin+1);
            }
        }

        if(peak_value != -1) {
            graph->SetPoint(idx, peak_pos, peak_value); idx++;
        }
    }

    graph->SetMarkerStyle(24);
    graph->SetMarkerColor(44);
    graph->Draw("AP");



}
