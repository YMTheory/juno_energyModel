#include "data_valid.h"

void datad()
{
    data_valid dv;
    int N = 5;
    int mom[5] = {1,2,4,6,8};
    double corr = 1;
    dv.read_data(N, mom);
    TGraphErrors* ge  = dv.fit_data();

    TF1* func = new TF1("func", "TMath::Sqrt([0]/x+[1]+[2]/x/x)", 0, 12);
    ge->Fit(func);
    ge->Draw("APL");
}
