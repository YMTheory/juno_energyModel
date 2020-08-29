void hist2d_text()
{

   auto c01 = new TCanvas("c01","c01",700,400);
   c01->Divide(2,1);
   auto htext1 = new TH1F("htext1","Option TEXT on 1D histograms ",10,-4,4);
   auto htext2 = new TH2F("htext2","Option TEXT on 2D histograms ",10,-4,4,10,-20,20);
   float px, py;
   for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      htext1->Fill(px,0.1);
      htext2->Fill(px,5*py,0.1);
   }
   gStyle->SetPaintTextFormat("4.1f m");
   htext2->SetMarkerSize(1.8);
   c01->cd(1);
   htext2->Draw("TEXT45");
   c01->cd(2);
   htext1->Draw();
   htext1->Draw("HIST TEXT0 SAME");


}
