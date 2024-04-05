// {
//   TH1F *h_narrow = new TH1F("h_narrow", "h_narrow", 200, -1., 1.);
//   h_narrow->FillRandom("gaus");
//   make sure that the new histogram has "compatible" bins' widths
//   TH1F *h_wide = new TH1F("h_wide", "h_wide", 1000, -5., 5.);
//   h_wide->Sumw2(kTRUE); // ensure that the errors are properly "transfered"
//   h_wide->Add(h_narrow);
//   h_wide->ResetStats(); // may be needed sometimes
//   h_wide->Draw();
//   h_narrow->SetLineColor(kRed);
//   h_narrow->Draw("SAME");
// }
int proba1(){
TH1F *h = new TH1F("h", "h", 200, -1., 1.); //your original histogram 
h->FillRandom("gaus");
int nbins = h->GetXaxis()->GetNbins(); 
TH1F *hnew = new TH1F("hnew","title",nbins,-3,5); 
for (int i=1;i<=nbins;i++) { 
    double y = h->GetBinContent(i); 
    double x = h->GetXaxis()->GetBinCenter(i); 
    double xnew = 2*x + 1; //your transformation 
    hnew->Fill(xnew,y); 
    hnew->Draw();
    h->SetLineColor(kRed);
    h->Draw("SAME");
}
return 0;
}
