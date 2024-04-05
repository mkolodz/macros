const int nf = 7;
const double peakratio = 1275./511;
TString fnames[nf] = {"/scratch1/gccb/data/TOFPET2/root/raw00454_single.root",
		      "/scratch1/gccb/data/TOFPET2/root/raw00455_single.root", 
		      "/scratch1/gccb/data/TOFPET2/root/raw00456_single.root", 
		      "/scratch1/gccb/data/TOFPET2/root/raw00457_single.root", 
		      "/scratch1/gccb/data/TOFPET2/root/raw00458_single.root", 
		      "/scratch1/gccb/data/TOFPET2/root/raw00459_single.root", 
		      "/scratch1/gccb/data/TOFPET2/root/raw00460_single.root"};

Int_t ov[nf]={14, 12, 10, 8, 6, 4, 2};
Double_t mu511[nf] = {25.6, 19.6, 15.2, 11.7, 8.5, 5.56, 2.45}; //change here!
// Double_t mu511[nf] = {29.5, 19.42, 13.7, 10.0, 6.86, 4.25, 1.81}; 
Double_t nlc(Double_t x){
Float_t p0, p1, p2, p3;
p0=8.0;
p1=1.04676;
p2=1.02734;
p3=0.31909;

// auto f_nlc = new TF1("f_nlc","[0]*pow([1], pow(x, [2]))+[3]*x-[0]",0,500);
// f_nlc->SetParameter(0,p0);
// f_nlc->SetParameter(1,p1);
// f_nlc->SetParameter(2,p2);
// f_nlc->SetParameter(3,p3);
// f_nlc->Draw();
return (p0*pow(p1, pow(x, p2))+p3*x-p0);
}  
  
TH1F* GetHisto(double ulimit, TString fname, int ov){
  TH1F* h= new TH1F(  Form("ov%d",ov),   Form("OV = %d V",ov), 200, 0, ulimit);
  h->GetXaxis()->SetTitle("qdc [ch]");
  h->GetXaxis()->SetLabelSize(0.1);
  h->GetYaxis()->SetTitle("counts");
  TFile* file = new TFile(fname, "READ");
  TTree *t = (TTree *)file->Get("events");
  Float_t energy;
  UInt_t channelID;
  t->SetBranchAddress("energy", &energy);
  t->SetBranchAddress("channelID", &channelID);
  //int N=t->GetEntries();;
  int N=1e7;
  for(Int_t i=0; i < N; ++i) {
    t->GetEntry(i);
    UShort_t ch = channelID;
    if(channelID >= 256) ch = channelID - 128;
    if(ch==139)
      h->Fill(energy);
//         h->Fill(nlc(energy));
  }
  h->SetMaximum(0.0006*N);
  
  file->Close();
  return h;
}



TCanvas* DrawAll(void){
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.15);
  TH1F* h[nf];
  TLine line;
  line.SetLineColor(kRed);
  TLatex tex;
  tex.SetTextColor(kRed);
  tex.SetTextSize(0.1);
  TCanvas* can = new TCanvas("can","can",500, 1500);
  can->Divide(1,nf);
  for(int i=0; i<nf; i++){
    h[i] = GetHisto(mu511[i]*3, fnames[i], ov[i]);
    can->cd(i+1);
    gPad->SetLogy(1);
    TF1* f = new TF1("f","gaus(0)",0.8*mu511[i], 1.2*mu511[i]);
    f->SetParameters(300,mu511[i],4);
    h[i]->Fit("f","R");
    h[i]->DrawClone("");
    double pp = h[i]->GetFunction("f")->GetParameter(1);
    double sigma = h[i]->GetFunction("f")->GetParameter(2);
    line.DrawLine(peakratio*pp,0,peakratio*pp,h[i]->GetMaximum());
    tex.DrawLatex(3.5*pp, h[i]->GetMaximum()*0.8, Form("#sigma_{511}/#mu_{511}=%.2f per cent",sigma/pp*100));
  }
  return can;
}
  
  

