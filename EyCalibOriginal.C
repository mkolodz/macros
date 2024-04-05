void EyCalib(int run=524){
   gROOT->Reset();
   gStyle->SetOptStat(1111111);
   TString fname = Form("/scratch1/gccb/data/TOFPET2/root/run00%i_single.root",run);
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
   if (!f) 
     f = new TFile(fname,"READ");
   
   TTree* events = (TTree*)f->Get("events");
   
//Declaration of leaves types
   Float_t         step1;
   Float_t         step2;
   Long64_t        time;
   UInt_t          channelID;
   Float_t         tot;
   Float_t         energy;
   UShort_t        tacID;
   Int_t           xi;
   Int_t           yi;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         tqT;
   Float_t         tqE;

//    const int sampleChID = 262403;
//    const int RFchannel = 4128+2*131072;
   const int RFchannel = 4128;
   Long64_t lastRFtime = 0;
   Long64_t lasttime = 0;
   Long64_t deltat = 0;
   Long64_t deltatref = 0;

   const int nfib=4;
   
   int channel[nfib][2];
   channel[0][0]=60;
   channel[0][1]=323;
   channel[1][0]=21;
   channel[1][1]=373;
   channel[2][0]=109;
   channel[2][1]=286;
   channel[3][0]=34;
   channel[3][1]=350;
   
//    for(int i=0; i<nfib; i++){
//      channel[i][0]+=(2*131072);
//      channel[i][1]+=(2*131072);
//    }
   
   TString htitle[nfib]={"M0L2F17","M0L0F18","M0L1F4","M0L3F26"};
   TString side[2]={"left","right"};
   TH1F* hEnergy[nfib][2];
   TH1F* hHitTimeInterval[nfib][2];
   for(int f=0; f<nfib; f++){
       hEnergy[f][0] = new TH1F(Form("hEnergy_%i_00",f),Form("%s_left",htitle[f].Data()),90,0,30);
       hEnergy[f][0]->SetLineColor(kMagenta);
       hEnergy[f][0]->SetFillColor(kMagenta);
       hEnergy[f][1] = new TH1F(Form("hEnergy_%i_01",f),Form("%s_right",htitle[f].Data()),90,0,30);
       hEnergy[f][1]->SetLineColor(kBlack);
       hEnergy[f][1]->SetLineWidth(2);
   }
   for(int f=0; f<nfib; f++){
       hHitTimeInterval[f][0] = new TH1F(Form("hHitTimeInterval_%i_00",f),Form("%s_left Time from last ref-hit;#Delta T [ns];",htitle[f].Data()),4000,-2000,2000);
       hHitTimeInterval[f][0]->SetLineColor(kMagenta);
       hHitTimeInterval[f][0]->SetLineWidth(2);
       hHitTimeInterval[f][1] = new TH1F(Form("hHitTimeInterval_%i_01",f),Form("%s_right Time from last ref-hit;#Delta T [ns];",htitle[f].Data()),4000,-2000,2000);
       hHitTimeInterval[f][1]->SetLineColor(kBlack);
       hHitTimeInterval[f][1]->SetLineWidth(2);
   }

//    TH1F* hHitTimeInterval = new TH1F("hHitTimeInterval","Time from last ref-hit;#Delta T [ns];",
// 				     2000,-2000,0);
//    hHitTimeInterval->SetLineColor(kBlue);
//    hHitTimeInterval->SetFillColor(kBlue);
//    hHitTimeInterval->SetLineWidth(2);


       // Set branch addresses.
   events->SetBranchAddress("step1",&step1);
   events->SetBranchAddress("step2",&step2);
   events->SetBranchAddress("time",&time);
   events->SetBranchAddress("channelID",&channelID);
   events->SetBranchAddress("tot",&tot);
   events->SetBranchAddress("energy",&energy);
   events->SetBranchAddress("tacID",&tacID);
   events->SetBranchAddress("xi",&xi);
   events->SetBranchAddress("yi",&yi);
   events->SetBranchAddress("x",&x);
   events->SetBranchAddress("y",&y);
   events->SetBranchAddress("z",&z);
   events->SetBranchAddress("tqT",&tqT);
   events->SetBranchAddress("tqE",&tqE);

   int nhitAfterRF = 0;
   
   Long64_t nentries = events->GetEntries();
   Long64_t prev_time=0;
   Long64_t curr_time=0;
   int t_diff=0;
   int negative_t_diff_counter=0;
   int positive_t_diff_counter=0;
//    events->BuildIndex("0","time");
//    TTreeIndex *I=(TTreeIndex*)events->GetTreeIndex(); // get the tree index
//    Long64_t* index=I->GetIndex(); //create an array of entries in sorted order

   std::cout <<"entering loop" << std::endl;
   std::cout << nentries << std::endl;
//       nentries=100000000;
   for (Long64_t i=0; i<nentries; i++) {

     //     if(!i%(nentries/10))
     //cout<<i*10<<"% processed"<<endl; 
     //for (Long64_t i=0; i<1e6; i++) {
//     events->GetEntry(nentries-index[i]);
       events->GetEntry(nentries-i);
//     events->GetEntry(i);
    curr_time=time;
    t_diff = curr_time-prev_time;
    prev_time=curr_time;
    if(t_diff<0) {
        negative_t_diff_counter++;
//         continue;
    }
    else {
        positive_t_diff_counter++;
        continue;
        
    }
//      events->GetEntry(index[i]);
//      std::cout << "chan: " << channelID << std::endl;
     if(channelID== RFchannel){
       nhitAfterRF = 0;
       lastRFtime=time;
     }
     else nhitAfterRF++;
     
     if(nhitAfterRF!=1) continue;
     for(int f=0; f<nfib; f++){
       for(int s=0; s<nfib; s++){
	 if(channelID == channel[f][s]){
	   deltat = (time-lastRFtime)/1e3;
	   hHitTimeInterval[f][s]->Fill(deltat);
	   if(deltat>-420 && deltat<-386)
	     hEnergy[f][s]->Fill(energy);
            }     
       }
     }
     
   }
   
   TCanvas* canTime = new TCanvas("canTime","canTime",2400,1200);
   canTime->DivideSquare(nfib);
   TCanvas* can = new TCanvas("can","can",2400,1200);
   can->DivideSquare(nfib);
   for(int i=0; i<nfib; i++){
     can->cd(i+1);
     if(hEnergy[i][0]->GetMaximum()>hEnergy[i][1]->GetMaximum()){
       hEnergy[i][0]->Draw("hist");
       hEnergy[i][1]->Draw("hist,same");
     }
     else{
       hEnergy[i][1]->Draw("hist");
       hEnergy[i][0]->Draw("hist,same");
     }
     gPad->SetGrid(1,0);
   }
   gPad->BuildLegend();
   
   can->cd(1);
   TLatex latex;
   latex.SetTextColor(kRed);
   latex.DrawLatexNDC(0.6, 0.8, Form("Run %i",run));
   


//    canTime->cd();
      for(int i=0; i<nfib; i++){
     canTime->cd(i+1);
     if(hHitTimeInterval[i][0]->GetMaximum()>hHitTimeInterval[i][1]->GetMaximum()){
       hHitTimeInterval[i][0]->Draw("hist");
       hHitTimeInterval[i][1]->Draw("hist,same");
     }
     else{
       hHitTimeInterval[i][1]->Draw("hist");
       hHitTimeInterval[i][0]->Draw("hist,same");
     }
     gPad->SetGrid(1,0);
   }
   gPad->BuildLegend();

   canTime->SaveAs(Form("EyCal_%i.root",run));
   canTime->SaveAs(Form("EyCal_%i.png",run));
  
   std::cout << "positive_t_diff_counter: " << positive_t_diff_counter << std::endl;
   std::cout << "negative_t_diff_counter: " << negative_t_diff_counter << std::endl;
}
