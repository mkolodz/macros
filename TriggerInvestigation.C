void TriggerInvestigation(int run=905){
   gROOT->Reset();
   gStyle->SetOptStat(111111);
   
   TString fname = Form("preselectedRun00%i_single.root",run);
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

   //   const int sampleChID = 262403;
//    const int RFchannel = 4128+2*131072;
   const int RFchannel = 4128;
   Long64_t lastRFtime = 0;
   Long64_t lasttime = 0;
   Long64_t deltat = 0;
   Long64_t deltatref = 0;

   const int nfib=4;
   
   int channel[nfib][2];
   channel[0][0]=122;
   channel[0][1]=284;
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
   
   TString htitle[nfib]={"M0L0S5","M0L0F18","M0L1F4","M0L3F26"};
   TString side[2]={"left","right"};
   TH1F* hEnergy[nfib][2];
   TH1F* hEnergyAll[nfib][2];
   for(int f=0; f<nfib; f++){
       hEnergy[f][0] = new TH1F(Form("hEnergy_%i_00",f),
				Form("%s_left",htitle[f].Data()),90,0,30);
       hEnergy[f][0]->SetLineColor(kMagenta);
       hEnergy[f][0]->SetLineWidth(2);
//        hEnergy[f][0]->SetFillColorAlpha(kMagenta,0.35);
       hEnergy[f][1] = new TH1F(Form("hEnergy_%i_01",f),
				Form("%s_right",htitle[f].Data()),90,0,30);
       hEnergy[f][1]->SetLineColorAlpha(kBlack,0.35);
       hEnergy[f][1]->SetLineWidth(2);

       hEnergyAll[f][0]=(TH1F*)hEnergy[f][0]->Clone();
       hEnergyAll[f][1]=(TH1F*)hEnergy[f][1]->Clone();
   }

   TH1F* hHitTimeInterval = new TH1F("hHitTimeInterval",
				     "Time from last ref-hit;#Delta T [ns];",
				     100000,-800000,200000);
   hHitTimeInterval->SetLineColor(kBlue);
   hHitTimeInterval->SetFillColorAlpha(kBlue,0.35);
   hHitTimeInterval->SetLineWidth(2);

   TH1F* hHitTimeIntervalFirst = (TH1F*)hHitTimeInterval->Clone("hHitTimeIntervalFirst");
   hHitTimeIntervalFirst->SetLineColor(kRed);
   hHitTimeIntervalFirst->SetFillColorAlpha(kRed,0.35);


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
   
   TStopwatch timer;
   timer.Start();
   events->BuildIndex("0","time");
   TTreeIndex *I=(TTreeIndex*)events->GetTreeIndex(); // get the tree index
   Long64_t* index=I->GetIndex(); //create an array of entries in sorted order
   timer.Stop();
   cout<<"Sorting took "<<endl;
   timer.Print();
   timer.Reset();

   timer.Start();
   for (Long64_t i=0; i<nentries; i++) {
     //     if(!i%(nentries/10))
     //cout<<i*10<<"% processed"<<endl; 
     //for (Long64_t i=0; i<1e6; i++) {
     events->GetEntry(index[nentries-i]);
     
     if(channelID== RFchannel){
       nhitAfterRF = 0;
       lastRFtime=time;
     }
     else nhitAfterRF++;
     
     // if(nhitAfterRF!=1) continue;
     for(int f=0; f<nfib; f++){
       for(int s=0; s<nfib; s++){
	 if(channelID == channel[f][s]){
	   deltat = (time-lastRFtime)/1e3;
	   hHitTimeInterval->Fill(deltat);
	   hEnergyAll[f][s]->Fill(energy);
	   if(nhitAfterRF==1)  hHitTimeIntervalFirst->Fill(deltat);
	   if(deltat>-420 && deltat<-386)
//        if(deltat>-480 && deltat<-400)
	     hEnergy[f][s]->Fill(energy);
	 }     
       }
     }
   }
   cout<<"Loooping over events took "<<endl;
   timer.Print();
   timer.Reset();
   
   TCanvas* can = new TCanvas("can","can",1200,600);
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

   TCanvas* canAll = new TCanvas("canAll","canAll",1200,600);
   canAll->DivideSquare(nfib);
   for(int i=0; i<nfib; i++){
     canAll->cd(i+1);
     if(hEnergyAll[i][0]->GetMaximum()>hEnergyAll[i][1]->GetMaximum()){
       hEnergyAll[i][0]->Draw("hist");
       hEnergyAll[i][1]->Draw("hist,same");
     }
     else{
       hEnergyAll[i][1]->Draw("hist");
       hEnergyAll[i][0]->Draw("hist,same");
     }
     gPad->SetGrid(1,0);
   }
   gPad->BuildLegend();
   
   can->cd(1);
   TLatex latex;
   latex.SetTextColor(kRed);
   latex.DrawLatexNDC(0.6, 0.8, Form("Run %i",run));
   
   TCanvas* can2 = new TCanvas("can2","can2",1200,600);
   gPad->SetLogy();
   hHitTimeInterval->Draw();
   hHitTimeIntervalFirst->Draw("same");
   

   can->SaveAs(Form("EyCal_%i.root",run));
   can->SaveAs(Form("EyCal_%i.png",run));
   canAll->SaveAs(Form("EyAll_%i.root",run));
   canAll->SaveAs(Form("EyAll_%i.png",run));
   can2->SaveAs(Form("TDiff_%i.root",run));
   can2->SaveAs(Form("TDiff_%i.png",run));
}
