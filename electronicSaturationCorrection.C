// The macro takes the original files, rescales the QDC value according to the equation provided in the TOFPET manual (electronic saturation/non-linearity correction) and copies the first nev events
// to a new file, identical tree structure 

Bool_t electronicSaturationCorrection(int run=491, Long64_t nev=-1){
   gROOT->Reset();

   //access input file
   TString finname = Form("run00%i_single.root",run);
   TString path = "/scratch1/gccb/data/Jan2023Beam/root/";
   TString fname = path+finname;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
   if (!f) 
     f = new TFile(fname,"READ");
   TTree* events = (TTree*)f->Get("events");
   
   Float_t p0, p1, p2, p3;
   p0=8.0;
   p1=1.04676;
   p2=1.02734;
   p3=0.31909;
   
   TH1F * hEnergy = new TH1F("hEnergy","hEnergy;QDC[a.u.];Counts",200,0,200);
   TH1F * hEnergyRescaled = new TH1F("hEnergyRescaled","hEnergyRescaled;QDC[a.u.];Counts",200,0,200);
   TH1F * hEnergyRescaledProcessed = new TH1F("hEnergyRescaledProcessed","hEnergyRescaledProcessed;QDC[a.u.];Counts",200,0,200);
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

   // Set branch addresses
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
    
   //check sanity of arguments
   Long64_t nentries = events->GetEntries();
   if(nev==-1){
    nev = nentries; //careful!
   }
   if(nev>nentries)
     { cout<<"You want more events than there are in the tree... "<<nev<<"\t"<<nentries<<endl;
       return kFALSE;
     }
   else {
     cout<<"Selecting first "<<nev<<" events from the tree of "<<nentries<<endl;
   }
   
   //create new tree
   TString foutname = finname;
   foutname.ReplaceAll("run","saturationCorrectedRun");
   TFile *fout = new TFile(foutname,"RECREATE");
   auto newtree = events->CloneTree(0);

   //sort events
   // events->BuildIndex("0","time");
   // TTreeIndex *I=(TTreeIndex*)events->GetTreeIndex(); // get the tree index
   // Long64_t* index=I->GetIndex(); //create an array of entries in sorted order

   TStopwatch timer;
   timer.Start();
   //loop over events
   for (Long64_t i=0; i<nev; i++) {
     //     events->GetEntry(index[i]);
     events->GetEntry(i);
     hEnergy->Fill(energy);
//      std::cout << "before rescaling: " << energy;
     energy = p0*pow(p1, pow(energy, p2))+p3*energy-p0;
    hEnergyRescaled->Fill(energy);
//      std::cout << "after rescaling: " << energy << std::endl;
     newtree->Fill();
   }
   newtree->Print();
   timer.Stop();
   timer.Print();
   fout->Write();
   fout->Close();
   
   TCanvas *c1= new TCanvas("c1", "c1", 300,300);
   c1->cd();
   hEnergy->Draw();
   hEnergyRescaled->SetLineColor(kRed);
   hEnergyRescaled->SetFillColor(kRed);
   hEnergyRescaled->Draw("same");
   TFile *f1 = new TFile("/scratch1/gccb/data/Jan2023Beam/root/run00491_single_nlc.root");
    if(!f1->IsOpen()) std::cerr << "file not open" << std::endl;
    TTree *t = (TTree*)f1->Get("events");
//     Long64_t time;
//     Float_t energy;
//     UInt_t channelID;
    t->Print();
    t->SetBranchAddress("time",&time);
    t->SetBranchAddress("channelID",&channelID);
    t->SetBranchAddress("energy",&energy); 

//    TFile *f1;
//     f1 = new TFile("/scratch1/gccb/data/Jan2023Beam/root/run00491_single_nlc.root","READ");
//    TTree* events1 = (TTree*)f1->Get("events");
   for (Long64_t i=0; i<nev; i++) {
     //     events->GetEntry(index[i]);
     t->GetEntry(i);
     hEnergyRescaledProcessed->Fill(energy);
     std::cout << energy << std::endl;
   }
      hEnergyRescaledProcessed->SetLineColor(kGreen);
      hEnergyRescaledProcessed->SetLineWidth(2);
    hEnergyRescaledProcessed->Draw("same");
   gPad->BuildLegend();
   
   return kTRUE;
}
