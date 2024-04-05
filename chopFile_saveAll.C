// The macro takes the original file and splits it into chunks of
// size "chunk" saved as separate files in current directory
// resulting trees are time-sorted

#include "TH1F.h"

using namespace std;
const Long64_t chunk = 25e6;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
Int_t chopFile(TString inpath = "/scratch1/gccb/data/Jan2023Beam/root/",
	       TString finname="run00506_single_nlc.root"){
   gROOT->Reset();
   auto pwd = gDirectory;

   //access input file
   TString fname = inpath+finname;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
   if (!f) 
     f = new TFile(fname,"READ");
   if (!f){
     cout<<"File "<<inpath+fname<<" does not exist"<<endl;
     return kFALSE;
   }
   TTree* events = (TTree*)f->Get("events");
   if (!events){
     cout<<"File "<<inpath+fname<<" does not contain a tree called events"<<endl;
     return kFALSE;
   }
   Long64_t nentries = events->GetEntries();
   const int ntrees = nentries/chunk +1;

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

   TStopwatch timer;
   timer.Start();
   //loop over events
   cout<<"The tree contains "<<nentries<<" entries, that will be split into "<<ntrees<<" chunks"<<endl;

   Long64_t i = 0;
   Long64_t nev = 0;
   for (Int_t n=0; n<ntrees; n++) {
	//for (Int_t n=0; n<2; n++) { // for testing purposes
     cout<<"\n\n*** Chunk #"<<n<<endl;
     TString outfname = inpath + "/chopping/" + finname;
     outfname.ReplaceAll(".root",Form("_all_chunk%03d.root",n+1));
     TFile* fout = new TFile(outfname,"RECREATE");
     auto chunktree = events->CloneTree(0);
     auto tmptree = events->CloneTree(0);
     nev = 0;
     while(nev<chunk){
       i++;
       if(i==nentries) break;
       if(i%25000000==0)
          cout<<i<<" events processed..."<<endl;
       events->GetEntry(i);
       tmptree->Fill();
       nev++;
     }
     //sorting output trees:
     tmptree->BuildIndex("0","time");
     TTreeIndex *I=(TTreeIndex*)tmptree->GetTreeIndex();
     Long64_t* index=I->GetIndex(); 
     Long64_t tmpnentries=tmptree->GetEntries();
     for (Long64_t i=0; i<tmpnentries; i++) {
       tmptree->GetEntry(index[i]);
       chunktree->Fill();
     }
     delete tmptree;
     fout->Write();
     fout->Close();
   }
   timer.Stop();
   timer.Print();
   return ntrees;
}
