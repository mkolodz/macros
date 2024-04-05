// The macro takes the original files and copies the first nev events
// to a new file, identical tree structure
// It also sorts the tree by ascending time of events

Bool_t subsetTOFPETTree(int run=544){
    gROOT->Reset();

    gStyle->SetOptStat(1111111);
    //access input file
//     TString finname = Form("run00%i_single.root",run);
    TString finname = Form("run00%i_single.root",run);
//     TString path = "/scratch1/gccb/data/TOFPET2/root/";
    TString path = "/scratch1/gccb/data/Jan2023Beam/root/";
    TString fname = path+finname;
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
    if (!f) 
        f = new TFile(fname,"READ");
    TTree* events = (TTree*)f->Get("events");
    TH1D * hTimeDiff = new TH1D("hTimeDiff", "hTimeDiff",200, -50, 50);
    TH1D * hTimeDiffSorted = new TH1D("hTimeDiffSorted", "hTimeDiffSorted",200, -50, 50);
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
    Int_t events_per_batch = 100000000;
    Int_t batchID = 0;
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
    std::cout << "nentries/events_per_batch" << nentries/events_per_batch << std::endl;
    Int_t n_batches = nentries/events_per_batch;
//     if(nev==-1){
//         nev=nentries;
//     }
//     if(nev>nentries){
//         cout<<"You want more events than there are in the tree... "<<nev<<"\t"<<nentries<<endl;
//         return kFALSE;
//     }
//     else {
//         cout<<"Selecting first "<<nev<<" events from the tree of "<<nentries<<endl;
//     }
    
    //create new tree
    TString foutname = finname;
    foutname.ReplaceAll("run",Form("run_b%i_", batchID));
    TFile *fout = new TFile(foutname,"RECREATE");
    auto newtree = events->CloneTree(0);

    TStopwatch timer;
    timer.Start();
    //loop over events
//     for (Long64_t i=100000000; i<nentries; i++) {
    for(int b=0;b<n_batches; b++){
            //sort events
//     events->BuildIndex("0","time");
//     TTreeIndex *I=(TTreeIndex*)events->GetTreeIndex(); // get the tree index
//     Long64_t* index=I->GetIndex(); //create an array of entries in sorted order

        for (Long64_t i=0; i<nentries; i++) {
            prev_time=time;
    //         events->GetEntry(index[i]);
            events->GetEntry(i);
            
            newtree->Fill();
            }
    }
    
    newtree->Print();
    timer.Stop();
    timer.Print();
    fout->Write();
    fout->Close();

    return kTRUE;
}
