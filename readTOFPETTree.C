#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
// #include "SiFiSeries.h"
// #include "SiFiDatabase.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include <RtypesCore.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <ctime>

int readTOFPETTree(){

    TFile *f; 
    f = new TFile("/scratch1/gccb/data/TOFPET2/root/run00423_single.root");
    if(!f->IsOpen()) std::cerr << "file not open" << std::endl;
    TTree *t = (TTree*)f->Get("data");
    Long64_t time,ltime;
    Float_t energy, lenergy;
    UInt_t channelID, lchannelID;
    t->Print();
    t->SetBranchAddress("time",&ltime);
    t->SetBranchAddress("channelID",&lchannelID);
    t->SetBranchAddress("energy",&lenergy); 

    TH1F *hTime   = new TH1F("hTime","hTime",100,-100,100);
    TH1F *hChannelID   = new TH1F("hChannelID","hChannelID",100,0,500);
    TH1F *hEnergy   = new TH1F("hEnergy","hEnergy",100,-100,100);

    Long64_t nentries = t->GetEntries();
//     for (Long64_t i=0;i<nentries;i++) {
//         t->GetEntry(i);
//         hTime->Fill(time);
//         hChannelID->Fill(channelID);
//         hEnergy->Fill(energy);
//     }
  
  TFile file_small("sifi_results_small_sample.root","recreate");
   TTree *data = new TTree("data","a friend Tree");
   data->Branch("time",&time);
    data->Branch("channelID",&channelID);
    data->Branch("energy",&energy); 

      for (Long64_t i=0;i<10000;i++) {
        t->GetEntry(i);
        time=ltime;
        energy=lenergy;
        channelID=lchannelID;
        data->Fill();
      }
      
         file_small.cd();
   data->Write();
//         hTime->Fill(time);
//         hChannelID->Fill(channelID);
//         hEnergy->Fill(energy);
//    
   
//     TCanvas * c = new TCanvas("c", "c", 900, 300);
//     c->Divide(3,1);
//     c->cd(1);
//     hTime->Draw();
//         c->cd(2);
//     hChannelID->Draw();
//         c->cd(3);
//     hEnergy->Draw();
return 0;
}
