#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TObjString.h"
#include "TKey.h"
#include "TMemFile.h"

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
#include <algorithm>

//----- scatterer geometry 4to1 (here fibers are actually SiPMs)
#define N_MODULES 1
#define N_LAYERS_PER_MODULE 4
#define N_FIBERS_PER_LAYER 28
#define N_SIDES 2

//----- canvas size
#define XCANVAS 1600
#define YCANVAS 1000

// //TOFPET 
#define NBINS_QDC 100
#define XLOW_QDC 0
#define XUP_QDC 35

#define NBINS_T 1000
#define XLOW_T 0
#define XUP_T 1000

#define NBINS_TDIFF 100
#define XLOW_TDIFF 0
#define XUP_TDIFF 10000
#define XUP_TDIFF_ZOOM 100

#define NBINS_MLR 100
#define XLOW_MLR -15
#define XUP_MLR 15

//----- cuts
#define TDIFF_MAX 10
#define QDC_MIN 0
#define T_MIN 0

#define max_n_files 100
int readTofpetTrees(std::vector<std::string> path){
    
    auto start = std::chrono::system_clock::now();
    gStyle->SetPalette(kBird);   
//     gStyle->SetOptStat(0);
    TFile * input_file[max_n_files]={nullptr};    ///< data input file
    TTree *t;
    Long64_t nentries;
    Long64_t time=-100;
    Float_t energy=-100;
    UInt_t channelID=0;
    size_t n_files;
    n_files = path.size();
//     const double deltaT = 10; // 10 ns window
    const double ps_to_ns = 1E3;
    TH1D *hQ[max_n_files*2] = {nullptr}; //for left and right, that's why multiplied by 2
//     hQ[0]=new TH1D("ch_m0l1elem8l_10mm","ch_m0l1elem8l; QDC[a.u.]; Counts", NBINS_QDC, XLOW_QDC, XUP_QDC);
//     hQ[1]=new TH1D("ch_m0l1elem8r_10mm","ch_m0l1elem8r; QDC[a.u.]; Counts", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[0]=new TH1D("ch_m0l0elem5l_10mm","ch_m0l0elem5l; QDC[a.u.]; Counts", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[1]=new TH1D("ch_m0l0elem5r_10mm","ch_m0l0elem5r; QDC[a.u.]; Counts", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[2]=new TH1D("ch115_90mm","ch115_90mm", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[3]=new TH1D("ch273_90mm","ch273_90mm", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[4]=new TH1D("ch115_50mm","ch115_50mm", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[5]=new TH1D("ch273_50mm","ch273_50mm", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[6]=new TH1D("ch115_bg","ch115_bg", NBINS_QDC, XLOW_QDC, XUP_QDC);
    hQ[7]=new TH1D("ch273_bg","ch273_bg", NBINS_QDC, XLOW_QDC, XUP_QDC);
    for(int k=0;k<n_files;k++){
        input_file[k] = new TFile(path[k].c_str());
    //     if(!input_file->IsOpen()) {
    //         std::cerr << "##### Error! Could not open input file!" << std::endl;
    //         std::cerr << path << std::endl;
    //         return false;
    //     }
        t = (TTree*)input_file[k]->Get("events");

        nentries=t->GetEntries();
    //     std::cout << "nentries " << nentries << std::endl;
        t->SetBranchAddress("time",&time);
        t->SetBranchAddress("channelID",&channelID);
        t->SetBranchAddress("energy",&energy);
    //         nentries=4000;
    //     int trg_counter = 0;
    //     int regular_channel_counter = 0;
        for(int i=0; i<nentries; i++){
            t->GetEntry(i);
            time /= ps_to_ns;        
    //         std::cout << channelID << " " << time << std::endl;
            if(channelID == 115) hQ[2*k]->Fill(energy);
//             if(channelID == 60) hQ[2*k+2]->Fill(energy);
//             if(channelID == 42) hQ[2*k+4]->Fill(energy);
//             if(channelID == 122) hQ[2*k+6]->Fill(energy);
            if(channelID == 273) hQ[2*k+1]->Fill(energy);
//             if(channelID == 323) hQ[2*k+3]->Fill(energy);
//             if(channelID == 339) hQ[2*k+5]->Fill(energy);
//             if(channelID == 284) hQ[2*k+7]->Fill(energy);
    //             regular_channel_counter++;
    //         if(channelID == 4128) trg_counter++;
            }
    }      
    TCanvas * c1 = new TCanvas("c1", "c1", 1200,1200);
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetGrid(1,1);
    hQ[0]->Draw();
    hQ[2]->SetLineColor(2);
    hQ[2]->Draw("same");
    hQ[4]->SetLineColor(3);
    hQ[4]->Draw("same");
    hQ[6]->SetLineColor(35);
    hQ[6]->Draw("same");
    auto legend = new TLegend(0.7,0.75,0.9,0.9); // option "C" allows to center the header
    legend->AddEntry(hQ[0],"10mm (closer to r)","l"); //p stands for "points"
    legend->AddEntry(hQ[2],"90mm (closer to l)","l");
//     legend->AddEntry(hQ[4],"50mm","l");
//     legend->AddEntry(hQ[6],"bg","l");
    legend->Draw();

    c1->cd(2);
    gPad->SetGrid(1,1);
    hQ[1]->Draw();
    hQ[3]->SetLineColor(2);
    hQ[3]->Draw("same");
    hQ[5]->SetLineColor(3);
    hQ[5]->Draw("same");
    hQ[7]->SetLineColor(35);
    hQ[7]->Draw("same");
    auto legend1 = new TLegend(0.7,0.75,0.9,0.9); // option "C" allows to center the header
    legend1->AddEntry(hQ[1],"10mm (closer to r)","l"); //p stands for "points"
    legend1->AddEntry(hQ[3],"90mm (closer to l)","l");
//     legend1->AddEntry(hQ[5],"50mm","l");
//     legend1->AddEntry(hQ[7],"bg","l");
    legend1->Draw();
    
//     std::cout << trg_counter << std::endl;    
    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> runtime = end-start;
    std::cout << "Time: " << runtime.count() << " s\n";
    std::cout << "Entries: " << nentries << "\n";
    std::cout << "Processing time per loop entry: " << runtime.count()/nentries << " s/entry"<< std::endl;
    return 0;
}
