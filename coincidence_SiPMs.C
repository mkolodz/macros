//#include "/home/gabriel/sifi-framework/lib/base/core/SCategoryManager.h"
#include "SCategoryManager.h"
// #include "/home/gabriel/sifi-framework/lib/base/util/SLoop.h"
#include "SLoop.h"
#include "SDDSamples.h"
// #include "/home/gabriel/sifi-framework/lib/fibers/SFibersRaw.h"
#include "SFibersRaw.h"
// #include "SiFiSeries.h"
// #include "SiFiDatabase.h"
R__ADD_INCLUDE_PATH(/home/magda/project/sifi-framework/sifi-framework-install/include)
#include "SFibersIdentification.h"
#include "SSiPMHit.h"
// #include "SSiPMCluster.h"
#include "SiFiSeries.h"

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
#include <algorithm>

// //----- prototype geometry
// #define N_MODULES 1
// #define N_LAYERS_PER_MODULE 7
// #define N_FIBERS_PER_LAYER 55
// #define N_SIDES 2

//----- prototype geometry TOFPET4to1 (SiPMs treated as fibers)
#define N_MODULES 1
#define N_LAYERS_PER_MODULE 4
#define N_FIBERS_PER_LAYER 28
#define N_SIDES 2

//to run this macro, load the following line first: source /scratch/gccb/macros/profile.sh

// This macro takes all the events at each time stamp and provides a heat map (right channels)x(left channels)
void coincidence_SiPMs(TString path /*= "/scratch1/gccb/data/Citiroc/results/Run624_sifi_results.root"*/)
{
    SLoop * loop = new SLoop();
    loop->addFile(std::string(path));
    loop->setInput({});
    
    
    
//     SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);
    SCategory * pCatSiPM = SCategoryManager::getCategory(SCategory::CatSiPMHit); 
    Int_t mod, lay, fi;
    char side;
//     SFibersRaw * pRaw;
    SSiPMHit * pHit;
    int nLoop = loop->getEntries();
    
    TCanvas *c = new TCanvas("c", "c");
    gPad->SetGrid(1,1);
    TH2F *h = new TH2F("CorrelationSiPMs", "CorrelationSiPMs", N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER, 0, N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER, 0, N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER); // x-> fibers (0,31); y->fibers (32,63)
    TH2F *h = new TH2F("HeatMapSiPMs", "HeatMapSiPMs", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    float channelQ[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]; //charges at every channel 
    for(int c = 0; c < N_LAYERS_PER_MODULE; c++){
    for( int a =0; a < N_FIBERS_PER_LAYER; a++) // filling with -100
        {
            for (int b = 0; b < 2; b++)
            {
                channelQ[0][c][a][b] = -100;
            }
        }
    }
    
//     int channelMap[64][2] =  {{39,3},{38,2},{37,1},{36,0}}; // # of channel on each side of a fibre i.e. fiber[0] has (ch39,ch3) on left and right respectively
    int treshold = 0; // for charges

    for (int i = 0; i < nLoop; ++i) // over each time-stamp
    {
        //cout << "time stamp#" << i << endl; //test
        size_t nCat = pCatSiPM->getEntries();
        int events =0;
        for (uint j = 0; j < nCat; ++j) // over each event in time stamp
        {
            pHit = (SSiPMHit*)pCatSiPM->getObject(j);
            pHit->getAddress(mod, lay, fi, side);
            channelQ[mod][lay][fi][0]=pHit->getQDC();
            channelQ[mod][lay][fi][1]=pHit->getQDC();
            /*cout << "event# = " << j <<", at (mod,lay,fib) = ("<< mod << "," << lay <<"," << fi << "), chargeL =" << pHit-> getQDC() <<" , chargeR =" << pHit-> getQDC() << "--> ("<< channelQ[mod][lay][fi][0] << "," << channelQ[mod][lay][fi][1] << ")" << endl;  
         */   
            //cout <<j << " "  << i << " "  << mod << " " << lay << " " << fi << " " << pHit->getQDCL() << " " << pHit->getTimeL() << " " << pHit->getQDCR() << " " << pHit->getTimeR() << std::endl;
        }
//         cout << "{{" << channelQ[mod][lay][0][0] << "," << channelQ[mod][lay][0][1] << "}, " <<
//         "{" << channelQ[mod][lay][1][0] << "," << channelQ[mod][lay][1][1] << "}," <<
//         "{" << channelQ[mod][lay][2][0] << "," << channelQ[mod][lay][2][1] << "}," <<
//         "{" << channelQ[mod][lay][3][0] << "," << channelQ[mod][lay][3][1] << "}}" << endl;
        for(int kk = 0; kk < N_LAYERS_PER_MODULE; kk++)
        {
        for(int k = 0; k < N_FIBERS_PER_LAYER; k++) // N_LAYERS_PER_MODULE = number of fibers, this loop looks for coincidences
        {
            for( int ll = 0; ll<N_LAYERS_PER_MODULE; ll++) 
            {
            for( int l = 0; l<N_FIBERS_PER_LAYER; l++) // for every ch VS every ch
            {
                if( 
//                     //channelQ[0][0][k][0]!= -100. &&  
//                     //channelQ[0][0][l][1] != -100. && 
                    channelQ[0][kk][k][0] > treshold && 
                    channelQ[0][ll][l][1] > treshold)
                {
                    h->Fill(ll*N_FIBERS_PER_LAYER+l,kk*N_FIBERS_PER_LAYER+k);
                    cout << ll*N_FIBERS_PER_LAYER+l << " " << kk*N_FIBERS_PER_LAYER+k << endl;
                    cout   << "i =" << i << ", filling with charges (" <<channelQ[0][ll][l][0] << "," << channelQ[0][ll][l][1] << ")" << endl;
                }
            }
            }
        }
        }
        
        for(int m = 0; m < N_MODULES; m++)
        {
            for(int l = 0; l < N_LAYERS_PER_MODULE; l++)
            {
                for(int f = 0; f<N_FIBERS_PER_LAYER; f++) 
                {
                    for(int s = 0; s<N_SIDES; s++) 
                        {
                            channelQ[m][l][f][s] = -100;
                        }
                }
            }
        }
//         cout << endl;
        //if(i == 10) break;
        loop->nextEvent();
    }
    h->GetXaxis()->SetTitle("Right ch");
    h->GetYaxis()->SetTitle("Left ch");
    h->GetXaxis()->SetNdivisions(004);
    h->GetYaxis()->SetNdivisions(004);
    gStyle->SetPalette(kBird);
    h->Draw("COLZ");
    TLine *l1 = new TLine(0,0,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER);
    TLine *l2 = new TLine(0,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER,N_LAYERS_PER_MODULE*N_FIBERS_PER_LAYER,0);
    l1->Draw();
    l2->Draw();
    
}
