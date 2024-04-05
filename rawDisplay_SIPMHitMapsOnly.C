#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
#include "SFibersRawCluster.h"
#include "SFibersRawClusterFinder.h"
#include "SSiPMHit.h"
#include "SSiPMCluster.h"
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

//----- prototype geometry 4to1
#define N_MODULES 1
#define N_LAYERS_PER_MODULE 4
#define N_SIPMS_PER_LAYER 28
#define N_SIDES 2

//----- canvas size
#define XCANVAS 1600
#define YCANVAS 1000

// //TOFPET 
#define NBINS_QDC 200
#define XLOW_QDC 0
#define XUP_QDC 100

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


int rawDisplay_SiPMHitMapsOnly(std::vector<TString> path)
{
    auto start = std::chrono::system_clock::now();  
    int ii=0;
    if(!path[ii].Contains("/") || !path[ii].BeginsWith("/"))
    {
        std::cout << "##### Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }
    
    gStyle->SetPalette(kBird);
    
    std::vector<std::string> str_path;
    for(int i=0; i < path.size(); i++){
        str_path.push_back(std::string(path[i]));
        printf("HERE: %s", str_path[i].c_str());
    }
    
	SLoop * loop = new SLoop();
	loop->addFiles(str_path);
	loop->setInput({});

    SCategory * pCatSiPMHit = SCategoryManager::getCategory(SCategory::CatSiPMHit);
	SCategory * pCatRawClus = SCategoryManager::getCategory(SCategory::CatFibersRawClus);
    
    TH2D * hSiPMHitMapTop = new TH2D("hSiPMHitMapTop", "hSiPMHitMapTop", N_SIPMS_PER_LAYER, 0, N_SIPMS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hSiPMHitMapBot = new TH2D("hSiPMHitMapBot", "hSiPMHitMapBot", N_SIPMS_PER_LAYER, 0, N_SIPMS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMapTopSemiUnique = new TH2D("hFiberHitMapTopSemiUnique", "hFiberHitMapTopSemiUnique", N_SIPMS_PER_LAYER, 0, N_SIPMS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMapBotSemiUnique = new TH2D("hFiberHitMapBotSemiUnique", "hFiberHitMapBotSemiUnique", N_SIPMS_PER_LAYER, 0, N_SIPMS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMapBotWhenTopSemiUnique = new TH2D("hFiberHitMapBotWhenTopSemiUnique", "hFiberHitMapBotWhenTopSemiUnique", N_SIPMS_PER_LAYER, 0, N_SIPMS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMapTopWhenBotSemiUnique = new TH2D("hFiberHitMapTopWhenBotSemiUnique", "hFiberHitMapTopWhenBotSemiUnique", N_SIPMS_PER_LAYER, 0, N_SIPMS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
     
	Int_t mod, lay, fi;
    char side;
    side = ' ';

	SSiPMHit * pHit;
    SFibersRawCluster * pRawClus;

 	int nLoop = loop->getEntries();
    std::vector<bool> isTopSemiUnique;
    std::vector<bool> isBotSemiUnique;
    int counter_rawClus = 0;
    int counter_SiPMHit = 0;
    loop->getEvent(0);
    for(int i=0; i<nLoop; i++){
        size_t nCat = pCatRawClus->getEntries();
        isTopSemiUnique.push_back(false);
        isBotSemiUnique.push_back(false);
		for (uint j = 0; j < nCat; ++j)
		{
            pRawClus = (SFibersRawCluster *)pCatRawClus->getObject(j);
            if (pRawClus->getFiberClusterLabel()==3) {
                isTopSemiUnique[isTopSemiUnique.size()-1]=true;
                break;
            }
            if (pRawClus->getFiberClusterLabel()==4) {
                isBotSemiUnique[isBotSemiUnique.size()-1]=true;
                break;
            }
            
        }
        loop->nextEvent();
    }
    loop->getEvent(0);
    for(int i=0; i<nLoop; i++){
        size_t nCat = pCatSiPMHit->getEntries();
		for (uint j = 0; j < nCat; ++j)
		{
            pHit = (SSiPMHit *)pCatSiPMHit->getObject(j);
			pHit->getAddress(mod, lay, fi, side);
            if(pHit->getQDC() > QDC_MIN)
            {
                if(side=='r') hSiPMHitMapTop->Fill(fi, lay);
                if(side=='l') hSiPMHitMapBot->Fill(fi, lay);
                if(isTopSemiUnique[i]==true) {
                    if(side=='r') hFiberHitMapTopSemiUnique->Fill(fi, lay);
                    if(side=='l') hFiberHitMapBotWhenTopSemiUnique->Fill(fi, lay);
                }
                if(isBotSemiUnique[i]==true){
                    if(side=='l') hFiberHitMapBotSemiUnique->Fill(fi, lay);
                    if(side=='r') hFiberHitMapTopWhenBotSemiUnique->Fill(fi, lay);
                }
            }
        }
        loop->nextEvent();
    }
    
    
    
// //----- histogram definition and filling    
// 	for (int i = 0; i < nLoop; ++i)
// 	{
// 		size_t nCat = pCatRawClus->getEntries();
//         
// 		for (uint j = 0; j < nCat; ++j)
// 		{
// 			pHit = (SSiPMHit *)pCatSiPMHit->getObject(j);
// 			pHit->getAddress(mod, lay, fi, side);
//             
// //             pRawClus = (SFibersRawCluster *)pCatRawClus->getObject(j);
// // 			pRawClus->getAddress(mod, lay, fi);
//             
// //             if(pHit->getQDCL() > QDC_MIN &&
// //             pHit->getQDCR() > QDC_MIN /*&&
// //             pRawClus->getTimeL() > T_MIN &&
// //             pRawClus->getTimeR() > T_MIN*/)
// //             {
//                 if(side=='r') hSiPMHitMapTop->Fill(fi, lay);
//                 if(side=='l') hSiPMHitMapBot->Fill(fi, lay);
// //                 if(side=='r' || pRawClus->getFiberClusterLabel()==3)hFiberHitMapTopSemiUnique->Fill(fi, lay);
// //                 if(side=='r' || pRawClus->getFiberClusterLabel()==4)hFiberHitMapTopSemiUnique->Fill(fi, lay);
// //                 if(side=='l' || pRawClus->getFiberClusterLabel()==3)hFiberHitMapBotSemiUnique->Fill(fi, lay);
// //                 if(side=='l' || pRawClus->getFiberClusterLabel()==4)hFiberHitMapBotSemiUnique->Fill(fi, lay);
// //             }	
// 		}
// 		
// 		loop->nextEvent();
// 	}
// 	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path[ii].ReplaceAll(".root", "_HITMAPS.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();
    hSiPMHitMapTop->Write();
    hSiPMHitMapBot->Write();
    
    TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(2,3);
    c1->cd(1);
    hSiPMHitMapTop->Draw("colz");
    c1->cd(2);
    hSiPMHitMapBot->Draw("colz");
    c1->cd(3);
    hFiberHitMapTopSemiUnique->Draw("colz");
    c1->cd(4);
    hFiberHitMapBotSemiUnique->Draw("colz");
    c1->cd(5);
    hFiberHitMapBotWhenTopSemiUnique->Draw("colz");
    c1->cd(6);
    hFiberHitMapTopWhenBotSemiUnique->Draw("colz");
    
    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
