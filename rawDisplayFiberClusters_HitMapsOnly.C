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
#define N_LAYERS_PER_MODULE 7
#define N_FIBERS_PER_LAYER 55
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



int rawDisplayFiberClusters_HitMapsOnly(std::vector<TString> path)
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

	SCategory * pCatRawClus = SCategoryManager::getCategory(SCategory::CatFibersRawClus);
//     SCategory * pCatRawClus = SCategoryManager::getCategory(SCategory::CatFibersRaw);
    
    TH2D * hFiberHitMap = new TH2D("hFiberHitMap", "hFiberHitMap", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH1D * hQDCL1 = new TH1D("hQDCL1", "hQDCL1", 200, 0, 50000);
    TH1D * hQDCL2 = new TH1D("hQDCL2", "hQDCL2", 200, 0, 50000);
    TH1D * hQDCL3 = new TH1D("hQDCL3", "hQDCL3", 200, 0, 50000);
    TH1D * hQDCL4 = new TH1D("hQDCL4", "hQDCL4", 200, 0, 50000);
    
	Int_t mod, lay, fi;
    Int_t mod1, lay1, fi1, mod2, lay2, fi2, mod3, lay3, fi3, mod4, lay4, fi4;
    mod1=0;
    lay1=0;
    fi1=10;
    
    mod2=0;
    lay2=2;
    fi2=15;
    
    mod3=0;
    lay3=3;
    fi3=40;
    
    mod4=0;
    lay4=6;
    fi4=35;
    
	SFibersRawCluster * pRawClus;
// 	SFibersRaw * pRawClus;

 	int nLoop = loop->getEntries();

//----- histogram definition and filling    
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCat = pCatRawClus->getEntries();
        
		for (uint j = 0; j < nCat; ++j)
		{
			pRawClus = (SFibersRawCluster *)pCatRawClus->getObject(j);
//             pRawClus = (SFibersRaw *)pCatRawClus->getObject(j);
			pRawClus->getAddress(mod, lay, fi);
            
            if(pRawClus->getQDCL() > QDC_MIN &&
            pRawClus->getQDCR() > QDC_MIN /*&&
            pRawClus->getTimeL() > T_MIN &&
            pRawClus->getTimeR() > T_MIN*/)
            {
                hFiberHitMap->Fill(fi, lay);
                if(mod==mod1 && lay ==lay1 && fi ==fi1) hQDCL1->Fill(pRawClus->getQDCL());
                if(mod==mod2 && lay ==lay2 && fi ==fi2) hQDCL2->Fill(pRawClus->getQDCL());
                if(mod==mod3 && lay ==lay3 && fi ==fi3) hQDCL3->Fill(pRawClus->getQDCL());
                if(mod==mod4 && lay ==lay4 && fi ==fi4) hQDCL4->Fill(pRawClus->getQDCL());
            }	

		}
		
		loop->nextEvent();
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path[ii].ReplaceAll(".root", "_HITMAPS.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();
    hFiberHitMap->Write();

    TCanvas * c = new TCanvas("c", "c", 800, 500);
    c->Divide(3,2);
    c->cd(1);
    hFiberHitMap->Draw("colz");
    c->cd(2);
    hQDCL1->Draw();
    c->cd(3);
    hQDCL2->Draw();
    c->cd(4);
    hQDCL3->Draw();
    c->cd(5);
    hQDCL4->Draw();
    
    TCanvas * c1 = new TCanvas("c1", "c1", 800, 500);
    c1->cd();
    hFiberHitMap->Draw("colz");
    
    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
