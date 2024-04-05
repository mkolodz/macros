R__ADD_INCLUDE_PATH(/scratch/gccb/software/framework/20230109-4to1-ubuntu20/include)
#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
#include "SFibersIdentification.h"
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
#define NBINS_QDC 100000
#define XLOW_QDC 0
#define XUP_QDC 100

#define NBINS_T 20
#define XLOW_T 0
#define XUP_T 20

#define NBINS_TDIFF 20
#define XLOW_TDIFF 0
#define XUP_TDIFF 20
// #define XUP_TDIFF_ZOOM 100

#define NBINS_MLR 100
#define XLOW_MLR -15
#define XUP_MLR 15

//----- cuts
#define TDIFF_MAX 10
#define QDC_MIN 10E-6
#define T_MIN 0


int nlcInvestigation(std::vector<TString> path)
{
    auto start = std::chrono::system_clock::now();  
    int ii=0;
    if(!path[ii].Contains("/") || !path[ii].BeginsWith("/"))
    {
        std::cout << "##### Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }
    
    gStyle->SetPalette(kBird);
    double rescaledQDC=0;
    double corrQDC=0;
    double calibFactor=13.9238;
    double corrCalibFactor=0;
    std::vector<std::string> str_path;
    for(int i=0; i < path.size(); i++){
        str_path.push_back(std::string(path[i]));
        printf("HERE: %s", str_path[i].c_str());
    }
    
	SLoop * loop = new SLoop();
	loop->addFiles(str_path);
	loop->setInput({});

	SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);
//     vector<Address> activeAddresses;
    
//     TH1D * hA[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
// 	TH1D * hQ[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
//     TH1D * hT[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
//     
// 	TH1D * hQSum[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
//     TH1D * hQAve[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
// 	TH1D * hTDif[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
//     TH1D * hMLR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
//     
// 	TH2D * hQLvsQR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
//     TH2D * hTLvsTR[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER]={nullptr};
//     
// 	TH1D * hFiberMult = new TH1D("hFiberMult", "hSiPMMult", N_FIBERS_PER_LAYER+1, -0.5, N_FIBERS_PER_LAYER+0.5); 

    TH1D * hQ;
    hQ = new TH1D("hQ", "hQ", NBINS_QDC, XLOW_QDC, XUP_QDC);
	Int_t mod, lay, fi, side;
	SFibersRaw * pRaw;
	UShort_t mult = 0;
    Int_t first_time_l = 0;
    Int_t first_time_r = 0;
    Int_t * activeSiPMs[N_MODULES][N_LAYERS_PER_MODULE][N_FIBERS_PER_LAYER][N_SIDES]={nullptr};
    double p0, p1, p2, p3;
    p0=8.0;
    p1=1.04676;
    p2=1.02734;
    p3=0.31909;	
    corrCalibFactor = (p0*pow(p1, pow(calibFactor, p2))+p3*calibFactor-p0);
    std::cout << corrCalibFactor << std::endl;
 	int nLoop = loop->getEntries();

//----- histogram definition and filling    
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCat = pCatRaw->getEntries();
        
		for (uint j = 0; j < nCat; ++j)
		{
			pRaw = (SFibersRaw *)pCatRaw->getObject(j);
			pRaw->getAddress(mod, lay, fi);
            

//             if(pRaw->getQDCL() > QDC_MIN &&
//             pRaw->getQDCR() > QDC_MIN)
//             {
//                 std::cout << "tL" <<pRaw->getTimeL()-first_time_l << std::endl;
//                 std::cout  << "tR"<< pRaw->getTimeR()- first_time_r<< std::endl;
//                 std::cout  << "tDiff"<< fabs(pRaw->getTimeL() - pRaw->getTimeR() ) << std::endl;
if(mod==0 && lay==0 && fi ==0){
//     rescaledQDC=pRaw->getQDCL()/13.9238*511.0;
//     hQ->Fill(rescaledQDC); //works (no NLC file, no NLC correction)
    
//     corrQDC = (p0*pow(p1, pow(pRaw->getQDCL(), p2))+p3*pRaw->getQDCL()-p0); //non-linearity correction
//     rescaledQDC=corrQDC/corrCalibFactor*511.0;
//     hQ->Fill(rescaledQDC); //works (no NLC file, factor and QDC are nlc-corrected)

/*    rescaledQDC=pRaw->getQDCL()/corrCalibFactor*511.0;
    hQ->Fill(rescaledQDC); //works (NLC file, factor is nlc-corrected, QDC is not nlc-corrected)  */ 
    
    hQ->Fill(pRaw->getQDCL()); 
    
}

                
//             }	
			
		}
		loop->nextEvent();
	}
	
 TCanvas * canQDC = new TCanvas("canQDC", "canQDC", 600,600); 	
	canQDC->cd();
    hQ->Draw();
return 0;
}
