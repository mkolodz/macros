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
#define NBINS_QDC 200
#define XLOW_QDC 0
// #define XUP_QDC 800
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



int fibersQDCsingleChannelsComparison(std::vector<TString> path)
{
    auto start = std::chrono::system_clock::now();  
    int ii=0;
    int fiberA=54;
    int layerA=0;
    int fiberB=46;
    int layerB=0;
    int entriesA=0;
    int entriesB=0;
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

	SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);

    TH2D * hFiberHeatMap = new TH2D("hFiberHeatMap", "hFiberHeatMap", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH1D * hQDCLChanA = new TH1D("hQDCLChanA", Form("QDCL_M0L%dF%d;QDC (a.u.);Counts", layerA, fiberA), NBINS_QDC, XLOW_QDC, XUP_QDC);
    TH1D * hQDCLChanB = new TH1D("hQDCLChanB", Form("QDCL_M0L%dF%d;QDC (a.u.);Counts", layerB, fiberB), NBINS_QDC, XLOW_QDC, XUP_QDC);
    TH1D * hQDCRChanA = new TH1D("hQDCRChanA", Form("QDCR_M0L%dF%d;QDC (a.u.);Counts", layerA, fiberA), NBINS_QDC, XLOW_QDC, XUP_QDC);
    TH1D * hQDCRChanB = new TH1D("hQDCRChanB", Form("QDCR_M0L%dF%d;QDC (a.u.);Counts", layerB, fiberB), NBINS_QDC, XLOW_QDC, XUP_QDC);
    TH1D * hQavChanA = new TH1D("hQavChanA", Form("Qav_M0L%dF%d;QDC (a.u.);Counts", layerA, fiberA), NBINS_QDC, XLOW_QDC, XUP_QDC);
    TH1D * hQavChanB = new TH1D("hQavChanB", Form("Qav_M0L%dF%d;QDC (a.u.);Counts", layerB, fiberB), NBINS_QDC, XLOW_QDC, XUP_QDC);
	Int_t mod, lay, fi;
	SFibersRaw * pRaw;

 	int nLoop = loop->getEntries();

//----- histogram definition and filling    
	for (int i = 0; i < nLoop; ++i)
	{
		size_t nCat = pCatRaw->getEntries();
        
		for (uint j = 0; j < nCat; ++j)
		{
			pRaw = (SFibersRaw *)pCatRaw->getObject(j);
			pRaw->getAddress(mod, lay, fi);
            
            if(pRaw->getQDCL() > QDC_MIN &&
            pRaw->getQDCR() > QDC_MIN &&
            pRaw->getTimeL() > T_MIN &&
            pRaw->getTimeR() > T_MIN)
            {
                hFiberHeatMap->Fill(fi, lay);
                if(mod==0 && lay==layerA && fi==fiberA){
                    hQDCLChanA->Fill(pRaw->getQDCL());
                    hQDCRChanA->Fill(pRaw->getQDCR());
                    hQavChanA->Fill(sqrt(pRaw->getQDCR()*pRaw->getQDCL()));
                    entriesA++;
                }
                if(mod==0 && lay==layerB && fi==fiberB){
                    hQDCLChanB->Fill(pRaw->getQDCL());
                    hQDCRChanB->Fill(pRaw->getQDCR());
                    hQavChanB->Fill(sqrt(pRaw->getQDCR()*pRaw->getQDCL()));
                    entriesB++;
                }
            }	

		}
		
		loop->nextEvent();
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;
    TCanvas * c = new TCanvas("c", "c", 1800, 600);
    c->Divide(3);
    c->cd(1);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    hQDCLChanA->Scale(1./entriesA);
    hQDCLChanA->Draw("HIST");
    hQDCLChanB->Scale(1./entriesB);
    hQDCLChanB->SetLineColor(kRed);
    hQDCLChanB->Draw("same,HIST");
    gPad->BuildLegend();
    c->cd(2);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    hQDCRChanA->Scale(1./entriesA);
    hQDCRChanA->Draw("HIST");
    hQDCRChanB->SetLineColor(kRed);
    hQDCRChanB->Scale(1./entriesB);
    hQDCRChanB->Draw("same,HIST");
    gPad->BuildLegend();
    c->cd(3);
    gPad->SetGrid(1,1);
    gPad->SetLogy();
    hQavChanA->Scale(1./entriesA);
    hQavChanA->Draw("HIST");
    hQavChanB->SetLineColor(kRed);
    hQavChanB->Scale(1./entriesB);
    hQavChanB->Draw("same,HIST");
    gPad->BuildLegend();
    
	path[ii].ReplaceAll(".root", "_HEATMAPS_ONLY_TEST.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();
    hFiberHeatMap->Write();


    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
