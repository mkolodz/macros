#include "SCategoryManager.h"
#include "SLoop.h"
#include "SDDSamples.h"
#include "SFibersRaw.h"
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
#include <string> 
#include <stdlib.h>
#include <cstring>
#include <chrono>
#include <ctime> 

//----- prototype geometry
#define N_MODULES 1
#define N_LAYERS 4
#define N_FIBERS 16
#define N_CROSSED_FIBERS 8
#define N_FIBERS_PER_LAYER 16
#define N_FIBERS_ALL 64
#define N_SIDES 2
#define N_READ_OUT_SIPMS 16
#define N_FIBERS_IN_CONTROL_LAYER 4
#define N_FIBERS_IN_EVENT_LAYER 4

//----- binning and histogram ranges
#define NBINS_ADC 200
#define XLOW_ADC 0
#define XUP_ADC 100

#define NBINS_QDC 100
#define XLOW_QDC 0
#define XUP_QDC 3000

#define NBINS_T0 500
#define XLOW_T0 0
#define XUP_T0 500

#define NBINS_T0DIFF 200
#define XLOW_T0DIFF -500
#define XUP_T0DIFF 500

#define NBINS_TOT 100
#define XLOW_TOT 0
#define XUP_TOT 200

#define NBINS_MLR 100
#define XLOW_MLR -5
#define XUP_MLR 5

//----- cuts
#define ADC_MAX 600
#define T0DIFF_MAX 10
#define QDC_MIN 0
#define TOT_MIN 0
#define T0_MIN 0
#define cutT0diff 50
#define HIGH_QDC_DEPOSITE 700

//----- canvas size
#define XCANVAS 1600
#define YCANVAS 1000

using namespace std;

void format_h_1D(TH1 * h)
{
// 	h->GetXaxis()->SetLabelSize(0.1);
// 	h->GetXaxis()->SetTitleSize(0.1);
// 	h->GetXaxis()->SetTitleOffset(1.1);
// 	h->GetXaxis()->SetNdivisions(505);
// 	h->GetYaxis()->SetLabelSize(0.1);
// 	h->GetYaxis()->SetTitleSize(0.1);
// 	h->GetYaxis()->SetTitleOffset(0.8);
 	h->GetYaxis()->SetMaxDigits(3);
// 	h->GetZaxis()->SetLabelSize(0.06);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);
// 	gPad->SetFrameLineWidth(1); 
// 	gPad->SetTicky(1);
// 	gStyle->SetTitleFontSize(0.1);
}



TString ConstructPdfName(TString pathname){

  TString out = pathname(0,26);
  out.Append("pdf/");
  TObjArray* Strings = pathname.Tokenize("/");
  TIter iString(Strings);
  TObjString* os = 0;
  TString tmp = "";
  
  
  while((os=(TObjString*)iString()))
  {
    tmp = os->GetString();
    if((tmp.Contains("2021") || tmp.Contains("2022") ) && tmp.Contains("_"))
      out.Append(tmp);
    else if(tmp.Contains(".root"))
    {
      out.Append("_");
      out.Append(tmp);
    }
  }
  
  out.ReplaceAll(".root",".pdf");
  delete Strings;
  
  return out;
}


TString ConstructFileNameOnly(TString pathname){

  TString out = pathname(26,43);
//   out.Append("pdf/");
//   TObjArray* Strings = pathname.Tokenize("/");
//   TIter iString(Strings);
//   TObjString* os = 0;
//   TString tmp = "";
//   
//   
//   while((os=(TObjString*)iString()))
//   {
//     tmp = os->GetString();
//     if((tmp.Contains("2021") || tmp.Contains("2022") ) && tmp.Contains("_"))
//       out.Append(tmp);
//     else if(tmp.Contains(".root"))
//     {
//       out.Append("_");
//       out.Append(tmp);
//     }
//   }
//   
//   out.ReplaceAll(".root",".pdf");
//   delete Strings;
  
  return out;
}




int rawDisplayCrossedFibers(TString path)
{
    
    if(!path.Contains("/") || !path.BeginsWith("/"))
    {
        cout << "##### Error! The functions needs the filename including the full absolute path..." << endl;
        std::abort();
    }
    
    gStyle->SetPalette(kBird);
//     gStyle->SetOptStat(0);
	SLoop * loop = new SLoop();
	loop->addFile(std::string(path));
	loop->setInput({});

	SCategory * pCatRaw = SCategoryManager::getCategory(SCategory::CatFibersRaw);
	SCategory * pCatSample = SCategoryManager::getCategory(SCategory::CatDDSamples);

    TH1D * hA[N_MODULES][N_LAYERS][N_FIBERS][N_SIDES];
	TH1D * hQ[N_MODULES][N_LAYERS][N_FIBERS][N_SIDES];
    TH1D * hT0[N_MODULES][N_LAYERS][N_FIBERS][N_SIDES];
    TH1D * hTOT[N_MODULES][N_LAYERS][N_FIBERS][N_SIDES];
    
	TH1D * hQSum[N_MODULES][N_LAYERS][N_FIBERS];
    TH1D * hQAve[N_MODULES][N_LAYERS][N_FIBERS];
	TH1D * hTDif[N_MODULES][N_LAYERS][N_FIBERS];
    TH1D * hMLR[N_MODULES][N_LAYERS][N_FIBERS];
    
    TH2D * hQLvsAL[N_MODULES][N_LAYERS][N_FIBERS];
    TH2D * hQRvsAR[N_MODULES][N_LAYERS][N_FIBERS];
	TH2D * hQLvsQR[N_MODULES][N_LAYERS][N_FIBERS];

	TH1D * hFiberMult = new TH1D("hFiberMult", "hFiberMult", N_FIBERS+1, -0.5, N_FIBERS+0.5); 
//     TH2D * hFiberHits = new TH2D("hFiberHits", "hFiberHits", N_FIBERS, -0.5, N_FIBERS-0.5, N_FIBERS, -0.5, N_FIBERS-0.5);
    TH2D * hSiPMMapping = new TH2D("hSiPMMapping", "hSiPMMapping", N_READ_OUT_SIPMS, -0.5, N_READ_OUT_SIPMS-0.5, N_READ_OUT_SIPMS, -0.5, N_READ_OUT_SIPMS-0.5);
    TH1D * hTeventTcontrolL[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER]; //change name to TcontrolTevent_diff or similar
    TH1D * hTeventTcontrolR[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER]; // -||-
    TH1D * hQeventQcontrolL[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER];
    TH1D * hQeventQcontrolR[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER];
    TH2D * h2QeventQcontrolL[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER];
    TH2D * h2QeventQcontrolR[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER];
    TH1D * hOnlyQcontrolL[N_FIBERS_IN_CONTROL_LAYER];
    TH1D * hOnlyQcontrolR[N_FIBERS_IN_CONTROL_LAYER];
    TH1D * hQDCHighDepL[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER];
    TH1D * hQDCHighDepR[N_FIBERS_IN_EVENT_LAYER][N_FIBERS_IN_CONTROL_LAYER];
    Int_t aSiPMHit[N_READ_OUT_SIPMS] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<Double_t> vTimeL[N_CROSSED_FIBERS],vTimeR[N_CROSSED_FIBERS];
    std::vector<Double_t> vQDCL[N_CROSSED_FIBERS],vQDCR[N_CROSSED_FIBERS];
	std::vector<Double_t> vQDCHighDepL,vQDCHighDepR, vFibNHighDepL, vFibNHighDepR;
    auto hProfL  = new TProfile("hProfL",";Fiber number;QDC [a.u.]",4,3.5,7.5,0,3000); //nbins,xlow,xup,ylow,yup
	auto hProfR  = new TProfile("hProfR",";Fiber number;QDC [a.u.]",4,3.5,7.5,0,3000); //nbins,xlow,xup,ylow,yup
    auto hProfNormL = new TProfile("hProfNormL",";Fiber number;QDC/QDC_E [a.u.]",4,3.5,7.5,0,3000);
    auto hProfNormR = new TProfile("hProfNormR",";Fiber number;QDC/QDC_E [a.u.]",4,3.5,7.5,0,3000);
    Int_t mod, lay, fi, side;
    
    for(int i = 0; i< N_FIBERS_IN_EVENT_LAYER; i++){
            std::string hnameQOnlyL = std::string("OnlyQcontrol_c") + std::to_string(i) + std::string("_L");
            std::string hnameQOnlyR = std::string("OnlyQcontrol_c") + std::to_string(i) + std::string("_R");
            hOnlyQcontrolL[i] = new TH1D(hnameQOnlyL.c_str(), hnameQOnlyL.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
            hOnlyQcontrolR[i] = new TH1D(hnameQOnlyR.c_str(), hnameQOnlyR.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);

        for(int j = 0; j< N_FIBERS_IN_CONTROL_LAYER; j++){
            std::string hnameTL = std::string("TeventTcontrol_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_L");
            std::string hnameTR = std::string("TeventTcontrol_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_R");
            std::string hnameQL = std::string("QeventQcontrol_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_L");
            std::string hnameQR = std::string("QeventQcontrol_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_R");
            std::string hnameQL2D = std::string("QeventQcontrol2D_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_L");
            std::string hnameQR2D = std::string("QeventQcontrol2D_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_R");
            std::string hnameQDepL = std::string("HighDepQeventQcontrol_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_L");
            std::string hnameQDepR = std::string("HighDepQeventQcontrol_e") + std::to_string(i) + std::string("_c") + std::to_string(j) + std::string("_R");
            
            hTeventTcontrolL[i][j] = new TH1D(hnameTL.c_str(), hnameTL.c_str(), NBINS_T0DIFF, XLOW_T0DIFF, XUP_T0DIFF);
            hTeventTcontrolR[i][j] = new TH1D(hnameTR.c_str(), hnameTR.c_str(), NBINS_T0DIFF, XLOW_T0DIFF, XUP_T0DIFF);
            hQeventQcontrolL[i][j] = new TH1D(hnameQL.c_str(), hnameQL.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
            hQeventQcontrolR[i][j] = new TH1D(hnameQR.c_str(), hnameQR.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
            h2QeventQcontrolL[i][j] = new TH2D(hnameQL2D.c_str(), hnameQL2D.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC, NBINS_QDC, XLOW_QDC, XUP_QDC);
            h2QeventQcontrolR[i][j] = new TH2D(hnameQR2D.c_str(), hnameQR2D.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC, NBINS_QDC, XLOW_QDC, XUP_QDC);
            hQDCHighDepL[i][j] = new TH1D(hnameQDepL.c_str(), hnameQDepL.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
            hQDCHighDepR[i][j] = new TH1D(hnameQDepR.c_str(), hnameQDepR.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);

        }
    }
    
	for(mod = 0; mod < N_MODULES; mod++)
	{
		for(lay = 0; lay < N_LAYERS; lay++)
		{
			for (fi = 0; fi < N_FIBERS; ++fi)
			{
				std::string suffix = std::string("_M") + std::to_string(mod) + std::string("L") + std::to_string(lay) + std::string("F") + std::to_string(fi);
                std::string A_L = std::string("A") + suffix + std::string("L");
                std::string A_R = std::string("A") + suffix + std::string("R");
				std::string Q_L = std::string("Q") + suffix + std::string("L");
				std::string Q_R = std::string("Q") + suffix + std::string("R");
				std::string T0_L = std::string("T0") + suffix + std::string("L");
				std::string T0_R = std::string("T0") + suffix + std::string("R");
				std::string TOT_L = std::string("TOT") + suffix + std::string("L");
				std::string TOT_R = std::string("TOT") + suffix + std::string("R");
				std::string TDif = std::string("TDif") + suffix;
                std::string QAve = std::string("QAve") + suffix;
				std::string QSum = std::string("QSum") + suffix;
				std::string MLR = std::string("MLR") + suffix;
				std::string QLvsQR = std::string("QLvsQR") + suffix;
				std::string QLvsAL = std::string("QLvsAL") + suffix;
				std::string QRvsAR = std::string("QRvsAR") + suffix;

				std::string title_A_L = A_L + std::string("; A_{L} [mV]; Counts");
				std::string title_A_R = A_R + std::string("; A_{R} [mV]; Counts");
				std::string title_Q_L = Q_L + std::string("; Q_{L} [a.u.]; Counts");
				std::string title_Q_R = Q_R + std::string("; Q_{R} [a.u.]; Counts");
				std::string title_T0_L = T0_L + std::string("; T0_{L} [ns]; Counts");
				std::string title_T0_R = T0_R + std::string("; T0_{R} [ns]; Counts");
				std::string title_TOT_L = TOT_L + std::string("; TOT_{L} [ns]; Counts");
				std::string title_TOT_R = TOT_R + std::string("; TOT_{R} [ns]; Counts");
				std::string title_TDif = TDif + std::string("; T_{dif} [ns]; Counts");
                std::string title_QAve = QAve + std::string("; Q_{ave} [a.u.]; Counts");
				std::string title_QSum = QSum + std::string("; Q_{L} + Q_{R} [a.u.]; Counts");
                std::string title_MLR = MLR + std::string("; M_{LR}; Counts");
				std::string title_QLvsQR = QLvsQR + std::string("; Q_{L} [a.u.]; Q_{R} [a.u.]");
                std::string title_QRvsAR = QRvsAR + std::string("; A_{R} [mV]; Q_{R} [a.u.]");
                std::string title_QLvsAL = QLvsAL + std::string("; A_{L} [mV]; Q_{L} [a.u.]");

				hA[mod][lay][fi][0] = new TH1D(A_L.c_str(), title_A_L.c_str(), NBINS_ADC, XLOW_ADC, XUP_ADC);
				hA[mod][lay][fi][1] = new TH1D(A_R.c_str(), title_A_R.c_str(), NBINS_ADC, XLOW_ADC, XUP_ADC);
				hQ[mod][lay][fi][0] = new TH1D(Q_L.c_str(), title_Q_L.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hQ[mod][lay][fi][1] = new TH1D(Q_R.c_str(), title_Q_R.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hT0[mod][lay][fi][0] = new TH1D(T0_L.c_str(), title_T0_L.c_str(), NBINS_T0, XLOW_T0, XUP_T0);
				hT0[mod][lay][fi][1] = new TH1D(T0_R.c_str(), title_T0_R.c_str(), NBINS_T0, XLOW_T0, XUP_T0);
				hTOT[mod][lay][fi][0] = new TH1D(TOT_L.c_str(), title_TOT_L.c_str(), NBINS_TOT, XLOW_TOT, XUP_TOT);
				hTOT[mod][lay][fi][1] = new TH1D(TOT_R.c_str(), title_TOT_R.c_str(), NBINS_TOT, XLOW_TOT, XUP_TOT);
				
                hTDif[mod][lay][fi] = new TH1D(TDif.c_str(), title_TDif.c_str(), NBINS_T0DIFF, XLOW_T0DIFF, XUP_T0DIFF);
                hQAve[mod][lay][fi] = new TH1D(QAve.c_str(), title_QAve.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC);
				hQSum[mod][lay][fi] = new TH1D(QSum.c_str(), title_QSum.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC*2);
                hMLR[mod][lay][fi] = new TH1D(MLR.c_str(), title_MLR.c_str(), NBINS_MLR, XLOW_MLR, XUP_MLR);
				
                hQLvsQR[mod][lay][fi] = new TH2D(QLvsQR.c_str(), title_QLvsQR.c_str(), NBINS_QDC, XLOW_QDC, XUP_QDC, NBINS_QDC, XLOW_QDC, XUP_QDC);
                hQRvsAR[mod][lay][fi] = new TH2D(QRvsAR.c_str(), title_QRvsAR.c_str(), NBINS_ADC, XLOW_ADC, XUP_ADC, NBINS_QDC, XLOW_QDC, XUP_QDC);
                hQLvsAL[mod][lay][fi] = new TH2D(QLvsAL.c_str(), title_QLvsAL.c_str(), NBINS_ADC, XLOW_ADC, XUP_ADC, NBINS_QDC, XLOW_QDC, XUP_QDC);
			}
		}
	}

	SDDSamples * samples;
	SFibersRaw * pRaw;
	SDDSignal  * sigL;
	SDDSignal  * sigR;
	UShort_t mult = 0;
    Int_t binCount = 0;
    Int_t rVetoCounter = 0;
 	int nLoop = loop->getEntries();
    
    char * outfilename="controllayer_50mm.txt";
    std::ofstream outfile(outfilename);

	for (int i = 0; i < nLoop; ++i) 
	{
		size_t nCat = pCatRaw->getEntries();
        for(int l=0;l<N_READ_OUT_SIPMS;l++) {aSiPMHit[l] = 0;}
        
		for (uint j = 0; j < nCat; ++j)
		{
			samples = (SDDSamples*)pCatSample->getObject(j);
			pRaw = (SFibersRaw *)pCatRaw->getObject(j);
			sigL = (SDDSignal*)samples->getSignalL();
			sigR = (SDDSignal*)samples->getSignalR();
			pRaw->getAddress(mod, lay, fi);
            for (int k = 0; k < N_CROSSED_FIBERS; k++){
                if(mod == 0 && lay == 2 && fi == k){
                    vTimeL[k].push_back(pRaw->getTimeL());
                    vTimeR[k].push_back(pRaw->getTimeR());
                    vQDCL[k].push_back(pRaw->getQDCL());
                    vQDCR[k].push_back(pRaw->getQDCR());
                }
            }
            if(pRaw->getQDCL() > QDC_MIN &&
               sigL->GetAmplitude() < ADC_MAX &&
               sigL->GetTOT() > TOT_MIN &&
               pRaw->getTimeL() > T0_MIN /*&&
               sigL->GetVeto() == 0*/)
            {
                aSiPMHit[fi] = 1;
            }
            if(pRaw->getQDCR() > QDC_MIN &&
               sigR->GetAmplitude() < ADC_MAX &&
               sigR->GetTOT() > TOT_MIN &&
               pRaw->getTimeR() > T0_MIN /*&&
               sigR->GetVeto() == 0*/)
            {
                aSiPMHit[fi+8] = 1;
            }
            if(sigR->GetVeto() == 0) rVetoCounter++;
			if(pRaw->getQDCL() > QDC_MIN &&
               pRaw->getQDCR() > QDC_MIN &&
               sigL->GetAmplitude() < ADC_MAX &&
               sigR->GetAmplitude() < ADC_MAX &&
               sigL->GetTOT() > TOT_MIN &&
               sigR->GetTOT() > TOT_MIN &&
               pRaw->getTimeL() > T0_MIN &&
               pRaw->getTimeR() > T0_MIN &&
               fabs(pRaw->getTimeL() - pRaw->getTimeR()) < T0DIFF_MAX/* &&
               sigL->GetVeto() == 0 &&
               sigR->GetVeto() == 0*/)
            {
                hA[mod][lay][fi][0]->Fill(sigL->GetAmplitude());
                hA[mod][lay][fi][1]->Fill(sigR->GetAmplitude());
				hQ[mod][lay][fi][0]->Fill(pRaw->getQDCL());
				hQ[mod][lay][fi][1]->Fill(pRaw->getQDCR());
                hT0[mod][lay][fi][0]->Fill(pRaw->getTimeL());
                hT0[mod][lay][fi][1]->Fill(pRaw->getTimeR());
                hTOT[mod][lay][fi][0]->Fill(sigL->GetTOT());
                hTOT[mod][lay][fi][1]->Fill(sigR->GetTOT());
                
                hQAve[mod][lay][fi]->Fill(sqrt(pRaw->getQDCL()*pRaw->getQDCR()));
                hQSum[mod][lay][fi]->Fill(pRaw->getQDCL() + pRaw->getQDCR());
                hMLR[mod][lay][fi]->Fill(log(sqrt(pRaw->getQDCR()/pRaw->getQDCL())));
                
                hQLvsQR[mod][lay][fi]->Fill(pRaw->getQDCL(), pRaw->getQDCR());
                hQRvsAR[mod][lay][fi]->Fill(sigR->GetAmplitude(), pRaw->getQDCR());
                hQLvsAL[mod][lay][fi]->Fill(sigL->GetAmplitude(), pRaw->getQDCL());
                
				if(abs(pRaw->getTimeL() - pRaw->getTimeR() ) < T0DIFF_MAX)
                {
					++mult;
				}
				
			}	

			if(pRaw->getQDCL() > QDC_MIN &&
               pRaw->getQDCR() > QDC_MIN &&
               sigL->GetAmplitude() < ADC_MAX &&
               sigR->GetAmplitude() < ADC_MAX &&
               sigL->GetTOT() > TOT_MIN &&
               sigR->GetTOT() > TOT_MIN &&
               pRaw->getTimeL() > T0_MIN &&
               pRaw->getTimeR() > T0_MIN &&
               sigL->GetVeto() == 0 &&
               sigR->GetVeto() == 0)
            {
				hTDif[mod][lay][fi]->Fill(pRaw->getTimeL() - pRaw->getTimeR() );
			}
		}
// 		for (int p = 0; p<N_READ_OUT_SIPMS; p++){std::cout << aSiPMHit[p] << " " ;}
		
		for(int m=0;m<N_READ_OUT_SIPMS;m++){
            for(int n=0;n<N_READ_OUT_SIPMS;n++){
                if(aSiPMHit[m]==1 && aSiPMHit[n]==1){
                    binCount = hSiPMMapping->GetBinContent(m+1,n+1);
                    hSiPMMapping->SetBinContent(m+1,n+1, binCount+1);
//                     std::cout << m <<" "<< n << " " << binCount<< std::endl;
                }
            }
        }
		hFiberMult->Fill(mult);

		mult = 0;
		loop->nextEvent();
	}
//     for(int m=0;m<N_CROSSED_FIBERS;m++){
//     std::cout << m << " " << vQDCL[m].size() << std::endl;
//     }
//     for(int m=0;m<N_CROSSED_FIBERS;m++){
//     std::cout << m << " " << vTimeR[m].size() << std::endl;
//     }
    
    for(int i = 4; i<N_CROSSED_FIBERS; i++){
        for(int j = 0; j<vTimeL[0].size(); j++){
            if(vQDCL[0][j] == -100  && vQDCR[0][j] == -100 && vQDCL[1][j] == -100  && vQDCR[1][j] == -100 && vQDCL[2][j] == -100  && vQDCR[2][j] == -100 && vQDCL[3][j] == -100  && vQDCR[3][j] == -100 && vQDCL[i][j] != -100  && vQDCR[i][j]!=-100){
                    hOnlyQcontrolL[i-4]->Fill(vQDCL[i][j]);
                    hOnlyQcontrolR[i-4]->Fill(vQDCR[i][j]);
//                 if(vQDCL[4][j] == -100  && vQDCR[4][j] == -100 && vQDCL[5][j] == -100  && vQDCR[5][j] == -100 && vQDCL[6][j] == -100  && vQDCR[6][j] == -100 && vQDCL[7][j] == -100  && vQDCR[7][j] == -100 && vQDCL[i-4][j] != -100  && vQDCR[i-4][j]!=-100){
//                 hOnlyQcontrolL[i-4]->Fill(vQDCL[i-4][j]);
//                 hOnlyQcontrolR[i-4]->Fill(vQDCR[i-4][j]);
            }
            
            
            for(int k=0; k<4; k++){
                if(vTimeL[k][j] != -100  && vTimeL[i][j] != -100)
                    hTeventTcontrolL[k][i-4]->Fill(vTimeL[k][j]-vTimeL[i][j]);
                if(vTimeR[k][j] != -100 && vTimeR[i][j] != -100)
                    hTeventTcontrolR[k][i-4]->Fill(vTimeR[k][j]-vTimeR[i][j]);
                
                if(fabs(vTimeL[k][j]-vTimeL[i][j]) < cutT0diff && vQDCL[i][j] !=-100){
                    hQeventQcontrolL[k][i-4]->Fill(vQDCL[i][j]);
                    h2QeventQcontrolL[k][i-4]->Fill(vQDCL[k][j], vQDCL[i][j]); //event, control
                }
                if(fabs(vTimeR[k][j]-vTimeR[i][j])<cutT0diff && vQDCR[i][j] !=-100){
                    hQeventQcontrolR[k][i-4]->Fill(vQDCR[i][j]);
                    h2QeventQcontrolR[k][i-4]->Fill(vQDCR[k][j], vQDCR[i][j]); //event, control
                }
                if(fabs(vTimeL[k][j]-vTimeL[i][j]) < cutT0diff && vQDCL[i][j] !=-100 && vQDCL[k][j]>HIGH_QDC_DEPOSITE){
                    outfile << vQDCL[4][j] << " " << vQDCL[5][j] << " "<< vQDCL[6][j] << " "<< vQDCL[7][j] << " " << std::endl; 
                    //hQDCHighDepL[k][i-4]->Fill(vQDCL[i][j]);
                    hProfL->Fill(i,vQDCL[i][j]);
                    hProfNormL->Fill(i,vQDCL[i][j]/vQDCL[k][j]);
                    //vQDCHighDepL.push_back(vQDCL[i][j]);
                    //vFibNHighDepL.push_back(i);
                }
                if(fabs(vTimeR[k][j]-vTimeR[i][j])<cutT0diff && vQDCR[i][j] !=-100 && vQDCR[k][j]>HIGH_QDC_DEPOSITE){
                    hQDCHighDepR[k][i-4]->Fill(vQDCR[i][j]);
                    hProfR->Fill(i,vQDCR[i][j]);
                    hProfNormR->Fill(i,vQDCR[i][j]/vQDCR[k][j]);
					//vQDCHighDepR.push_back(vQDCR[i][j]);
					//vFibNHighDepR.push_back(i);
                }
                

            }
        }
    }

	std::cout << "\n\nLoop entries: " << nLoop << std::endl;
    
    std::cout << "\n\nrVetoCounter: " << rVetoCounter << std::endl;

/*    
    //char command_4[500] = "rm ";
    strcat(command_4, meas_path.c_str());
    strcat(command_4, " .");
    std::cout << command_4 << std::endl;
    puts(command_4);
    status = system(command_4);
    std::cout << "Status: " << status << std::endl;*/

    
	path.ReplaceAll(".root", "_HISTOS.root");
	TFile *output = new TFile(path,"RECREATE");
	output->cd();
   
	for(Int_t m = 0; m < N_MODULES; m++)
	{
		output->mkdir(Form("Module%d/", m));
		output->cd(Form("Module%d/", m));
		for(Int_t l = 0; l < N_LAYERS; l++)
		{
			output->mkdir(Form("Module%d/Layer%d", m, l));
            output->cd(Form("Module%d/Layer%d", m, l));
			for(Int_t f = 0; f < N_FIBERS; f++)
			{
				output->mkdir(Form("Module%d/Layer%d/Fiber%d",m, l, f));
				output->cd(Form("Module%d/Layer%d/Fiber%d",m, l, f));
				hQ[m][l][f][0]->SetDirectory(gDirectory);
				hQ[m][l][f][0]->Write();
				hQ[m][l][f][1]->SetDirectory(gDirectory);
				hQ[m][l][f][1]->Write();
                hA[m][l][f][0]->SetDirectory(gDirectory);
                hA[m][l][f][0]->Write();
                hA[m][l][f][1]->SetDirectory(gDirectory);
                hA[m][l][f][1]->Write();
                hT0[m][l][f][0]->SetDirectory(gDirectory);
                hT0[m][l][f][0]->Write();
                hT0[m][l][f][1]->SetDirectory(gDirectory);
                hT0[m][l][f][1]->Write();
                hTOT[m][l][f][0]->SetDirectory(gDirectory);
                hTOT[m][l][f][0]->Write();
                hTOT[m][l][f][1]->SetDirectory(gDirectory);
                hTOT[m][l][f][1]->Write();
                
				hTDif[m][l][f]->SetDirectory(gDirectory);
				hTDif[m][l][f]->Write();
				hQSum[m][l][f]->SetDirectory(gDirectory);
				hQSum[m][l][f]->Write();
				hQAve[m][l][f]->SetDirectory(gDirectory);
				hQAve[m][l][f]->Write();
				hMLR[m][l][f]->SetDirectory(gDirectory);
				hMLR[m][l][f]->Write();
				
                hQLvsQR[m][l][f]->SetDirectory(gDirectory);
				hQLvsQR[m][l][f]->Write();
                hQLvsAL[m][l][f]->SetDirectory(gDirectory);
				hQLvsAL[m][l][f]->Write();
                hQRvsAR[m][l][f]->SetDirectory(gDirectory);
				hQRvsAR[m][l][f]->Write();
			}
		}
	}    
	
	TCanvas * canQDC = new TCanvas("QDC", "QDC", XCANVAS, YCANVAS); 
    canQDC->SetLeftMargin(0.2);
	canQDC->DivideSquare(N_FIBERS);

    TCanvas * canADC = new TCanvas("ADC", "ADC", XCANVAS, YCANVAS); 
    canADC->SetLeftMargin(0.2);
    canADC->DivideSquare(N_FIBERS);

    TCanvas * canT0 = new TCanvas("T0", "T0", XCANVAS, YCANVAS); 
    canT0->SetLeftMargin(0.2); 
    canT0->DivideSquare(N_FIBERS);

    TCanvas * canTOT = new TCanvas("TOT", "TOT", XCANVAS, YCANVAS);
    canTOT->SetLeftMargin(0.2); 
    canTOT->DivideSquare(N_FIBERS);
    
	TCanvas * canTDif = new TCanvas("TDif", "TDif", XCANVAS, YCANVAS);
    canTDif->SetLeftMargin(0.2);
	canTDif->DivideSquare(N_FIBERS);
    
    TCanvas * canQAve = new TCanvas("QAve", "QAve", XCANVAS, YCANVAS); 
    canQAve->SetLeftMargin(0.2); 
    canQAve->DivideSquare(N_FIBERS);
    
	TCanvas * canQSum = new TCanvas("QSum", "QSum", XCANVAS, YCANVAS); 
    canQSum->SetLeftMargin(0.2); 
	canQSum->DivideSquare(N_FIBERS);
    
	TCanvas * canMLR = new TCanvas("MLR", "MLR", XCANVAS, YCANVAS); 
    canMLR->SetLeftMargin(0.2); 
	canMLR->DivideSquare(N_FIBERS);
    
	TCanvas * canQLvsQR = new TCanvas("QLvsQR", "Q_{L} vs Q_{R}", XCANVAS, YCANVAS); 
    canQLvsQR->SetLeftMargin(0.2);
	canQLvsQR->DivideSquare(N_FIBERS);
    
	TCanvas * canQRvsAR = new TCanvas("QRvsAR", "Q_{R} vs A_{R}", XCANVAS, YCANVAS); 
    canQRvsAR->SetLeftMargin(0.2);
	canQRvsAR->DivideSquare(N_FIBERS);
    
	TCanvas * canQLvsAL = new TCanvas("QLvsAL", "Q_{L} vs A_{L}", XCANVAS, YCANVAS); 
    canQLvsAL->SetLeftMargin(0.2);
	canQLvsAL->DivideSquare(N_FIBERS);
    
	TCanvas * canMult = new TCanvas("Mult", "Mult", 600, 400); 
    canMult->SetLeftMargin(0.2);
    
    TCanvas * canMapping = new TCanvas("Mapping", "Mapping", 600, 600); 
    canMapping->SetLeftMargin(0.2);
    
    TCanvas * canTimeDiffL = new TCanvas("TimeDiffL", "TimeDiffL", XCANVAS, YCANVAS); 
    canTimeDiffL->SetLeftMargin(0.2);
    canTimeDiffL->DivideSquare(16);
    
    TCanvas * canTimeDiffR = new TCanvas("TimeDiffR", "TimeDiffR", XCANVAS, YCANVAS); 
    canTimeDiffR->SetLeftMargin(0.2);
    canTimeDiffR->DivideSquare(16);

    TCanvas * canQDCcutL = new TCanvas("QDCcutL", "QDCcutL", XCANVAS, YCANVAS); 
    canQDCcutL->SetLeftMargin(0.2);
    canQDCcutL->DivideSquare(16);
    
    TCanvas * canQDCcutR = new TCanvas("QDCcutR", "QDCcutR", XCANVAS, YCANVAS); 
    canQDCcutR->SetLeftMargin(0.2);
    canQDCcutR->DivideSquare(16);
    
    TCanvas * canQeQc2DL = new TCanvas("QeQc2DL", "QeQc2DL", XCANVAS, YCANVAS); 
    canQeQc2DL->SetLeftMargin(0.2);
    canQeQc2DL->DivideSquare(16);
    
    TCanvas * canQeQc2DR = new TCanvas("QeQc2DR", "QeQc2DR", XCANVAS, YCANVAS); 
    canQeQc2DR->SetLeftMargin(0.2);
    canQeQc2DR->DivideSquare(16);
    
    TCanvas * canQDCControlOnlyL = new TCanvas("QDCControlOnlyL", "QDCControlOnlyL", XCANVAS, YCANVAS); 
    canQDCControlOnlyL->SetLeftMargin(0.2);
    canQDCControlOnlyL->DivideSquare(4);
    
    TCanvas * canQDCControlOnlyR = new TCanvas("QDCControlOnlyR", "QDCControlOnlyR", XCANVAS, YCANVAS); 
    canQDCControlOnlyR->SetLeftMargin(0.2);
    canQDCControlOnlyR->DivideSquare(4);
    
    TCanvas * canQDCHighDepL = new TCanvas("QDCHighDepL", "QDCHighDepL", XCANVAS, YCANVAS); 
    canQDCHighDepL->SetLeftMargin(0.2);
    canQDCHighDepL->DivideSquare(16);
    
    TCanvas * canQDCHighDepR = new TCanvas("QDCHighDepR", "QDCHighDepR", XCANVAS, YCANVAS); 
    canQDCHighDepR->SetLeftMargin(0.2);
    canQDCHighDepR->DivideSquare(16);
    
//     TCanvas * canHits = new TCanvas("Hits", "Hits", 600, 400);
//     canHits->SetLeftMargin(0.2);

    Int_t colL = kBlack;
    Int_t colR = kRed;
    Int_t j = 2;
    Int_t jj = 0;
    	for(Int_t i = 0; i < N_MODULES; i++)
	{
// 		for(Int_t j = 0; j < N_LAYERS; j++)
// 		{
        
			for(Int_t k = 0; k < N_FIBERS; k++)
			{
				canQDC->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);
                gPad->SetGrid(1, 1);
                hQ[i][j][k][0]->SetLineColor(colL);
                hQ[i][j][k][1]->SetLineColor(colR);
				hQ[i][j][k][0]->Draw();
                hQ[i][j][k][1]->Draw("same");
                hQ[i][j][k][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hQ[i][j][k][0]->GetBinContent(hQ[i][j][k][0]->GetMaximumBin()),
                                                                    hQ[i][j][k][1]->GetBinContent(hQ[i][j][k][1]->GetMaximumBin())) + 50);
 				format_h_1D(&hQ[i][j][k][0][0]);
 				format_h_1D(&hQ[i][j][k][1][0]);

                canADC->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);  
                gPad->SetGrid(1, 1);
                hA[i][j][k][0]->SetLineColor(colL);
                hA[i][j][k][1]->SetLineColor(colR);
                hA[i][j][k][0]->Draw();
                hA[i][j][k][1]->Draw("same");
                hA[i][j][k][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hA[i][j][k][0]->GetBinContent(hA[i][j][k][0]->GetMaximumBin()),
                                                                    hA[i][j][k][1]->GetBinContent(hA[i][j][k][1]->GetMaximumBin())) + 50);
                format_h_1D(hA[i][j][k][0]);
                format_h_1D(hA[i][j][k][1]);
                
                canT0->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);
                gPad->SetGrid(1, 1);  
                hT0[i][j][k][0]->SetLineColor(colL);
                hT0[i][j][k][1]->SetLineColor(colR);
                hT0[i][j][k][0]->Draw();
                hT0[i][j][k][1]->Draw("same");
                hT0[i][j][k][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hT0[i][j][k][0]->GetBinContent(hT0[i][j][k][0]->GetMaximumBin()),
                                                                     hT0[i][j][k][1]->GetBinContent(hT0[i][j][k][1]->GetMaximumBin())) + 50);
                format_h_1D(hT0[i][j][k][0]);
                format_h_1D(hT0[i][j][k][1]);
                
                canTOT->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);  
                gPad->SetGrid(1, 1);
                hTOT[i][j][k][0]->SetLineColor(colL);
                hTOT[i][j][k][1]->SetLineColor(colR);
                hTOT[i][j][k][0]->Draw();
                hTOT[i][j][k][1]->Draw("same");
                hTOT[i][j][k][0]->GetYaxis()->SetRangeUser(0, TMath::Max(hTOT[i][j][k][0]->GetBinContent(hTOT[i][j][k][0]->GetMaximumBin()),
                                                                      hTOT[i][j][k][1]->GetBinContent(hTOT[i][j][k][1]->GetMaximumBin())) + 50);
                format_h_1D(hTOT[i][j][k][0]);
                format_h_1D(hTOT[i][j][k][1]);

				canTDif->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);  
                gPad->SetGrid(1, 1);
                hTDif[i][j][k]->Draw();
 				format_h_1D(hTDif[i][j][k]);

                canQAve->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);  
                gPad->SetGrid(1, 1);
                hQAve[i][j][k]->Draw();
                format_h_1D(hQAve[i][j][k]);

				canQSum->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);  
                gPad->SetGrid(1, 1);
				hQSum[i][j][k]->Draw();
 				format_h_1D(hQSum[i][j][k]);

				canMLR->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);  
                gPad->SetGrid(1, 1);
				hMLR[i][j][k]->Draw();
 				format_h_1D(hMLR[i][j][k]);

				canQLvsQR->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);  
                gPad->SetGrid(1, 1);
				hQLvsQR[i][j][k]->Draw("COLZ");
 				format_h_1D(hQLvsQR[i][j][k]);
                
                canQRvsAR->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);
                gPad->SetGrid(1, 1);
                hQRvsAR[i][j][k]->Draw("COLZ");
                format_h_1D(hQRvsAR[i][j][k]);
                
                canQLvsAL->cd(i*N_LAYERS*N_FIBERS + jj*N_FIBERS + k + 1);
                gPad->SetGrid(1, 1);
                hQLvsAL[i][j][k]->Draw("COLZ");
                format_h_1D(hQLvsAL[i][j][k]);
			}
// 		}
	}

	canMult->cd();  
    gPad->SetGrid(1, 1);
	hFiberMult->Draw();
 	format_h_1D(hFiberMult);
	hFiberMult->GetXaxis()->SetTitle("# fibers");
	hFiberMult->GetYaxis()->SetTitle("Counts");
    
    canMapping->cd();  
    gPad->SetGrid(1, 1);
	hSiPMMapping->Draw("colz");
 	format_h_1D(hSiPMMapping);
	hSiPMMapping->GetXaxis()->SetTitle("SiPM #");
	hSiPMMapping->GetYaxis()->SetTitle("SiPM #");
    
    for(Int_t k = 0; k < N_FIBERS_IN_EVENT_LAYER; k++)
    {
            canQDCControlOnlyL->cd(k+1);  
            gPad->SetGrid(1, 1);
            hOnlyQcontrolL[k]->Draw();
            format_h_1D(hOnlyQcontrolL[k]);
            hOnlyQcontrolL[k]->GetXaxis()->SetTitle("QDC");
            hOnlyQcontrolL[k]->GetYaxis()->SetTitle("Counts");   
            
            canQDCControlOnlyR->cd(k+1);  
            gPad->SetGrid(1, 1);
            hOnlyQcontrolR[k]->Draw();
            format_h_1D(hOnlyQcontrolR[k]);
            hOnlyQcontrolR[k]->GetXaxis()->SetTitle("QDC");
            hOnlyQcontrolR[k]->GetYaxis()->SetTitle("Counts");

        for(Int_t j = 0; j < N_FIBERS_IN_CONTROL_LAYER; j++)
        {
            canTimeDiffL->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            hTeventTcontrolL[k][j]->Draw();
            format_h_1D(hTeventTcontrolL[k][j]);
            hTeventTcontrolL[k][j]->GetXaxis()->SetTitle("Time diff [ns]");
            hTeventTcontrolL[k][j]->GetYaxis()->SetTitle("Counts");   
            
            canTimeDiffR->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            hTeventTcontrolR[k][j]->Draw();
            format_h_1D(hTeventTcontrolR[k][j]);
            hTeventTcontrolR[k][j]->GetXaxis()->SetTitle("Time diff [ns]");
            hTeventTcontrolR[k][j]->GetYaxis()->SetTitle("Counts");   
            
            canQDCcutL->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            hQeventQcontrolL[k][j]->Draw();
            format_h_1D(hQeventQcontrolL[k][j]);
            hQeventQcontrolL[k][j]->GetXaxis()->SetTitle("QDC");
            hQeventQcontrolL[k][j]->GetYaxis()->SetTitle("Counts");   
            
            canQDCcutR->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            hQeventQcontrolR[k][j]->Draw();
            format_h_1D(hQeventQcontrolR[k][j]);
            hQeventQcontrolR[k][j]->GetXaxis()->SetTitle("QDC");
            hQeventQcontrolR[k][j]->GetYaxis()->SetTitle("Counts");   
            
            canQeQc2DL->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            h2QeventQcontrolL[k][j]->Draw("COLZ");
            format_h_1D(hQeventQcontrolL[k][j]);
            h2QeventQcontrolL[k][j]->GetXaxis()->SetTitle("Event layer QDC");
            h2QeventQcontrolL[k][j]->GetYaxis()->SetTitle("Control layer QDC");   
            
            canQeQc2DR->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            h2QeventQcontrolR[k][j]->Draw("COLZ");
            format_h_1D(hQeventQcontrolL[k][j]);
            h2QeventQcontrolR[k][j]->GetXaxis()->SetTitle("Event layer QDC");
            h2QeventQcontrolR[k][j]->GetYaxis()->SetTitle("Control layer QDC");  
            
            canQDCHighDepL->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            hQDCHighDepL[k][j]->Draw();
            format_h_1D(hQDCHighDepL[k][j]);
            hQDCHighDepL[k][j]->GetXaxis()->SetTitle("QDC");
            hQDCHighDepL[k][j]->GetYaxis()->SetTitle("Counts");   
            
            canQDCHighDepR->cd(N_FIBERS_IN_CONTROL_LAYER*k+(j+1));  
            gPad->SetGrid(1, 1);
            hQDCHighDepR[k][j]->Draw();
            format_h_1D(hQDCHighDepR[k][j]);
            hQDCHighDepR[k][j]->GetXaxis()->SetTitle("QDC");
            hQDCHighDepR[k][j]->GetYaxis()->SetTitle("Counts");  
        }
    }
    
	output->mkdir("Canvases");
	output->cd("Canvases");
    
	canQDC->Write();
    canADC->Write();   
    canT0->Write();
    canTOT->Write();
	canTDif->Write();
    canQAve->Write();
	canQSum->Write();
    canMLR->Write();
	canQLvsQR->Write();    
    canQRvsAR->Write();   
    canQLvsAL->Write();
    canMult->Write();
    canMapping->Write();
    canQDCControlOnlyL->Write();
    canQDCControlOnlyR->Write();
    canTimeDiffL->Write();
    canTimeDiffR->Write();
    canQDCcutL->Write();
    canQeQc2DL->Write();
    canQeQc2DR->Write();
    canQDCcutR->Write();
    canQDCHighDepL->Write();
    canQDCHighDepR->Write();
    
//     output->Close();
    
    TString pdfname = ConstructPdfName(path);
    TString filenameonly = ConstructFileNameOnly(path);
    canQDC->Print(Form("%s(",pdfname.Data()), "pdf");
    canADC->Print(pdfname, "pdf");
    canT0->Print(pdfname,"pdf");
    canTOT->Print(pdfname,"pdf");
    canTDif->Print(pdfname,"pdf");
    canQAve->Print(pdfname,"pdf");
    canQSum->Print(pdfname,"pdf");
    canMLR->Print(pdfname,"pdf");
    canQLvsQR->Print(pdfname,"pdf");
    canQRvsAR->Print(pdfname,"pdf");
    canQLvsAL->Print(pdfname,"pdf");
    canMult->Print(pdfname,"pdf");
    canMapping->Print(pdfname,"pdf");
    canQDCControlOnlyL->Print(pdfname,"pdf");
    canQDCControlOnlyR->Print(pdfname,"pdf");
    canTimeDiffL->Print(pdfname,"pdf");
    canTimeDiffR->Print(pdfname,"pdf");
    canQDCcutL->Print(pdfname,"pdf");
    canQDCcutR->Print(pdfname,"pdf");
    canQeQc2DL->Print(pdfname,"pdf");
    canQeQc2DR->Print(pdfname,"pdf");
    canQDCHighDepL->Print(pdfname,"pdf");
    canQDCHighDepR->Print(pdfname,"pdf");

    
    
    
    auto canProfile = new TCanvas("Profile","Profile",1000,500);
    canProfile->SetLeftMargin(0.4);
    canProfile->Divide(2);
    canProfile->cd(1);
    gPad->SetLeftMargin(0.15);

    hProfL->Draw();
    hProfL->SetLineWidth(4);
    hProfL->SetLineColor(kRed);
    hProfL->SetAxisRange(3,8,"X");
    
    std::cout << "hProfL # of entries: " << hProfL->GetBinEntries(0) << std::endl;
    std::cout << "hProfL # of entries: " << hProfL->GetBinEntries(1) << std::endl;
    std::cout << "hProfL # of entries: " << hProfL->GetBinEntries(2) << std::endl;
    std::cout << "hProfL # of entries: " << hProfL->GetBinEntries(3) << std::endl;
    std::cout << "hProfL # of entries: " << hProfL->GetBinEntries(4) << std::endl;
    std::cout << "hProfL # of entries: " << hProfL->GetBinEntries(5) << std::endl;    
    
    std::cout << "hProfL bin content: " << hProfL->GetBinContent(0) << std::endl;
    std::cout << "hProfL bin content: " << hProfL->GetBinContent(1) << std::endl;
    std::cout << "hProfL bin content: " << hProfL->GetBinContent(2) << std::endl;
    std::cout << "hProfL bin content: " << hProfL->GetBinContent(3) << std::endl;
    std::cout << "hProfL bin content: " << hProfL->GetBinContent(4) << std::endl;
    std::cout << "hProfL bin content: " << hProfL->GetBinContent(5) << std::endl;
    
        std::cout << "hProfL # of entries: " << hProfR->GetBinEntries(0) << std::endl;
    std::cout << "hProfL # of entries: " << hProfR->GetBinEntries(1) << std::endl;
    std::cout << "hProfL # of entries: " << hProfR->GetBinEntries(2) << std::endl;
    std::cout << "hProfL # of entries: " << hProfR->GetBinEntries(3) << std::endl;
    std::cout << "hProfL # of entries: " << hProfR->GetBinEntries(4) << std::endl;
    std::cout << "hProfL # of entries: " << hProfR->GetBinEntries(5) << std::endl;    
    
    std::cout << "hProfL bin content: " << hProfR->GetBinContent(0) << std::endl;
    std::cout << "hProfL bin content: " << hProfR->GetBinContent(1) << std::endl;
    std::cout << "hProfL bin content: " << hProfR->GetBinContent(2) << std::endl;
    std::cout << "hProfL bin content: " << hProfR->GetBinContent(3) << std::endl;
    std::cout << "hProfL bin content: " << hProfR->GetBinContent(4) << std::endl;
    std::cout << "hProfL bin content: " << hProfR->GetBinContent(5) << std::endl;

    canProfile->cd(2);
        gPad->SetLeftMargin(0.15);

    hProfR->Draw();
    hProfR->SetLineWidth(4);
    hProfR->SetLineColor(kRed);
    hProfR->SetAxisRange(3,8,"X");

    
    auto canProfileNorm = new TCanvas("ProfileNorm","ProfileNorm",1000,500);
    
    canProfileNorm->SetLeftMargin(0.4);
    canProfileNorm->Divide(2);
    canProfileNorm->cd(1);
    gPad->SetLeftMargin(0.15);
    hProfNormL->Draw();
    hProfNormL->SetLineWidth(4);
    hProfNormL->SetLineColor(kRed);
    hProfNormL->SetAxisRange(3,8,"X");
    
    canProfileNorm->cd(2);
    gPad->SetLeftMargin(0.15);
    hProfNormR->Draw();
    hProfNormR->SetLineWidth(4);
    hProfNormR->SetLineColor(kRed);
    hProfNormR->SetAxisRange(3,8,"X");
    
    
   // canProfile->SaveAs("/scratch/gccb/data/202111/pdf/profile_2021_11_28_05_42_50cm.pdf");
    canProfile->Print(pdfname,"pdf");
    canProfileNorm->Print(pdfname,"pdf");
    canMult->Print(Form("%s)",pdfname.Data()),"pdf");
	return 0;
}
