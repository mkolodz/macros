#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include "FitterFactory.h"
#include "TString.h"

Bool_t Init(FitterFactory & ff, TH1D * hL, TH1D * hR, TH1D* hAv)
{
    TString functions = "gaus(0) pol0(3)+[4]*TMath::Exp((x-[5])*[6])";
   
    TF1 *fun_exp_tmp = new TF1("fun_exp_tmp", "[0]*TMath::Exp((x-[1])*[2])");
    
    TSpectrum *spec= new TSpectrum(10);
    
    //----- left
    Int_t npeaksL = spec->Search(hL, 10, "", 0.5);
    Double_t *peaksL = spec->GetPositionX();
    Double_t peakL = TMath::MaxElement(npeaksL, peaksL);

    Double_t xmin = peakL - 200;
    Double_t xmax = peakL +300;
    Double_t p_const = hL->GetBinContent(hL->FindBin(peakL));
    Double_t mean = peakL;
    Double_t sigma = 30;
    
    Double_t exp_fit_min = xmin;
    Double_t exp_fit_max = xmin + 100;

    hL->Fit("fun_exp_tmp", "", "", exp_fit_min, exp_fit_max);

    Double_t par3 = 20;
    Double_t par4 = fun_exp_tmp->GetParameter(0);
    Double_t par5 = fun_exp_tmp->GetParameter(1);   
    Double_t par6 = fun_exp_tmp->GetParameter(2);
    
    auto pL = std::make_unique<HistogramFitParams>(hL->GetName(), "gaus(0)", "pol0(3)+[4]*TMath::Exp((x-[5])*[6])", xmin, xmax);
    pL->setParam(0, p_const, ParamValue::FitMode::Free);
    pL->setParam(1, mean, ParamValue::FitMode::Free);
    pL->setParam(2, sigma, ParamValue::FitMode::Free);
    pL->setParam(3, par3, ParamValue::FitMode::Free);
    pL->setParam(4, par4, ParamValue::FitMode::Free);
    pL->setParam(5, par5, ParamValue::FitMode::Free);
    pL->setParam(5, par6, ParamValue::FitMode::Free);
        
    //----- right
    Int_t npeaksR = spec->Search(hR, 12, "", 0.5);
    Double_t *peaksR = spec->GetPositionX();
    Double_t peakR = TMath::MaxElement(npeaksR, peaksR);
    
    xmin = peakR - 200;
    xmax = peakR + 300;
    p_const = hR->GetBinContent(hR->FindBin(peakR));
    mean = peakR;
    sigma = 30;
    
    exp_fit_min = xmin;
    exp_fit_max = xmin + 100;
    
    hR->Fit("fun_exp_tmp", "", "", exp_fit_min, exp_fit_max);

    par3 = 20;
    par4 = fun_exp_tmp->GetParameter(0);
    par5 = fun_exp_tmp->GetParameter(1);   
    par6 = fun_exp_tmp->GetParameter(2);

    auto pR = std::make_unique<HistogramFitParams>(hR->GetName(), "gaus(0)", "pol0(3)+[4]*TMath::Exp((x-[5])*[6])", xmin, xmax);
    pR->setParam(0, p_const, ParamValue::FitMode::Free);
    pR->setParam(1, mean, ParamValue::FitMode::Free);
    pR->setParam(2, sigma, ParamValue::FitMode::Free);
    pR->setParam(3, par3, ParamValue::FitMode::Free);
    pR->setParam(4, par4, ParamValue::FitMode::Free);
    pR->setParam(5, par5, ParamValue::FitMode::Free);
    pR->setParam(5, par6, ParamValue::FitMode::Free);
           
    //----- average
    Int_t npeaksAv = spec->Search(hAv, 12, "", 0.5);
    Double_t *peaksAv = spec->GetPositionX();
    Double_t peakAv = TMath::MaxElement(npeaksAv, peaksAv);
    
    xmin = peakAv - 200;
    xmax = peakAv + 300;
    p_const = hAv->GetBinContent(hR->FindBin(peakAv));
    mean = peakAv;
    sigma = 30;
    
    exp_fit_min = xmin;
    exp_fit_max = xmin + 100;
    
    hAv->Fit("fun_exp_tmp", "", "", exp_fit_min, exp_fit_max);

    par3 = 20;
    par4 = fun_exp_tmp->GetParameter(0);
    par5 = fun_exp_tmp->GetParameter(1);   
    par6 = fun_exp_tmp->GetParameter(2);

    auto pAV = std::make_unique<HistogramFitParams>(hAv->GetName(), "gaus(0)", "pol0(3)+[4]*TMath::Exp((x-[5])*[6])", xmin, xmax);
    pAV->setParam(0, p_const, ParamValue::FitMode::Free);
    pAV->setParam(1, mean, ParamValue::FitMode::Free);
    pAV->setParam(2, sigma, ParamValue::FitMode::Free);
    pAV->setParam(3, par3, ParamValue::FitMode::Free);
    pAV->setParam(4, par4, ParamValue::FitMode::Free);
    pAV->setParam(5, par5, ParamValue::FitMode::Free);
    pAV->setParam(5, par6, ParamValue::FitMode::Free);

    ff.insertParameters(std::move(pL));
    ff.insertParameters(std::move(pR));
    ff.insertParameters(std::move(pAV));

    return kTRUE;
}

void AddEntries(FitterFactory & ff, TString hname_l, TString hname_r, TString hname_av)
{
    auto pL = ff.findParams("Q_M0L0F0L");
    if (pL)
        ff.insertParameters(pL->clone(hname_l));

    auto pR = ff.findParams("Q_M0L0F0R");
    if (pR)
        ff.insertParameters(pR->clone(hname_r));

    auto pA = ff.findParams("QAve_M0L0F0");
    if (pA)
         ff.insertParameters(pA->clone(hname_av));
}

Bool_t attenuation(TString logname, Int_t mod, Int_t lay, Int_t fib)
{
    
    Int_t colL = kBlack;
    Int_t colR = kRed;
    Int_t colAv = kBlue + 2;
    
    gStyle->SetPadLeftMargin(0.15);
    
    //parse logname
    std::string logname_str = std::string(logname);
    int         nletters    = logname_str.length();
    char        letters[nletters];
    strcpy(letters, logname_str.c_str());
  
    int istop = -1;
  
    for (int i = nletters; i > 0; i--)
    {
        if (letters[i] == '/')
        {
            istop = i;
            break;
        }
    }

  int istart = 0;
   
  if (istop == -1)
  {
      std::cerr << "##### Error in main()!" << std::endl;
      std::cerr << "Cannot interpret log name!" << std::endl;
      return 0;
  }

  std::string data_path = std::string(&letters[istart], &letters[istop + 1]) + "../";
  std::string file_name = std::string(&letters[istop+1], &letters[nletters-4]);
  std::cout << "data path: " << data_path << std::endl;
  std::cout << "file name: " << file_name << std::endl;
    
    //----- reading logfile
    std::ifstream logfile(logname);
    
    if(!logfile.is_open())
    {
        std::cout << "Couldn't open log file " << logname << std::endl;
        return kFALSE;
    }
    
    std::string dummy;
    std::string line;
    TString  dir_names[50];
    Double_t positions[50];
    Int_t    npoints = 0;
    
    while(logfile.good())
    {
     getline(logfile, line);
     getline(logfile, line);
     getline(logfile, line);
     logfile >> dir_names[npoints];
     logfile >> dummy >> dummy >> positions[npoints];
     std::cout << dir_names[npoints] << "\t" << positions[npoints] << std::endl;
     getline(logfile, line);
     getline(logfile, line);
     npoints++;
   }
   
   std::cout << "npoints: " << npoints << std::endl;
   
   //----- accessing data & filling histograms
   TFile *file = nullptr;
//    TTree *tree = nullptr;
   
   std::vector <TH1D*> hL;
   std::vector <TH1D*> hR;
   std::vector <TH1D*> hQAve;
   std::vector <TH1D*> hRat;
   
   for(Int_t i=0; i<npoints; i++)
   {
       TString full_root_name = TString(data_path) + dir_names[i] + "/sifi_results_HISTOS.root";
       std::cout << full_root_name << std::endl;
       file = new TFile(full_root_name, "READ");
       
//        tree = (TTree*)file->Get("S");
       
//        tree->Draw(Form("SDDSamples.data.signal_l.fPE>>htemp_l%i(250, 0, 2500)", i),
//                   Form("SDDSamples.data.module == %i && "
//                        "SDDSamples.data.layer == %i && "
//                        "SDDSamples.data.fiber == %i", mod, lay, fib), "");
//        hL.push_back((TH1D*)gROOT->FindObjectAny(Form("htemp_l%i", i))->Clone(Form("hL_%.1fmm", positions[i])));
       
       TString dir = Form("/Module%i/Layer%i/Fiber%i/", mod, lay, fib);
       TString hname = Form("Q_M%iL%iF%iL", mod, lay, fib);
       std::cout << dir+hname << std::endl;
       file->cd(dir);
       hL.push_back((TH1D*)gDirectory->FindObjectAny(hname));
       if(hL[i] == nullptr)
       {
           std::cerr << "Could not get histogram: " << dir+hname << std::endl;
           std::abort();
       }
       
//        tree->Draw(Form("SDDSamples.data.signal_r.fPE>>htemp_r%i(250, 0, 2500)", i),
//                   Form("SDDSamples.data.module == %i && "
//                        "SDDSamples.data.layer == %i && "
//                        "SDDSamples.data.fiber == %i", mod, lay, fib), "");
//        hR.push_back((TH1D*)gROOT->FindObjectAny(Form("htemp_r%i", i))->Clone(Form("hR_%.1fmm", positions[i])));
       
       dir = Form("/Module%i/Layer%i/Fiber%i/", mod, lay, fib);
       hname = Form("Q_M%iL%iF%iR", mod, lay, fib);
       std::cout << dir+hname << std::endl;
       file->cd(dir);
       hR.push_back((TH1D*)gDirectory->FindObjectAny(hname));
       if(hR[i] == nullptr)
       {
           std::cerr << "Could not get histogram: " << dir+hname << std::endl;
           std::abort();
       }
       
//        tree->Draw(Form("sqrt(SDDSamples.data.signal_r.fPE*SDDSamples.data.signal_l.fPE)>>htemp_ave%i(250, 0, 2500)", i),
//                   Form("SDDSamples.data.module == %i && "
//                        "SDDSamples.data.layer == %i && "
//                        "SDDSamples.data.fiber == %i", mod, lay, fib), "");
//        hQAve.push_back((TH1D*)gROOT->FindObjectAny(Form("htemp_ave%i", i))->Clone(Form("hAve_%.1fmm", positions[i])));
       
       dir = Form("/Module%i/Layer%i/Fiber%i/", mod, lay, fib);
       hname = Form("QAve_M%iL%iF%i", mod, lay, fib);
       std::cout << dir+hname << std::endl;
       file->cd(dir);
       hQAve.push_back((TH1D*)gDirectory->FindObjectAny(hname));
       if(hQAve[i] == nullptr)
       {
           std::cerr << "Could not get histogram: " << dir+hname << std::endl;
           std::abort();
       }
       
//        tree->Draw(Form("log(sqrt(SDDSamples.data.signal_l.fPE/SDDSamples.data.signal_r.fPE))"
//                         ">>htemp_rat%i(200, -5, 5)", i),
//                   Form("SDDSamples.data.module == %i && "
//                        "SDDSamples.data.layer == %i && "
//                        "SDDSamples.data.fiber == %i", mod, lay, fib), "");
//        hRat.push_back((TH1D*)gROOT->FindObjectAny(Form("htemp_rat%i", i))->Clone(Form("hrat_%.1fmm", positions[i])));
       
       dir = Form("/Module%i/Layer%i/Fiber%i/", mod, lay, fib);
       hname = Form("MLR_M%iL%iF%i", mod, lay, fib);
       std::cout << dir+hname << std::endl;
       file->cd(dir);
       hRat.push_back((TH1D*)gDirectory->FindObjectAny(hname));
       if(hRat[i] == nullptr)
       {
           std::cerr << "Could not get histogram: " << dir+hname << std::endl;
           std::abort();
       }
   }
   
   //----- fitting MLR and calculating attenuation length
   TF1 *fgaus = new TF1("fgaus", "gaus" ,-5, 5);
   
   TGraphErrors *gMLR = new TGraphErrors(npoints);
   gMLR->SetMarkerStyle(4);
   gMLR->SetName("gMLR");
   gMLR->SetTitle("M_{LR} vs. source position");
   gMLR->GetXaxis()->SetTitle("source position [mm]");
   gMLR->GetYaxis()->SetTitle("M_{LR}");
  
   for(Int_t i=0; i<npoints; i++)
   {
        hRat[i]->Fit("fgaus","Q");
        gMLR->SetPoint(i, positions[i], fgaus->GetParameter(1));
        gMLR->SetPointError(i, 1.5, fgaus->GetParError(1));
   }
   
   TF1 *fpol1 = new TF1("fpol1", "pol1", positions[0], positions[npoints-1]);
//    fpol1->SetParameters(0.1,1E-3);
   fpol1->SetParameters(0.1,-1E-3);
   gMLR->Fit("fpol1", "R");
   
   Double_t att =  1./fpol1->GetParameter(1);
   Double_t att_err = fpol1->GetParError(1)/pow(fpol1->GetParameter(1),2);
   
   std::cout << "\n\nAttenuation length (MLR method): " << att << " mm +/- "
             << att_err << " mm" << std::endl;

   //----- attenuation for separate channels 
   TGraphErrors *gL = new TGraphErrors(npoints);
   gL->SetMarkerStyle(4);
   gL->SetMarkerColor(colL);
   gL->SetLineColor(colL);
   gL->SetName("gL");
   gL->SetTitle("L and R vs. source position");
   gL->GetXaxis()->SetTitle("source position [mm]");
   gL->GetYaxis()->SetTitle("511 keV peak position");
  
   TGraphErrors *gR = new TGraphErrors(npoints);
   gR->SetMarkerStyle(4);
   gR->SetMarkerColor(colR);
   gR->SetLineColor(colR);
   gR->SetName("gR");
   
   TGraphErrors *gAv = new TGraphErrors(npoints);
   gAv->SetMarkerStyle(4);
   gAv->SetMarkerColor(colAv);
   gAv->SetName("gAv");
   
   TGraphErrors *gERL = new TGraphErrors(npoints);
   gERL->SetMarkerStyle(4);
   gERL->SetMarkerColor(colL);
   gERL->SetLineColor(colL);
   gERL->SetName("gERL");
   gERL->SetTitle("energy resolution vs. source position");
   gERL->GetXaxis()->SetTitle("source position [mm]");
   gERL->GetYaxis()->SetTitle("energy resolution [%]");
   
   TGraphErrors *gERR = new TGraphErrors(npoints);
   gERR->SetMarkerStyle(4);
   gERR->SetMarkerColor(colR);
   gERR->SetLineColor(colR);
   gERR->SetName("gERR");
   gERR->SetTitle("energy resolution vs. source position");
   gERR->GetXaxis()->SetTitle("source position [mm]");
   gERR->GetYaxis()->SetTitle("energy resolution [%]");
   
   TGraphErrors *gERAve = new TGraphErrors(npoints);
   gERAve->SetMarkerStyle(4);
   gERAve->SetMarkerColor(colAv);
   gERAve->SetName("gERAve");
   gERAve->SetTitle("energy resolution (from average spectra) vs. source position");
   gERAve->GetXaxis()->SetTitle("source position [mm]");
   gERAve->GetYaxis()->SetTitle("energy resolution [%]");

   for(Int_t i=0; i<npoints; i++)
   {
          FitterFactory fitter;
          fitter.initFactoryFromFile((TString(data_path) + dir_names[i] + "/fitconfig.txt").Data(),
                              (TString(data_path) + dir_names[i] + "/fitparams.out").Data());

//        HistFitParams* hfp = fitter.findParams(hL[i]->GetName());

       if (!fitter.count())
       {
            Init(fitter, hL[i], hR[i], hQAve[i]);
       }
       else
       {
            HistogramFitParams* hfp2 = fitter.findParams(hL[i]->GetName());

            if (!hfp2)
                AddEntries(fitter, hL[i]->GetName(), hR[i]->GetName(), hQAve[i]->GetName());
       }
        
        std::cout << hL[i]->GetName() << std::endl;
        std::cout << hR[i]->GetName() << std::endl;
        std::cout << hQAve[i]->GetName() << std::endl;
        
        //----- fit left
        auto histFPL = fitter.findParams(hL[i]->GetName());
//         printf("fl = %d for %s\n", fl, hL[i]->GetName());
        fitter.fit(histFPL, hL[i]);
//         histFPL->update();
        TF1 * fittedL = (TF1*)histFPL->function_sum.Clone();
        
        gL->SetPoint(i, positions[i], fittedL->GetParameter(1));
        gL->SetPointError(i, 1.5, fittedL->GetParError(1));
        
        Double_t ER = (fabs(fittedL->GetParameter(2))/fittedL->GetParameter(1))*100;
        Double_t ER_er = ER * sqrt(pow(fittedL->GetParError(1), 2) /
                                  pow(fittedL->GetParameter(1), 2) +
                                  pow(fittedL->GetParError(2), 2) /
                                  pow(fittedL->GetParameter(2), 2));
        gERL->SetPoint(i, positions[i], ER);
        gERL->SetPointError(i, 1.5,  ER_er);
        
        std::cout << "ER (left) = " << ER << " +/- " << ER_er << std::endl;
        
        //----- fit right
        auto histFPR = fitter.findParams(hR[i]->GetName());
//         printf("fr = %d for %s\n", fr, hR[i]->GetName());
        fitter.fit(histFPR, hR[i]);
//         histFPR->update();
        TF1 * fittedR = (TF1*)histFPR->function_sum.Clone();
                
        gR->SetPoint(i, positions[i], fittedR->GetParameter(1));
        gR->SetPointError(i, 1.5, fittedR->GetParError(1));
        
        ER = (fabs(fittedR->GetParameter(2))/fittedR->GetParameter(1))*100;
        ER_er = ER * sqrt(pow(fittedR->GetParError(1), 2) /
                          pow(fittedR->GetParameter(1), 2) +
                          pow(fittedR->GetParError(2), 2) /
                          pow(fittedR->GetParameter(2), 2));
        
        std::cout << "ER (right) = " << ER << " +/- " << ER_er << std::endl;
        
        gERR->SetPoint(i, positions[i], ER);
        gERR->SetPointError(i, 1.5, ER_er);
        
        //----- fit average
        auto histFPAv = fitter.findParams(hQAve[i]->GetName());
//         printf("fav = %d for %s\n", fav, hQAve[i]->GetName());
        fitter.fit(histFPAv, hQAve[i]);
//         histFPAv->update();
        
        TF1 * fittedAv = (TF1*)histFPAv->function_sum.Clone();
                
        gAv->SetPoint(i, positions[i], fittedAv->GetParameter(1));
        gAv->SetPointError(i, 1.5, fittedAv->GetParError(1));
        
        ER = (fabs(fittedAv->GetParameter(2))/fittedAv->GetParameter(1))*100;
        ER_er = ER * sqrt(pow(fittedAv->GetParError(1), 2) /
                          pow(fittedAv->GetParameter(1), 2) +
                          pow(fittedAv->GetParError(2), 2) /
                          pow(fittedAv->GetParameter(2), 2));
        
        std::cout << "ER (averaged) = " << ER << " +/- " << ER_er << std::endl;
        
        gERAve->SetPoint(i, positions[i], ER);
        gERAve->SetPointError(i, 1.5, ER_er);
        
        std::cout << "position: " << positions[i] << std::endl;
        std::cout << "\n\n" << std::endl;

        fitter.exportFactoryToFile();
   }
   
   //----- drawing histograms & graphs
   TCanvas *canSpec = new TCanvas("canSpec", "canSpec", 1500, 1000);
   canSpec->DivideSquare(npoints);
   
   TCanvas *canRat = new TCanvas("canRat", "canrat", 1500, 1000);
   canRat->DivideSquare(npoints);
   
   TCanvas *canAve = new TCanvas("canAve", "canave", 1500, 1000);
   canAve->DivideSquare(npoints);
      
   TCanvas *tmp_can_l = new TCanvas("tmp_can_l", "tmp_can_l", 1500, 1000);
   tmp_can_l->DivideSquare(npoints);
   
   TCanvas *tmp_can_r = new TCanvas("tmp_can_r", "tmp_can_r", 1500, 1000);
   tmp_can_r->DivideSquare(npoints);
   
   for(Int_t i=0; i<npoints; i++)
   {
       tmp_can_l->cd(i+1);
       hL[i]->SetTitle(Form("Q_{L} %.1f mm", positions[i]));
       hL[i]->GetXaxis()->SetTitle("charge [a.u]");
       hL[i]->GetYaxis()->SetTitle("counts");
       gPad->SetGrid(1, 1);
       hL[i]->DrawClone();
    
       auto fun_sum_l = std::unique_ptr<TF1>(hL[i]->GetFunction("f_" + TString(hL[i]->GetName())));
       auto fun_gauss = std::make_unique<TF1>("fun_gauss", "gaus", 0, 6000);
       fun_gauss->SetLineColor(kBlue);
       fun_gauss->SetParameter(0, fun_sum_l->GetParameter(0));
       fun_gauss->SetParameter(1, fun_sum_l->GetParameter(1));
       fun_gauss->SetParameter(2, fun_sum_l->GetParameter(2));
       
       auto fun_pol0 = std::make_unique<TF1>("fun_pol0", "pol0", 0, 6000);
       fun_pol0->SetLineColor(kMagenta);
       fun_pol0->SetParameter(0, fun_sum_l->GetParameter(3));
       
       auto fun_exp = std::make_unique<TF1>("fun_exp", "[0]*TMath::Exp((x-[1])*[2])", 0, 6000);
       fun_exp->SetLineColor(kGreen+2);
       fun_exp->SetParameter(0, fun_sum_l->GetParameter(4));
       fun_exp->SetParameter(1, fun_sum_l->GetParameter(5));
       fun_exp->SetParameter(2, fun_sum_l->GetParameter(6));
       
       fun_gauss->DrawClone("same");
       fun_pol0->DrawClone("same");
       fun_exp->DrawClone("same");
       
       tmp_can_r->cd(i+1);
       hR[i]->SetTitle(Form("Q_{R} %.1f mm", positions[i]));
       hR[i]->GetXaxis()->SetTitle("charge [a.u]");
       hR[i]->GetYaxis()->SetTitle("counts");
       gPad->SetGrid(1, 1);
       hR[i]->DrawClone();
    
       auto fun_sum_r = std::unique_ptr<TF1>(hR[i]->GetFunction("f_" + TString(hR[i]->GetName())));

       fun_gauss->SetParameter(0, fun_sum_r->GetParameter(0));
       fun_gauss->SetParameter(1, fun_sum_r->GetParameter(1));
       fun_gauss->SetParameter(2, fun_sum_r->GetParameter(2));
       
       fun_pol0->SetParameter(0, fun_sum_r->GetParameter(3));
       
       fun_exp->SetParameter(0, fun_sum_r->GetParameter(4));
       fun_exp->SetParameter(1, fun_sum_r->GetParameter(5));
       fun_exp->SetParameter(2, fun_sum_r->GetParameter(6));
       
       fun_gauss->DrawClone("same");
       fun_pol0->DrawClone("same");
       fun_exp->DrawClone("same");
       
       canSpec->cd(i+1);
       gPad->SetGrid(1, 1);
       hL[i]->SetStats(0);
       hL[i]->SetTitle(Form("Q_{L} and Q_{R} %.1f mm", positions[i]));
       hL[i]->GetYaxis()->SetRangeUser(0, 500);
       hL[i]->SetLineColor(colL);
       hL[i]->Draw();
       hR[i]->SetStats(0);
       hR[i]->SetLineColor(colR);
       hR[i]->Draw("same");
       fun_sum_l->DrawClone("same");
       fun_sum_r->DrawClone("same");
       
       canRat->cd(i+1);
       gPad->SetGrid(1, 1);
       hRat[i]->SetTitle(Form("log(sqrt(L/R)) %.1f mm", positions[i]));
       hRat[i]->GetXaxis()->SetTitle("log(sqrt(R/L))");
       hRat[i]->GetYaxis()->SetTitle("counts");
       hRat[i]->Draw();
       
       canAve->cd(i+1);
       gPad->SetGrid(1,1);
       hQAve[i]->SetTitle(Form("#sqrt(Q_{L} Q_{R}) %.1f mm", positions[i]));
       hQAve[i]->GetXaxis()->SetTitle("sqrt(RL) [a.u.]");
       hQAve[i]->GetYaxis()->SetTitle("counts");
       hQAve[i]->Draw();
       
       auto fun_sum_ave = std::unique_ptr<TF1>(hQAve[i]->GetFunction("f_" + TString(hQAve[i]->GetName())));

       fun_gauss->SetParameter(0, fun_sum_ave->GetParameter(0));
       fun_gauss->SetParameter(1, fun_sum_ave->GetParameter(1));
       fun_gauss->SetParameter(2, fun_sum_ave->GetParameter(2));
       
       fun_pol0->SetParameter(0, fun_sum_ave->GetParameter(3));
       
       fun_exp->SetParameter(0, fun_sum_ave->GetParameter(4));
       fun_exp->SetParameter(1, fun_sum_ave->GetParameter(5));
       fun_exp->SetParameter(2, fun_sum_ave->GetParameter(6));
       
       fun_sum_ave->DrawClone("same");
       fun_gauss->DrawClone("same");
       fun_pol0->DrawClone("same");
       fun_exp->DrawClone("same");
   }
   
   canSpec->cd(1);
   
   TLatex legend;
   legend.SetNDC(kTRUE);
   legend.SetTextColor(colL);
   legend.DrawLatex(0.6, 0.80, "LEFT");
   legend.SetTextColor(colR);
   legend.DrawLatex(0.6, 0.75, "RIGHT");
   
   TCanvas *canAtt = new TCanvas("canAtt", "canAtt", 1500, 500);
   canAtt->Divide(3, 1);
   
   canAtt->cd(1);
   gPad->SetGrid(1, 1);
   gMLR->Draw("AP");
   
   TLatex text;
   text.SetNDC(kTRUE);
   text.DrawLatex(0.2, 0.8, Form("#lambda = %.2f #pm %.2f mm", att, att_err));
   
   canAtt->cd(2);
   gPad->SetGrid(1, 1);
   
   Double_t x0, y0, x1, y1;
   gL->GetPoint(npoints-1, x0, y0);
   gR->GetPoint(0, x1, y1);
   Double_t gmin = y0 < y1 ? y0 : y1;
   gL->GetPoint(0, x0, y0);
   gR->GetPoint(npoints-1, x1, y1);
   Double_t gmax = y0 > y1 ? y0 : y1;
   
   gL->GetYaxis()->SetRangeUser(500,1400);
   //gL->GetYaxis()->SetRangeUser(0.8*gmin, 1.2*gmax);
//    gL->SetMinimum(0.8*gmin);
//    gL->SetMaximum(1.2*gmax);
   gL->Draw("AP");
   gR->Draw("P");
   gAv->Draw("P");
   
   legend.SetTextColor(colL);
   legend.DrawLatex(0.6, 0.80, "LEFT");
   legend.SetTextColor(colR);
   legend.DrawLatex(0.6, 0.75, "RIGHT");
   legend.SetTextColor(colAv);
   legend.DrawLatex(0.6, 0.70, "AVERAGED");
   
   canAtt->cd(3);
   gPad->SetGrid(1, 1);
   gERL->GetPoint(npoints-1, x0, y0);
   gERR->GetPoint(0, x1, y1);
   gmin = y0 < y1 ? y0 : y1;
   gERL->GetPoint(0, x0, y0);
   gERR->GetPoint(npoints-1, x1, y1);
   gmax = y0 > y1 ? y0 : y1;
   
//    gERL->GetYaxis()->SetRangeUser(0.8*gmin, 1.2*gmax);
   gERL->GetYaxis()->SetRangeUser(5, 20);
   gERL->Draw("AP");
   gERR->Draw("P");
   gERAve->Draw("P");
   
   
   //----- saving canvases
   TString file_name_root = file_name + Form("_M%iL%iF%i_attenuation.root", mod, lay, fib);
   TString full_out_root_name = TString(data_path) + "hist/" + file_name_root;
   std::cout << "output ROOT file name: " << full_out_root_name << std::endl;
   
   TFile *outfile = new TFile(full_out_root_name, "RECREATE");
   canSpec->Write();
   tmp_can_l->Write();
   tmp_can_r->Write();
   canRat->Write();
   canAtt->Write();
   canAve->Write();
   outfile->Close();
   
    //output lambda to file
    //TFile *outlambdafile = new TFile("attenuation_lengths.txt", "UPDATE");
    char * log_name="lambdas.txt";
    std::string a,b;
    std::ifstream log(log_name);
    while(log.good())
    {
        getline(log,a);
        log >> b;
        cout << "a: " << a << "b: " << b << endl;
    }
    //log << mod << " " << lay << " " << fib << " " << att << std::endl;
    
   
   //-----exporting
   
   TString pdf_name = TString(data_path) + "pdf/" + TString(file_name) + Form("_M%iL%iF%i", mod, lay, fib) + "_attenuation.pdf";
   canSpec->Print(Form("%s(", pdf_name.Data()), "pdf");
   canAve->Print(pdf_name, "pdf");
   canRat->Print(pdf_name, "pdf");
   canAtt->Print(Form("%s)", pdf_name.Data()), "pdf");
   
   return kTRUE;
}
