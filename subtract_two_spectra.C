double aminusb(TString name="single_fibel_hamamatsu_1.root", TString name2="single_fibel_no_source_hamamatsu_1.root" ){
  
TFile *f1 = new TFile(name,"READ");
TFile *f2 = new TFile(name2,"READ");

TH1D *h0 = new TH1D("h0", ";tot [ns];Counts", 100, 340, 1100); //single fiber with source

TH1D *h2 = new TH1D("h2", ";tot [ns];Counts", 100, 340, 1100);//single fiber no source


TH1D *t1 = new TH1D("t1", "", 10800, 0, 10800);
TH1D *t2 = new TH1D("t2", "", 10800, 0, 10800);


/*if (f1==0)
   
    {std::cout<<"Error opening file "<<name << std::endl;
      return -100;
    }*/
    
TTree *tree1 = (TTree *)f1->Get("UserObjects/GeneralInfo/tree");
TTree *tree2 = (TTree *)f2->Get("UserObjects/GeneralInfo/tree");

/*if (tree1==0)
    
    {std::cout<<"Tree is not found in: "<<name<< std::endl;
      return -100;
    }
   */
gStyle->SetOptStat(0);
gStyle->SetOptFit(0010);
gStyle-> SetStatY(0.75); 
gStyle-> SetStatX(0.97);
gStyle-> SetStatW(0.11);
gStyle-> SetStatH(0.1);
gStyle->SetOptDate(1);

tree1->Draw("rising[20]/10e9>>t1");
tree2->Draw("rising[20]/10e9>>t2");  
    
//SCALE    
    
double t0=-1;
double tk=-1;
int bin=1;

while(t1->GetBinContent(bin)==0)
    {bin++;
    }
    t0=t1->GetBinCenter(bin);
    bin=t1->GetXaxis()->GetNbins();
    
while(t1->GetBinContent(bin)==0)
    {bin--;
    }
    tk=t1->GetBinCenter(bin);
	std::cout<<"tk-to "<<tk-t0<< std::endl;
	
//double cps=t1->GetEntries()/(tk-t0);

double t02=-1;
double tk2=-1;
int bin2=1;

while(t2->GetBinContent(bin2)==0)
    {bin2++;
    }
    t02=t2->GetBinCenter(bin2);
    bin2=t2->GetXaxis()->GetNbins();
    
while(t2->GetBinContent(bin2)==0)
    {bin2--;
    }
    tk2=t2->GetBinCenter(bin2);

//std::cout<<"cps "<<cps<< std::endl;


tree1->Draw("falling[20]-rising[20]>>h0");
h0->Scale(1./(tk-t0));
h0->Sumw2(kFALSE);

tree2->Draw("falling[20]-rising[20]>>h2");
h2->Scale(1./(tk2-t02));
h2->Sumw2(kFALSE);



if(name.Contains("no_source"))
    return 0;
    
//HO->H1
    
TH1F* h1 = (TH1F*) h0->Clone("clone");
h1->Add(h2,-1);

h0->Draw("");
h0->SetLineColor(kViolet);
h2->Draw("same");
h2->SetLineColor(kGreen);
h1->Draw("same");


  
TSpectrum* spec = new TSpectrum(10,25);
    TH1F* hBg = (TH1F*)spec->Background(h1,15,"same nosmoothing");
    
    TH1F* hPeaksOnly = (TH1F*) h1->Clone("hpeaks");
    hPeaksOnly->Add(hBg,-1);
    hPeaksOnly->SetLineColor(kBlack);
    hPeaksOnly->Draw("same");
  

  
Int_t fNuOfPeaks = spec->Search(hPeaksOnly,5,"background,new,nodraw",0.1);
Double_t* fxPeaks = spec->GetPositionX();
Double_t* fyPeaks = spec->GetPositionY();
cout<<"Found "<<fNuOfPeaks<<" peaks at:"<<endl;
  
for(int i=0; i<fNuOfPeaks; i++)
    cout<<fxPeaks[i]<<",\t "<<fyPeaks[i]<<endl;
    hBg->SetLineColor(kGray);
    hBg->Draw("same");
  
  
if(fNuOfPeaks==0)
    return -100;
    
/*if(fNuOfPeaks!=0)
    return ; */
  
double x511 = fxPeaks[fNuOfPeaks-1];
double y511 = fyPeaks[fNuOfPeaks-1];
double a= 0.85*x511;
double b= 1.3*x511;
hBg->Fit("expo","Q","",a, b);
  
double lambda = hBg->GetFunction("expo")->GetParameter(1);
double A0 = hBg->GetFunction("expo")->GetParameter(0);
hBg->GetFunction("expo")->Delete();
  
  
//FITTING
TF1 *fun = new TF1("linefit","expo(1)+[0]+gaus(3)",0,1200);
    fun->SetLineColor(kRed);
    fun->SetLineWidth(2);
    fun->SetParNames("Yoffset","A_{0}","#lambda","A","Mean","#sigma");
    fun->SetParameter(0,0);
    fun->SetParameter(1, A0);
    fun->SetParameter(2,lambda);
    fun->SetParameter(3,y511);
    fun->SetParameter(4,x511);
    fun->SetParameter(5,25);
  
   h1->Fit(fun,"","same",a,b);
   h1->SetLineColor(kBlue);
  
//CALIBRATION
 
h1->Fit("gaus","+","same",1020,1100);

double x1274 = h1->GetFunction("gaus")->GetParameter(1);

TF1* fCalibration = new TF1("calibration","pol1",0,1500);
      fCalibration->SetParameter(1,(1274.0-511)/(x1274-x511));
      fCalibration->SetParameter(0,511-fCalibration->GetParameter(1)*x511);

      cout<<"Calibration parameters: p0 = "<<
      fCalibration->GetParameter(0)<<" keV, p1 = "<<
      fCalibration->GetParameter(1)<<" keV/ch"<<endl;

      cout<<"Pedestal at "<<-fCalibration->GetParameter(0)/fCalibration->GetParameter(1)<<endl;


double mean = fCalibration->Eval(x511);
double sigma = fCalibration->GetParameter(1)*fun->GetParameter(5);
double resolution= 100*sigma/mean; 

/*
 std::cout<<"res: "<<resolution<< std::endl;
cout<<"sigma(E)/E_cal = "<<Form("%.2f",100*sigma/mean)<<"%"<<endl;*/


//INTEGRAL 	 
TF1 *fgaus = new TF1("gaus","gaus");
hPeaksOnly->Fit(fgaus,"","same",a,b);
double sig=fgaus->Integral(-10e6,10e6)/h1->GetXaxis()->GetBinWidth(1);
cout<<"integral = "<<sig<<endl; 	 
double x0=fCalibration->GetX(0,0,1100);
double nrbin_start=h1->GetXaxis()->FindBin(x0);
double nrbin_stop=h1->GetNbinsX();
double all=h1->Integral(nrbin_start,nrbin_stop);  

cout<<"nrbinstart = "<<nrbin_start<<endl;
cout<<"nrbinstop = "<<nrbin_stop<<endl;
cout<<"all = "<<all<<endl;
hPeaksOnly->GetFunction("gaus")->Delete(); 


//TLATEX
TLatex wynik;
      wynik.DrawLatexNDC(0.85,0.30,TString::Format("#sigma_{E}/E  = %1.2f %%",resolution)); 
TLatex wynik2;  	
      wynik2.DrawLatexNDC(0.85,0.20,TString::Format("#sigma/all  = %1.2f ",sig/all)); 
/*TLatex wynikcps;
      wynikcps.DrawLatexNDC(0.85,0.10,TString::Format("cps  = %4.2f",cps)); 
      */
      
//SECOND XAXIS -energetic (dodatkowa oś -skalibrowana energia)
TGaxis *axis2 = new TGaxis(h1->GetXaxis()->GetXmin(),
     1.1*h1->GetMaximum(),
     h1->GetXaxis()->GetXmax(),
     1.10*h1->GetMaximum(),
     fCalibration->Eval(h1->GetXaxis()->GetXmin()),
     fCalibration->Eval(h1->GetXaxis()->GetXmax()),
     510,"");
   
  
     axis2->SetName("axis1");
     axis2->SetTitle("E [keV]");
     axis2->SetTitleOffset(-0.20);
     axis2->SetLabelColor(kRed);
     axis2->SetTitleColor(kRed);
     axis2->SetLineColor(kRed);
     axis2->Draw();
     gPad->SetGrid(1,0);
  
  
     h1->GetXaxis()->SetTitle("tot [ns]");
     h1->GetXaxis()->SetLabelFont(42);
     h1->GetXaxis()->SetTitleOffset(1);
     h1->GetXaxis()->SetTitleFont(42);
     h1->GetYaxis()->SetTitle("Counts");
     h1->GetYaxis()->SetLabelFont(42);
     h1->GetYaxis()->SetTitleFont(42);
     h1->GetZaxis()->SetLabelFont(42);
     h1->GetZaxis()->SetTitleOffset(1);
     h1->GetZaxis()->SetTitleFont(42);
     
//H TITTLE
name.ReplaceAll(".root"," ");
h1->SetTitle(name);
  
  
gPad->SetGrid(1,1);
gPad->SetLogy(1);
gPad->SetRightMargin(0.25);
   
      return resolution;
   

   
 
}

