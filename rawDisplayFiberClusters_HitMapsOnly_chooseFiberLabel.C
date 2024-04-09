#include <iostream>

#include "SCategoryManager.h"
#include "SLoop.h"
#include "SFibersRawCluster.h"

//----- prototype geometry 4to1
#define N_LAYERS_PER_MODULE 7
#define N_FIBERS_PER_LAYER 55

//----- cuts
#define QDC_MIN 0

//how to  run: 
//root -l
//.x rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel.C({"/path/to/sifi_tree1.root","/path/to/sifi_tree2.root"})

int rawDisplayFiberClusters_HitMapsOnly_chooseFiberLabel(std::vector<TString> path)
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
    
    //if needed, change the category below to any other category valid in sifi-framework (and add appropriate header)
	SCategory * pCatRawClus = SCategoryManager::getCategory(SCategory::CatFibersRawClus);
    
    TH2D * hFiberHitMap = new TH2D("hFiberHitMap", "hFiberHitMap", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    TH2D * hFiberHitMapNoSemiUnique = new TH2D("hFiberHitMapNoSemiUnique", "hFiberHitMapNoSemiUnique", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
     
	Int_t mod, lay, fi;
    
    SFibersRawCluster * pRawClus;
 	int nLoop = loop->getEntries();
    
    for(int i=0; i<nLoop; i++)
    {
        loop->getEvent(i);
        size_t nCat = pCatRawClus->getEntries();
		for (uint j = 0; j < nCat; j++)
        {
            pRawClus = (SFibersRawCluster *)pCatRawClus->getObject(j);
            pRawClus->getAddress(mod, lay, fi);
            if(pRawClus->getQDCL() > QDC_MIN &&
            pRawClus->getQDCR() > QDC_MIN)
            {
                hFiberHitMap->Fill(fi, lay);
                if(pRawClus->getFiberClusterLabel()>0 && pRawClus->getFiberClusterLabel()<3) hFiberHitMapNoSemiUnique->Fill(fi, lay);
            }
        }
    }

	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path[ii].ReplaceAll(".root", "_HITMAPS_chooseFiberLabel.root");
	TFile *output = new TFile(path[ii],"RECREATE");
	output->cd();
    hFiberHitMapNoSemiUnique->Write();
    hFiberHitMap->Write();
    
    TCanvas * c1 = new TCanvas("c1", "c1", 1600, 800);
    c1->Divide(2,1);
    c1->cd(1);
    hFiberHitMapNoSemiUnique->Draw("colz");
    c1->cd(2);
    hFiberHitMap->Draw("colz");

    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
