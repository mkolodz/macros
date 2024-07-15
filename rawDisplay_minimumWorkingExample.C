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
//.x rawDisplayFiberClusters_HitMapsOnly.C("/path/to/sifi_tree.root")

int rawDisplayFiberClusters_HitMapsOnly(TString path)
{
    auto start = std::chrono::system_clock::now();  
    if(!path.Contains("/") || !path.BeginsWith("/"))
    {
        std::cout << "##### Error! The functions needs the filename including the full absolute path..." << std::endl;
        std::abort();
    }
    
    gStyle->SetPalette(kBird);
    
//     std::vector<std::string> str_path;
//     for(int i=0; i < path.size(); i++){
//         str_path.push_back(std::string(path[i]));
//         printf("HERE: %s", str_path[i].c_str());
//     }
    
	SLoop * loop = new SLoop();
	loop->addFile(std::string(path));
	loop->setInput({});

    //if needed, change the category below to any other category valid in sifi-framework (and add appropriate header)
	SCategory * pCatRawClus = SCategoryManager::getCategory(SCategory::CatFibersRawClus);
    
    TH2D * hFiberHitMap = new TH2D("hFiberHitMap", "hFiberHitMap", N_FIBERS_PER_LAYER, 0, N_FIBERS_PER_LAYER, N_LAYERS_PER_MODULE, 0, N_LAYERS_PER_MODULE);
    
	Int_t mod, lay, fi;
    
	SFibersRawCluster * pRawClus;
 	int nLoop = loop->getEntries();

	for (int i = 0; i < nLoop; i++)
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
            }	
		}
	}
	
	std::cout << "\n\nLoop entries: " << nLoop << std::endl;

	path.ReplaceAll(".root", "_HITMAPS.root");
	TFile *output = new TFile(path,"RECREATE");
	output->cd();
    hFiberHitMap->Write();
    
    TCanvas * c1 = new TCanvas("c1", "c1", 800, 500);
    c1->cd();
    hFiberHitMap->Draw("colz");
    
    auto end = std::chrono::system_clock::now();   
    std::chrono::duration<double> time = end-start;
    std::cout << "Time: " << time.count() << " s\n";
    
	return 0;
}
