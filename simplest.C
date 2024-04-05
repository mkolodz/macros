#include <SLoop.h>
#include <SCategoryManager.h>
#include <SCategory.h>
// #include "/home/magda/project/sifi-framework/sifi-framework-install/include/SCategoryManager.h"
R__ADD_INCLUDE_PATH(/home/magda/project/sifi-framework/sifi-framework-install/include)
// #include </home/magda/project/sifi-framework/sifi-framework-install/include/SSiPMHit.h>
#include <SSiPMHit.h>
#include <iostream>


int simplest(TString path)
{
    
TFile * input_file;    ///< data input file
TTree *t;
input_file = new TFile(path); //here problem with libraries
t = (TTree*)input_file->Get("S");
t->SetBranchStatus("*", false);
t->SetBranchStatus("SSiPMHit.data.swSiPMID", true);
Int_t variable;
t->SetBranchAddress("SSiPMHit.data.swSiPMID", &variable);
for (int iEntry = 0; t->LoadTree(iEntry) >= 0; ++iEntry) {
//    t->GetEntry(iEntry); //this line causes segmentation viol.
   printf("%d\n", variable); //always the same, irrelevant and very high number is printed here
}

return 0;
}
