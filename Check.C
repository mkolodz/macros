#include <SLoop.h>
#include <SCategoryManager.h>
#include <SCategory.h>
#include <SSiPMHit.h>
#include <SFibersRaw.h>
void Check() {
    SLoop *loop = new SLoop();
    loop->addFile("/scratch1/gccb/data/TOFPET2/results/00449sifi.root");
    loop->setInput({});
    SCategory * pCatSiPM = SCategoryManager::getCategory(SCategory::CatSiPMHit);
    SSiPMHit* pHit = nullptr;
    int nLoop = loop->getEntries();
    for (int i = 0; i < nLoop; ++i)
    {
        loop->getEvent(i);
        size_t nCat = pCatSiPM->getEntries();

        printf("%zu\n", nCat);
        for (uint j = 0; j < nCat; ++j)
        {
            pHit = dynamic_cast<SSiPMHit*>(pCatSiPM->getObject(j));
            printf("%f\n", pHit->getQDC());
        }
    }
}

