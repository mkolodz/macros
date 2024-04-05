#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
TCanvas * DrawCanvas(UShort_t from, UShort_t to, const char * label, std::vector<TH1F *> &vec) {
    TCanvas *c = new TCanvas(label, label, 1000, 1000);
    c->Divide(4, 4);
    for(Int_t p=from; p < to; ++p) {
        c->cd(p-from+1);
        vec[p]->SetTitle(Form("%d", p) );
        vec[p]->SetLabelSize(.1, "x");
        vec[p]->SetLineWidth(2);
        vec[p]->SetNdivisions(505, "x");
        vec[p]->Draw();
    }
    return c;
}
void Fill(const char *filename) {
    TFile *f = new TFile(filename, "READ");
    TTree *t = (TTree *)f->Get("events");
    std::vector<TH1F *> vec;
    Float_t tot;
    UInt_t channelID;
    t->SetBranchAddress("tot", &tot);
    t->SetBranchAddress("channelID", &channelID);
    TFile *fOut = new TFile("output.root", "RECREATE");
    for(Int_t j=0; j < 256; ++j) {
        vec.push_back(new TH1F(Form("h%d", j), Form("%d;tot[ns]", j), 100, 0, 500) );
    }
    Int_t N = t->GetEntries();
    for(Int_t i=0; i < N; ++i) {
        t->GetEntry(i);
        printProgress(1.0*i/N );
//        if(channelID > 340) {
//            printf("%d %f\n", channelID, energy);
//            continue;
//        }
        UShort_t ch = channelID - 131200;
        vec[ch]->Fill(tot*1e-3);
    }
//    for(Int_t j=0; j < 256; ++j) {
//        if(vec[j]->GetEntries() > 0)
//        vec[j]->Write();
//    }

    gStyle->SetTitleSize(.2, "t");
    gStyle->SetOptStat(111110);
//    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.4);
    gStyle->SetStatH(0.4);
    TLatex l;
    for(UShort_t i=0; i < 16; ++i) {
        TCanvas *c0 = DrawCanvas(i*16, i*16+16, Form("ch%d to %d", i*16, i*16+16-1), vec);
        c0->SaveAs(Form("ch%d_%d.pdf", i*16, i*16+16-1) );
        c0->Write();
    }
    fOut->Write();
}
