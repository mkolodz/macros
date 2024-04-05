/**
 * plot_dynamic_range.C: ROOT Macro for plotting the collection of ADC vs injected pulse height.
 * 
 * This macro reads a map of vectors filled manually and creates a PDF file for the publication.
 *
 * usage: root -b plot_dynamic_range.C 
 *
 */
void plot_dynamic_range() {
    std::vector<float> pulseheight = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 1000, 1500, 1750, 2000, 2250}; // mV
    std::vector<float> reference(15, 0);
    for(int i=0; i<pulseheight.size(); ++i) {
        reference[i] = 1E+12 * (33E-12 * pulseheight[i]*1E-3 ); // Q = C x V in pC
    }
    std::map<std::string, std::vector<float> > _measured = {
        {"klaus6b10bitLG", {729, 738, 745, 755, 762, 769, 777, 785, 793, 801, 861, 930, 962, 993, 1023} },
        {"twinpeaks", {1049, 1200, 1264, 1313, 1353, 1389, 0, 1449, 0, 0, 0, 0, 0, 0, 0} }, //tot ns
        {"citirocLG", {882, 1288, 1672, 2069, 2470, 2861, 3266, 3668, 0, 0, 0, 0, 0, 0, 0} },
        {"tofpet2b", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} },
        {"drs4", {0, 0, 0, 279, 326, 374, 424, 467, 512, 0, 0, 0, 0, 0, 0} }
    };
    // fit Gauss and use error of the mean wherever possible
    std::map<std::string, std::vector<float> > _err = {
        {"klaus6b11bitLG", {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} },
        {"twinpeaks", {1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0} },
        {"citirocLG", {27, 28, 15, 14, 11, 9, 9, 9, 0, 0, 0, 0, 0, 0, 0} },
        {"tofpet2b", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} },
        {"drs4", {0, 0, 0, 5, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0} }
    };

    std::map<std::string, std::vector<float> > _data;
    for(std::map<std::string, std::vector<float> >::iterator it = _measured.begin(); it != _measured.end(); ++it) {
        for(Int_t i=0; i<it->second.size(); ++i) {
            if(it->first.compare("klaus6b10bitLG") == 0) {
                _data[it->first].push_back((_measured[it->first][i]-707)/reference[i]);
            } else if(it->first.compare("twinpeaks") == 0) {
                _data[it->first].push_back((_measured[it->first][i]-340)/reference[i]);
            } else {
                _data[it->first].push_back(_measured[it->first][i]/reference[i]);
            }
        }
    }

    TGraphErrors *g1 = new TGraphErrors(reference.size(), &reference[0], &_data["klaus6b10bitLG"][0], 0, &_err["klaus6b10bitLG"][0]);
    g1->SetTitle("Klaus6b 10bit LG");
    g1->SetLineWidth(2);
    g1->SetLineColor(kRed);
    TGraphErrors *g2 = new TGraphErrors(reference.size(), &reference[0], &_data["twinpeaks"][0], 0, &_err["twinpeaks"][0]);
    g2->SetTitle("TwinPeaks");
    g2->SetLineWidth(2);
    g2->SetLineColor(kBlue);
    TGraphErrors *g3 = new TGraphErrors(reference.size(), &reference[0], &_data["citirocLG"][0], 0, &_err["citirocLG"][0]);
    g3->SetTitle("Citiroc LG");
    g3->SetLineWidth(2);
    g3->SetLineColor(kGreen+3);
    TGraphErrors *g4 = new TGraphErrors(reference.size(), &reference[0], &_data["tofpet2b"][0], 0, &_err["tofpet2b"][0]);
    g4->SetTitle("TOFPET2b");
    g4->SetLineWidth(2);
    g4->SetLineColor(kMagenta);
    TGraphErrors *g5 = new TGraphErrors(reference.size(), &reference[0], &_data["drs4"][0], 0, &_err["drs4"][0]);
    g5->SetTitle("DRS4");
    g5->SetLineWidth(2);
    g5->SetLineColor(kBlack);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";Charge[pC];#frac{ADC}{Charge[pC]}");
    mg->Add(g1, "ALP");
    mg->Add(g2, "ALP");
    mg->Add(g3, "ALP");
    mg->Add(g4, "ALP");
    mg->Add(g5, "ALP");
    TCanvas *c = new TCanvas("c", "c");
    mg->Draw("A");
    c->Draw();
    TLegend *legend = c->BuildLegend(0.13, 0.66, 0.43, 0.87);
    legend->SetLineWidth(0);
    c->SaveAs("plot_dynamic_range.pdf");
} 
