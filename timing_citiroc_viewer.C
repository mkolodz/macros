

Bool_t timing_citiroc_viewer(TString path_name, Int_t ch)
{
    int nHeaderLines = 9;
    double duration = 0;
    //----- checking channel number
    if (ch > 63 || ch < 0)
    {
        std::cerr << "##### Error! Incorrect channel number!" << std::endl;
        std::cerr << "Channel numbers from 0 to 63 are correct" << std::endl;
        return kFALSE;
    }

    //----- opening CSV file
    std::ifstream input(path_name);

    if (!input.is_open())
    {
        std::cerr << "##### Error! Cannot open binary file!" << std::endl;
        std::cerr << path_name << std::endl;
        return kFALSE;
    }

    std::string csvLine = " ";
    std::string tmp = " ";

//     input >> tmp >> tmp >> gIPOINTS;
    for(int i=0; i<nHeaderLines; i++) getline(input, csvLine);

//     std::cout << "\n\n-----------------------------------------" << std::endl;
//     std::cout << "Reading oscilloscope file: " << path_name << std::endl;
//     std::cout << "Channel number: " << ch << std::endl;
//     std::cout << "Number of samples: " << gIPOINTS << std::endl;
//     std::cout << "-----------------------------------------\n" << std::endl;

    //----- setting histogram & canvas
    TCanvas* can = new TCanvas(Form("Channel_%i", ch), Form("Channel_%i", ch), 800, 800);
    TH1F* h = new TH1F("h", "h", gIPOINTS, 0, gIPOINTS);
    h->GetXaxis()->SetTitle("time [ns]");
    h->GetYaxis()->SetTitle("amplitude [mV]");

    //----- setting lines
    TLine line_t0;
    TLine line_thr;
    line_t0.SetLineColor(kRed);
    line_thr.SetLineColor(kRed);

    Int_t counter = 0;
    Float_t baseline = 0.;
    Float_t true_thr = 0.;
    Float_t t0 = 0.;
    Float_t t0_ns = 0.;

    Float_t time_start = 0;
    Float_t time_stop = 0;
    Float_t bin_wdth = 0;
    Float_t xmin = 0;
    Float_t xmax = 0;

    std::vector<std::string> csvRow;
    std::string csvElement;

    //----- loop
    while (input.good())
    {
        counter++;
        baseline = 0.;
        true_thr = 0.;

        //----- filling non-base-line-subtracted histogram
        for (Int_t i = 0; i < gIPOINTS; i++)
        {
            csvRow.clear();
            std::getline(input, csvLine);
            std::stringstream stream(csvLine);

            while (std::getline(stream, csvElement, ','))
            {
                csvRow.push_back(csvElement);
            }

            h->SetBinContent(i + 1, std::stof(csvRow[ch + 1]) * 1E3); // mV

            if (i == 0) time_start = std::stof(csvRow[0]);
            if (i == gIPOINTS - 1) time_stop = std::stof(csvRow[0]) - time_start;
        }

        time_start = 0;

        //----- calculating and subtracting base line
        for (Int_t i = 1; i < IBL + 1; i++)
            baseline += h->GetBinContent(i);

        baseline = baseline / IBL; // mV

        if (BL_flag)
        {
            for (Int_t i = 0; i < gIPOINTS; i++)
                h->SetBinContent(i + 1, h->GetBinContent(i + 1) - baseline);

            true_thr = thr; // mV
        }
        else
        {
            true_thr = baseline + thr; // mV
        }

        //----- xaxis range calculation
        bin_wdth = (time_stop - time_start) / (gIPOINTS - 1);

        xmin = (time_start - bin_wdth / 2.) * 1E9; // ns
        xmax = (time_stop + bin_wdth / 2.) * 1E9;  // ns

        h->GetXaxis()->SetLimits(xmin, xmax);
        t0 = GetT0(h, true_thr, BL_flag); // bins
        t0_ns = t0 * bin_wdth * 1E9;      // ns

        //----- drawing
        gPad->SetGrid(1, 1);
        h->SetTitle(Form("Signal number %i", counter));
        h->SetStats(kFALSE);
        h->Draw();

        Float_t ymin = h->GetBinContent(h->GetMinimumBin());
        Float_t ymax = h->GetBinContent(h->GetMaximumBin());
        line_t0.DrawLine(t0_ns, ymin, t0_ns, ymax);
        line_thr.DrawLine(xmin, true_thr, xmax, true_thr);

        can->Update();
        can->WaitPrimitive();
    }

    input.close();
    return kTRUE;
}

