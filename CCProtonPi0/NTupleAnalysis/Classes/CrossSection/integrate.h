#ifndef integrate_h
#define integrate_h

double integrate(TH1D* h, double lower_bound, double upper_bound)
{
    
    int lower_bin = h->FindBin(lower_bound);
    int upper_bin = h->FindBin(upper_bound);
    
    double nbkg = h->Integral(lower_bin, upper_bin);
    
        // assume uniform event distribution within the bin,
        // calculate the over counting for the lower bin
    double low_edge = h->GetBinLowEdge(lower_bin);
        //std::cout << "lower bound " << lower_bound << std::endl;
        //std::cout << "low_edge    " << low_edge    << std::endl;
    double dm1 = lower_bound - low_edge;
    assert(dm1 >= 0.0);
    
    double bin_width1   = h->GetBinWidth(lower_bin);
    double bin_content1 = h->GetBinContent(lower_bin);
    double overcount1   = (dm1/bin_width1) * bin_content1;
    
        // assume uniform event distribution within the bin,
        // calculate the over counting for the upper bin
    double high_edge = h->GetBinLowEdge(upper_bin) + h->GetBinWidth(upper_bin);
        //std::cout << "upper bound " << upper_bound << std::endl;
        //std::cout << "high_edge    " << high_edge    << std::endl;
    double dm2 = high_edge - upper_bound;
    assert(dm2 >= 0.0);
    
    double bin_width2   = h->GetBinWidth(upper_bin);
    double bin_content2 = h->GetBinContent(upper_bin);
    double overcount2   = (dm2/bin_width2) * bin_content2;
    
        // correct for over-counting
    double corrected_nbkg = nbkg - overcount1 - overcount2;

        /*
    std::cout << " uncorrected integral, corrected integral: "
              << nbkg << "," << corrected_nbkg
              << std::endl;
        */
    
    return corrected_nbkg;
            
}

#endif
