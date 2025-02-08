#include <iostream>
#include <memory>
#include <string>

// ROOT Includes
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include <TLatex.h>

#pragma once

/**
* @file
* @author Matthew McEneaney
* @date 16/Dec./24
* @version 0.0.0
* @brief Find ideal bin limits and compute bin migration fractions.
*/

namespace saga {

namespace bins {

/**
* @brief Find bin limits for equal bin statistics
*
* Given the desired number of bins in a distribution, find the bin limits
* that will ensure all bins have roughly equal statistics.
*
* @param frame ROOT RDataFrame with which to find bin limits
* @param varname Bin variable name
* @param nbins Number of bins
*
* @return Bin limits
*/
std::vector<double> findBinLimits(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        std::string varname,
        const int nbins
    ) {

    // Initialize vectors
    std::vector<double> binmeans;
    std::vector<double> binmeans_err;
    std::vector<double> binlims;

    // Get overall limits and count and the necessary bin count with nbins
    double varmin = (double)*frame.Min(varname);
    double varmax = (double)*frame.Max(varname);
    int    count  = (int)*frame.Count();
    double targetbincount = (double)count/nbins;

    // Set initial step size and add initial bin limit
    double step = (double)(varmax-varmin)/nbins;
    binlims.push_back(varmin);

    // Loop each bin and find its upper limit
    int nsteps;
    double threshold = 0.01*targetbincount; //NOTE: This is a threshold in counts, which will be an int.
    for (int i=0; i<nbins; i++) {

        // Set the bin max cut starting at the bin min 
        double binmin = (i==0 ? (double)varmin : binlims.at(binlims.size()-1)); //NOTE: NEED TO SET OUTSIDE OR ADD.....
        double binmax = (i==nbins-1 ? (double)varmax: binmin);
        std::string bin_cut = Form("%s>=%.16f && %s<%.16f",varname.c_str(),binmin,varname.c_str(),binmax);
        int bincount = (int)*frame.Filter(bin_cut).Count();
        bool pass_flag = false;

        // Set the initial adjustment step
        step = (double)(varmax-varmin)/nbins;//IMPORTANT! NEED TO RESET IN CASE IT'S NEGATIVE BUT ALSO SO IT DOESN'T START REALLY SMALL.
        double delta = TMath::Abs(bincount-targetbincount);

        // Adjust the bin max cut until the statistics match within the threshold value
        if (i<nbins-1) {
            while (delta>threshold) {

                // Set the flags for whether or not you surpassed the bin count
                if (bincount>targetbincount)  { pass_flag = true;  }
                if (bincount<=targetbincount) { pass_flag = false; }

                // Reset upper bin limit
                binmax += step;

                // Reset bin cut with new bin lims
                bin_cut = Form("%s>=%.16f && %s<%.16f",varname.c_str(),binmin,varname.c_str(),binmax);

                // Reset bin count with new bin lims
                bincount = (int)*frame.Filter(bin_cut.c_str()).Count();

                // Reset step size if passed targetbincount and switch step direction
                if (pass_flag && bincount<=targetbincount) step *= -0.3; //NOTE: IF YOU DID PASS THE BIN COUNT ALREADY AND ARE NOW LESS 
                if (!pass_flag && bincount>targetbincount) step *= -0.3; //NOTE: IF YOU DIDN'T 

                // Reset delta of bin count to target
                delta = TMath::Abs(bincount-targetbincount);

            } // while(delta<threshold)
        } // if (i<nbins-1)

        // Compute the bin statistics
        double binmean = (double)*frame.Filter(bin_cut).Mean(varname.c_str());
        double binstd  = (double)*frame.Filter(bin_cut).StdDev(varname.c_str());

        // Add to bin lims vector
        binlims.push_back(binmax);

        // Show results message
        std::cout<<"------------------------------------------------------------"<<std::endl;
        std::cout<<" i        = "<<i<<std::endl;
        std::cout<<" varname  = "<<varname.c_str()<<std::endl;
        std::cout<<" bin_cut   = "<<bin_cut.c_str()<<std::endl;
        std::cout<<" varmax   = "<<varmax<<std::endl;
        std::cout<<" varmin   = "<<varmin<<std::endl;
        std::cout<<" mean     = "<<binmean<<std::endl;
        std::cout<<" stddev   = "<<binstd<<std::endl;
        std::cout<<" bincount = "<<bincount<<std::endl;
        std::cout<<" step                = "<<step<<std::endl;
        std::cout<<" bincount            = "<<bincount<<std::endl;
        std::cout<<" |bincount - target| = "<<delta<<std::endl;
        std::cout<<" \% diff             = "<<100*(delta/targetbincount)<<"\%"<<std::endl;
        std::cout<<" bin_cut             = "<<bin_cut<<std::endl;

    } // for (int i=0; i<nbins; i++) {

    // Print out bin limits
    std::cout<<varname<<" = [";
    for (int i=0; i<nbins; i++) {
        double binmin = binlims.at(i);
        std::string limstring = Form(" %.4f,",binmin);
        std::cout<<limstring.c_str();
    }
    std::cout<<" ]"<<std::endl;

    return binlims;

} // std::vector<double> findBinLimits()

/**
* @brief Compute bin limits on a regular interval between a minimum and maximum.
*
* @param nbins Number of bins
* @param xmin Minimum bound of bin variable
* @param xmax Maximum bound of bin variable
*
* @return Bin limits
*/
std::vector<double> getBinLims(
        const int nbins,
        double xmax,
        double xmin
    ) {

    std::vector<double> binlims;
    double step = (xmax-xmin)/nbins;
    for (int i=0; i<nbins+1; i++) {
        double lim = xmin + i*step; 
        binlims.push_back(lim);
    }

    return binlims;

} // std::vector<double> getBinLims(

/**
* @brief Compute 1D bin migration fractions and store in a histogram.
*
* @param frame ROOT RDataframe from which to compute bin migration fraction
* @param varname Bin variable name
* @param vartitle Bin variable title
* @param mcvarname MC bin variable name
* @param mcvartitle MC bin variable title
* @param bins_ Bin limits
* @param drawopt ROOT histogram draw option
* @param f ROOT file in which to save histogram
*
*/
void getBinMigrationHistograms1D(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    std::string varname,
    std::string vartitle,
    std::string mcvarname,
    std::string mcvartitle,
    std::vector<double> bins_,
    const char *drawopt,
    TFile *f
    ) {

    // Convert bin limits vector to array
    const int nbins = bins_.size()-1;
    double bins[nbins+1]; for (int i=0; i<nbins+1; i++) { bins[i] = bins_.at(i); }
    
    // Create histogram name and title bases
    std::string name  = Form("h2d_bin_migration_%s",varname.c_str());
    std::string title = Form("Bin Migration in %s",vartitle.c_str());

    // Create bin migration histograms
    TH1D h1mc_ = (TH1D)*frame.Histo1D({"h1mc_",title.c_str(),nbins,bins},mcvarname.c_str());
    TH1D *h1mc = (TH1D*)h1mc_.Clone(Form("h1mc_%s",name.c_str()));
    TH2D h2_    = (TH2D)*frame.Histo2D({"h2_",title.c_str(),nbins,bins,nbins,bins},mcvarname.c_str(),varname.c_str());
    TH2D *h2   = (TH2D*)h2_.Clone(name.c_str());
    h2->GetXaxis()->SetTitle(mcvartitle.c_str());
    h2->GetYaxis()->SetTitle(vartitle.c_str());

    // Normalize bin migration histogram
    for (int i=1; i<=nbins; i++) { // Loop generated //NOTE: Bin indices begin at 1 in ROOT!
        double divisor = h1mc->GetBinContent(i); //NOTE: Divide by number of generated events!
        for (int j=1; j<=nbins; j++) { // Loop reconstructed
            double bincontent = (divisor==0) ? 0.0 : (double)h2->GetBinContent(i,j)/divisor; // f[i->j] = [# generated in i AND reconstructed in j] / [# generated in bin i]
            h2->SetBinContent(i,j,bincontent); //NOTE: Now this should give the correct matrix to read into numpy
        }
    }

    // Create canvas and draw histogram
    TCanvas *c1 = new TCanvas(Form("c2d_bin_migration_%s",varname.c_str()));
    c1->SetBottomMargin(0.125);
    c1->cd();
    h2->Draw(drawopt);

    // Save canvas
    c1->Write();
    c1->SaveAs(Form("%s.pdf",c1->GetName()));

    // Save to file for future use
    h2->Write();

} // void getBinMigrationHistograms1D()

} // namespace bins {

} // namespace saga {
