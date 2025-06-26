#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <map>

// ROOT Includes
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <ROOT/RDataFrame.hxx>

// RooFit Includes
#include <RooCategory.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooGenericPdf.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>

// Local includes
#include <data.h>

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 16/Jun./2025
* @version 0.0.0
* @brief Fit resolution distributions using RooFit unbinned Maximum Likelihood methods.
*/

namespace saga {

namespace resolution {

using namespace RooFit;
using RNode = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;

/**
* @brief Fit a resolution distribution.
*
* Fit a resolution distribution, that is, the difference in reconstructed and true values \f$\Delta X = X_{Rec} - X_{True}\f$
* with a generic PDF, although this will default to Gaussian.  Starting parameter values and limits may be loaded from a yaml file for each bin.
*
* The returned vector will contain, in order:
* - The total bin count
* - The average bin variable means and corresponding standard deviations
* - The \f$\chi^2/NDF\f$ of the fit from a 1D histogram in each fit variable
* - The parameter value and error for each fit parameter
*
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param resfitvars List of resolution fit variables
* @param yamlfile Path to YAML file specifying the remaining fit arguments
* @param pdf_name Base name of PDF.  Note, the actual PDF name will be: `<pdf_name>_<binid>`.
* @param fitformula The PDF formula in ROOT TFormula format
* @param parnames List of PDF parameter names
* @param partitles List of PDF parameter titles
* @param parunits List of PDF parameter unit titles
* @param parinits List of PDF parameter initial values
* @param parlims List of PDF parameter minimum and maximum bounds
* @param plot_title Title of fit plot
* @param lg_text_size Size of TLegend text
* @param lg_margin Margin of TLegend
* @param lg_ncols Number of columns in TLegend
* @param use_sumw2error Option to use `RooFit::SumW2Error(true)` option when fitting to dataset which is necessary if using a weighted dataset
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_binned_fit Option to use a binned fit to the data
* @param out Output stream
*
* @return List containing fit results
*/
std::vector<double> fitResolution(
        RooWorkspace                    *w,
        std::string                      dataset_name,
        std::string                      binid,
        std::string                      bincut,
        std::vector<std::string>         binvars,
        std::vector<std::string>         resfitvars,
        std::string                      yamlfile,
        std::string                      pdf_name         = "gauss",
        std::string                      fitformula       = "gaus(x[0],x[1],x[2],x[3])",
        std::vector<std::string>         parnames         = {"constant","mu","sigma"},
        std::vector<std::string>         partitles        = {"C","#mu","#sigma"},
        std::vector<std::string>         parunits         = {"","",""},
        std::vector<double>              parinits         = {1.0, 0.0, 0.1},
        std::vector<std::vector<double>> parlims          = {{1.0, 1.0}, {-1.0, 1.0}, {0.0, 1.0}},
        std::string                      plot_title       = "Fit Resolution",
        double                           lg_text_size     = 0.04,
        double                           lg_margin        = 0.1,
        int                              lg_ncols         = 1,
        bool                             use_sumw2error   = false,
        bool                             use_extended_nll = true,
        bool                             use_binned_fit   = false,
        std::ostream                    &out              = std::cout
    ) {

    std::string method_name = "fitResolution";

    // Load YAML
    YAML::Node node;
    bool loaded_yaml = false;
    if (yamlfile!="") {
        try {
            node = YAML::LoadFile(yamlfile.c_str());
            loaded_yaml = true;
        } catch (std::exception& e) {
            std::cerr<<"WARNING: "<<method_name.c_str()<<": Could not load yaml: "<<yamlfile.c_str()<<std::endl;
            std::cerr << e.what() << std::endl;
        }
    }

    // Parse YAML
    if (loaded_yaml) {

        // Set parsing parameters
        std::string message_prefix = "INFO: ";
        bool verbose = true;
        std::ostream &yamlargout = out;

        // Parse arguments
        pdf_name = saga::util::getYamlArg<std::string>(node, "pdf_name", pdf_name, message_prefix, verbose, yamlargout); //NOTE: This must be non-empty!
        fitformula = saga::util::getYamlArg<std::string>(node, "fitformula", fitformula, message_prefix, verbose, yamlargout); //NOTE: This is parsed by RooGenericPdf using TFormula
        parnames = saga::util::getYamlArg<std::vector<std::string>>(node, "parnames", parnames, message_prefix, verbose, yamlargout);
        partitles = saga::util::getYamlArg<std::vector<std::string>>(node, "partitles", partitles, message_prefix, verbose, yamlargout);
        parunits = saga::util::getYamlArg<std::vector<std::string>>(node, "parunits", parunits, message_prefix, verbose, yamlargout);
        parinits = saga::util::getYamlArg<std::vector<double>>(node, "parinits", parinits, message_prefix, verbose, yamlargout);
        parlims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "parlims", parlims, message_prefix, verbose, yamlargout);
        lg_text_size = saga::util::getYamlArg<double>(node, "lg_text_size", lg_text_size, message_prefix, verbose, yamlargout);
        lg_margin = saga::util::getYamlArg<double>(node, "lg_margin", lg_margin, message_prefix, verbose, yamlargout);
        lg_ncols = saga::util::getYamlArg<double>(node, "lg_ncols", lg_ncols, message_prefix, verbose, yamlargout);
        use_sumw2error = saga::util::getYamlArg<bool>(node, "use_sumw2error", use_sumw2error, message_prefix, verbose, yamlargout);
        use_extended_nll = saga::util::getYamlArg<bool>(node, "use_extended_nll", use_extended_nll, message_prefix, verbose, yamlargout);
        use_binned_fit = saga::util::getYamlArg<bool>(node, "use_binned_fit", use_binned_fit, message_prefix, verbose, yamlargout);
    }

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Load dataset
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Cut dataset
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(Form("%s", bincut.c_str()));

    // Get count
    auto count = (int)bin_ds->sumEntries();

    // Get bin variable means and errors
    std::vector<double> binvarmeans;
    std::vector<double> binvarerrs;
    RooRealVar * b[(const int)binvars.size()];
    for (int i=0; i<binvars.size(); i++) {
        b[i] = w->var(binvars[i].c_str());
        double mean   = bin_ds->mean(*b[i]);
        double stddev = TMath::Sqrt(bin_ds->moment(*b[i],2.0));
        binvarmeans.push_back(mean);
        binvarerrs.push_back(stddev);
    }

    // Load fit variables from workspace
    RooRealVar * f[(const int)resfitvars.size()];
    for (int i=0; i<resfitvars.size(); i++) {
        f[i] = w->var(resfitvars[i].c_str());
    }

    // Create mass fit signal parameters
    int nparams = parinits.size();
    RooRealVar *pars[nparams];
    for (int aa=0; aa<nparams; aa++) {
        pars[aa] = new RooRealVar(parnames[aa].c_str(),partitles[aa].c_str(),parinits[aa],parlims[aa][0],parlims[aa][1]);
    }

    // Add parameters to argument list in order
    RooArgSet *argset = new RooArgSet();
    for (int ff=0; ff<resfitvars.size(); ff++) { // Fit independent variables
        argset->add(*f[ff]);
    }
    for (int aa=0; aa<nparams; aa++) { // Fit parameters
        argset->add(*pars[aa]);
    }

    // Create fit PDF
    std::string model_name = Form("%s_%s",pdf_name.c_str(),binid.c_str());
    RooGenericPdf _model(Form("_%s",model_name.c_str()),Form("_%s",model_name.c_str()),fitformula.c_str(),*argset);

    // Create extended fit PDF
    RooRealVar nsig(Form("nsig_%s",model_name.c_str()), "number of events", count, 0.0, 2.0*count);
    RooExtendPdf model(model_name.c_str(), model_name.c_str(), _model, nsig);

    // Fit the PDF to data
    std::unique_ptr<RooFitResult> r;
    if (use_binned_fit) {

        // Create binned data
        std::unique_ptr<RooDataHist> dh = (std::unique_ptr<RooDataHist>)bin_ds->binnedClone();

        // Fit PDF
        if (use_extended_nll) {
            r = (std::unique_ptr<RooFitResult>)model.fitTo(*dh, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        } else {
            r = (std::unique_ptr<RooFitResult>)_model.fitTo(*dh, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        }

    } else {

        // Fit PDF
        if (use_extended_nll) {
            r = (std::unique_ptr<RooFitResult>)model.fitTo(*bin_ds, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        } else {
            r = (std::unique_ptr<RooFitResult>)_model.fitTo(*bin_ds, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        }
    }
    
    // Print fit result
    r->Print("v");

    // Extract covariance and correlation matrix as TMatrixDSym
    const TMatrixDSym &corMat = r->correlationMatrix();
    const TMatrixDSym &covMat = r->covarianceMatrix();

    // Print correlation, covariance matrix
    std::cout << "correlation matrix" << std::endl;
    corMat.Print();
    std::cout << "covariance matrix" << std::endl;
    covMat.Print();

    // Get signal fit parameter values and errors
    std::vector<double> fitpars;
    std::vector<double> fitparerrs;
    for (int aa=0; aa<nparams; aa++) {
        fitpars.push_back((double)pars[aa]->getVal());
        fitparerrs.push_back((double)pars[aa]->getError());
    }

    // Compute chi2 from 1D histograms
    std::vector<double> chi2ndfs;
    RooDataHist *rdhs_1d[(const int)resfitvars.size()];
    for (int i=0; i<resfitvars.size(); i++) {

        // Import TH1 histogram into RooDataHist //NOTE: For now just use the first fit variable.
        std::string dh_title = Form("%s",f[i]->GetTitle());
        rdhs_1d[i] = new RooDataHist(dh_title.c_str(), dh_title.c_str(), *f[i], *bin_ds);

        // Compute chi2/NDF value
        OwningPtr<RooAbsReal> chi2;
        if (use_extended_nll) {
            model.createChi2(*rdhs_1d[i], Range("fullRange"),
                    Extended(use_extended_nll), DataError(RooAbsData::Poisson));
        } else {
            _model.createChi2(*rdhs_1d[i], Range("fullRange"),
                    Extended(use_extended_nll), DataError(RooAbsData::Poisson));
        }
        int nparameters;
        if (use_extended_nll) {
            nparameters = (int)model.getParameters(RooArgSet(*f[i]))->size();
        } else {
            nparameters = (int)_model.getParameters(RooArgSet(*f[i]))->size();
        }
        int ndf = f[i]->getBins() - nparameters; //NOTE: ASSUME ALL BINS NONZERO
        double chi2ndf = (double) chi2->getVal()/ndf;
        chi2ndfs.push_back(chi2ndf);
    }

    //---------------------------------------- Plot projections ----------------------------------------//
    for (int i=0; i<resfitvars.size(); i++) {

        // Plot dataset and PDF
        RooPlot *frame = f[i]->frame();
        frame->SetTitle(plot_title.c_str());
        bin_ds->plotOn(frame);
        if (use_extended_nll) {
            model.plotOn(frame);
        } else {
            _model.plotOn(frame);
        }

        // Create legend
        TLegend *legend = new TLegend();
        legend->SetTextSize(lg_text_size);
        legend->SetMargin(lg_margin);
        if (lg_ncols>1) legend->SetNColumns(lg_ncols);

        // Create legend entries
        std::string str_chi2 = Form("#chi^{2}/NDF = %.3g",chi2ndfs[i]);

        // Add legend entries
        legend->AddEntry((TObject*)0, str_chi2.c_str(), Form(" %g ",0.0));

        // Create and add legend entries for PDF parameter values and errors
        for (int i=0; i<nparams; i++) {
            std::string par_str = Form("%s = %.3g #pm %.3g %s", pars[i]->GetTitle(), pars[i]->getVal(), pars[i]->getError(), parunits[i].c_str());
            legend->AddEntry((TObject*)0, par_str.c_str(), Form(" %g ",chi2ndfs[i]));
        }

        // Create canvas and draw
        std::string cname = Form("c_%s_%s_%s", method_name.c_str(), binid.c_str(), f[i]->GetName());
        TCanvas *c = new TCanvas(cname.c_str());
        c->cd();
        gPad->SetLeftMargin(0.15);
        frame->GetYaxis()->SetTitleOffset(1.6);
        frame->Draw();
        legend->Draw();

        // Save to PDF
        c->Print(Form("%s.pdf", cname.c_str()));
    }

    // Show fit info
    out << "------------------------------------------------------------" <<std::endl;
    out << " "<<method_name.c_str()<<"():" << std::endl;
    out << " binid:       = " << binid << std::endl;
    out << " bincut:      = " << bincut << std::endl;
    out << " bincount     = " << count << std::endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvarmeans[idx] << "±" << binvarerrs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << " ]" << std::endl;
    out << " resfitvars      = [" ;
    for (int idx=0; idx<resfitvars.size(); idx++) {
        out << resfitvars[idx];
        if (idx<resfitvars.size()-1) { out << " , "; }
    }
    out << " ]" << std::endl;
    out << " chi2/ndfs  = [" ;
    for (int idx=0; idx<resfitvars.size(); idx++) {
        out << resfitvars[idx] << " : " << chi2ndfs[idx];
        if (idx<resfitvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitformula   = " << fitformula.c_str() << std::endl;
    out << " nparams      = " << nparams <<std::endl;
    out << " parinits     = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << parinits[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << " ]" << std::endl;
    out << " pars         = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << fitpars[idx] << "±" << fitparerrs[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << " ]" << std::endl;
    out << "------------------------------------------------------------" <<std::endl;

    // Return bin info and fit info
    std::vector<double> arr;
    arr.push_back(count);
    for (int idx=0; idx<binvars.size(); idx++) {
        arr.push_back(binvarmeans[idx]);
        arr.push_back(binvarerrs[idx]);
    }
    for (int idx=0; idx<chi2ndfs.size(); idx++) {
        arr.push_back(chi2ndfs[idx]);
    }
    for (int i = 0; i<nparams; i++) {
        arr.push_back(fitpars[i]);
        arr.push_back(fitparerrs[i]);
    }
    return arr;

} // std::vector<double> fitResolution()

/**
* @brief Loop kinematic bins and fit a resolution distribution.
*
* Loop bins cuts and fit a resolution distribution with the `saga::resolution::fitResolution()` method.
*
* Results will be saved in a csv file with the following columns:
*
* - `bin_id`: The unique bin id
*
* - `count`: The total number of counts in the bin
*
* - For each bin variable `binvar`
*
*   - `<binvar>`: Mean value
*
*   - `<binvar>_err`: Standard deviation
*
* - For each independent fit variable:
*
*   - `<chi2ndf>`: \f$\chi^2/NDF\f$ value of the 1D projection of the full PDF in that variable
*
* - For each resolution fit PDF parameter `fitpar`
*
*   - `<fitpar>`: Final parameter value
*
*   - `<fitpar>_err`: Final parameter error
*
* @param scheme_name Name bin scheme and basename of output csv file
* @param frame ROOT RDataframe from which to create RooFit datasets
* @param workspace_name Name of workspace in which to work
* @param workspace_title Title of workspace in which to work
* @param dataset_name Dataset name
* @param dataset_title Dataset title
* @param categories_as_float List of category variables to include as asymmetry fit variables in dataset
* @param helicity Name of helicity variable
* @param helicity_states Map of state names to helicity values
* @param tspin Name of target spin variable
* @param tspin_states Map of state names to target spin values
* @param htspin Name of helicity times target spin variable
* @param htspin_states Map of state names to helicity times target spin values
* @param combined_spin_state Name of combined spin state variable
* @param bincuts Map of unique bin id ints to bin variable cuts for bin
* @param binvars List of kinematic binning variables names
* @param binvar_titles List of kinematic binning variables titles
* @param binvar_lims List kinematic binning variable minimum and maximum bounds 
* @param binvar_bins List of kinematic binning variables bins
* @param depolvars List of depolarization variables names
* @param depolvar_titles List of depolarization variables titles
* @param depolvar_lims List depolarization variable minimum and maximum bounds 
* @param depolvar_bins List of depolarization variables bins
* @param resfitvars List of resolution fit variables names
* @param resfitvar_titles List of resolution fit variables titles
* @param resfitvar_lims List resolution fit variable minimum and maximum bounds 
* @param resfitvar_bins List of resolution fit variables bins
* @param massfitvars List of invariant mass fit variables names
* @param massfitvar_titles List of invariant mass fit variables titles
* @param massfitvar_lims List invariant mass fit variable minimum and maximum bounds 
* @param massfitvar_bins List of invariant mass fit variables bins
*
* @param yamlfile_map Map of bin ids to the paths of yaml files specifying the remaining resolution fit arguments.  Note that the values specified here will function as the defaults.
* @param pdf_name Base name of the resolution PDF.  Note, the actual PDF name will be: `<pdf_name>_<binid>`.
* @param fitformula The resolution PDF formula in ROOT TFormula format
* @param parnames List of resolution PDF parameter names
* @param partitles List of resolution PDF parameter titles
* @param parunits List of resolution PDF parameter unit titles
* @param parinits List of resolution PDF parameter initial values
* @param parlims List of resolution PDF parameter minimum and maximum bounds
*
* @param lg_text_size Size of TLegend text for the signal and background mass fit
* @param lg_margin Margin of TLegend for the signal and background mass fit
* @param lg_ncols Number of columns in TLegend for the signal and background mass fit
* @param use_sumw2error Option to use `RooFit::SumW2Error(true)` option for the signal and background mass fit which is necessary if using a weighted dataset 
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization for the signal and background mass fit
* @param use_binned_fit Option to use a binned fit to the data for the signal and background mass fit
*
* @param out Output stream
*/
void getKinBinnedResolutions(
        std::string                      scheme_name,
        RNode                            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        std::string                      workspace_name,
        std::string                      workspace_title,

        // parameters passed to data::createDataset()
        std::string                      dataset_name,
        std::string                      dataset_title,
        std::vector<std::string>         categories_as_float,
        std::string                      helicity,
        std::map<std::string,int>        helicity_states,
        std::string                      tspin,
        std::map<std::string,int>        tspin_states,
        std::string                      htspin,
        std::map<std::string,int>        htspin_states,
        std::string                      combined_spin_state,
        std::map<int,std::string>        bincuts,
        std::vector<std::string>         binvars,
        std::vector<std::string>         binvar_titles,
        std::vector<std::vector<double>> binvar_lims,
        std::vector<int>                 binvar_bins,
        std::vector<std::string>         depolvars,
        std::vector<std::string>         depolvar_titles,
        std::vector<std::vector<double>> depolvar_lims,
        std::vector<int>                 depolvar_bins,
        std::vector<std::string>         resfitvars,
        std::vector<std::string>         resfitvar_titles,
        std::vector<std::vector<double>> resfitvar_lims,
        std::vector<int>                 resfitvar_bins,
        std::vector<std::string>         massfitvars,
        std::vector<std::string>         massfitvar_titles,
        std::vector<std::vector<double>> massfitvar_lims,
        std::vector<int>                 massfitvar_bins,

        // parameters passed to saga::signal::fitResolution() //TODO: Add init fit parameter value and limits arguments here...assuming you always want a chebychev polynomial background...
        std::map<std::string,std::string> yamlfile_map,
        std::string                       pdf_name,
        std::string                       fitformula,
        std::vector<std::string>          parnames,
        std::vector<std::string>          partitles,
        std::vector<std::string>          parunits,
        std::vector<double>               parinits,
        std::vector<std::vector<double>>  parlims,

        // Parameters passed to signal::fitResolution()
        double                           lg_text_size     = 0.04,
        double                           lg_margin        = 0.1,
        int                              lg_ncols         = 1,
        bool                             use_sumw2error   = false,
        bool                             use_extended_nll = true,
        bool                             use_binned_fit   = false,

        // Ouput stream
        std::ostream                    &out              = std::cout
    ) {

    // Check arguments
    if (binvars.size()<1) {std::cerr<<"ERROR: Number of bin variables is <1.  Exiting...\n"; return;}

    // Starting message
    std::string method_name = "getKinBinnedResolutions";
    out << "----------------------- "<<method_name.c_str()<<" ----------------------\n";
    out << "bincuts = { ";
    for (auto it = bincuts.begin(); it != bincuts.end(); it++) {
        out << it->first << " : " << it->second.c_str() << " , ";
    }
    out << " }\n";

    // Open output CSV
    std::string csvpath = Form("%s.csv",scheme_name.c_str());
    std::ofstream csvoutf; csvoutf.open(csvpath.c_str());
    std::ostream &csvout = csvoutf;
    std::string csv_separator = ",";

    // Set CSV column headers
    // COLS: bin_id,count,{binvarmean,binvarerr},{chi2ndf},{resfitvar,resfitvarerr}
    csvout << "bin_id" << csv_separator.c_str();
    csvout << "count" << csv_separator.c_str();
    for (int bb=0; bb<binvars.size(); bb++) {
        csvout << binvars[bb].c_str() << csv_separator.c_str();
        csvout << binvars[bb].c_str() << "_err" << csv_separator.c_str();
    }
    for (int bb=0; bb<binvars.size(); bb++) {
        csvout << "chi2ndf_" << binvars[bb].c_str() << csv_separator.c_str();
    }
    for (int aa=0; aa<parinits.size(); aa++) {
        csvout << parnames[aa].c_str() << csv_separator.c_str();
        csvout << parnames[aa].c_str() << "_err";
        if (aa<parinits.size()-1) csvout << csv_separator.c_str();
        else csvout << std::endl;//NOTE: IMPORTANT!
    }

    // Loop bins and get data
    for (auto it = bincuts.begin(); it != bincuts.end(); it++) {

        // Get bin id and cut
        int         bin_id  = it->first;
        std::string bin_cut = it->second;

        // Set bin id string
        std::string scheme_binid = Form("scheme_%s_bin_%d",scheme_name.c_str(),bin_id);

        // Create workspace
        RooWorkspace *ws    = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());

        // Make bin cut on frame
        auto binframe = frame.Filter(bin_cut.c_str());

        // Create bin dataset
        data::createDataset(
            binframe,
            ws,
            dataset_name,
            dataset_title,
            categories_as_float,
            helicity,
            helicity_states,
            tspin,
            tspin_states,
            htspin,
            htspin_states,
            combined_spin_state,
            binvars,
            binvar_titles,
            binvar_lims,
            binvar_bins,
            depolvars,
            depolvar_titles,
            depolvar_lims,
            depolvar_bins,
            resfitvars,
            resfitvar_titles,
            resfitvar_lims,
            resfitvar_bins,
            massfitvars,
            massfitvar_titles,
            massfitvar_lims,
            massfitvar_bins
        );

        // Set plot title
        std::string plot_title = "";
        std::string bincut_title = bin_cut;
        for (int idx=0; idx<binvars.size(); idx++) {
            saga::util::replaceAll(bincut_title,binvars[idx],binvar_titles[idx]);
        }
        saga::util::replaceAll(bincut_title,"&&","&");
        for (int idx=0; idx<resfitvars.size(); idx++) { plot_title = (plot_title.size()==0) ? Form("%s",resfitvar_titles[idx].c_str()) : Form("%s & %s",plot_title.c_str(),resfitvar_titles[idx].c_str()); }
        plot_title = Form("%s : %s",plot_title.c_str(),bincut_title.c_str());

        // Get resolution fit results
        std::string yamlfile = yamlfile_map[scheme_binid];
        std::vector<double> resfitresult = fitResolution(
                                ws,
                                dataset_name,
                                scheme_binid,
                                bin_cut,
                                binvars,
                                resfitvars,
                                yamlfile,
                                pdf_name,
                                fitformula,
                                parnames,
                                partitles,
                                parunits,
                                parinits,
                                parlims,
                                plot_title,
                                lg_text_size,
                                lg_margin,
                                lg_ncols,
                                use_sumw2error,
                                use_extended_nll,
                                use_binned_fit,
                                out
                            );

        // Initialize data
        int    nbinvars = binvars.size();
        int    nparams  = parinits.size();
        double xs[nbinvars];
        double exs[nbinvars];
        double chi2ndfs[nbinvars];
        double ys[nparams];
        double eys[nparams];

        // Get resolution fit bin data
        int k = 0;
        int count = (int)resfitresult[k++];
        for (int idx=0; idx<binvars.size(); idx++) {
            xs[idx]     = resfitresult[k++];
            exs[idx]    = resfitresult[k++];
        }
        for (int idx=0; idx<binvars.size(); idx++) {
            chi2ndfs[idx]     = resfitresult[k++];
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx] = resfitresult[k++];
            eys[idx] = resfitresult[k++];
        }

        // Write out a row of data to csv
        // COLS: bin_id,count,{binvarmean,binvarerr},{chi2ndf},{resfitvar,resfitvarerr}
        csvout << bin_id << csv_separator.c_str();
        csvout << count << csv_separator.c_str();
        for (int bb=0; bb<binvars.size(); bb++) {
            csvout << xs[bb] << csv_separator.c_str();
            csvout << exs[bb] << csv_separator.c_str();
        }
        for (int bb=0; bb<binvars.size(); bb++) {
            csvout << chi2ndfs[bb] << csv_separator.c_str();
        }
        for (int aa=0; aa<nparams; aa++) {
            csvout << ys[aa] << csv_separator.c_str();
            csvout << eys[aa];
            if (aa<nparams-1) csvout << csv_separator.c_str();
            else csvout << std::endl;//NOTE: IMPORTANT!
        }
    }

    csvoutf.close();
    out << " Saved resolution fit results to " << csvpath.c_str() << std::endl;

    // Ending message
    out << "------------------- END of "<<method_name.c_str()<<" -------------------\n";

} // getKinBinnedResolutions()

} // namespace saga {

} // namespace resolution {