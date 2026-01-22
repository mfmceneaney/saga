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

// Local Includes
#include <log.h>
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

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::is_same;
using std::map;
using std::ofstream;
using std::ostream;
using std::runtime_error;
using std::string;
using std::unique_ptr;
using std::vector;
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
vector<double> fitResolution(
        RooWorkspace                    *w,
        string                      dataset_name,
        string                      binid,
        string                      bincut,
        vector<string>         binvars,
        vector<string>         resfitvars,
        string                      yamlfile,
        string                      pdf_name         = "gauss",
        string                      fitformula       = "gaus(x[0],x[1],x[2],x[3])",
        vector<string>         parnames         = {"constant","mu","sigma"},
        vector<string>         partitles        = {"C","#mu","#sigma"},
        vector<string>         parunits         = {"","",""},
        vector<double>              parinits         = {1.0, 0.0, 0.1},
        vector<vector<double>> parlims          = {{1.0, 1.0}, {-1.0, 1.0}, {0.0, 1.0}},
        string                      plot_title       = "Fit Resolution",
        double                           lg_text_size     = 0.04,
        double                           lg_margin        = 0.1,
        int                              lg_ncols         = 1,
        bool                             use_sumw2error   = false,
        bool                             use_extended_nll = true,
        bool                             use_binned_fit   = false,
        ostream                    &out              = cout
    ) {

    string method_name = "fitResolution";

    // Load YAML
    YAML::Node node;
    bool loaded_yaml = false;
    if (yamlfile!="") {
        try {
            LOG_DEBUG(Form("[%s]: Loading yaml file %s", method_name.c_str(), yamlfile.c_str()));
            node = YAML::LoadFile(yamlfile.c_str());
            loaded_yaml = true;
        } catch (exception& e) {
            LOG_ERROR(Form("[%s]: Could not load yaml file %s", method_name.c_str(), yamlfile.c_str()));
            throw runtime_error(e.what());
        }
    }

    // Parse YAML
    if (loaded_yaml) {

        // Set parsing parameters
        string message_prefix = "["+method_name+"]: ";
        bool verbose = true;
        ostream &yamlargout = out;

        // Parse arguments
        pdf_name = saga::util::getYamlArg<string>(node, "pdf_name", pdf_name, message_prefix, verbose); //NOTE: This must be non-empty!
        fitformula = saga::util::getYamlArg<string>(node, "fitformula", fitformula, message_prefix, verbose); //NOTE: This is parsed by RooGenericPdf using TFormula
        parnames = saga::util::getYamlArg<vector<string>>(node, "parnames", parnames, message_prefix, verbose);
        partitles = saga::util::getYamlArg<vector<string>>(node, "partitles", partitles, message_prefix, verbose);
        parunits = saga::util::getYamlArg<vector<string>>(node, "parunits", parunits, message_prefix, verbose);
        parinits = saga::util::getYamlArg<vector<double>>(node, "parinits", parinits, message_prefix, verbose);
        parlims = saga::util::getYamlArg<vector<vector<double>>>(node, "parlims", parlims, message_prefix, verbose);
        lg_text_size = saga::util::getYamlArg<double>(node, "lg_text_size", lg_text_size, message_prefix, verbose);
        lg_margin = saga::util::getYamlArg<double>(node, "lg_margin", lg_margin, message_prefix, verbose);
        lg_ncols = saga::util::getYamlArg<double>(node, "lg_ncols", lg_ncols, message_prefix, verbose);
        use_sumw2error = saga::util::getYamlArg<bool>(node, "use_sumw2error", use_sumw2error, message_prefix, verbose);
        use_extended_nll = saga::util::getYamlArg<bool>(node, "use_extended_nll", use_extended_nll, message_prefix, verbose);
        use_binned_fit = saga::util::getYamlArg<bool>(node, "use_binned_fit", use_binned_fit, message_prefix, verbose);
    }

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Load dataset
    LOG_DEBUG(Form("[%s]: Loading dataset from workspace...", method_name.c_str()));
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Cut dataset
    LOG_DEBUG(Form("[%s]: Applying bin cut...", method_name.c_str()));
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(Form("%s", bincut.c_str()));

    // Get count
    auto count = (int)bin_ds->sumEntries();

    // Get bin variable means and errors
    LOG_DEBUG(Form("[%s]: Getting bin variable means and errors...", method_name.c_str()));
    vector<double> binvarmeans;
    vector<double> binvarerrs;
    RooRealVar * b[(const int)binvars.size()];
    for (int i=0; i<binvars.size(); i++) {
        b[i] = w->var(binvars[i].c_str());
        double mean   = bin_ds->mean(*b[i]);
        double stddev = TMath::Sqrt(bin_ds->moment(*b[i],2.0));
        binvarmeans.push_back(mean);
        binvarerrs.push_back(stddev);
    }

    // Load fit variables from workspace
    LOG_DEBUG(Form("[%s]: Loading fit variables from workspace...", method_name.c_str()));
    RooRealVar * f[(const int)resfitvars.size()];
    for (int i=0; i<resfitvars.size(); i++) {
        f[i] = w->var(resfitvars[i].c_str());
    }

    // Create fit parameters
    LOG_DEBUG(Form("[%s]: Creating fit parameters...", method_name.c_str()));
    int nparams = parinits.size();
    RooRealVar *pars[nparams];
    for (int aa=0; aa<nparams; aa++) {
        pars[aa] = new RooRealVar(parnames[aa].c_str(),partitles[aa].c_str(),parinits[aa],parlims[aa][0],parlims[aa][1]);
    }

    // Add parameters to argument list in order
    LOG_DEBUG(Form("[%s]: Adding fit parameters to RooArgSet...", method_name.c_str()));
    RooArgSet *argset = new RooArgSet();
    for (int ff=0; ff<resfitvars.size(); ff++) { // Fit independent variables
        argset->add(*f[ff]);
    }
    for (int aa=0; aa<nparams; aa++) { // Fit parameters
        argset->add(*pars[aa]);
    }

    // Create fit PDF
    LOG_DEBUG(Form("[%s]: Creating pdf %s from formula %s", method_name.c_str(), pdf_name.c_str(), fitformula.c_str()));
    string model_name = Form("%s_%s",pdf_name.c_str(),binid.c_str());
    RooGenericPdf _model(Form("_%s",model_name.c_str()),Form("_%s",model_name.c_str()),fitformula.c_str(),*argset);

    // Create extended fit PDF
    LOG_DEBUG(Form("[%s]: Creating extended pdf...", method_name.c_str()));
    RooRealVar nsig(Form("nsig_%s",model_name.c_str()), "number of events", count, 0.0, 2.0*count);
    RooExtendPdf model(model_name.c_str(), model_name.c_str(), _model, nsig);

    // Fit the PDF to data
    unique_ptr<RooFitResult> r;
    if (use_binned_fit) {

        // Create binned data
        LOG_DEBUG(Form("[%s]: Creating binned dataset...", method_name.c_str()));
        unique_ptr<RooDataHist> dh = (unique_ptr<RooDataHist>)bin_ds->binnedClone();

        // Fit PDF
        LOG_DEBUG(Form("[%s]: Fitting pdf to dataset...", method_name.c_str()));
        if (use_extended_nll) {
            r = (unique_ptr<RooFitResult>)model.fitTo(*dh, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        } else {
            r = (unique_ptr<RooFitResult>)_model.fitTo(*dh, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        }

    } else {

        // Fit PDF
        LOG_DEBUG(Form("[%s]: Fitting pdf to dataset...", method_name.c_str()));
        if (use_extended_nll) {
            r = (unique_ptr<RooFitResult>)model.fitTo(*bin_ds, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        } else {
            r = (unique_ptr<RooFitResult>)_model.fitTo(*bin_ds, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
        }
    }
    
    // Print fit result
    r->Print("v");

    // Extract covariance and correlation matrix as TMatrixDSym
    const TMatrixDSym &corMat = r->correlationMatrix();
    const TMatrixDSym &covMat = r->covarianceMatrix();

    // Print correlation, covariance matrix
    cout << "correlation matrix" << endl;
    corMat.Print();
    cout << "covariance matrix" << endl;
    covMat.Print();

    // Get signal fit parameter values and errors
    LOG_DEBUG(Form("[%s]: Getting fit parameter values and errors...", method_name.c_str()));
    vector<double> fitpars;
    vector<double> fitparerrs;
    for (int aa=0; aa<nparams; aa++) {
        RooRealVar *fitpar = (RooRealVar*)w->var(pars[aa]->GetName()); //NOTE: Load from workspace since parameters are copied to work space when you import the PDF.
        fitpars.push_back((double)fitpar->getVal());
        fitparerrs.push_back((double)fitpar->getError());
    }

    // Compute chi2 from 1D histograms
    LOG_DEBUG(Form("[%s]: Getting chi2 values from 1d histograms...", method_name.c_str()));
    vector<double> chi2ndfs;
    RooDataHist *rdhs_1d[(const int)resfitvars.size()];
    for (int i=0; i<resfitvars.size(); i++) {

        // Import TH1 histogram into RooDataHist //NOTE: For now just use the first fit variable.
        string dh_name = Form("h_%s",f[i]->GetName());
        string dh_title = Form("%s",f[i]->GetTitle());
        rdhs_1d[i] = new RooDataHist(dh_name.c_str(), dh_title.c_str(), *f[i], *bin_ds);

        // Compute chi2/NDF value
        OwningPtr<RooAbsReal> chi2;
        if (use_extended_nll) {
            chi2 = model.createChi2(*rdhs_1d[i], Range("fullRange"),
                    Extended(use_extended_nll), DataError(RooAbsData::Poisson));
        } else {
            chi2 = _model.createChi2(*rdhs_1d[i], Range("fullRange"),
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

        LOG_DEBUG(Form("[%s]: Plotting 1d fit projection in %s", method_name.c_str(), f[i]->GetName()));

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
        string str_chi2 = Form("#chi^{2}/NDF = %.3g",chi2ndfs[i]);

        // Add legend entries
        legend->AddEntry((TObject*)0, str_chi2.c_str(), Form(" %g ",0.0));

        // Create and add legend entries for PDF parameter values and errors
        for (int i=0; i<nparams; i++) {
            string par_str = Form("%s = %.3g #pm %.3g %s", pars[i]->GetTitle(), fitpars[i], fitparerrs[i], parunits[i].c_str());
            legend->AddEntry((TObject*)0, par_str.c_str(), Form(" %g ",chi2ndfs[i]));
        }

        // Create canvas and draw
        string cname = Form("c_%s_%s_%s", method_name.c_str(), binid.c_str(), f[i]->GetName());
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
    out << "------------------------------------------------------------" <<endl;
    out << " "<<method_name.c_str()<<"():" << endl;
    out << " binid:       = " << binid << endl;
    out << " bincut:      = " << bincut << endl;
    out << " bincount     = " << count << endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvarmeans[idx] << "±" << binvarerrs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << " ]" << endl;
    out << " resfitvars      = [" ;
    for (int idx=0; idx<resfitvars.size(); idx++) {
        out << resfitvars[idx];
        if (idx<resfitvars.size()-1) { out << " , "; }
    }
    out << " ]" << endl;
    out << " chi2/ndfs  = [" ;
    for (int idx=0; idx<resfitvars.size(); idx++) {
        out << resfitvars[idx] << " : " << chi2ndfs[idx];
        if (idx<resfitvars.size()-1) { out << " , "; }
    }
    out << "]" << endl;
    out << " fitformula   = " << fitformula.c_str() << endl;
    out << " nparams      = " << nparams <<endl;
    out << " parinits     = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << parinits[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << " ]" << endl;
    out << " pars         = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << fitpars[idx] << "±" << fitparerrs[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << " ]" << endl;
    out << "------------------------------------------------------------" <<endl;

    // Return bin info and fit info
    vector<double> arr;
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

} // vector<double> fitResolution()

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
* @param weight_name Name of weight variable, ignored if empty
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
*
* @throws Runtime error
*/
void getKinBinnedResolutions(
        string                      scheme_name,
        RNode                            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        string                      workspace_name,
        string                      workspace_title,

        // parameters passed to data::createDataset()
        string                      dataset_name,
        string                      dataset_title,
        string                      weight_name,
        vector<string>         categories_as_float,
        string                      helicity,
        map<string,int>        helicity_states,
        string                      tspin,
        map<string,int>        tspin_states,
        string                      htspin,
        map<string,int>        htspin_states,
        string                      combined_spin_state,
        map<int,string>        bincuts,
        vector<string>         binvars,
        vector<string>         binvar_titles,
        vector<vector<double>> binvar_lims,
        vector<int>                 binvar_bins,
        vector<string>         depolvars,
        vector<string>         depolvar_titles,
        vector<vector<double>> depolvar_lims,
        vector<int>                 depolvar_bins,
        vector<string>         resfitvars,
        vector<string>         resfitvar_titles,
        vector<vector<double>> resfitvar_lims,
        vector<int>                 resfitvar_bins,
        vector<string>         massfitvars,
        vector<string>         massfitvar_titles,
        vector<vector<double>> massfitvar_lims,
        vector<int>                 massfitvar_bins,

        // parameters passed to saga::signal::fitResolution()
        map<string,string> yamlfile_map,
        string                       pdf_name,
        string                       fitformula,
        vector<string>          parnames,
        vector<string>          partitles,
        vector<string>          parunits,
        vector<double>               parinits,
        vector<vector<double>>  parlims,

        // Parameters passed to signal::fitResolution()
        double                           lg_text_size     = 0.04,
        double                           lg_margin        = 0.1,
        int                              lg_ncols         = 1,
        bool                             use_sumw2error   = false,
        bool                             use_extended_nll = true,
        bool                             use_binned_fit   = false,

        // Ouput stream
        ostream                    &out              = cout
    ) {

    // Set method name
    string method_name = "getKinBinnedResolutions";

    // Check arguments
    if (binvars.size()<1) {
        string msg = Form("[%s]: Number of bin variables is <1", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }

    // Starting message
    out << "----------------------- " << method_name.c_str() << " ----------------------\n";
    out << "bincuts = { ";
    for (auto it = bincuts.begin(); it != bincuts.end(); it++) {
        out << it->first << " : " << it->second.c_str() << " , ";
    }
    out << " }\n";

    // Open output CSV
    string csvpath = Form("%s.csv",scheme_name.c_str());
    LOG_DEBUG(Form("[%s]: Opening csv file %s", method_name.c_str(), csvpath.c_str()));
    ofstream csvoutf; csvoutf.open(csvpath.c_str());
    ostream &csvout = csvoutf;
    string csv_separator = ",";

    // Set CSV column headers
    // COLS: bin_id,count,{binvarmean,binvarerr},{chi2ndf},{resfitvar,resfitvarerr}
    LOG_DEBUG(Form("[%s]: Writing csv headers...", method_name.c_str()));
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
        else csvout << endl;//NOTE: IMPORTANT!
    }

    // Loop bins and get data
    for (auto it = bincuts.begin(); it != bincuts.end(); it++) {

        // Get bin id and cut
        int         bin_id  = it->first;
        string bin_cut = it->second;

        // Set bin id string
        string scheme_binid = Form("scheme_%s_bin_%d",scheme_name.c_str(),bin_id);

        // Create workspace
        LOG_DEBUG(Form("[%s]: Creating workspace %s", method_name.c_str(), workspace_name.c_str()));
        RooWorkspace *ws    = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());

        // Make bin cut on frame
        LOG_DEBUG(Form("[%s]: Filtering frame with bin cut %s", method_name.c_str(), bin_cut.c_str()));
        auto binframe = frame.Filter(bin_cut.c_str());

        // Create bin dataset
        LOG_DEBUG(Form("[%s]: Creating dataset", method_name.c_str()));
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
        string plot_title = "";
        string bincut_title = bin_cut;
        for (int idx=0; idx<binvars.size(); idx++) {
            saga::util::replaceAll(bincut_title,binvars[idx],binvar_titles[idx]);
        }
        saga::util::replaceAll(bincut_title,"&&","&");
        for (int idx=0; idx<resfitvars.size(); idx++) { plot_title = (plot_title.size()==0) ? Form("%s",resfitvar_titles[idx].c_str()) : Form("%s & %s",plot_title.c_str(),resfitvar_titles[idx].c_str()); }
        plot_title = Form("%s : %s",plot_title.c_str(),bincut_title.c_str());

        // Get resolution fit results
        string yamlfile = yamlfile_map[scheme_binid];
        LOG_DEBUG(Form("[%s]: Fitting resolution for bin scheme %s", method_name.c_str(), scheme_binid.c_str()));
        vector<double> resfitresult = fitResolution(
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
        LOG_DEBUG(Form("[%s]: Initializing arrays...", method_name.c_str()));
        int    nbinvars = binvars.size();
        int    nparams  = parinits.size();
        double xs[nbinvars];
        double exs[nbinvars];
        double chi2ndfs[nbinvars];
        double ys[nparams];
        double eys[nparams];

        // Get resolution fit bin data
        LOG_DEBUG(Form("[%s]: Grabbing data...", method_name.c_str()));
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
        LOG_DEBUG(Form("[%s]: Writing csv data...", method_name.c_str()));
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
            else csvout << endl;//NOTE: IMPORTANT!
        }
    }

    csvoutf.close();
    out << " Saved resolution fit results to " << csvpath.c_str() << endl;

    // Ending message
    out << "------------------- END of "<<method_name.c_str()<<" -------------------\n";

} // getKinBinnedResolutions()

} // namespace saga {

} // namespace resolution {