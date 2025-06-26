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
// #include <TH1.h>
#include <ROOT/RDataFrame.hxx>
// #include <Fit/Fitter.h>
// #include <Fit/BinData.h>
// #include <Fit/Chi2FCN.h>
// #include <Math/WrappedMultiTF1.h>
// #include <HFitInterface.h>
// #include <TGraphErrors.h>
// #include <TRandom.h>
// #include <TF2.h>
// #include <TLatex.h>

// RooFit Includes
#include <RooCategory.h>
#include <RooRealVar.h>
// #include <RooFormulaVar.h>
#include <RooProduct.h>
#include <RooDataSet.h>
#include <RooPlot.h>
// #include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <RooExtendPdf.h>
#include <RooSimultaneous.h>
// #include <RooFFTConvPdf.h>
// #include <RooCrystalBall.h>
// #include <RooLandau.h>
// #include <RooGaussian.h>
// #include <RooChebychev.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>

// RooStats includes
#include <RooStats/SPlot.h>

// Yaml includes
#include <yaml-cpp/yaml.h>

// Local includes
#include <util.h>

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 28/May/2025
* @version 0.0.0
* @brief Fit invariant mass distributions using RooFit unbinned Maximum Likelihood methods.
*/

namespace saga {

namespace signal {

using namespace RooFit;
using RNode = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;

// TODO: Write method to load mass fit parameter initial values and limits and titles and names from a yaml file path and bin id -> Could just leave this for a bigger change to load all parameters from YAML...

/**
* @brief Create a PDF for fitting a combined signal and background invariant mass distribution.
* 
* Create a PDF given the formulas for the signal and background distributions.
* The PDF will be constructed internally using <a href="https://root.cern.ch/doc/master/classRooGenericPdf.html">RooGenericPdf</a>
* in the form:
*
* @f[
* \begin{aligned}
* PDF(x_0, x_1, ..., &a_{SG,0}, a_{SG,1}, ..., a_{BG,0}, a_{BG,1}, ...) = \\
* & u\cdot f_{SG}(\vec{x}, \vec{a}_{SG}) \\
* & + (1-u)\cdot f_{BG}(\vec{x}, \vec{a}_{BG}), \\
* \end{aligned}
* @f]
* where \f$u\f$ is the ratio of signal events to total events in the dataset and
* \f$f_{C}(x_0, x_1, ..., a_{C,0}, a_{C,1}, ...)\f$ for \f$C\in(SG,BG)\f$ are the given signal and background PDFs.
* The `x_<int>` denote the independent fit variables, i.e., the mass variables,
* and the `a_<int>` denote the dependent fit variables, i.e., the signal and background PDF parameters.
*
* The variable names in each fit formula should (separately) follow the <a href="https://root.cern.ch/doc/master/classTFormula.html">TFormula</a> notation, e.g.,
* `x_0`\f$\rightarrow\f$`x[0]`, `x_1`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[N_x]`, `a_1`\f$\rightarrow\f$`x[N_x+1]`, etc.
*
* The indexing of parameters in the fit formula is separate for the signal and background PDFs.
* 
* @param w RooWorkspace in which to work
* @param argset_sg Argument set for the signal PDF
* @param argset_bg Argument set for the background PDF
* @param pdf_name Base name of full PDF
* @param sgYield_name Base name of the signal yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param bgYield_name Base name of the background yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param binid Unique bin id, used to name the PDF
* @param fitformula_sg The signal PDF formula in ROOT TFormula format
* @param fitformula_bg The background PDF formula in ROOT TFormula format
* @param initsgfrac Initial value of the ratio of signal events to total events \f$u\f$ in the dataset
* @param count Bin count
* @param use_extended_nll Option to use an extended likelihood term
* 
* @return List of the names of the combined, signal, and background pdfs and the signal and background yield variabless in that order
*/
std::vector<std::string> getGenMassPdf(
    RooWorkspace *w,
    RooArgSet *argset_sg,
    RooArgSet *argset_bg,
    std::string pdf_name,
    std::string sgYield_name,
    std::string bgYield_name,
    std::string binid,
    std::string fitformula_sg,
    std::string fitformula_bg,
    double initsgfrac,
    int count,
    bool use_extended_nll
) {

    // Create the summed PDF
    RooAddPdf * model;
    std::string model_name = Form("%s_%s",pdf_name.c_str(),binid.c_str());

    //----- Signal PDF -----//

    // Create signal PDF
    RooGenericPdf sg(Form("%s_sg",model_name.c_str()), fitformula_sg.c_str(), *argset_sg);

    // Create extended signal PDF
    RooRealVar sgYield(Form("%s_%s",sgYield_name.c_str(),binid.c_str()), "number of signal events", count*initsgfrac, 0.0, count);

    //----- Background PDF -----//

    // Create signal PDF
    RooGenericPdf bg(Form("%s_bg",model_name.c_str()), fitformula_bg.c_str(), *argset_bg);

    // Create extended signal PDF
    RooRealVar bgYield(Form("%s_%s",bgYield_name.c_str(),binid.c_str()), "number of background events", count*(1.0-initsgfrac), 0.0, count);

    //----- Summed PDF -----//

    // Create the signal fraction parameter
    RooRealVar sgfrac(Form("sgfrac_%s",binid.c_str()), "signal fraction", initsgfrac, 0.0, 1.0);

    // Create the summed PDF
    if (use_extended_nll) {
        model = new RooAddPdf(model_name.c_str(), model_name.c_str(), RooArgList(sg,bg), RooArgList(sgYield,bgYield)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit
    }
    else {
        model = new RooAddPdf(model_name.c_str(), model_name.c_str(), RooArgList(sg,bg), RooArgList(sgfrac)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit
    }

    // Import yield variables and model to workspace
    w->import(sgYield);
    w->import(bgYield);
    w->import(*model);

    // Create and fill return array
    std::vector<std::string> arr;
    arr.push_back(model->GetName());
    arr.push_back(sg.GetName());
    arr.push_back(bg.GetName());
    arr.push_back(sgYield.GetName());
    arr.push_back(bgYield.GetName());

    return arr;

} // std::string getSimGenAsymPdf()


/**
* @brief Fit a combined signal and background mass distribution.
*
* Compute the bin count, bin variable mean values and variances,
* and fit the combined signal and background mass distribution
* with a binned or unbinned dataset using a maximum likelihood fit method with an optional extended likelihood term.
* Note that for the maximum likelihood fit, the given PDF formulas \f$f_{(SG,BG)}(x_0, x_1, ..., a_{SG,0}, a_{SG,1}, ..., a_{BG,0}, a_{BG,1}, ...)\f$
* will be used internally by `getGenMassPdf()` to construct a PDF of the form:
* @f[
* \begin{aligned}
* PDF(x_0, x_1, ..., &a_{SG,0}, a_{SG,1}, ..., a_{BG,0}, a_{BG,1}, ...) = \\
* & u\cdot f_{SG}(\vec{x}, \vec{a}_{SG}) \\
* & + (1-u)\cdot f_{BG}(\vec{x}, \vec{a}_{BG}), \\
* \end{aligned}
* @f]
* where \f$u\f$ is the ratio of signal events to total events in the dataset.
* The `x_<int>` denote the independent fit variables, i.e., the mass variables,
* and the `a_<int>` denote the dependent fit variables, i.e., the signal and background PDF parameters.
*
* The variable names in each fit formula should (separately) follow the <a href="https://root.cern.ch/doc/master/classTFormula.html">TFormula</a> notation, e.g.,
* `x_0`\f$\rightarrow\f$`x[0]`, `x_1`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[N_x]`, `a_1`\f$\rightarrow\f$`x[N_x+1]`, etc.
*
* The mass fit parameters (`massfit_*`) may optionally be loaded from a yaml file containing all these parameters by name
* to allow one to specify bin-dependent starting parameters, limits, and even PDF formulas.
*
* The returned vector will contain, in order:
* - Total bin count
*
* - Signal PDF integral \f$N_{SG}^{PDF}\f$  and error in the signal region
*
* - Background PDF integral \f$N_{BG}^{PDF}\f$ and error in the signal region
*
* - Full PDF integral \f$N^{PDF}\f$ and error in the signal region
*
* - Full dataset sum \f$N^{DS}\f$  and Poissonian error \f$\sqrt{N^{DS}}\f$ in the signal region
*
* - Background fraction \f$\varepsilon_{1} = \frac{N_{BG}^{PDF}}{N^{DS}}\f$ and error
*
* - Background fraction \f$\varepsilon_{2} = 1 - \frac{N_{SG}^{PDF}}{N^{DS}}\f$ and error
*
* - Background fraction \f$\varepsilon_{3} = 1 - \frac{N_{SG}^{PDF}}{N^{PDF}}\f$ and error
*
* - \f$\chi^2/NDF\f$ computed from the 1D histogram in each fit variable
*
* - Bin variable mean values and standard deviations
*
* - Signal PDF parameters and errors
*
* - Background PDF parameters and errors
*
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param fitvars List of mass fit variables
* @param yamlfile Path to YAML file specifying the remaining mass fit arguments
* @param massfit_pdf_name Base name of the combined signal and background pdf.  Note, the actual PDF name will be: `<massfit_pdf_name>_<binid>`.
* @param massfit_formula_sg The signal PDF formula in ROOT TFormula format
* @param massfit_formula_bg The background PDF formula in ROOT TFormula format
* @param massfit_sgYield_name Base name of the signal yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param massfit_bgYield_name Base name of the background yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param massfit_initsgfrac Initial value of the ratio of signal events to total events \f$u\f$ in the dataset
* @param massfit_parinits_sg List of signal PDF parameter initial values
* @param massfit_parnames_sg List of signal PDF parameter names
* @param massfit_partitles_sg List of signal PDF parameter titles
* @param massfit_parunits_sg List of signal PDF parameter unit titles
* @param massfit_parlims_sg List of signal PDF parameter minimum and maximum bounds
* @param massfit_parinits_bg List of background PDF parameter initial values
* @param massfit_parnames_bg List of background PDF parameter names
* @param massfit_partitles_bg List of background PDF parameter titles
* @param massfit_parunits_bg List of background PDF parameter unit titles
* @param massfit_parlims_bg List of background PDF parameter minimum and maximum bounds
* @param massfit_sgregion_lims List of signal region minimum and maximum bounds for each fit variable
* @param massfit_plot_bg_pars Option to plot background pdf parameters on TLegend
* @param massfit_lg_text_size Size of TLegend text
* @param massfit_lg_margin Margin of TLegend
* @param massfit_lg_ncols Number of columns in TLegend
* @param massfit_use_sumw2error Option to use `RooFit::SumW2Error(true)` option when fitting to dataset which is necessary if using a weighted dataset
* @param massfit_use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param massfit_use_binned_fit Option to use a binned fit to the data
* @param out Output stream
*
* @return List of bin count, bin variable means and errors, fit parameters and errors
*/
std::vector<double> fitMass(
        RooWorkspace                    *w,
        std::string                      dataset_name,
        std::string                      binid,
        std::string                      bincut,
        std::vector<std::string>         binvars,
        std::vector<std::string>         fitvars,
        std::string                      yamlfile,
        std::string                      massfit_pdf_name,
        std::string                      massfit_formula_sg,
        std::string                      massfit_formula_bg,
        std::string                      massfit_sgYield_name,
        std::string                      massfit_bgYield_name,
        double                           massfit_initsgfrac,
        std::vector<double>              massfit_parinits_sg,
        std::vector<std::string>         massfit_parnames_sg,
        std::vector<std::string>         massfit_partitles_sg,
        std::vector<std::string>         massfit_parunits_sg,
        std::vector<std::vector<double>> massfit_parlims_sg,
        std::vector<double>              massfit_parinits_bg,
        std::vector<std::string>         massfit_parnames_bg,
        std::vector<std::string>         massfit_partitles_bg,
        std::vector<std::string>         massfit_parunits_bg,
        std::vector<std::vector<double>> massfit_parlims_bg,
        std::vector<std::vector<double>> massfit_sgregion_lims,
        double                           massfit_lg_text_size     = 0.04,
        double                           massfit_lg_margin        = 0.1,
        int                              massfit_lg_ncols         = 1,
        bool                             massfit_plot_bg_pars     = false,
        bool                             massfit_use_sumw2error   = false,
        bool                             massfit_use_extended_nll = true,
        bool                             massfit_use_binned_fit   = false,
        std::ostream                    &out              = std::cout
    ) {

    // Set method name
    std::string method_name = "fitMass";

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
        massfit_pdf_name = saga::util::getYamlArg<std::string>(node, "massfit_pdf_name", massfit_pdf_name, message_prefix, verbose, yamlargout); //NOTE: This must be non-empty!
        massfit_formula_sg = saga::util::getYamlArg<std::string>(node, "massfit_formula_sg", massfit_formula_sg, message_prefix, verbose, yamlargout); //NOTE: This is parsed by RooGenericPdf using TFormula
        massfit_formula_bg = saga::util::getYamlArg<std::string>(node, "massfit_formula_bg", massfit_formula_bg, message_prefix, verbose, yamlargout); //NOTE: This is parsed by RooGenericPdf using TFormula
        massfit_sgYield_name = saga::util::getYamlArg<std::string>(node, "massfit_sgYield_name", massfit_sgYield_name, message_prefix, verbose, yamlargout);
        massfit_bgYield_name = saga::util::getYamlArg<std::string>(node, "massfit_bgYield_name", massfit_bgYield_name, message_prefix, verbose, yamlargout);
        massfit_initsgfrac = saga::util::getYamlArg<double>(node, "massfit_initsgfrac", massfit_initsgfrac, message_prefix, verbose, yamlargout);
        massfit_parnames_sg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parnames_sg", massfit_parnames_sg, message_prefix, verbose, yamlargout);
        massfit_partitles_sg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_partitles_sg", massfit_partitles_sg, message_prefix, verbose, yamlargout);
        massfit_parunits_sg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parunits_sg", massfit_parunits_sg, message_prefix, verbose, yamlargout);
        massfit_parinits_sg = saga::util::getYamlArg<std::vector<double>>(node, "massfit_parinits_sg", massfit_parinits_sg, message_prefix, verbose, yamlargout);
        massfit_parlims_sg = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfit_parlims_sg", massfit_parlims_sg, message_prefix, verbose, yamlargout);
        massfit_parnames_bg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parnames_bg", massfit_parnames_bg, message_prefix, verbose, yamlargout);
        massfit_partitles_bg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_partitles_bg", massfit_partitles_bg, message_prefix, verbose, yamlargout);
        massfit_parunits_bg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parunits_bg", massfit_parunits_bg, message_prefix, verbose, yamlargout);
        massfit_parinits_bg = saga::util::getYamlArg<std::vector<double>>(node, "massfit_parinits_bg", massfit_parinits_bg, message_prefix, verbose, yamlargout);
        massfit_parlims_bg = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfit_parlims_bg", massfit_parlims_bg, message_prefix, verbose, yamlargout);
        massfit_sgregion_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfit_sgregion_lims", massfit_sgregion_lims, message_prefix, verbose, yamlargout);
        massfit_lg_text_size = saga::util::getYamlArg<double>(node, "massfit_lg_text_size", massfit_lg_text_size, message_prefix, verbose, yamlargout);
        massfit_lg_margin = saga::util::getYamlArg<double>(node, "massfit_lg_margin", massfit_lg_margin, message_prefix, verbose, yamlargout);
        massfit_lg_ncols = saga::util::getYamlArg<double>(node, "massfit_lg_ncols", massfit_lg_ncols, message_prefix, verbose, yamlargout);
        massfit_plot_bg_pars = saga::util::getYamlArg<bool>(node, "massfit_plot_bg_pars", massfit_plot_bg_pars, message_prefix, verbose, yamlargout);
        massfit_use_sumw2error = saga::util::getYamlArg<bool>(node, "massfit_use_sumw2error", massfit_use_sumw2error, message_prefix, verbose, yamlargout);
        massfit_use_extended_nll = saga::util::getYamlArg<bool>(node, "massfit_use_extended_nll", massfit_use_extended_nll, message_prefix, verbose, yamlargout);
        massfit_use_binned_fit = saga::util::getYamlArg<bool>(node, "massfit_use_binned_fit", massfit_use_binned_fit, message_prefix, verbose, yamlargout);
    }

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Load fit variables from workspace
    RooRealVar * f[(const int)fitvars.size()];
    for (int i=0; i<fitvars.size(); i++) {
        f[i] = w->var(fitvars[i].c_str());
        f[i]->setRange("fullRange", f[i]->getMin(), f[i]->getMax());//NOTE: DEFINE FULL RANGE FOR COMPUTING CHI2 VARIABLE.
    }

    // Load dataset from workspace
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Apply bin cuts
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(bincut.c_str());

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

    // Create mass fit signal parameters
    int nparams_sg = massfit_parinits_sg.size();
    RooRealVar *a_sg[nparams_sg];
    for (int aa=0; aa<nparams_sg; aa++) {
        a_sg[aa] = new RooRealVar(massfit_parnames_sg[aa].c_str(),massfit_partitles_sg[aa].c_str(),massfit_parinits_sg[aa],massfit_parlims_sg[aa][0],massfit_parlims_sg[aa][1]);
    }

    // Add parameters to signal argument list in order
    RooArgSet *argset_sg = new RooArgSet();
    for (int ff=0; ff<fitvars.size(); ff++) { // Fit independent variables
        argset_sg->add(*f[ff]);
    }
    for (int aa=0; aa<nparams_sg; aa++) { // Fit parameters
        argset_sg->add(*a_sg[aa]);
    }

    // Create mass fit background parameters
    int nparams_bg = massfit_parinits_bg.size();
    RooRealVar *a_bg[nparams_bg];
    for (int aa=0; aa<nparams_bg; aa++) {
        a_bg[aa] = new RooRealVar(massfit_parnames_bg[aa].c_str(),massfit_partitles_bg[aa].c_str(),massfit_parinits_bg[aa],massfit_parlims_bg[aa][0],massfit_parlims_bg[aa][1]);
    }

    // Add parameters to background argument list in order
    RooArgSet *argset_bg = new RooArgSet();
    for (int ff=0; ff<fitvars.size(); ff++) { // Fit independent variables
        argset_bg->add(*f[ff]);
    }
    for (int aa=0; aa<nparams_bg; aa++) { // Fit parameters
        argset_bg->add(*a_bg[aa]);
    }

    // Create and load combined signal and background mass PDF names and yield variables to workspace
    std::vector<std::string> ws_obj_names = getGenMassPdf(
        w,
        argset_sg,
        argset_bg,
        massfit_pdf_name,
        massfit_sgYield_name,
        massfit_bgYield_name,
        binid,
        massfit_formula_sg,
        massfit_formula_bg,
        massfit_initsgfrac,
        count,
        massfit_use_extended_nll
    );

    // Set pdf and variable names for signal+background fit
    std::string model_name    = ws_obj_names[0];
    std::string sg_name       = ws_obj_names[1];
    std::string bg_name       = ws_obj_names[2];
    std::string _sgYield_name = ws_obj_names[3];
    std::string _bgYield_name = ws_obj_names[4];

    // Load pdfs and yield variables from workspace
    RooAbsPdf *model    = w->pdf(model_name.c_str());
    RooAbsPdf *sg       = w->pdf(sg_name.c_str());
    RooAbsPdf *bg       = w->pdf(bg_name.c_str());
    RooRealVar *sgYield = w->var(_sgYield_name.c_str());
    RooRealVar *bgYield = w->var(_bgYield_name.c_str());

    // Fit the PDF to data
    std::unique_ptr<RooFitResult> r;
    if (massfit_use_binned_fit) {

        // Create binned data
        std::unique_ptr<RooDataHist> dh = (std::unique_ptr<RooDataHist>)bin_ds->binnedClone();

        // Fit PDF
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*dh, Save(), SumW2Error(massfit_use_sumw2error), PrintLevel(-1));

    } else {

        // Fit PDF
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*bin_ds, Save(), SumW2Error(massfit_use_sumw2error), PrintLevel(-1));
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
    std::vector<double> massfit_pars_sg;
    std::vector<double> massfit_parerrs_sg;
    for (int aa=0; aa<nparams_sg; aa++) {
        RooRealVar *avar_sg = (RooRealVar*)w->var(a_sg[aa]->GetName()); //NOTE: Load from workspace since parameters are copied to work space when you import the PDF.
        massfit_pars_sg.push_back((double)avar_sg->getVal());
        massfit_parerrs_sg.push_back((double)avar_sg->getError());
    }

    // Get background fit parameter values and errors
    std::vector<double> massfit_pars_bg;
    std::vector<double> massfit_parerrs_bg;
    for (int aa=0; aa<nparams_bg; aa++) {
        RooRealVar *avar_bg = (RooRealVar*)w->var(a_bg[aa]->GetName()); //NOTE: Load from workspace since parameters are copied to work space when you import the PDF.
        massfit_pars_bg.push_back((double)avar_bg->getVal());
        massfit_parerrs_bg.push_back((double)avar_bg->getError());
    }

    // Compute chi2 from 1D histograms
    std::vector<double> chi2ndfs;
    RooDataHist *rdhs_1d[(const int)fitvars.size()];
    for (int i=0; i<fitvars.size(); i++) {

        // Import TH1 histogram into RooDataHist
        std::string dh_title = Form("%s",f[i]->GetTitle());
        rdhs_1d[i] = new RooDataHist(dh_title.c_str(), dh_title.c_str(), *f[i], *bin_ds);

        // Compute chi2 value
        OwningPtr<RooAbsReal> chi2 = model->createChi2(*rdhs_1d[i], Range("fullRange"),
                    Extended(massfit_use_extended_nll), DataError(RooAbsData::Poisson));
        int nparameters = (int)model->getParameters(RooArgSet(*f[i]))->size();
        int ndf = f[i]->getBins() - nparameters; //NOTE: ASSUME ALL BINS NONZERO
        double chi2ndf = (double) chi2->getVal()/ndf;
        chi2ndfs.push_back(chi2ndf);
    }

    //---------------------------------------- Compute integrals and background fractions ----------------------------------------//

    // Define signal regions in fit variables and add fit variables to integration and normalizaation listss
    RooArgSet *iset = new RooArgSet();
    RooArgSet *_nset = new RooArgSet();
    for (int i=0; i<fitvars.size(); i++) {
        f[i]->setRange("signal", massfit_sgregion_lims[i][0], massfit_sgregion_lims[i][1]);
        iset->add(*f[i]);
        _nset->add(*f[i]);
    }
    auto nset = NormSet(*_nset);

    // Integrate the signal PDF
    std::unique_ptr<RooAbsReal> int_sg{sg->createIntegral(*iset, nset, Range("signal"))};
    RooProduct int_sg_pdf{"int_sg_pdf", "int_sg_pdf", {*int_sg, *sgYield}};
    double int_sg_pdf_val = (double)int_sg_pdf.getVal();
    double int_sg_pdf_err = (double)int_sg_pdf.getPropagatedError(*r, *_nset); //NOTE: Need fit result to be saved for this to work!

    // Integrate the background PDF
    std::unique_ptr<RooAbsReal> int_bg{bg->createIntegral(*iset, nset, Range("signal"))};
    RooProduct int_bg_pdf{"int_bg_pdf", "int_bg_pdf", {*int_bg, *bgYield}};
    double int_bg_pdf_val = (double)int_bg_pdf.getVal();
    double int_bg_pdf_err = (double)int_bg_pdf.getPropagatedError(*r, *_nset); //NOTE: Need fit result to be saved for this to work!

    // Integrate the full PDF
    std::unique_ptr<RooAbsReal> int_model{model->createIntegral(*iset, nset, Range("signal"))};
    RooRealVar int_model_pdf("int_model_pdf", "int_model_pdf",sgYield->getVal()+bgYield->getVal());
    double int_model_pdf_val = (double)int_model->getVal() * int_model_pdf.getVal(); //NOTE: Not sure why you need to use the actual integral here instead of the RooRealVar...
    double int_model_pdf_err = (double)int_model->getPropagatedError(*r, *_nset) * int_model_pdf.getVal(); //NOTE: Need fit result to be saved for this to work!

    // Sum on the dataset
    std::string sgcut = saga::util::addLimitCuts("",fitvars,massfit_sgregion_lims);
    double int_ds_val = bin_ds->sumEntries(sgcut.c_str());
    double int_ds_err = TMath::Sqrt(int_ds_val);

    // Compute epsilon in a couple different ways: from the background or signal PDF and normalized by the dataset integral or full PDF integral
    double eps_bg_pdf = int_bg_pdf_val / int_ds_val;
    double eps_sg_pdf = 1.0 - int_sg_pdf_val / int_ds_val;
    double eps_pdf    = 1.0 - int_sg_pdf_val / int_model_pdf_val;

    // Compute epsilon errors in a couple different ways: from the background or signal PDF and normalized by the dataset integral or full PDF integral
    double eps_bg_pdf_err = (double)TMath::Sqrt( TMath::Power(int_bg_pdf_err / int_ds_val,2) + TMath::Power(int_bg_pdf_val / (int_ds_val*int_ds_val) * int_ds_err,2) );
    double eps_sg_pdf_err = (double)TMath::Sqrt( TMath::Power(int_sg_pdf_err / int_ds_val,2) + TMath::Power(int_sg_pdf_val / (int_ds_val*int_ds_val) * int_ds_err,2) );
    double eps_pdf_err    = (double)TMath::Sqrt( TMath::Power(int_sg_pdf_err / int_model_pdf_val,2) + TMath::Power(int_sg_pdf_val / (int_model_pdf_val*int_model_pdf_val) * int_model_pdf_err,2) );

    //---------------------------------------- Plot projections ----------------------------------------//
    for (int i=0; i<fitvars.size(); i++) {

        // Create a histogram from the background PDF
        OwningPtr<RooDataHist> h_bg = bg->generateBinned(RooArgSet(*f[i]), (double)bgYield->getVal());

        // Create the signal histogram by subtracting the background histogram from the full histogram
        RooDataHist *h = (RooDataHist*)rdhs_1d[i]->Clone(Form("%s_clone",rdhs_1d[i]->GetName()));
        std::string mass_cut = Form("%.8f<%s && %s<%.8f",f[i]->getMin(),f[i]->GetName(),f[i]->GetName(),f[i]->getMax());
        h->add(*h_bg, mass_cut.c_str(), -1.);

        // Plot invariant mass fit from RooFit
        RooPlot *mframe_1d = f[i]->frame(Title(Form("1D PDF fit in %s",f[i]->GetTitle())));
        h->plotOn(mframe_1d, LineStyle(kDashed), LineColor(kBlack));
        bin_ds->plotOn(mframe_1d);//NOTE: PDF Normalization correspond to the last data object plotted.
        model->plotOn(mframe_1d);
        model->plotOn(mframe_1d, Components(*sg), LineStyle(kDashed), LineColor(kRed));
        model->plotOn(mframe_1d, Components(*bg), LineStyle(kDashed), LineColor(kBlue));

        // Plot on a TCanvas
        TCanvas *c_massfit = new TCanvas(Form("c_%s_%s_%s",method_name.c_str(),binid.c_str(),f[i]->GetName()));
        c_massfit->cd();
        gPad->SetLeftMargin(0.15);
        mframe_1d->GetYaxis()->SetTitleOffset(1.6);
        mframe_1d->Draw();

        // Create Legend
        TLegend *legend = new TLegend();
        legend->SetTextSize(massfit_lg_text_size);
        legend->SetMargin(massfit_lg_margin);
        if (massfit_lg_ncols>1) legend->SetNColumns(massfit_lg_ncols);

        // Create legend entries
        std::string str_chi2 = Form("#chi^{2}/NDF = %.3g",chi2ndfs[i]);
        std::string str_ntot = Form("N_{Tot} = %.2e #pm %.0f",int_ds_val,int_ds_err);
        std::string str_nbg  = Form("N_{bg} = %.2e #pm %.0f",int_bg_pdf_val,int_bg_pdf_err);
        std::string str_eps  = Form("#varepsilon = %.3f #pm %.3f",eps_bg_pdf,eps_bg_pdf_err); //NOTE: Default to computing background fraction from bg PDF integral / sum(dataset)

        // Add legend entries
        legend->AddEntry((TObject*)0, str_chi2.c_str(), Form(" %g ",int_ds_val));
        legend->AddEntry((TObject*)0, str_ntot.c_str(), Form(" %g ",int_ds_val));
        legend->AddEntry((TObject*)0, str_nbg.c_str(),  Form(" %g ",int_ds_val));
        legend->AddEntry((TObject*)0, str_eps.c_str(),  Form(" %g ",int_ds_val));

        // Create and add legend entries for signal PDF parameter values and errors
        for (int i=0; i<nparams_sg; i++) {
            std::string par_str = Form("%s = %.3g #pm %.3g %s", a_sg[i]->GetTitle(), a_sg[i]->getVal(), a_sg[i]->getError(), massfit_parunits_sg[i].c_str());
            legend->AddEntry((TObject*)0, par_str.c_str(), Form(" %g ",int_ds_val));
        }

        // Optionally create and add legend entries for background PDF parameter values and errors
        if (massfit_plot_bg_pars) {
            for (int i=0; i<nparams_bg; i++) {
                std::string par_str = Form("%s = %.3g #pm %.3g %s", a_bg[i]->GetTitle(), a_bg[i]->getVal(), a_bg[i]->getError(), massfit_parunits_bg[i].c_str());
                legend->AddEntry((TObject*)0, par_str.c_str(), Form(" %g ",int_ds_val));
            }
        }

        // Draw the legend
        legend->Draw();

        // Save the canvas
        c_massfit->SaveAs(Form("%s_%s.pdf",c_massfit->GetName(),binid.c_str()));

    } // for (int i=0; i<fitvars.size(); i++) {

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " "<<method_name.c_str()<<"():" << std::endl;
    out << " bincut     = " << bincut.c_str() << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvarmeans[idx] << "±" << binvarerrs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitvars  = [" ;
    for (int idx=0; idx<fitvars.size(); idx++) {
        out << fitvars[idx];
        if (idx<fitvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " chi2/ndfs  = [" ;
    for (int idx=0; idx<fitvars.size(); idx++) {
        out << fitvars[idx] << " : " << chi2ndfs[idx];
        if (idx<fitvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " massfit_formula_sg  = " << massfit_formula_sg.c_str() << std::endl;
    out << " nparams_sg     = " << nparams_sg <<std::endl;
    out << " massfit_parinits_sg = [" ;
    for (int idx=0; idx<nparams_sg; idx++) {
        out << massfit_parinits_sg[idx];
        if (idx<nparams_sg-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " massfit_pars_sg = [" ;
    for (int idx=0; idx<nparams_sg; idx++) {
        out << massfit_pars_sg[idx] << "±" << massfit_parerrs_sg[idx];
        if (idx<nparams_sg-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " massfit_formula_bg  = " << massfit_formula_bg.c_str() << std::endl;
    out << " nparams_bg     = " << nparams_bg <<std::endl;
    out << " massfit_parinits_bg = [" ;
    for (int idx=0; idx<nparams_bg; idx++) {
        out << massfit_parinits_bg[idx];
        if (idx<nparams_bg-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " massfit_pars_bg = [" ;
    for (int idx=0; idx<nparams_bg; idx++) {
        out << massfit_pars_bg[idx] << "±" << massfit_parerrs_bg[idx];
        if (idx<nparams_bg-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " sgYield       = " << sgYield->getVal() << " ± " << sgYield->getError() << std::endl;
    out << " bgYield       = " << bgYield->getVal() << " ± " << bgYield->getError() << std::endl;
    out << " sgcut         = " << sgcut.c_str() << std::endl;
    out << " int_sg_pdf    = " << int_sg_pdf_val << " ± " << int_sg_pdf_err << std::endl;
    out << " int_bg_pdf    = " << int_bg_pdf_val << " ± " << int_bg_pdf_err << std::endl;
    out << " int_model_pdf = " << int_model_pdf_val << " ± " << int_model_pdf_err << std::endl;
    out << " int_ds        = " << int_ds_val << " ± " << int_ds_err << std::endl;
    out << " eps_bg_pdf    = " << eps_bg_pdf << " ± " << eps_bg_pdf_err << std::endl;
    out << " eps_sg_pdf    = " << eps_sg_pdf << " ± " << eps_sg_pdf_err << std::endl;
    out << " eps_pdf       = " << eps_pdf << " ± " << eps_pdf_err << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Create return array
    std::vector<double> arr; //NOTE: Dimension = 1(bin count)+8(signal region integration)+6(background fractions)+fitvars.size()+2*binvars.size()+2*massfit_parinits_sg.size()+2*massfit_parinits_bg.size()

    // Add total bin count
    arr.push_back(count);

    // Add signal region integration values
    arr.push_back(int_sg_pdf_val);
    arr.push_back(int_sg_pdf_err);
    arr.push_back(int_bg_pdf_val);
    arr.push_back(int_bg_pdf_err);
    arr.push_back(int_model_pdf_val);
    arr.push_back(int_model_pdf_err);
    arr.push_back(int_ds_val);
    arr.push_back(int_ds_err);

    // Add background fractions
    arr.push_back(eps_bg_pdf);
    arr.push_back(eps_bg_pdf_err);
    arr.push_back(eps_sg_pdf);
    arr.push_back(eps_sg_pdf_err);
    arr.push_back(eps_pdf);
    arr.push_back(eps_pdf_err);
    for (int idx=0; idx<chi2ndfs.size(); idx++) {
        arr.push_back(chi2ndfs[idx]);
    }
    for (int idx=0; idx<binvars.size(); idx++) {
        arr.push_back(binvarmeans[idx]);
        arr.push_back(binvarerrs[idx]);
    }
    for (int idx=0; idx<nparams_sg; idx++) {
        arr.push_back(massfit_pars_sg[idx]);
        arr.push_back(massfit_parerrs_sg[idx]);
    }
    for (int idx=0; idx<nparams_bg; idx++) {
        arr.push_back(massfit_pars_bg[idx]);
        arr.push_back(massfit_parerrs_bg[idx]);
    }

    return arr;

} // std::vector<double> fitMass()

/**
* @brief Apply the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
*
* Apply sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a> given a dataset, yield variables, and a PDF model and 
* add the sWeighted datasets to the workspace.
* 
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param model_name Full PDF name
* @param dataset_sg_name Name of dataset with signal sweights
* @param dataset_bg_name Name of dataset with background sweights
*/
void applySPlot(
        RooWorkspace *w,
        std::string dataset_name,
        std::string sgYield_name,
        std::string bgYield_name,
        std::string model_name,
        std::string dataset_sg_name,
        std::string dataset_bg_name
    ) {

    // Get variables from workspace
    RooRealVar *sgYield = w->var(sgYield_name.c_str());
    RooRealVar *bgYield = w->var(bgYield_name.c_str());

    // Get pdf from workspace
    RooAbsPdf *model = w->pdf(model_name.c_str());

    // Get dataset from workspace
    RooDataSet *rooDataSetResult = (RooDataSet*)w->data(dataset_name.c_str());

    // Run sPlot and create weighted datasets
    RooStats::SPlot sData{"sData", "sPlot Data", *rooDataSetResult, model, RooArgList(*sgYield, *bgYield)};
    auto& data = static_cast<RooDataSet&>(*rooDataSetResult);
    RooDataSet data_sg_sw{dataset_sg_name.c_str(), data.GetTitle(), &data, *data.get(), nullptr, Form("%s_sw",sgYield_name.c_str())};
    RooDataSet data_bg_sw{dataset_bg_name.c_str(), data.GetTitle(), &data, *data.get(), nullptr, Form("%s_sw",bgYield_name.c_str())};

    // Import sweighted datasets into workspace
    w->import(data_sg_sw);
    w->import(data_bg_sw);

}

/**
* @brief Set the background fraction of a dataset bin by bin from a map of bin ids to background fraction values.
*
* Set the background fractions of a dataset bin by bin from a map of unique integer bin identifiers to vectors of background fraction values.
* The background fractions are created first from the RDataFrame since defining a conditional variable
* for a RooDataSet is nigh impossible. Then, they are merged into the existing dataset `rds` and a new
* dataset containing the background fraction columns is created and uploaded to the workspace.
*
* @param w RooWorkspace in which to work
* @param rds RooDataSet to which to add the background fraction columns
* @param frame ROOT RDataframe from which to define background fraction columns
* @param rds_out_name Name of the new RooDataSet containing the background fraction columns to import into the RooWorkspace
* @param bincuts Map of unique bin id ints to bin variable cuts for bin
* @param bgfracvar Background fraction variable name
* @param bgfracs_map Map of unique integer bin identifiers to background fraction values
* @param bgfrac_idx Index of the background fraction of interest in the vector of background fraction values
* @param bgfracs_default Weight variable default value for events outside provided cuts
*/
void getBinnedBGFractionsDataset(
        RooWorkspace                      *w,
        RooAbsData                        *rds,
        RNode                             frame,
        std::string                       rds_out_name,
        std::map<int,std::string>         bincuts,
        std::string                       bgfracvar,
        std::map<int,std::vector<double>> bgfracs_map,
        int                               bgfrac_idx      = 0,
        double                            bgfracs_default = 0.0
    ) {

    // Create the conditional bgfrac variable formula
    std::string bgfracvar_formula = "";
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get bin id and cut
        int idx = it->first;
        std::string bincut = it->second;

        // Get background fraction by index
        double bgfrac = bgfracs_map[idx][bgfrac_idx];
        
        // Add to the background fraction variable formula
        if (bgfracvar_formula.size()==0) {
            bgfracvar_formula = Form("if (%s) return %.8f;",bincut.c_str(),bgfrac);
        } else {
            std::string bin_condition = Form("else if (%s) return %.8f;",bincut.c_str(),bgfrac);
            bgfracvar_formula = Form("%s %s",bgfracvar_formula.c_str(),bin_condition.c_str());
        }
    }

    // Set the default background fraction condition
    bgfracvar_formula = Form("%s else return %.8f;",bgfracvar_formula.c_str(),bgfracs_default);

    // Define background fraction in dataframe
    auto frame_bgfracs = frame.Define(bgfracvar.c_str(),bgfracvar_formula.c_str());

    // Create background fraction dataset from dataframe
    RooRealVar * bgf = w->var(bgfracvar.c_str());
    ROOT::RDF::RResultPtr<RooDataSet> rds_bgfracs = frame_bgfracs.Book<double>(
                RooDataSetHelper(rds_out_name.c_str(),rds_out_name.c_str(),RooArgSet(*bgf)),
                {bgfracvar.c_str()}
            );

    // Merge dataset into dataset with bgfracs
    static_cast<RooDataSet&>(*rds_bgfracs).merge(&static_cast<RooDataSet&>(*rds));

    // Import dataset with background fraction columns into workspace
    w->import(*rds_bgfracs);

} // void getBinnedBGFractionsDataset()

/**
* @brief Apply a generic mass fit and set the background fraction for a dataset bin by bin.
*
* Apply a generic mass fit in each asymmetry fit variable bin using `fitMass()`
* and set the background fraction column for the given dataset, which should only contain events from either the signal or sideband region.
* Background fractions \f$\varepsilon\f$ will be taken from one of three choices specified by index:
*
* - `0`: Background fraction \f$\varepsilon_{1} = \frac{N_{BG}^{PDF}}{N^{DS}}\f$ and error
*
* - `1`: Background fraction \f$\varepsilon_{2} = 1 - \frac{N_{SG}^{PDF}}{N^{DS}}\f$ and error
*
* - `2`: Background fraction \f$\varepsilon_{3} = 1 - \frac{N_{SG}^{PDF}}{N^{PDF}}\f$ and error
*
* The naming scheme for the bins will be `<binid>__<asymfitvar_binid>`.
*
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param fitvars List of mass fit variables
* @param yamlfile_map Map of bin ids to paths of YAML files specifying the remaining mass fit arguments
* @param massfit_pdf_name Base name of the combined signal and background pdf.  Note, the actual PDF name will be: `<massfit_pdf_name>_<binid>`.
* @param massfit_formula_sg The signal PDF formula in ROOT TFormula format
* @param massfit_formula_bg The background PDF formula in ROOT TFormula format
* @param massfit_sgYield_name Base name of the signal yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param massfit_bgYield_name Base name of the background yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param massfit_initsgfrac Initial value of the ratio of signal events to total events \f$u\f$ in the dataset
* @param massfit_parinits_sg List of signal PDF parameter initial values
* @param massfit_parnames_sg List of signal PDF parameter names
* @param massfit_partitles_sg List of signal PDF parameter titles
* @param massfit_parunits_sg List of signal PDF parameter unit titles
* @param massfit_parlims_sg List of signal PDF parameter minimum and maximum bounds
* @param massfit_parinits_bg List of background PDF parameter initial values
* @param massfit_parnames_bg List of background PDF parameter names
* @param massfit_partitles_bg List of background PDF parameter titles
* @param massfit_parunits_bg List of background PDF parameter unit titles
* @param massfit_parlims_bg List of background PDF parameter minimum and maximum bounds
* @param massfit_sgregion_lims List of signal region minimum and maximum bounds for each fit variable
* @param frame ROOT RDataframe in which to define the background fraction variable
* @param bgcut Background invariant mass region cut
* @param asymfitvars List of asymmetry fit variables names
* @param asymfitvar_bincuts Map of unique bin id ints to bin variable cuts for bin
* @param rds_out_name Name of signal region RooDataSet under which to import it into the RooWorkspace
* @param sb_rds_out_name Name of sideband region RooDataSet under which to import it into the RooWorkspace
* @param bgfracvar Background fraction variable name
* @param bgfracvar_lims Background fraction variable limits
* @param massfit_plot_bg_pars Option to plot background pdf parameters on TLegend
* @param massfit_lg_text_size Size of TLegend text
* @param massfit_lg_margin Margin of TLegend
* @param massfit_lg_ncols Number of columns in TLegend
* @param massfit_use_sumw2error Option to use `RooFit::SumW2Error(true)` option when fitting to dataset which is necessary if using a weighted dataset
* @param massfit_use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param massfit_use_binned_fit Option to use a binned fit to the data
* @param bgfrac_idx Index of the background fraction of interest from the available choices
* @param bgfracs_default Weight variable default value for events outside provided cuts
* @param out Output stream
*/
void setBinnedBGFractions(
        RooWorkspace                     *w, // fitMass() arguments
        std::string                       dataset_name,
        std::string                       binid,
        std::string                       bincut,
        std::vector<std::string>          binvars,
        std::vector<std::string>          fitvars,
        std::map<std::string,std::string> yamlfile_map,
        std::string                       massfit_pdf_name,
        std::string                       massfit_formula_sg,
        std::string                       massfit_formula_bg,
        std::string                       massfit_sgYield_name,
        std::string                       massfit_bgYield_name,
        double                            massfit_initsgfrac,
        std::vector<double>               massfit_parinits_sg,
        std::vector<std::string>          massfit_parnames_sg,
        std::vector<std::string>          massfit_partitles_sg,
        std::vector<std::string>          massfit_parunits_sg,
        std::vector<std::vector<double>>  massfit_parlims_sg,
        std::vector<double>               massfit_parinits_bg,
        std::vector<std::string>          massfit_parnames_bg,
        std::vector<std::string>          massfit_partitles_bg,
        std::vector<std::string>          massfit_parunits_bg,
        std::vector<std::vector<double>>  massfit_parlims_bg,
        std::vector<std::vector<double>>  massfit_sgregion_lims,
        RNode                             frame, // arguments for this method
        std::string                       bgcut, 
        std::vector<std::string>          asymfitvars,
        std::map<int,std::string>         asymfitvar_bincuts,
        std::string                       rds_out_name,
        std::string                       sb_rds_out_name,
        std::string                       bgfracvar,
        std::vector<double>               bgfracvar_lims           = {0., 1.0},
        double                            massfit_lg_text_size     = 0.04, // fitMass() arguments
        double                            massfit_lg_margin        = 0.1,
        int                               massfit_lg_ncols         = 1,
        bool                              massfit_plot_bg_pars     = false,
        bool                              massfit_use_sumw2error   = false,
        bool                              massfit_use_extended_nll = true,
        bool                              massfit_use_binned_fit   = false,
        int                               bgfrac_idx               = 0,
        double                            bgfracs_default          = 0.0, // arguments for this method
        std::ostream                     &out                      = std::cout
    ) {

    // Load asymmetry fit variables from workspace
    RooRealVar *rrvars[asymfitvars.size()];
    for (int idx=0; idx<asymfitvars.size(); idx++) {
        rrvars[idx] = w->var(asymfitvars[idx].c_str());
    }

    // Create binning scheme if not provided
    if (asymfitvar_bincuts.size()==0) {

        // Create binning scheme
        std::map<std::string,std::vector<double>> binscheme;
        for (int idx=0; idx<asymfitvars.size(); idx++) {
            std::string asymfitvar = asymfitvars[idx];
            std::vector<double> binlims = saga::bins::getBinLims(
                                                            rrvars[idx]->getBins(),
                                                            rrvars[idx]->getMin(),
                                                            rrvars[idx]->getMax()
                                                        );
            binscheme[asymfitvar] = binlims;
        }

        // Set bin cuts from bin scheme
        asymfitvar_bincuts = saga::bins::getBinCuts(binscheme,0);
    }

    // Load dataset from workspace
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Apply bin cuts
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(bincut.c_str());
    RNode bin_frame = frame.Filter(bincut.c_str());

    // Get the signal region cut
    std::string sgcut = saga::util::addLimitCuts("",fitvars,massfit_sgregion_lims);

    // Apply signal and background cuts
    RooDataSet *sg_bin_ds = (RooDataSet*)bin_ds->reduce(sgcut.c_str());
    RNode sg_bin_frame = bin_frame.Filter(sgcut.c_str());
    RooDataSet *bg_bin_ds = (RooDataSet*)bin_ds->reduce(bgcut.c_str());
    RNode bg_bin_frame = bin_frame.Filter(bgcut.c_str());

    // Loop bins and apply fits recording background fractions and errors
    std::map<int,std::vector<double>> bgfracs_map;
    for (auto it = asymfitvar_bincuts.begin(); it != asymfitvar_bincuts.end(); ++it) {

        // Get bin id and cut
        int id = it->first;
        std::string asymfitvar_bincut = it->second;

        // Get bin unique id
        std::string binid_asymfitvars = Form("%s__%d",binid.c_str(),id);

        // Create full bin cut
        std::string bincut_full = Form("%s && %s", bincut.c_str(), asymfitvar_bincut.c_str());

        // Set yaml path for mass fit parameters
        std::string yamlfile = yamlfile_map[binid_asymfitvars]; //NOTE: This should just return an empty string if not found which will use the default parameters

        // Fit the mass spectrum
        std::vector<double> massfit_result = fitMass(
                w, // RooWorkspace                    *w,
                dataset_name, // std::string                      dataset_name,
                binid_asymfitvars, // std::string                      binid,
                bincut_full, // std::string                      bincut,
                binvars, // std::vector<std::string>         binvars,
                fitvars, // std::vector<std::string>         fitvars,
                yamlfile, // std::string                      yamlfile,
                massfit_pdf_name, // std::string                      massfit_pdf_name,
                massfit_formula_sg, // std::string                      massfit_formula_sg,
                massfit_formula_bg, // std::string                      massfit_formula_bg,
                massfit_sgYield_name, // std::string                      massfit_sgYield_name,
                massfit_bgYield_name, // std::string                      massfit_bgYield_name,
                massfit_initsgfrac, // double                           massfit_initsgfrac,
                massfit_parinits_sg, // std::vector<double>              massfit_parinits_sg,
                massfit_parnames_sg, // std::vector<std::string>         massfit_parnames_sg,
                massfit_partitles_sg, // std::vector<std::string>         massfit_partitles_sg,
                massfit_parunits_sg, // std::vector<std::string>         massfit_parunits_sg,
                massfit_parlims_sg, // std::vector<std::vector<double>> massfit_parlims_sg,
                massfit_parinits_bg, // std::vector<double>              massfit_parinits_bg,
                massfit_parnames_bg, // std::vector<std::string>         massfit_parnames_bg,
                massfit_partitles_bg, // std::vector<std::string>         massfit_partitles_bg,
                massfit_parunits_bg, // std::vector<std::string>         massfit_parunits_bg,
                massfit_parlims_bg, // std::vector<std::vector<double>> massfit_parlims_bg,
                massfit_sgregion_lims, // std::vector<std::vector<double>> massfit_sgregion_lims,
                massfit_lg_text_size, // double                           massfit_lg_text_size     = 0.04,
                massfit_lg_margin, // double                           massfit_lg_margin        = 0.1,
                massfit_lg_ncols, // int                              massfit_lg_ncols         = 1,
                massfit_plot_bg_pars, // bool                             massfit_plot_bg_pars     = false,
                massfit_use_sumw2error, // bool                             massfit_use_sumw2error   = true,
                massfit_use_extended_nll, // bool                             massfit_use_extended_nll = false,
                massfit_use_binned_fit, // bool                             massfit_use_binned_fit   = false,
                out // std::ostream                    &out              = std::cout
        );

        // Compute the background fractions
        double eps1 = massfit_result[9];
        double eps2 = massfit_result[11];
        double eps3 = massfit_result[13];
        bgfracs_map[id] = {eps1, eps2, eps3};
    }

    // Create and import background fraction variable so that it is visible to both datasets
    RooRealVar *bgf = new RooRealVar(bgfracvar.c_str(),bgfracvar.c_str(),bgfracvar_lims[0],bgfracvar_lims[1]);
    w->import(*bgf);

    // Set data set background fractions from the binned mass fits for the signal region
    getBinnedBGFractionsDataset(
        w,
        sg_bin_ds,
        sg_bin_frame,
        rds_out_name,
        asymfitvar_bincuts,
        bgfracvar,
        bgfracs_map,
        bgfrac_idx,
        bgfracs_default
    );

    // Set data set background fractions from the binned mass fits for the sideband region
    getBinnedBGFractionsDataset(
        w,
        bg_bin_ds,
        bg_bin_frame,
        sb_rds_out_name,
        asymfitvar_bincuts,
        bgfracvar,
        bgfracs_map,
        bgfrac_idx,
        bgfracs_default
    );

} // void setBinnedBGFractions(

} // namespace signal {

} // namespace saga {
