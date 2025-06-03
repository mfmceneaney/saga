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
* @param pdf_name Base name of the combined signal and background pdf.  Note, the actual PDF name will be: `<pdf_name>_<binid>`.
* @param fitformula_sg The signal PDF formula in ROOT TFormula format
* @param fitformula_bg The background PDF formula in ROOT TFormula format
* @param sgYield_name Base name of the signal yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param bgYield_name Base name of the background yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param initsgfrac Initial value of the ratio of signal events to total events \f$u\f$ in the dataset
* @param fitparinits_sg List of signal PDF parameter initial values
* @param fitparnames_sg List of signal PDF parameter names
* @param fitpartitles_sg List of signal PDF parameter titles
* @param fitparunits_sg List of signal PDF parameter unit titles
* @param fitparlims_sg List of signal PDF parameter minimum and maximum bounds
* @param fitparinits_bg List of background PDF parameter initial values
* @param fitparnames_bg List of background PDF parameter names
* @param fitpartitles_bg List of background PDF parameter titles
* @param fitparunits_bg List of background PDF parameter unit titles
* @param fitparlims_bg List of background PDF parameter minimum and maximum bounds
* @param sgregion_lims List of signal region minimum and maximum bounds for each fit variable
* @param plot_bg_pars Option to plot background pdf parameters on TLegend
* @param lg_text_size Size of TLegend text
* @param lg_margin Margin of TLegend
* @param lg_ncols Number of columns in TLegend
* @param use_sumw2error Option to use `RooFit::SumW2Error(true)` option when fitting to dataset which is necessary if using a weighted dataset
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_binned_fit Option to use a binned fit to the data
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
        std::string                      pdf_name,
        std::string                      fitformula_sg,
        std::string                      fitformula_bg,
        std::string                      sgYield_name,
        std::string                      bgYield_name,
        double                           initsgfrac,
        std::vector<double>              fitparinits_sg,
        std::vector<std::string>         fitparnames_sg,
        std::vector<std::string>         fitpartitles_sg,
        std::vector<std::string>         fitparunits_sg,
        std::vector<std::vector<double>> fitparlims_sg,
        std::vector<double>              fitparinits_bg,
        std::vector<std::string>         fitparnames_bg,
        std::vector<std::string>         fitpartitles_bg,
        std::vector<std::string>         fitparunits_bg,
        std::vector<std::vector<double>> fitparlims_bg,
        std::vector<std::vector<double>> sgregion_lims,
        double                           lg_text_size     = 0.04,
        double                           lg_margin        = 0.1,
        int                              lg_ncols         = 1,
        bool                             plot_bg_pars     = false,
        bool                             use_sumw2error   = false,
        bool                             use_extended_nll = true,
        bool                             use_binned_fit   = false,
        std::ostream                    &out              = std::cout
    ) {

    // Set method name
    std::string method_name = "fitMass";

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
    int nparams_sg = fitparinits_sg.size();
    RooRealVar *a_sg[nparams_sg];
    for (int aa=0; aa<nparams_sg; aa++) {
        a_sg[aa] = new RooRealVar(fitparnames_sg[aa].c_str(),fitpartitles_sg[aa].c_str(),fitparinits_sg[aa],fitparlims_sg[aa][0],fitparlims_sg[aa][1]);
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
    int nparams_bg = fitparinits_bg.size();
    RooRealVar *a_bg[nparams_bg];
    for (int aa=0; aa<nparams_bg; aa++) {
        a_bg[aa] = new RooRealVar(fitparnames_bg[aa].c_str(),fitpartitles_bg[aa].c_str(),fitparinits_bg[aa],fitparlims_bg[aa][0],fitparlims_bg[aa][1]);
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
        pdf_name,
        sgYield_name,
        bgYield_name,
        binid,
        fitformula_sg,
        fitformula_bg,
        initsgfrac,
        count,
        use_extended_nll
    );

    // Set pdf and variable names for signal+background fit
    std::string model_name   = ws_obj_names[0];
    std::string sg_name      = ws_obj_names[1];
    std::string bg_name      = ws_obj_names[2];
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
    if (use_binned_fit) {

        // Create binned data
        std::unique_ptr<RooDataHist> dh = (std::unique_ptr<RooDataHist>)bin_ds->binnedClone();

        // Fit PDF
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*dh, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));

    } else {

        // Fit PDF
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*bin_ds, Save(), SumW2Error(use_sumw2error), PrintLevel(-1));
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
    std::vector<double> fitpars_sg;
    std::vector<double> fitparerrs_sg;
    for (int aa=0; aa<nparams_sg; aa++) {
        fitpars_sg.push_back((double)a_sg[aa]->getVal());
        fitparerrs_sg.push_back((double)a_sg[aa]->getError());
    }

    // Get background fit parameter values and errors
    std::vector<double> fitpars_bg;
    std::vector<double> fitparerrs_bg;
    for (int aa=0; aa<nparams_bg; aa++) {
        fitpars_bg.push_back((double)a_bg[aa]->getVal());
        fitparerrs_bg.push_back((double)a_bg[aa]->getError());
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
                    Extended(use_extended_nll), DataError(RooAbsData::Poisson));
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
        f[i]->setRange("signal", sgregion_lims[i][0], sgregion_lims[i][1]);
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
    std::string sgcut = saga::util::addLimitCuts("",fitvars,sgregion_lims);
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
        TCanvas *c_massfit = new TCanvas(Form("c_%s_%s",method_name.c_str(),binid.c_str()));
        c_massfit->cd();
        gPad->SetLeftMargin(0.15);
        mframe_1d->GetYaxis()->SetTitleOffset(1.6);
        mframe_1d->Draw();

        // Create Legend
        TLegend *legend = new TLegend();
        legend->SetTextSize(lg_text_size);
        legend->SetMargin(lg_margin);
        if (lg_ncols>1) legend->SetNColumns(lg_ncols);

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
            std::string par_str = Form("%s = %.3g #pm %.3g %s", a_sg[i]->GetTitle(), a_sg[i]->getVal(), a_sg[i]->getError(), fitparunits_sg[i].c_str());
            legend->AddEntry((TObject*)0, par_str.c_str(), Form(" %g ",int_ds_val));
        }

        // Optionally create and add legend entries for background PDF parameter values and errors
        if (plot_bg_pars) {
            for (int i=0; i<nparams_bg; i++) {
                std::string par_str = Form("%s = %.3g #pm %.3g %s", a_bg[i]->GetTitle(), a_bg[i]->getVal(), a_bg[i]->getError(), fitparunits_bg[i].c_str());
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
    out << " fitformula_sg  = " << fitformula_sg.c_str() << std::endl;
    out << " nparams_sg     = " << nparams_sg <<std::endl;
    out << " fitparinits_sg = [" ;
    for (int idx=0; idx<nparams_sg; idx++) {
        out << fitparinits_sg[idx];
        if (idx<nparams_sg-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitpars_sg = [" ;
    for (int idx=0; idx<nparams_sg; idx++) {
        out << fitpars_sg[idx] << "±" << fitparerrs_sg[idx];
        if (idx<nparams_sg-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitformula_bg  = " << fitformula_bg.c_str() << std::endl;
    out << " nparams_bg     = " << nparams_bg <<std::endl;
    out << " fitparinits_bg = [" ;
    for (int idx=0; idx<nparams_bg; idx++) {
        out << fitparinits_bg[idx];
        if (idx<nparams_bg-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitpars_bg = [" ;
    for (int idx=0; idx<nparams_bg; idx++) {
        out << fitpars_bg[idx] << "±" << fitparerrs_bg[idx];
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
    std::vector<double> arr; //NOTE: Dimension = 1(bin count)+8(signal region integration)+6(background fractions)+fitvars.size()+2*binvars.size()+2*fitparinits_sg.size()+2*fitparinits_bg.size()

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
        arr.push_back(fitpars_sg[idx]);
        arr.push_back(fitparerrs_sg[idx]);
    }
    for (int idx=0; idx<nparams_bg; idx++) {
        arr.push_back(fitpars_bg[idx]);
        arr.push_back(fitparerrs_bg[idx]);
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
* @brief Weight a dataset bin by bin from a map of bin ids to weight vectors.
*
* Weight a dataset bin by bin from a map of unique integer bin identifiers to weight vectors.
* The weights are created first from the RDataFrame since defining a conditional weight variable
* for a RooDataSet is nigh impossible. Then, they are merged into your dataset and a new
* weighted dataset is created and uploaded to the workspace.
*
* @param w RooWorkspace in which to work
* @param rds RooDataSet to weight
* @param frame ROOT RDataframe from which to create weights
* @param rds_weighted_name Name of weighted RooDataSet under which to import it into the RooWorkspace
* @param bincuts Map of unique bin id ints to bin variable cuts for bin
* @param sgcut Signal invariant mass region cut
* @param bgcut Background invariant mass region cut
* @param weightvar Weight variable name
* @param weightvar_lims Weight variable limits
* @param weights_map Map of unique integer bin identifiers to weight vectors
* @param weights_default Weight variable default value for events outside provided cuts
* @param use_raw_weights Option to use raw signal and background weights instead of interpreting first entry of weight vectors as the background fraction \f$\varepsilon\f$
*/
void getMassFitWeightedData(
    RooWorkspace                      *w,
    RooAbsData                        *rds,
    RNode                             frame,
    std::string                       rds_weighted_name,
    std::map<int,std::string>         bincuts,
    std::string                       sgcut,
    std::string                       bgcut,
    std::string                       weightvar,
    std::vector<double>               weightvar_lims,
    std::map<int,std::vector<double>> weights_map,
    double                            weights_default = 0.0,
    bool                              use_raw_weights = false
    ) {

    // Create the conditional weight variable formula
    std::string weightvar_formula = "";
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get bin id and cut
        int idx = it->first;
        std::string bincut = it->second;

        // Get signal and background weights directly OR assuming `epss` contains the background fractions eps = N_BG / N_TOT in the massfit signal region
        double sgweight, bgweight;
        if (use_raw_weights) {
            sgweight = weights_map[idx][0];
            bgweight = weights_map[idx][1];
        } else {
            double eps = weights_map[idx][0];
            sgweight = (eps != 0.0) ? 1.0/(1.0-eps) : 0.0;
            bgweight = (eps != 0.0) ? -eps/(1.0-eps) : 0.0;
        }
        
        // Add to the weight variable formula
        if (weightvar_formula.size()==0) {
            weightvar_formula = Form("if (%s && %s) return %.8f; else if (%s && %s) return %.8f;",bincut.c_str(),sgcut.c_str(),sgweight,bincut.c_str(),bgcut.c_str(),bgweight);
        } else {
            std::string bin_condition = Form("else if (%s && %s) return %.8f; else if (%s && %s) return %.8f;",bincut.c_str(),sgcut.c_str(),sgweight,bincut.c_str(),bgcut.c_str(),bgweight);
            weightvar_formula = Form("%s %s",weightvar_formula.c_str(),bin_condition.c_str());
        }
    }

    // Set the default weight condition
    weightvar_formula = Form("%s else return %.8f;",weightvar_formula.c_str(),weights_default);

    // Define weight in dataframe
    auto frame_weighted = frame.Define(weightvar.c_str(),weightvar_formula.c_str());

    // Create weights dataset from dataframe
    if (weightvar_lims.size()!=2) weightvar_lims = {-9999.,9999.};
    RooRealVar wv(weightvar.c_str(),weightvar.c_str(),weightvar_lims[0],weightvar_lims[1]);
    ROOT::RDF::RResultPtr<RooDataSet> rds_weights = frame_weighted.Book<double>(
                RooDataSetHelper(Form("%s_weights",rds->GetName()),"Data Set Weights",RooArgSet(wv)),
                {weightvar.c_str()}
            );

    // Merge dataset with weights
    static_cast<RooDataSet&>(*rds).merge(&static_cast<RooDataSet&>(*rds_weights));

    // Create new data set with weights variable
    auto& data = static_cast<RooDataSet&>(*rds);
    RooDataSet rds_weighted{rds_weighted_name.c_str(), data.GetTitle(), &data, *data.get(), nullptr, weightvar.c_str()};

    // Import weighted dataset into workspace
    w->import(rds_weighted);

} // void getMassFitWeightedData()

/**
* @brief Apply a generic mass fit and weight a dataset bin by bin.
*
* Apply a generic mass fit in each asymmetry fit variable bin
* and weight the given dataset in the signal and background regions.
* Only sideband events will be weighted and the weight \f$W_{sb}\f$ will be computed using the
* background fraction \f$\varepsilon\f$ taken from the integral of the background PDF
* and the sum of the dataset within the signal region from the mass fit.  The weights are scaled from sideband to
* signal region counts to ensure non-negative normalization:
* @f[
* W_{sb} = - \varepsilon \cdot \frac{N_{sg}}{N_{sb}}.
* @f]
*
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param fitvars List of mass fit variables
* @param pdf_name Base name of the combined signal and background pdf.  Note, the actual PDF name will be: `<pdf_name>_<binid>`.
* @param fitformula_sg The signal PDF formula in ROOT TFormula format
* @param fitformula_bg The background PDF formula in ROOT TFormula format
* @param sgYield_name Base name of the signal yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param bgYield_name Base name of the background yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param initsgfrac Initial value of the ratio of signal events to total events \f$u\f$ in the dataset
* @param fitparinits_sg List of signal PDF parameter initial values
* @param fitparnames_sg List of signal PDF parameter names
* @param fitpartitles_sg List of signal PDF parameter titles
* @param fitparunits_sg List of signal PDF parameter unit titles
* @param fitparlims_sg List of signal PDF parameter minimum and maximum bounds
* @param fitparinits_bg List of background PDF parameter initial values
* @param fitparnames_bg List of background PDF parameter names
* @param fitpartitles_bg List of background PDF parameter titles
* @param fitparunits_bg List of background PDF parameter unit titles
* @param fitparlims_bg List of background PDF parameter minimum and maximum bounds
* @param sgregion_lims List of signal region minimum and maximum bounds for each fit variable
* @param frame ROOT RDataframe in which to define the weight variable
* @param bgcut Background invariant mass region cut
* @param asymfitvars List of asymmetry fit variables names
* @param asymfitvar_bincuts Map of unique bin id ints to bin variable cuts for bin
* @param rds_weighted_name Name of weighted RooDataSet under which to import it into the RooWorkspace
* @param weightvar Weight variable name
* @param weightvar_lims Weight variable limits
* @param plot_bg_pars Option to plot background pdf parameters on TLegend
* @param lg_text_size Size of TLegend text
* @param lg_margin Margin of TLegend
* @param lg_ncols Number of columns in TLegend
* @param use_sumw2error Option to use `RooFit::SumW2Error(true)` option when fitting to dataset which is necessary if using a weighted dataset
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_binned_fit Option to use a binned fit to the data
* @param weights_default Weight variable default value for events outside provided cuts
* @param out Output stream
*/
void setWeightsFromMassFit(
        RooWorkspace                    *w, // fitMass() arguments
        std::string                      dataset_name,
        std::string                      binid,
        std::string                      bincut,
        std::vector<std::string>         binvars,
        std::vector<std::string>         fitvars,
        std::string                      pdf_name,
        std::string                      fitformula_sg,
        std::string                      fitformula_bg,
        std::string                      sgYield_name,
        std::string                      bgYield_name,
        double                           initsgfrac,
        std::vector<double>              fitparinits_sg,
        std::vector<std::string>         fitparnames_sg,
        std::vector<std::string>         fitpartitles_sg,
        std::vector<std::string>         fitparunits_sg,
        std::vector<std::vector<double>> fitparlims_sg,
        std::vector<double>              fitparinits_bg,
        std::vector<std::string>         fitparnames_bg,
        std::vector<std::string>         fitpartitles_bg,
        std::vector<std::string>         fitparunits_bg,
        std::vector<std::vector<double>> fitparlims_bg,
        std::vector<std::vector<double>> sgregion_lims,
        RNode                            frame, // arguments for this method
        std::string                      bgcut, 
        std::vector<std::string>         asymfitvars,
        std::map<int,std::string>        asymfitvar_bincuts,
        std::string                      rds_weighted_name,
        std::string                      weightvar,
        std::vector<double>              weightvar_lims,
        double                           lg_text_size     = 0.04, // fitMass() arguments
        double                           lg_margin        = 0.1,
        int                              lg_ncols         = 1,
        bool                             plot_bg_pars     = false,
        bool                             use_sumw2error   = false,
        bool                             use_extended_nll = true,
        bool                             use_binned_fit   = false,
        double                           weights_default  = 0.0, // arguments for this method
        std::ostream                    &out              = std::cout
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

    // Get the signal region cut
    std::string sgcut = saga::util::addLimitCuts("",fitvars,sgregion_lims);

    // Loop bins and apply fits recording background fractions and errors
    std::map<int,std::vector<double>> weights_map;
    for (auto it = asymfitvar_bincuts.begin(); it != asymfitvar_bincuts.end(); ++it) {

        // Get bin id and cut
        int id = it->first;
        std::string asymfitvar_bincut = it->second;

        // Get bin unique id
        std::string binid_asymfitvars = Form("%s__asymfitvars_bin_%d",binid.c_str(),id);

        // Get bin dataset
        RooDataSet *bin_ds_asymfitvar = (RooDataSet*)bin_ds->reduce(asymfitvar_bincut.c_str());

        // Get the signal and sideband counts
        int sg_count = (int)bin_ds_asymfitvar->sumEntries(sgcut.c_str());
        int bg_count = (int)bin_ds_asymfitvar->sumEntries(bgcut.c_str());

        // Create full bin cut
        std::string bincut_full = Form("%s && %s", bincut.c_str(), binid_asymfitvars.c_str());

        // Fit the mass spectrum
        std::vector<double> massfit_result = fitMass(
                w, // RooWorkspace                    *w,
                dataset_name, // std::string                      dataset_name,
                binid_asymfitvars, // std::string                      binid,
                bincut_full, // std::string                      bincut,
                binvars, // std::vector<std::string>         binvars,
                fitvars, // std::vector<std::string>         fitvars,
                pdf_name, // std::string                      pdf_name,
                fitformula_sg, // std::string                      fitformula_sg,
                fitformula_bg, // std::string                      fitformula_bg,
                sgYield_name, // std::string                      sgYield_name,
                bgYield_name, // std::string                      bgYield_name,
                initsgfrac, // double                           initsgfrac,
                fitparinits_sg, // std::vector<double>              fitparinits_sg,
                fitparnames_sg, // std::vector<std::string>         fitparnames_sg,
                fitpartitles_sg, // std::vector<std::string>         fitpartitles_sg,
                fitparunits_sg, // std::vector<std::string>         fitparunits_sg,
                fitparlims_sg, // std::vector<std::vector<double>> fitparlims_sg,
                fitparinits_bg, // std::vector<double>              fitparinits_bg,
                fitparnames_bg, // std::vector<std::string>         fitparnames_bg,
                fitpartitles_bg, // std::vector<std::string>         fitpartitles_bg,
                fitparunits_bg, // std::vector<std::string>         fitparunits_bg,
                fitparlims_bg, // std::vector<std::vector<double>> fitparlims_bg,
                sgregion_lims, // std::vector<std::vector<double>> sgregion_lims,
                lg_text_size, // double                           lg_text_size     = 0.04,
                lg_margin, // double                           lg_margin        = 0.1,
                lg_ncols, // int                              lg_ncols         = 1,
                plot_bg_pars, // bool                             plot_bg_pars     = false,
                use_sumw2error, // bool                             use_sumw2error   = true,
                use_extended_nll, // bool                             use_extended_nll = false,
                use_binned_fit, // bool                             use_binned_fit   = false,
                out // std::ostream                    &out              = std::cout
        );

        // Compute the sideband weights
        double eps = massfit_result[9];//TODO: Maybe add an option to grab the desired background fraction
        double sb_weight = - eps * sg_count / bg_count;
        weights_map[id] = {1.0, sb_weight};
    }

    // Set data set weights from the binned mass fits
    getMassFitWeightedData(
        w,
        ds,
        frame,
        rds_weighted_name,
        asymfitvar_bincuts,
        sgcut,
        bgcut,
        weightvar,
        weightvar_lims,
        weights_map,
        weights_default,
        true
    );

} // void setWeightsFromMassFit(

} // namespace signal {

} // namespace saga {
