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
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TF2.h>
#include <TLatex.h>

// RooFit Includes
#include <RooCategory.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooProduct.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <RooExtendPdf.h>
#include <RooSimultaneous.h>
#include <RooFFTConvPdf.h>
#include <RooCrystalBall.h>
#include <RooLandau.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>

// RooStats includes
#include <RooStats/SPlot.h>

// Local includes
#include <data.h>
#include <bins.h>
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

namespace massfit {

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
* @param sgYield_name Base name of the signal yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param bgYield_name Base name of the background yield variable for an extended fit, the bin id will be appended for uniqueness in the workspace
* @param method_name Method name, used to name the PDF
* @param binid Unique bin id, used to name the PDF
* @param fitformula_sg The signal PDF formula in ROOT TFormula format
* @param fitformula_bg The background PDF formula in ROOT TFormula format
* @param initsgfrac Initial value of the ratio of signal events to total events \f$u\f$ in the dataset
* @param count Bin count
* @param use_extended_nll Option to use an extended likelihood term
* 
* @return std::string
*/
std::string getGenMassPdf(
    RooWorkspace *w,
    RooArgSet *argset_sg,
    RooArgSet *argset_bg,
    std::string sgYield_name,
    std::string bgYield_name,
    std::string method_name,
    std::string binid,
    std::string fitformula_sg,
    std::string fitformula_bg,
    double initsgfrac,
    int count,
    bool use_extended_nll
) {

    // Create the summed PDF
    RooAddPdf * model;
    std::string model_name = Form("model_%s_%s",method_name.c_str(),binid.c_str());

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

    w->import(sgYield);
    w->import(bgYield);
    w->import(*model);
    return model_name;

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
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param fitvars List of mass fit variables
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
        bool use_sumw2error              = true,
        bool use_extended_nll            = false,
        bool use_binned_fit              = false,
        std::ostream &out                = std::cout
    ) {

    // Set method name
    std::string method_name = "fitMass";

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Load fit variables from workspace
    RooRealVar * f[(const int)fitvars.size()];
    for (int i=0; i<fitvars.size(); i++) {
        f[i] = w->var(fitvars[i].c_str());
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
    for (int aa=0; aa<nparams; aa++) { // Fit parameters
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
    for (int aa=0; aa<nparams; aa++) { // Fit parameters
        argset_bg->add(*a_bg[aa]);
    }

    // Create and load combined signal and background mass PDF
    std::string model_name = getGenMassPdf(
        w,
        argset_sg,
        argset_bg,
        sgYield_name,
        bgYield_name,
        method_name,
        binid,
        fitformula_sg,
        fitformula_bg,
        initsgfrac,
        count,
        use_extended_nll
    );
    RooAbsPdf *model = w->pdf(model_name.c_str());
    RooAbsPdf *sg    = model->pdfList()[0];
    RooAbsPdf *bg    = model->pdfList()[1];

    // Fit the PDF to data
    std::unique_ptr<RooFitResult> r;
    if (use_binned_fit) {

        // Create binned data
        std::unique_ptr<RooDataHist> dh = (std::unique_ptr<RooDataHist>)bin_ds->binnedClone();

        // Fit PDF
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*dh, RooFit::Save(), RooFit::SumW2Error(use_sumw2error), RooFit::PrintLevel(-1));

    } else {

        // Fit PDF
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*bin_ds, RooFit::Save(), RooFit::SumW2Error(use_sumw2error), RooFit::PrintLevel(-1));
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
        rdhs_1d[i] = new RooDataHist(dh_title.c_str(), dh_title.c_str(), *f[i], *bin_ds);

        // Compute chi2 value
        RooFit::OwningPtr<RooAbsReal> chi2 = model->createChi2(*rdhs_1d[i], Range("fullRange"),
                    Extended(use_extended_nll), DataError(RooAbsData::Poisson));
        int nparameters = (int)model->getParameters(RooArgSet(*f[i]))->size();
        int ndf = f[i]->getBins() - nparameters; //NOTE: ASSUME ALL BINS NONZERO
        double chi2ndf = (double) chi2->getVal()/ndf;
        chi2ndfs.push_back(chi2ndf);
    }

    //TODO: Integrate and calculate background fractions...

    //---------------------------------------- Compute epsilon ----------------------------------------//
    // Bg hist
    RooFit::OwningPtr<RooDataHist> bghist_roofit = bg.generateBinned(RooArgSet(*m), (double)bgYield.getVal());
    TH1F *bghist = (TH1F*)bghist_roofit->createHistogram("bghist",*m); //bg->GetHistogram();
    bghist->SetTitle("");
    bghist->SetBins(mass_nbins_hist,mass_min,mass_max);
    bghist->SetLineWidth(1);
    bghist->SetLineColor(kAzure);
    bghist->Draw("SAME PE"); // Rebinning is a pain and fn is automatically binned to 100, so just keep 100 bins.

    // Signal hist
    TH1D *hist = (TH1D*)h->Clone("hist");
    hist->Add(bghist,-1);
    hist->SetLineWidth(1);
    hist->SetLineColor(8);
    hist->Draw("SAME PE");

    // Integrate signal histogram
    auto i_sig_err = 0.0;
    auto i_sig = hist->IntegralAndError(hist->FindBin(sg_region_min),hist->FindBin(sg_region_max),i_sig_err);

    // Define a range named "signal" in x from -5,5
    m->setRange("signal", sg_region_min, sg_region_max);
 
    // Integrate the signal PDF
    std::unique_ptr<RooAbsReal> igm_sig{sig->createIntegral(*m, NormSet(*m), Range("signal"))};
    RooProduct i_sg_yield{"i_sg_yield", "i_sg_yield", {*igm_sig, sgYield}};
    Double_t integral_sg_value = i_sg_yield.getVal();
    Double_t integral_sg_value_error = i_sg_yield.getPropagatedError(*fit_result_data, *m); // Note Need fit result saved for this to work!

    // Integrate the background PDF
    std::unique_ptr<RooAbsReal> igm_bg{bg.createIntegral(*m, NormSet(*m), Range("signal"))};
    RooProduct i_bg_yield{"i_bg_yield", "i_bg_yield", {*igm_bg, bgYield}};
    Double_t integral_bg_value = i_bg_yield.getVal();
    Double_t integral_bg_value_error = i_bg_yield.getPropagatedError(*fit_result_data, *m); // Note Need fit result saved for this to work!                                                                                                                

    // Get Full Model PDF integral
    std::unique_ptr<RooAbsReal> igm_model{model.createIntegral(*m, NormSet(*m), Range("signal"))};
    RooRealVar model_yield("model_yield", "model_yield",sgYield.getVal()+bgYield.getVal());
    Double_t integral_model_value = igm_model->getVal() * model_yield.getVal();
    Double_t integral_model_value_error = igm_model->getPropagatedError(*fit_result_data, *m) * model_yield.getVal(); // Note Need fit result saved for this to work!                                                                                                                 

    // Sum on dataset
    std::string signal_cut = Form("%.8f<%s && %s<%.8f",(double)sg_region_min,m->GetName(),m->GetName(),(double)sg_region_max);
    double i_ds = rds->sumEntries(signal_cut.c_str());
    double i_ds_err = TMath::Sqrt(i_ds);

    // Compute epsilon in a couple different ways
    double eps_pdf_int = integral_bg_value / i_ds;
    double eps_sg_th1_int = 1.0 - i_sig / i_ds ;
    double eps_sg_pdf_int = 1.0 - integral_sg_value / i_ds ;

    // Compute epsilon errors in a couple different ways
    double eps_pdf_int_err = integral_bg_value_error / i_ds;
    double eps_sg_th1_int_err = i_sig_err / i_ds ;
    double eps_sg_pdf_int_err = integral_sg_value_error / i_ds ;

    // // Now compute true epsilon
    // double true_bg_count = (double)*frame.Filter(mccuts_true_bg.c_str()).Filter(signal_cut.c_str()).Count();
    // double true_full_count = (double)*frame.Filter(signal_cut.c_str()).Count();
    // double eps_true = (double) true_bg_count / true_full_count;

    //TODO: Plot mass fit

        // Plot invariant mass fit from RooFit
    RooPlot *mframe_1d = m->frame(Title("1D PDF fit mass_ppim."));
    rds->plotOn(mframe_1d);
    model.plotOn(mframe_1d);
    model.plotOn(mframe_1d, Components(*sig), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(mframe_1d, Components(bg), LineStyle(kDashed), LineColor(kBlue));
    TCanvas *c_massfit = new TCanvas("c_massfit");
    c_massfit->cd();
    gPad->SetLeftMargin(0.15);
    mframe_1d->GetYaxis()->SetTitleOffset(1.6);
    mframe_1d->Draw();

    // Plot sig and bg histograms
    bghist->Draw("SAME");
    hist->Draw("SAME");

    // Create Legend Entries
    TString s_chi2, s_ntot, s_nbg, s_epsilon;
    std::string str_chi2 = Form("#chi^{2}/NDF = %.3g",chi2ndf);
    std::string str_ntot = Form("N_{Tot} = %.2e #pm %.0f",i_ds,i_ds_err);
    std::string str_nbg  = Form("N_{bg} = %.2e #pm %.0f",integral_bg_value,integral_bg_value_error);
    std::string str_eps  = Form("#varepsilon = %.3f #pm %.3f",eps_pdf_int,eps_pdf_int_err);

    // Draw Legend
    //TODO: Set legend location from arguments
    TLegend *legend=new TLegend(0.45,0.2,0.875,0.625);
    legend->SetTextSize(0.04);
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, str_chi2.c_str(), Form(" %g ",chi2ndf));
    legend->AddEntry((TObject*)0, str_ntot.c_str(), Form(" %g ",chi2ndf));
    legend->AddEntry((TObject*)0, str_nbg.c_str(),  Form(" %g ",chi2ndf));
    legend->AddEntry((TObject*)0, str_eps.c_str(),  Form(" %g ",chi2ndf));

    //TODO: Add signal PDF parameter values and errors
    //TODO: Add histograms optionally
    //TODO: Add pdfs optionally
    for (int i=0; i<nparams_sg; i++) {
        RooRealVar *par = a_sg[i];
        std::string par_str = Form("%s = %.3g #pm %.3g %s", a_sg[i]->GetTitle(), a_sg[i]->getVal(), a_sg[i]->getError(), fitparunits_sg[i].c_str());
        legend->AddEntry((TObject*)0, par_str.c_str(), Form(" %g ",chi2ndf));
    }
    if (plot_bg_pars) {
        for (int i=0; i<nparams_bg; i++) {
            RooRealVar *par = a_bg[i];
            std::string par_str = Form("%s = %.3g #pm %.3g %s", a_bg[i]->GetTitle(), a_bg[i]->getVal(), a_bg[i]->getError(), fitparunits_bg[i].c_str());
            legend->AddEntry((TObject*)0, par_str.c_str(), Form(" %g ",chi2ndf));
        }
    }
    legend->Draw();

    // Save Canvas
    c_massfit->SaveAs(Form("%s_%s_%s.pdf",c_massfit->GetName(),sig_pdf_name.c_str(),bin_id.c_str()));

    //TODO: Output integration and background fraction values in messages below

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
    out << "--------------------------------------------------" << std::endl;

    // Fill return array
    std::vector<double> arr; //NOTE: Dimension = 5+chi2ndfs.size()+2*binvars.size()+2*fitparinits_sg.size()+2*fitparinits_bg.size()
    arr.push_back(count);
    arr.push_back(sgYield_sgregion);
    arr.push_back(bgYield_sgregion);
    arr.push_back(epsilon);
    arr.push_back(epsilon_err);
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
* @brief Apply a \f$\Lambda\f$ baryon mass fit
*
* Apply a \f$\Lambda\f$ baryon mass fit with signal function chosen from
* (`"gauss"`, `"landau"`, `"cb"`, `"landau_X_gauss"`, `"cb_X_gauss"`, `"cb_gauss"`) and
* Chebychev polynomial background function and save model and
* yield variables to workspace for use with sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* This will also return \f$\varepsilon\f$ which is the fraction of events
* in the signal region (`sig_region_min`,`sig_region_max`) which
* are background.  The background fraction is computed using the sums of the observed
* distribution and the histogrammed background function in the signal region respectively.
*
* @param w RooWorkspace in which to work
* @param massfitvars Invariant mass fit variable names
* @param rds RooDataSet to use for fit
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param frame ROOT RDataframe from which to create a histogram of invariant mass
* @param mass_nbins_conv Number of invariant mass bins for convolution of PDFs
* @param model_name Full PDF name
* @param sig_pdf_name Signal PDF name
* @param sg_region_min Invariant mass signal region lower bound
* @param sg_region_max Invariant mass signal region upper bound
* @param use_poly4_bg Use a 4th order Chebychev polynomial background instead of 2nd order
* @param bin_id Unique bin identifier string
*
* @return List containing background fraction \f$\varepsilon\f$ and its statistical error
*/
std::vector<double> applyLambdaMassFit(
    RooWorkspace *w,
    std::vector<std::string> massfitvars,
    RooAbsData *rds,
    std::string sgYield_name,
    std::string bgYield_name,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    int mass_nbins_conv,
    std::string model_name,
    std::string sig_pdf_name,
    double sg_region_min,
    double sg_region_max,
    int use_poly4_bg,
    std::string bin_id
    ) {

    using namespace RooFit;

    std::string massvar = massfitvars[0]; //NOTE: ASSUME ONE MASS FIT VARIABLE.

    // Get variables from workspace
    RooRealVar *m = w->var(massvar.c_str());
    double mass_min = m->getMin();
    double mass_max = m->getMax();
    int    mass_nbins_hist =  m->getBins();

    // Get dataset length
    int count = (int)rds->numEntries();

    // Create histogram
    auto h1 = (TH1D) *frame.Histo1D({"h1",massvar.c_str(),mass_nbins_hist,mass_min,mass_max},massvar.c_str());
    TH1D *h = (TH1D*)h1.Clone(massvar.c_str());
    h->SetTitle("");
    h->GetXaxis()->SetTitle(Form("%s (GeV)",m->GetTitle()));
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetYaxis()->SetTitle("Counts");
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.87);

    // Set y axes limits
    h->GetYaxis()->SetRangeUser(0.0,1.05*h->GetMaximum());

    // Construct signal gauss(t,mg,sg)
    RooRealVar mg("mg", "mg", 1.1157, mass_min, mass_max);
    RooRealVar sg("sg", "sg", 0.008, 0.0, 0.1);
    RooGaussian gauss(Form("gauss%s",bin_id.c_str()), "gauss", *m, mg, sg);

    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml", "mean landau", 1.1157, mass_min, mass_max);
    RooRealVar sl("sl", "sigma landau", 0.005, 0.0, 0.1);
    RooLandau landau(Form("landau%s",bin_id.c_str()), "landau", *m, ml, sl);

    // Construct signal parameters and function
    RooRealVar mu("mu","#mu",1.1157,0.0,2.0);
    RooRealVar s("s","#sigma",0.008,0.0,1.0);
    RooRealVar a_left("a_left","alpha_left",10.0);
    RooRealVar n_left("n_left","n_left",10.0);
    RooRealVar a("a","#alpha",1.0,0.0,3.0);
    RooRealVar n("n","n",2.0,2.0,10.0);
    RooCrystalBall cb(Form("cb%s",bin_id.c_str()), "crystal_ball", *m, mu, s, a_left, n_left, a, n); //NOTE: Include model name for uniqueness within bin.

    // Construct addition component pdfs
    RooRealVar sg_add("sg_add", "sg_add", 0.0001, 0.0, 0.001);
    RooGaussian gauss_add(Form("gauss_add%s",bin_id.c_str()), "gauss_add", *m, mu, sg_add);
    RooCrystalBall cb_add(Form("cb_add%s",bin_id.c_str()), "crystal_ball", *m, mu, s, a_left, n_left, a, n); //NOTE: Include model name for uniqueness within bin.
    double cbFrac_init = 0.1;
    RooRealVar cbFrac(Form("cbFrac%s",bin_id.c_str()), "fitted yield for signal", cbFrac_init, 0., 1.0);
    RooAddPdf cb_gauss(Form("cb_gauss%s",bin_id.c_str()), Form("cb_add%s+gauss_add%s",bin_id.c_str(),bin_id.c_str()), RooArgList(cb_add,gauss_add), RooArgList(cbFrac)); //NOTE: N-1 Coefficients! 

    // Construct convolution component pdfs
    RooRealVar mg_conv("mg_conv", "mg_conv", 0.0);
    RooRealVar sg_conv("sg_conv", "sg_conv", 0.008, 0.0, 0.1);
    RooGaussian gauss_landau_conv(Form("gauss_landau_conv%s",bin_id.c_str()), "gauss_landau_conv", *m, mg_conv, sg_conv);
    RooGaussian gauss_cb_conv(Form("gauss_cb_conv%s",bin_id.c_str()), "gauss_cb_conv", *m, mg_conv, sg_conv);
    RooLandau landau_conv(Form("landau_conv%s",bin_id.c_str()), "landau_conv", *m, ml, sl);
    RooCrystalBall cb_conv(Form("cb_conv%s",bin_id.c_str()), "crystal_ball_conv", *m, mu, s, a_left, n_left, a, n);
    
    // Set #bins to be used for FFT sampling to 10000
    m->setBins(mass_nbins_conv, "cache");
    
    // Construct Convolution PDFs
    RooFFTConvPdf landau_X_gauss(Form("landau_X_gauss%s",bin_id.c_str()), "CB (X) gauss_conv", *m, landau_conv, gauss_landau_conv);
    RooFFTConvPdf cb_X_gauss(Form("cb_X_gauss%s",bin_id.c_str()), "CB (X) gauss_conv", *m, cb_conv, gauss_cb_conv);

    // Import signal functions to workspace
    w->import(gauss);
    w->import(landau);
    w->import(cb);
    w->import(landau_X_gauss);
    w->import(cb_X_gauss);
    w->import(cb_gauss);

    // Pick out signal function based on preference
    std::string sig_pdf_name_unique = Form("%s%s",sig_pdf_name.c_str(),bin_id.c_str());
    RooAbsPdf *sig = w->pdf(sig_pdf_name_unique.c_str());

    // Consruct background parameters and function
    RooRealVar b1("b1","b_{1}",  0.72,-10.0,10.0);
    RooRealVar b2("b2","b_{2}", -0.17,-10.0,10.0);
    RooRealVar b3("b3","b_{3}",  0.05,-10.0,10.0);
    RooRealVar b4("b4","b_{4}", -0.01,-10.0,10.0);
    std::string bg_pdf_name_unique = Form("bg%s",bin_id.c_str());
    RooChebychev bg(bg_pdf_name_unique.c_str(),bg_pdf_name_unique.c_str(),*m,(use_poly4_bg==1 ? RooArgList(b1,b2,b3,b4) : RooArgList(b1,b2)));
    
    // Combine signal and background functions
    double sgfrac = 0.1;
    double sgYield_init = sgfrac * count;
    double bgYield_init = (1.0-sgfrac) * count;
    RooRealVar sgYield(Form("%s%s",sgYield_name.c_str(),bin_id.c_str()), "fitted yield for signal", sgYield_init, 0., 2.0*count);
    RooRealVar bgYield(Form("%s%s",bgYield_name.c_str(),bin_id.c_str()), "fitted yield for background", bgYield_init, 0., 2.0*count);
    RooAddPdf model(Form("%s%s",model_name.c_str(),bin_id.c_str()), Form("%s+%s",sig_pdf_name_unique.c_str(),bg_pdf_name_unique.c_str()), RooArgList(*sig,bg), RooArgList(sgYield,bgYield)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit

    // Fit invariant mass spectrum
    std::unique_ptr<RooFitResult> fit_result_data{model.fitTo(*rds, Save(), PrintLevel(-1))};

    //---------------------------------------- Compute chi2 ----------------------------------------//
    // Import TH1 histogram into RooDataHist
    RooDataHist dh("dh", "dh", *m, Import(*h));

    // Compute chi2 value
    RooFit::OwningPtr<RooAbsReal> chi2 = model.createChi2(dh, Range("fullRange"),
                 Extended(true), DataError(RooAbsData::Poisson));
    int nparameters = (int) model.getParameters(RooArgSet(*m))->size();
    int ndf = mass_nbins_hist - nparameters; //NOTE: ASSUME ALL BINS NONZERO
    double chi2ndf = (double) chi2->getVal()/ndf;

    //---------------------------------------- Compute epsilon ----------------------------------------//
    // Bg hist
    RooFit::OwningPtr<RooDataHist> bghist_roofit = bg.generateBinned(RooArgSet(*m), (double)bgYield.getVal());
    TH1F *bghist = (TH1F*)bghist_roofit->createHistogram("bghist",*m); //bg->GetHistogram();
    bghist->SetTitle("");
    bghist->SetBins(mass_nbins_hist,mass_min,mass_max);
    bghist->SetLineWidth(1);
    bghist->SetLineColor(kAzure);
    bghist->Draw("SAME PE"); // Rebinning is a pain and fn is automatically binned to 100, so just keep 100 bins.

    // Signal hist
    TH1D *hist = (TH1D*)h->Clone("hist");
    hist->Add(bghist,-1);
    hist->SetLineWidth(1);
    hist->SetLineColor(8);
    hist->Draw("SAME PE");

    // Integrate signal histogram
    auto i_sig_err = 0.0;
    auto i_sig = hist->IntegralAndError(hist->FindBin(sg_region_min),hist->FindBin(sg_region_max),i_sig_err);

    // Define a range named "signal" in x from -5,5
    m->setRange("signal", sg_region_min, sg_region_max);
 
    // Integrate the signal PDF
    std::unique_ptr<RooAbsReal> igm_sig{sig->createIntegral(*m, NormSet(*m), Range("signal"))};
    RooProduct i_sg_yield{"i_sg_yield", "i_sg_yield", {*igm_sig, sgYield}};
    Double_t integral_sg_value = i_sg_yield.getVal();
    Double_t integral_sg_value_error = i_sg_yield.getPropagatedError(*fit_result_data, *m); // Note Need fit result saved for this to work!

    // Integrate the background PDF
    std::unique_ptr<RooAbsReal> igm_bg{bg.createIntegral(*m, NormSet(*m), Range("signal"))};
    RooProduct i_bg_yield{"i_bg_yield", "i_bg_yield", {*igm_bg, bgYield}};
    Double_t integral_bg_value = i_bg_yield.getVal();
    Double_t integral_bg_value_error = i_bg_yield.getPropagatedError(*fit_result_data, *m); // Note Need fit result saved for this to work!                                                                                                                

    // Get Full Model PDF integral
    std::unique_ptr<RooAbsReal> igm_model{model.createIntegral(*m, NormSet(*m), Range("signal"))};
    RooRealVar model_yield("model_yield", "model_yield",sgYield.getVal()+bgYield.getVal());
    Double_t integral_model_value = igm_model->getVal() * model_yield.getVal();
    Double_t integral_model_value_error = igm_model->getPropagatedError(*fit_result_data, *m) * model_yield.getVal(); // Note Need fit result saved for this to work!                                                                                                                 

    // Sum on dataset
    std::string signal_cut = Form("%.8f<%s && %s<%.8f",(double)sg_region_min,m->GetName(),m->GetName(),(double)sg_region_max);
    double i_ds = rds->sumEntries(signal_cut.c_str());
    double i_ds_err = TMath::Sqrt(i_ds);

    // Compute epsilon in a couple different ways
    double eps_pdf_int = integral_bg_value / i_ds;
    double eps_sg_th1_int = 1.0 - i_sig / i_ds ;
    double eps_sg_pdf_int = 1.0 - integral_sg_value / i_ds ;

    // Compute epsilon errors in a couple different ways
    double eps_pdf_int_err = integral_bg_value_error / i_ds;
    double eps_sg_th1_int_err = i_sig_err / i_ds ;
    double eps_sg_pdf_int_err = integral_sg_value_error / i_ds ;

    // // Now compute true epsilon
    // double true_bg_count = (double)*frame.Filter(mccuts_true_bg.c_str()).Filter(signal_cut.c_str()).Count();
    // double true_full_count = (double)*frame.Filter(signal_cut.c_str()).Count();
    // double eps_true = (double) true_bg_count / true_full_count;

    // Plot invariant mass fit from RooFit
    RooPlot *mframe_1d = m->frame(Title("1D pdf fit mass_ppim."));
    rds->plotOn(mframe_1d);
    model.plotOn(mframe_1d);
    model.plotOn(mframe_1d, Components(*sig), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(mframe_1d, Components(bg), LineStyle(kDashed), LineColor(kBlue));
    TCanvas *c_massfit = new TCanvas("c_massfit");
    c_massfit->cd();
    gPad->SetLeftMargin(0.15);
    mframe_1d->GetYaxis()->SetTitleOffset(1.6);
    mframe_1d->Draw();

    // Plot sig and bg histograms
    bghist->Draw("SAME");
    hist->Draw("SAME");

    // Create Legend Entries
    TString s_chi2, s_ntot, s_nbg, s_epsilon;
    s_chi2.Form("#chi^{2}/NDF = %.2f",chi2ndf);
    s_ntot.Form("N_{Tot} = %.2e #pm %.0f",i_ds,i_ds_err);
    s_nbg.Form("N_{bg} = %.2e #pm %.0f",integral_bg_value,integral_bg_value_error);
    s_epsilon.Form("#varepsilon = %.3f #pm %.3f",eps_pdf_int,eps_pdf_int_err);
    TString s_mg, s_sg;
    s_mg.Form("#mu_{Gaus} = %.4f #pm %.4f GeV",mg.getVal(),mg.getError());
    s_sg.Form("#sigma_{Gaus} = %.4f #pm %.4f GeV",sg.getVal(),sg.getError());
    TString s_mg_conv, s_sg_conv;
    s_mg_conv.Form("#mu_{Gaus} = %.4f #pm %.4f GeV",mg_conv.getVal(),mg_conv.getError());
    s_sg_conv.Form("#sigma_{Gaus} = %.4f #pm %.4f GeV",sg_conv.getVal(),sg_conv.getError());
    TString s_mg_add, s_sg_add;
    s_sg_add.Form("#sigma_{Gaus} = %.4f #pm %.4f GeV",sg_add.getVal(),sg_add.getError());
    TString s_ml, s_sl;
    s_ml.Form("#mu_{Landau} = %.4f #pm %.4f GeV",ml.getVal(),ml.getError());
    s_sl.Form("#sigma_{Landau} = %.4f #pm %.4f GeV",sl.getVal(),sl.getError());
    TString s_alpha, s_n, s_sigma, s_mu, s_c1;
    s_alpha.Form("#alpha = %.3f #pm %.3f",a.getVal(),a.getError());
    s_n.Form("n = %.2f #pm %.2f",n.getVal(),n.getError());
    s_sigma.Form("#sigma = %.5f #pm %.5f GeV",s.getVal(),s.getError());
    s_mu.Form("#mu = %.5f #pm %.2f GeV",mu.getVal(),mu.getError());
    s_c1.Form("C = %.5f #pm %.5f GeV",sgYield.getVal(),sgYield.getError());

    // Draw Legend
    TLegend *legend=new TLegend(0.45,0.2,0.875,0.625); //NOTE: FOR WITHOUT MC DECOMP
    legend->SetTextSize(0.04);
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, s_chi2, Form(" %g ",chi2ndf));
    if (sig_pdf_name=="gauss") {
        legend->AddEntry((TObject*)0, s_mg, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sg, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="landau") {
        legend->AddEntry((TObject*)0, s_ml, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sl, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="cb") {
        legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_n, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mu, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="landau_X_gauss") {
        legend->AddEntry((TObject*)0, s_ml, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sl, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mg_conv, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sg_conv, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="cb_X_gauss") {
        legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_n, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mu, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mg_conv, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sg_conv, Form(" %g ",chi2ndf));
    }
    else if (sig_pdf_name=="cb_gauss") {
        legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_n, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mu, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sg_add, Form(" %g ",chi2ndf));
    }
    else {
        legend->AddEntry((TObject*)0, s_alpha, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_n, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_sigma, Form(" %g ",chi2ndf));
        legend->AddEntry((TObject*)0, s_mu, Form(" %g ",chi2ndf));
    }
    legend->AddEntry((TObject*)0, s_ntot, Form(" %g ",chi2ndf));
    legend->AddEntry((TObject*)0, s_nbg, Form(" %g ",chi2ndf));
    legend->AddEntry((TObject*)0, s_epsilon, Form(" %g ",chi2ndf));
    legend->Draw();

    // Save Canvas
    c_massfit->SaveAs(Form("%s_%s_%s.pdf",c_massfit->GetName(),sig_pdf_name.c_str(),bin_id.c_str()));

    // Add yield variables to workspace
    w->import(sgYield);
    w->import(bgYield);

    // Add model to workspace
    w->import(model);

    // Return background fraction and error
    std::vector<double> result;
    result.push_back(eps_sg_th1_int);
    result.push_back(eps_sg_th1_int_err);
    return result;
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
    RooWorkspace                                                 *w,
    RooAbsData                                                   *rds,
    RNode frame,
    std::string                                                   rds_weighted_name,
    std::map<int,std::string>                                     bincuts,
    std::string                                                   sgcut,
    std::string                                                   bgcut,
    std::string                                                   weightvar,
    std::vector<double>                                           weightvar_lims,
    std::map<int,std::vector<double>>                             weights_map,
    double                                                        weights_default = 0.0,
    bool                                                          use_raw_weights = false
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
* @brief Apply a \f$\Lambda\f$ mass fit and weight a dataset bin by bin.
*
* Apply a \f$\Lambda\f$ mass fit in each asymmetry fit variable bin
* and weight the given dataset in the signal and background regions.
* Only sideband events will be weighted and the weight will be computed using the
* background fraction \f$\varepsilon\f$ from the mass fit and scaled from sideband to
* signal region counts to ensure non-negative normalization:
* @f[
* W_{sb} = - \varepsilon \cdot \frac{N_{sg}}{N_{sb}}.
* @f]
*
* @param w RooWorkspace in which to work
* @param rds RooDataSet to fit and weight
* @param massfitvars Invariant mass fit variable names
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param frame ROOT RDataframe from which to create a histogram of invariant mass
* @param massfit_nbins_conv Number of invariant mass bins for convolution of PDFs
* @param massfit_model_name Full PDF name
* @param massfit_sig_pdf_name Signal PDF name
* @param massfit_sg_region_min Invariant mass signal region lower bound
* @param massfit_sg_region_max Invariant mass signal region upper bound
* @param use_poly4_bg Use a 4th order Chebychev polynomial background instead of 2nd order
* @param bin_id Unique bin identifier string
* @param sgcut Signal invariant mass region cut
* @param bgcut Background invariant mass region cut
* @param asymfitvars List of asymmetry fit variables names
* @param bincuts Map of unique bin id ints to bin variable cuts for bin
* @param rds_weighted_name Name of weighted RooDataSet under which to import it into the RooWorkspace
* @param weightvar Weight variable name
* @param weightvar_lims Weight variable limits
* @param weights_default Weight variable default value for events outside provided cuts
*/
void setWeightsFromLambdaMassFit(
    RooWorkspace                                                 *w,
    RooAbsData                                                   *rds,
    std::vector<std::string>                                      massfitvars,
    std::string                                                   sgYield_name,
    std::string                                                   bgYield_name,
    RNode frame,
    int                                                           massfit_nbins_conv,
    std::string                                                   massfit_model_name,
    std::string                                                   massfit_sig_pdf_name,
    double                                                        massfit_sg_region_min,
    double                                                        massfit_sg_region_max,
    int                                                           use_poly4_bg,
    std::string                                                   bin_id,
    std::string                                                   sgcut,
    std::string                                                   bgcut,
    std::vector<std::string>                                      asymfitvars,
    std::map<int,std::string>                                     bincuts,
    std::string                                                   rds_weighted_name,
    std::string                                                   weightvar,
    std::vector<double>                                           weightvar_lims,
    double                                                        weights_default = 0.0
    ) {

    // Load asymmetry fit variables from workspace
    RooRealVar *rrvars[asymfitvars.size()];
    for (int idx=0; idx<asymfitvars.size(); idx++) {
        rrvars[idx] = w->var(asymfitvars[idx].c_str());
    }

    // Create binning scheme if not provided
    if (bincuts.size()==0) {

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
        bincuts = saga::bins::getBinCuts(binscheme,0);
    }

    // Loop bins and apply fits recording background fractions and errors
    std::map<int,std::vector<double>> weights_map;
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get bin id and cut
        int id = it->first;
        std::string bincut = it->second;

        // Get bin dataset and frame
        RooDataSet *bin_ds = (RooDataSet*)rds->reduce(bincut.c_str());
        auto binframe = frame.Filter(bincut.c_str());

        // Get bin unique id
        std::string bin_unique_id = Form("%s__asymfitvars_bin_%d",bin_id.c_str(),id);

        // Get the signal and sideband counts
        int sg_count = (int)*binframe.Filter(sgcut.c_str()).Count();
        int bg_count = (int)*binframe.Filter(bgcut.c_str()).Count();

        // Get Lambda mass fit
        std::vector<double> massfit_result = applyLambdaMassFit(
                w,
                massfitvars,
                bin_ds,
                Form("%s_%s",sgYield_name.c_str(),bin_unique_id.c_str()),
                Form("%s_%s",bgYield_name.c_str(),bin_unique_id.c_str()),
                binframe,
                massfit_nbins_conv,
                massfit_model_name,
                massfit_sig_pdf_name,
                massfit_sg_region_min,
                massfit_sg_region_max,
                use_poly4_bg,
                bin_unique_id
        );

        // Compute the sideband weights
        double eps = massfit_result[0];
        double sb_weight = - eps * sg_count / bg_count;
        weights_map[id] = {1.0, sb_weight};
    }

    // Set data set weights from the binned mass fits
    getMassFitWeightedData(
        w,
        rds,
        frame,
        rds_weighted_name,
        bincuts,
        sgcut,
        bgcut,
        weightvar,
        weightvar_lims,
        weights_map,
        weights_default,
        true
    );

} // void setWeightsFromLambdaMassFit(

} // namespace analysis {

} // namespace saga {
