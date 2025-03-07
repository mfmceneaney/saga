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

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 12/Dec./2024
* @version 0.0.0
* @brief Fit asymmetries using RooFit unbinned extended Maximum Likelihood methods and sideband subtraction 
* or the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a> for background correction.
*/

namespace saga {

namespace analysis {

/**
* @brief Apply a \f$\Lambda\f$ baryon mass fit
*
* Apply a \f$\Lambda\f$ baryon mass fit with signal function chosen from
* (`"gauss"`, `"landau"`, `"cb"`, `"landau_X_gauss"`, `"cb_X_gauss"`, `"cb_gauss"`) and
* Chebychev polynomial background function and save model and
* yield variables to workspace for use with sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* This will also return \f$\varepsilon\f$ which is the fraction of events
* in the signal region (`sig_region_min`,`sig_region_max`) which
* are background based on the difference between the observed
* distribution and the histogrammed background function.
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
* @return List containing background fraction epsilon and its statistical error
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
    c_massfit->SaveAs(Form("%s_%s_%s_%s.pdf",c_massfit->GetName(),sig_pdf_name.c_str(),bin_id.c_str(),bin_id.c_str()));

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
* @param use_raw_weights Option to raw signal and background weights instead of interpretting first entry of weight vectors as the background fraction \f$\varepsilon\f$
*/
void getMassFitWeightedData(
    RooWorkspace                                                 *w,
    RooAbsData                                                   *rds,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
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
* background fraction $\varepsilon$ from the mass fit and scaled from sideband to
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
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
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
        massfit_results,
        weights_default,
        true
    );

} // void setWeightsFromLambdaMassFit(

/**
* @brief Apply the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
*
* Apply sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a> given a dataset, yield variables, and a model and 
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
* @brief Fit an asymmetry.
*
* Compute the bin count, bin variable mean values and variances, depolarization variable values and errors,
* and fit the asymmetry with a binned or unbinned dataset using a maximum likelihood fit method with an optional extended likelihood term.
* Note that for the maximum likelihood fit, the given asymmetry formula \f$ A(x_0, x_1, ..., a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) \f$ will be converted internally using
* <a href="https://root.cern.ch/doc/master/classRooGenericPdf.html">RooGenericPdf</a> to a PDF of the form:
*
* @f[
* \begin{aligned}
* PDF(h, x_0, x_1, ..., &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + h \cdot P \cdot A(x_0, x_1, ..., a_0, a_1, a_2, ..., d_0, d_1, d_2, ...),
* \end{aligned}
* @f]
* and a simultaneous fit will be applied over the data subsets distinguished by the helicity states.  The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the corresponding depolarization factors.
*
* The variable names in the fit formula should follow the <a href="https://root.cern.ch/doc/master/classTFormula.html">TFormula</a> notation, e.g.,
* `x_0`\f$\rightarrow\f$`x[0]`, `x_1`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[N_x]`, `a_1`\f$\rightarrow\f$`x[N_x+1]`, etc.
*
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param pol Luminosity averaged beam polarization
* @param helicity Name of helicity variable
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param depolvars List of depolarization variables
* @param fitvars List of asymmetry fit variables
* @param fitformula The asymmetry formula in ROOT TFormula format
* @param initparams List of initial values for asymmetry parameters
* @param initparamlims List of initial asymmetry parameter minimum and maximum bounds
* @param use_sumw2error Option to use RooFit::SumW2Error(true) option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_binned_fit Option to use a binned fit to the data
* @param out Output stream
*
* @return List of bin count, bin variable means and errors, depolarization variable means and errors, fit parameters and errors
*/
std::vector<double> fitAsym(
        RooWorkspace                    *w,
        std::string                      dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
        double                           pol,
        std::string                      helicity,
        std::string                      binid,
        std::string                      bincut,
        std::vector<std::string>         binvars,
        std::vector<std::string>         depolvars,
        std::vector<std::string>         fitvars,
        std::string                      fitformula,
        std::vector<double>              initparams,
        std::vector<std::vector<double>> initparamlims,
        bool use_sumw2error              = true,
        bool use_average_depol           = false,
        bool use_extended_nll            = false,
        bool use_binned_fit              = false,
        std::ostream &out                = std::cout
    ) {

    // Set method name
    std::string method_name = "fitAsym";

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Load helicity variable from workspace
    RooCategory * h = w->cat(helicity.c_str());

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

    // Get depolarization factor means and errors
    std::vector<double> depols;
    std::vector<double> depolerrs;
    RooRealVar * d[(const int)depolvars.size()];
    for (int i=0; i<depolvars.size(); i++) {
        d[i] = w->var(depolvars[i].c_str());
        double mean   = bin_ds->mean(*d[i]);
        double stddev = TMath::Sqrt(bin_ds->moment(*d[i],2.0));
        depols.push_back(mean);
        depolerrs.push_back(stddev);
    }

    // Create asymmetry amplitude parameters
    std::vector<std::string> anames;
    int nparams = initparams.size();
    RooRealVar *a[nparams];
    for (int aa=0; aa<nparams; aa++) {
        std::string aname = Form("a%d",aa);
        anames.push_back(aname);
        a[aa] = new RooRealVar(anames[aa].c_str(),anames[aa].c_str(),initparams[aa],initparamlims[aa][0],initparamlims[aa][1]);
    }

    // Add parameters to argument list in order
    RooArgSet *argset = new RooArgSet();
    for (int ff=0; ff<fitvars.size(); ff++) { // Fit independent variables
        argset->add(*f[ff]);
    }
    for (int aa=0; aa<nparams; aa++) { // Fit asymmetry amplitude parameters
        argset->add(*a[aa]);
    }
    if (!use_average_depol) {
        for (int dd=0; dd<depolvars.size(); dd++) { // Fit depolarization factor variables
            argset->add(*d[dd]);
        }
    }

    // Create pdf positive helicity
    std::string fitformula_pos = Form("1.0+%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf _model_pos("_model_pos", fitformula_pos.c_str(), *argset);

    // Create extended pdf positive helicity
    RooRealVar nsig_pos("nsig_pos", "number of signal events", count/2, 0.0, count);
    RooExtendPdf model_pos("model_pos", "extended signal pdf", _model_pos, nsig_pos);

    // Create pdf negative helicity
    std::string fitformula_neg = Form("1.0-%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf _model_neg("_model_neg", fitformula_neg.c_str(), *argset);

    // Create extended pdf negative helicity
    RooRealVar nsig_neg("nsig_neg", "number of signal events", count/2, 0.0, count);
    RooExtendPdf model_neg("model_neg", "extended signal pdf", _model_neg, nsig_neg);

    // Create simultaneous pdf
    RooSimultaneous * model;
    std::string model_name = Form("model_%s_%s",method_name.c_str(),binid.c_str()); //TODO: Make model names more specific above to avoid naming conflicts...
    if (use_extended_nll) { model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf", {{"plus", &model_pos}, {"minus", &model_neg}}, *h); } //TODO: Set these from helicity states...
    else { model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf", {{"plus", &_model_pos}, {"minus", &_model_neg}}, *h); }

    // Fit the pdf to data
    std::unique_ptr<RooFitResult> r;
    if (use_binned_fit) {

        // Create binned data
        std::unique_ptr<RooDataHist> dh = (std::unique_ptr<RooDataHist>)bin_ds->binnedClone();

        // Fit pdf
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*dh, RooFit::Save(), RooFit::SumW2Error(use_sumw2error), RooFit::PrintLevel(-1));

    } else {

        // Fit pdf
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

    // Define the asymmetry as a function
    RooFormulaVar f_asym("f_asym","Asymmetry function",Form("%.3f*%s",pol,fitformula.c_str()), *argset);//NOTE: NEED TO CORRECT FOR POLARIZATION FACTOR.

    // Loop fit variables and plot fit projections
    for (int idx=0; idx<fitvars.size(); idx++) {

        // Plot projection of fitted distribution in fit variable
        RooPlot *xframe = f[idx]->frame(RooFit::Title(Form("%s Projection, Bin: %s",f[idx]->GetTitle(),bincut.c_str())));
        bin_ds->plotOn(xframe, RooFit::Asymmetry(*h));
        f_asym.plotOn(xframe, RooFit::LineColor(kRed));

        // Draw the frame on the canvas
        std::string c1_x_name = Form("c1_%s__fitvar_%s",binid.c_str(),fitvars[idx].c_str());
        TCanvas *c1_x = new TCanvas(c1_x_name.c_str(), c1_x_name.c_str());
        gPad->SetLeftMargin(0.15);
        xframe->GetYaxis()->SetTitleOffset(1.4);
        xframe->Draw();
        c1_x->Print(Form("%s.pdf",c1_x_name.c_str()));
    }

    // Get fit parameter values and errors
    std::vector<double> params;
    std::vector<double> paramerrs;
    for (int aa=0; aa<nparams; aa++) {
        params.push_back((double)a[aa]->getVal());
        paramerrs.push_back((double)a[aa]->getError());
    }

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " "<<method_name.c_str()<<"():" << std::endl;
    out << " pol        = " << pol << std::endl;
    out << " bincut     = " << bincut.c_str() << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvarmeans[idx] << "±" << binvarerrs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " depolvars  = [" ;
    for (int idx=0; idx<depolvars.size(); idx++) {
        out << depolvars[idx];
        if (idx<depolvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " depols  = [" ;
    for (int idx=0; idx<depols.size(); idx++) {
        out << depols[idx];
        if (idx<depols.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitvars  = [" ;
    for (int idx=0; idx<fitvars.size(); idx++) {
        out << fitvars[idx];
        if (idx<fitvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " initial params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << initparams[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx] << "±" << paramerrs[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    if (use_extended_nll) {
        out << " nsig_pos = " << (double)nsig_pos.getVal() << "±" << (double)nsig_pos.getError() << std::endl;
        out << " nsig_neg = " << (double)nsig_neg.getVal() << "±" << (double)nsig_neg.getError() << std::endl;
    }
    out << "--------------------------------------------------" << std::endl;

    // Fill return array
    std::vector<double> arr; //NOTE: Dimension = 1+2*binvars.size()+2*depolvars.size()+2*nparams
    arr.push_back(count);
    for (int idx=0; idx<binvars.size(); idx++) {
        arr.push_back(binvarmeans[idx]);
        arr.push_back(binvarerrs[idx]);
    }
    for (int idx=0; idx<depolvars.size(); idx++) {
        arr.push_back(depols[idx]);
        arr.push_back(depolerrs[idx]);
    }
    for (int idx=0; idx<nparams; idx++) {
        arr.push_back(params[idx]);
        arr.push_back(paramerrs[idx]);
    }

    return arr;

} // std::vector<double> fitAsym()

/**
* @brief Loop kinematic bins and fit an asymmetry, correcting for background with sideband subtraction or <a href="http://arxiv.org/abs/physics/0402083">sPlots</a>.
*
* Loop bins cuts and fit an asymmetry with the analysis::fitAsym() method.  Optionally apply an invariant mass fit and background correction using the
* sideband subtraction method or the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* Results will be saved in a csv file.
*
* @param scheme_name Name bin scheme and basename of output csv file
* @param frame ROOT RDataframe from which to create RooFit datasets
* @param workspace_name Name of workspace in which to work
* @param workspace_title Title of workspace in which to work
* @param dataset_name Dataset name
* @param dataset_title Dataset title
* @param helicity Name of helicity variable
* @param helicity_states Map of state names to helicity values
* @param bincuts Map of unique bin id ints to bin variable cuts for bin
* @param binvars List of kinematic binning variables names
* @param binvar_titles List of kinematic binning variables titles
* @param binvar_lims List kinematic binning variable minimum and maximum bounds 
* @param binvar_bins List of kinematic binning variables bins
* @param depolvars List of depolarization variables names
* @param depolvar_titles List of depolarization variables titles
* @param depolvar_lims List depolarization variable minimum and maximum bounds 
* @param depolvar_bins List of depolarization variables bins
* @param asymfitvars List of asymmetry fit variables names
* @param asymfitvar_titles List of asymmetry fit variables titles
* @param asymfitvar_lims List asymmetry fit variable minimum and maximum bounds 
* @param asymfitvar_bins List of asymmetry fit variables bins
* @param massfitvars List of invariant mass fit variables names
* @param massfitvar_titles List of invariant mass fit variables titles
* @param massfitvar_lims List invariant mass fit variable minimum and maximum bounds 
* @param massfitvar_bins List of invariant mass fit variables bins
* @param pol Luminosity averaged beam polarization
* @param asymfit_formula The asymmetry formula in ROOT TFormula format
* @param asymfitpar_inits List of initial values for asymmetry fit variables
* @param asymfitpar_initlims List of initial asymmetry fit variables minimum and maximum bounds
* @param use_sumw2error Option to use RooFit::SumW2Error(true) option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_binned_fit Option to use a binned fit to the data
* @param massfit_model_name Full PDF name
* @param massfit_nbins_conv Number of invariant mass bins for convolution of PDFs
* @param massfit_sig_pdf_name Signal PDF name for invariant mass fit
* @param massfit_sg_region_min Invariant mass signal region lower bound
* @param massfit_sg_region_max Invariant mass signal region upper bound
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param use_splot Option to use sPlot method and perform fit with sWeighted dataset
* @param massfit_sgcut Signal region cut for sideband subtraction background correction
* @param massfit_bgcut Signal region cut for sideband subtraction background correction
* @param use_sb_subtraction Option to use sideband subtraction for background correction
* @param use_binned_sb_weights Option to use weights from invariant mass fits binned in the asymmetry fit variable for background correction
* @param asymfitvar_bincuts Map of unique bin id ints to bin variable cuts for asymmetry fit variable bins
* @param out Output stream
*/
void getKinBinnedAsym(
        std::string                      scheme_name,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        std::string                      workspace_name,
        std::string                      workspace_title,

        // parameters passed to data::createDataset()
        std::string                      dataset_name,
        std::string                      dataset_title,
        std::string                      helicity,
        std::map<std::string,int>        helicity_states,
        std::map<int,std::string>        bincuts,
        std::vector<std::string>         binvars,
        std::vector<std::string>         binvar_titles,
        std::vector<std::vector<double>> binvar_lims,
        std::vector<int>                 binvar_bins,
        std::vector<std::string>         depolvars,
        std::vector<std::string>         depolvar_titles,
        std::vector<std::vector<double>> depolvar_lims,
        std::vector<int>                 depolvar_bins,
        std::vector<std::string>         asymfitvars,
        std::vector<std::string>         asymfitvar_titles,
        std::vector<std::vector<double>> asymfitvar_lims,
        std::vector<int>                 asymfitvar_bins,
        std::vector<std::string>         massfitvars,
        std::vector<std::string>         massfitvar_titles,
        std::vector<std::vector<double>> massfitvar_lims,
        std::vector<int>                 massfitvar_bins,

        // parameterss passed to analysis::fitAsym()
        double                           pol,
        std::string                      asymfit_formula,
        std::vector<double>              asymfitpar_inits,
        std::vector<std::vector<double>> asymfitpar_initlims,
        bool                             use_sumw2error,
        bool                             use_average_depol,
        bool                             use_extended_nll,
        bool                             use_binned_fit,

        // parameters passed to data::createDataset() and analysis::applyLambdaMassFit() //TODO: Add init fit parameter value and limits arguments here...assuming you always want a chebychev polynomial background...
        std::string                      massfit_model_name,
        int                              massfit_nbins_conv,
        std::string                      massfit_sig_pdf_name, //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")
        double                           massfit_sg_region_min,
        double                           massfit_sg_region_max,

        // Parameters passed to analysis::applySPlots()
        std::string                      sgYield_name,
        std::string                      bgYield_name,
        bool                             use_splot,

        // Parameters used for sb subtraction
        std::string                      massfit_sgcut,
        std::string                      massfit_bgcut,
        bool                             use_sb_subtraction,
        bool                             use_binned_sb_weights,
        std::map<int,std::string>        asymfitvar_bincuts,

        // Ouput stream
        std::ostream &out                = std::cout
    ) {

    // Check arguments
    if (binvars.size()<1) {std::cerr<<" *** ERROR *** Number of bin variables is <1.  Exiting...\n"; return;}
    if (depolvars.size()!=asymfitvars.size()) {std::cerr<<" *** WARNING *** depolvars.size() does not match the number of parameters injected."<<std::endl;}
    if (use_sb_subtraction && use_splot) {std::cerr<<" *** ERROR *** Cannot simultaneously use sideband subtraction and sPlot.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedAsym ----------------------\n";
    out << "bincuts = { ";
    for (auto it = bincuts.begin(); it != bincuts.end(); it++) {
        out << it->first << " : " << it->second.c_str() << " , ";
    }
    out << " }\n";

    // Filter frames for signal and sideband
    auto frame_sg = frame.Filter(massfit_sgcut.c_str());
    auto frame_sb = frame.Filter(massfit_bgcut.c_str());

    // Open output CSV
    std::string csvpath = Form("%s.csv",scheme_name.c_str());
    std::ofstream csvoutf; csvoutf.open(csvpath.c_str());
    std::ostream &csvout = csvoutf;
    std::string csv_separator = ",";

    // Set CSV column headers
    // COLS: bin_id,count,{binvarmean,binvarerr},{depolvarmean,depolvarerr},{asymfitvar,asymfitvarerr}
    csvout << "bin_id" << csv_separator.c_str();
    csvout << "count" << csv_separator.c_str();
    for (int bb=0; bb<binvars.size(); bb++) {
        csvout << binvars[bb].c_str() << csv_separator.c_str();
        csvout << binvars[bb].c_str() << "_err" << csv_separator.c_str();
    }
    for (int dd=0; dd<depolvars.size(); dd++) {
        csvout << depolvars[dd].c_str() << csv_separator.c_str();
        csvout << depolvars[dd].c_str() << "_err" << csv_separator.c_str();
    }
    for (int aa=0; aa<asymfitpar_inits.size(); aa++) {
        csvout << Form("a%d",aa) << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
        csvout << Form("a%d",aa) << "_err";
        if (aa<asymfitpar_inits.size()-1) csvout << csv_separator.c_str();
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
        RooWorkspace *ws_sg = new RooWorkspace(Form("%s_sg",workspace_name.c_str()),Form("%s_signal",workspace_title.c_str()));
        RooWorkspace *ws_sb = new RooWorkspace(Form("%s_sb",workspace_name.c_str()),Form("%s_sideband",workspace_title.c_str())); //NOTE: Use separate signal and sideband workspaces for dataset, variable, and pdf name uniqueness.

        // Make bin cut on frame
        auto binframe = frame.Filter(bin_cut.c_str());
        auto binframe_sg = frame_sg.Filter(bin_cut.c_str());

        // Create bin dataset
        data::createDataset(
            binframe,
            ws,
            dataset_name,
            dataset_title,
            helicity,
            helicity_states,
            binvars,
            binvar_titles,
            binvar_lims,
            binvar_bins,
            depolvars,
            depolvar_titles,
            depolvar_lims,
            depolvar_bins,
            asymfitvars,
            asymfitvar_titles,
            asymfitvar_lims,
            asymfitvar_bins,
            massfitvars,
            massfitvar_titles,
            massfitvar_lims,
            massfitvar_bins
        );

        // Apply Lambda mass fit to FULL bin frame
        RooAbsData *rooDataSetResult = ws->data(dataset_name.c_str());
        std::vector<double> epss = {0.0, 0.0};
        if (massfit_sig_pdf_name.size()>0) {
            epss = applyLambdaMassFit(
                    ws,
                    massfitvars,
                    rooDataSetResult,
                    sgYield_name,
                    bgYield_name,
                    binframe,
                    massfit_nbins_conv,
                    massfit_model_name,
                    massfit_sig_pdf_name,
                    massfit_sg_region_min,
                    massfit_sg_region_max,
                    1,//use_poly4_bg //NEWTODO: Add argument map for polynomial order
                    scheme_binid
                );
        }

        // Apply sPlot
        std::string fit_dataset_name = dataset_name; // -> Use this for sPlot
        if (use_splot) {
            std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
            std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
            applySPlot(
                ws,
                dataset_name,
                Form("%s%s",sgYield_name.c_str(),scheme_binid.c_str()),//NOTE: applyLambdaMassFit renames these variables to ensure workspace uniqueness
                Form("%s%s",bgYield_name.c_str(),scheme_binid.c_str()),
                Form("%s%s",massfit_model_name.c_str(),scheme_binid.c_str()),
                dataset_sg_name,
                dataset_bg_name
            );
            fit_dataset_name = dataset_sg_name;
        }

        // Weight dataset from binned mass fits
        if (use_binned_sb_weights) {
            std::string rds_weighted_name = (std::string)Form("%s_binned_sb_ws",dataset_name.c_str());
            setWeightsFromLambdaMassFit(
                ws,
                rooDataSetResult,
                massfitvars,
                sgYield_name,
                bgYield_name,
                binframe,
                massfit_nbins_conv,
                massfit_model_name,
                massfit_sig_pdf_name,
                massfit_sg_region_min,
                massfit_sg_region_max,
                1, //use_poly4_bg
                scheme_binid,
                massfit_sgcut,
                massfit_bgcut,
                asymfitvars,
                asymfitvar_bincuts,
                rds_weighted_name,
                "binned_sb_w", //weightvar
                {-999.,999.}, //weightvar_lims,
                0.0 //weights_default
            );
            fit_dataset_name = rds_weighted_name;
        }

        // Create signal region dataset for sideband subtraction
        if (use_sb_subtraction) {
            data::createDataset(
                binframe_sg,
                ws_sg, //NOTE: Use separate sideband workspace for dataset, variable, and pdf name uniqueness.
                dataset_name,
                dataset_title,
                helicity,
                helicity_states,
                binvars,
                binvar_titles,
                binvar_lims,
                binvar_bins,
                depolvars,
                depolvar_titles,
                depolvar_lims,
                depolvar_bins,
                asymfitvars,
                asymfitvar_titles,
                asymfitvar_lims,
                asymfitvar_bins,
                massfitvars,
                massfitvar_titles,
                massfitvar_lims,
                massfitvar_bins
            );
        }

        // Compute signal region bin results
        std::vector<double> asymfit_result = fitAsym(
                                (use_sb_subtraction ? ws_sg : ws),
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                pol,
                                helicity,
                                scheme_binid,
                                bin_cut,
                                binvars,
                                depolvars,
                                asymfitvars,
                                asymfit_formula,
                                asymfitpar_inits,
                                asymfitpar_initlims,
                                use_sumw2error,
                                use_average_depol,
                                use_extended_nll,
                                use_binned_fit,
                                out
                            );

        // Compute sideband region bin results
        std::string sb_scheme_binid = Form("sb_%s",scheme_binid.c_str());
        std::vector<double> asymfit_result_sb;
        if (use_sb_subtraction) {

            // Make bin cut on sideband frame
            auto binframe_sb = frame_sb.Filter(bin_cut.c_str());

            // Create sideband dataset
            data::createDataset(
                binframe_sb,
                ws_sb, //NOTE: Use separate sideband workspace for dataset, variable, and pdf name uniqueness.
                dataset_name,
                dataset_title,
                helicity,
                helicity_states,
                binvars,
                binvar_titles,
                binvar_lims,
                binvar_bins,
                depolvars,
                depolvar_titles,
                depolvar_lims,
                depolvar_bins,
                asymfitvars,
                asymfitvar_titles,
                asymfitvar_lims,
                asymfitvar_bins,
                massfitvars,
                massfitvar_titles,
                massfitvar_lims,
                massfitvar_bins
            );

            // Compute sideband bin results
            asymfit_result_sb = fitAsym(
                                ws_sb,
                                dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                pol,
                                helicity,
                                sb_scheme_binid,
                                bin_cut,
                                binvars,
                                depolvars,
                                asymfitvars,
                                asymfit_formula,
                                asymfitpar_inits,
                                asymfitpar_initlims,
                                use_sumw2error,
                                use_average_depol,
                                use_extended_nll,
                                use_binned_fit,
                                out
                            );
        }

        // Initialize data
        int nbinvars = binvars.size();
        int nparams  = asymfitpar_inits.size();
        double xs[nbinvars];
        double exs[nbinvars];
        int    count;

        double ys[nparams];
        double eys[nparams];
        double ys_sb[nparams];
        double eys_sb[nparams];
        double depols[nparams];
        double edepols[nparams];
        double ys_corrected[nparams];
        double eys_corrected[nparams];

        // Get bin data
        int k = 0;
        count = (int)asymfit_result[k++];
        for (int idx=0; idx<binvars.size(); idx++) {
            xs[idx]     = asymfit_result[k++]; //NOTE: ASSUME THIS IS 1D BINNING FOR NOW
            exs[idx]    = asymfit_result[k++];
        }
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx] = asymfit_result[k++];
            edepols[idx] = asymfit_result[k++];
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx] = asymfit_result[k++];
            eys[idx] = asymfit_result[k++];
        }

        // Apply sideband subtraction to asymmetries assuming that depolarization factors do not vary much
        double epsilon, epsilon_err;
        if (use_sb_subtraction) {
            int k2 = 1 + + depolvars.size();
            epsilon = epss[0];
            epsilon_err = epss[1];
            for (int idx=0; idx<nparams; idx++) {
                ys_sb[idx] = asymfit_result_sb[k2++];
                eys_sb[idx] = asymfit_result_sb[k2++];
                ys[idx]  = (ys[idx] - epsilon * ys_sb[idx]) / (1.0 - epsilon);
                eys[idx] = TMath::Sqrt(eys[idx]*eys[idx] + epsilon * epsilon * eys_sb[idx]*eys_sb[idx]) / (1.0 - epsilon);
            }
        }

        // Divide out depolarization factors
        if (use_average_depol) {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx] = ys[idx] / depols[idx];
                eys_corrected[idx] = eys[idx] / depols[idx];
            }
        } else {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx] = ys[idx];
                eys_corrected[idx] = eys[idx];
            }
        }

        // Output message
        out << "--- Accpetance, background, and depolarization corrected signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            if (use_average_depol) {
                out << " ys["<< idx <<"]             = " << ys[idx] << "\n";
                out << " eys["<< idx <<"]            = " << eys[idx] << "\n";
                out << " depols["<< idx <<"]         = " << depols[idx] << "\n";
                out << " edepols["<< idx <<"]        = " << edepols[idx] << "\n";
            } if (use_sb_subtraction) {
                out << " ys_sb["<< idx <<"]       = " << ys_sb[idx] << "\n";
                out << " eys_sb["<< idx <<"]      = " << eys_sb[idx] << "\n";
            }
            out << " ys_corrected["<< idx <<"]   = " << ys_corrected[idx] << "\n";
            out << " eys_corrected["<< idx <<"]  = " << eys_corrected[idx] << "\n";
        }
        out << "---------------------------\n";

        // Write out a row of data to csv
        // COLS: bin_id,count,{binvarmean,binvarerr},{depolvarmean,depolvarerr},{asymfitvar,asymfitvarerr}
        csvout << bin_id << csv_separator.c_str();
        csvout << count << csv_separator.c_str();
        for (int bb=0; bb<binvars.size(); bb++) {
            csvout << xs[bb] << csv_separator.c_str();
            csvout << exs[bb] << csv_separator.c_str();
        }
        for (int dd=0; dd<depolvars.size(); dd++) {
            csvout << depols[dd] << csv_separator.c_str();
            csvout << edepols[dd] << csv_separator.c_str();
        }
        for (int aa=0; aa<nparams; aa++) {
            csvout << ys_corrected[aa] << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
            csvout << eys_corrected[aa];
            if (aa<nparams-1) csvout << csv_separator.c_str();
            else csvout << std::endl;//NOTE: IMPORTANT!
        }
    }

    csvoutf.close();
    out << " Saved asymmetry fit results to " << csvpath.c_str() << std::endl;

    // Ending message
    out << "------------------- END of getKinBinnedAsym -------------------\n";

} // getKinBinnedAsym()

} // namespace analysis {

} // namespace saga {
