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

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 12/Dec./2024
* @version 0.0.0
* @brief Fit asymmetries using RooFit unbinned extended Maximum Likelihood methods and sideband subtraction 
* or the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a> for background correction.
*/

namespace analysis {

/**
* @brief Apply a \f$\Lambda\f$ baryon mass fit
*
* Apply a \f$\Lambda\f$ baryon mass fit with crystal ball signal and Chebychev polynomial
* background and save model and yield variables to workspace.
*
* @param w RooWorkspace in which to work
* @param massvar Invariant mass variable name
* @param dataset_name Dataset name
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param model_name Full PDF name
*/
void applyLambdaMassFit(
    RooWorkspace *w,
    std::string massvar,
    std::string dataset_name,
    std::string sgYield_name,
    std::string bgYield_name,
    std::string model_name
    ) {

    using namespace RooFit;

    // Get variables from workspace
    RooRealVar *m = w->var(massvar.c_str());

    // Get dataset from workspace
    RooAbsData *rooDataSetResult = w->data(dataset_name.c_str());

    // Get dataset length
    int count = (int)rooDataSetResult->numEntries();

    // Construct signal parameters and function
    RooRealVar mu("mu","#mu",1.1157,0.0,2.0);
    RooRealVar s("s","#sigma",0.008,0.0,1.0);
    RooRealVar a_left("a_left","alpha_left",10.0);
    RooRealVar n_left("n_left","n_left",10.0);
    RooRealVar a("a","#alpha",1.0,0.0,3.0);
    RooRealVar n("n","n",2.0,2.0,10.0);
    RooCrystalBall sig("sig", "sig", *m, mu, s, a_left, n_left, a, n);

    // Consruct background parameters and function
    RooRealVar b1("b1","b_{1}",  0.72,-10.0,10.0);
    RooRealVar b2("b2","b_{2}", -0.17,-10.0,10.0);
    RooRealVar b3("b3","b_{3}",  0.05,-10.0,10.0);
    RooRealVar b4("b4","b_{4}", -0.01,-10.0,10.0);
    RooChebychev bg("bg","bg",*m,RooArgList(b1,b2,b3,b4));
    
    // Combine signal and background functions
    double sgfrac = 0.1;
    double sgYield_init = sgfrac * count;
    double bgYield_init = (1.0-sgfrac) * count;
    RooRealVar sgYield(sgYield_name.c_str(), "fitted yield for signal", sgYield_init, 0., 2.0*count);
    RooRealVar bgYield(bgYield_name.c_str(), "fitted yield for background", bgYield_init, 0., 2.0*count);
    RooAddPdf model(model_name.c_str(), "sig+bg", RooArgList(sig,bg), RooArgList(sgYield,bgYield)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit

    // Fit invariant mass spectrum
    model.fitTo(*rooDataSetResult, Save(), PrintLevel(-1));

    // Plot invariant mass fit
    RooPlot *mframe_1d = m->frame(Title("1D pdf fit mass_ppim."));
    rooDataSetResult->plotOn(mframe_1d);
    model.plotOn(mframe_1d);
    model.plotOn(mframe_1d, Components(sig), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(mframe_1d, Components(bg), LineStyle(kDashed), LineColor(kBlue));
    TCanvas *c_massfit = new TCanvas("c_massfit");
    c_massfit->cd();
    gPad->SetLeftMargin(0.15);
    mframe_1d->GetYaxis()->SetTitleOffset(1.6);
    mframe_1d->Draw();
    c_massfit->SaveAs(Form("%s.pdf",c_massfit->GetName()));

    // Add yield variables to workspace
    w->import(sgYield);
    w->import(bgYield);

    // Add model to workspace
    w->import(model);

    return;
}

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
* @param massvar Invariant mass variable name
* @param dataset_name Dataset name
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param frame ROOT RDataframe from which to create a histogram of invariant mass
* @param mass_nbins_hist Number of bins for the invariant mass histogram
* @param mass_min Minimum bound for invariant mass variable
* @param mass_max Maximum bound for invariant mass variable
* @param mass_nbins_conv Number of invariant mass bins for convolution of PDFs
* @param model_name Full PDF name
* @param sig_pdf_name Signal PDF name
* @param sg_region_min Invariant mass signal region lower bound
* @param sg_region_max Invariant mass signal region upper bound
* @param ws_unique_id Identifier string to ensure PDFs uniqueness in workspace
* @param use_poly4_bg Use a 4th order Chebychev polynomial background instead 2nd order
* @param bin_id Unique bin identifier string
*
* @return List containing background fraction epsilon and its statistical error
*/
std::vector<double> applyLambdaMassFit(
    RooWorkspace *w,
    std::string massvar,
    std::string dataset_name,
    std::string sgYield_name,
    std::string bgYield_name,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    int mass_nbins_hist,
    double mass_min,
    double mass_max,
    int mass_nbins_conv,
    std::string model_name,
    std::string sig_pdf_name,
    double sg_region_min,
    double sg_region_max,
    std::string ws_unique_id,
    int use_poly4_bg,
    std::string bin_id
    ) {

    using namespace RooFit;

    // Get variables from workspace
    RooRealVar *m = w->var(massvar.c_str());

    // Get dataset from workspace
    RooAbsData *rooDataSetResult = w->data(dataset_name.c_str());

    // Get dataset length
    int count = (int)rooDataSetResult->numEntries();

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
    RooGaussian gauss(Form("gauss%s",ws_unique_id.c_str()), "gauss", *m, mg, sg);

    // Construct landau(t,ml,sl) ;
    RooRealVar ml("ml", "mean landau", 1.1157, mass_min, mass_max);
    RooRealVar sl("sl", "sigma landau", 0.005, 0.0, 0.1);
    RooLandau landau(Form("landau%s",ws_unique_id.c_str()), "landau", *m, ml, sl);

    // Construct signal parameters and function
    RooRealVar mu("mu","#mu",1.1157,0.0,2.0);
    RooRealVar s("s","#sigma",0.008,0.0,1.0);
    RooRealVar a_left("a_left","alpha_left",10.0);
    RooRealVar n_left("n_left","n_left",10.0);
    RooRealVar a("a","#alpha",1.0,0.0,3.0);
    RooRealVar n("n","n",2.0,2.0,10.0);
    RooCrystalBall cb(Form("cb%s",ws_unique_id.c_str()), "crystal_ball", *m, mu, s, a_left, n_left, a, n); //NOTE: Include model name for uniqueness within bin.

    // Construct addition component pdfs
    RooRealVar sg_add("sg_add", "sg_add", 0.0001, 0.0, 0.001);
    RooGaussian gauss_add(Form("gauss_add%s",ws_unique_id.c_str()), "gauss_add", *m, mu, sg_add);
    RooCrystalBall cb_add(Form("cb_add%s",ws_unique_id.c_str()), "crystal_ball", *m, mu, s, a_left, n_left, a, n); //NOTE: Include model name for uniqueness within bin.
    double cbFrac_init = 0.1;
    RooRealVar cbFrac(Form("cbFrac%s",ws_unique_id.c_str()), "fitted yield for signal", cbFrac_init, 0., 1.0);
    RooAddPdf cb_gauss(Form("cb_gauss%s",ws_unique_id.c_str()), Form("cb_add%s+gauss_add%s",ws_unique_id.c_str(),ws_unique_id.c_str()), RooArgList(cb_add,gauss_add), RooArgList(cbFrac)); //NOTE: N-1 Coefficients! 

    // Construct convolution component pdfs
    RooRealVar mg_conv("mg_conv", "mg_conv", 0.0);
    RooRealVar sg_conv("sg_conv", "sg_conv", 0.008, 0.0, 0.1);
    RooGaussian gauss_landau_conv(Form("gauss_landau_conv%s",ws_unique_id.c_str()), "gauss_landau_conv", *m, mg_conv, sg_conv);
    RooGaussian gauss_cb_conv(Form("gauss_cb_conv%s",ws_unique_id.c_str()), "gauss_cb_conv", *m, mg_conv, sg_conv);
    RooLandau landau_conv(Form("landau_conv%s",ws_unique_id.c_str()), "landau_conv", *m, ml, sl);
    RooCrystalBall cb_conv(Form("cb_conv%s",ws_unique_id.c_str()), "crystal_ball_conv", *m, mu, s, a_left, n_left, a, n);
    
    // Set #bins to be used for FFT sampling to 10000
    m->setBins(mass_nbins_conv, "cache");
    
    // Construct Convolution PDFs
    RooFFTConvPdf landau_X_gauss(Form("landau_X_gauss%s",ws_unique_id.c_str()), "CB (X) gauss_conv", *m, landau_conv, gauss_landau_conv);
    RooFFTConvPdf cb_X_gauss(Form("cb_X_gauss%s",ws_unique_id.c_str()), "CB (X) gauss_conv", *m, cb_conv, gauss_cb_conv);

    // Import signal functions to workspace
    w->import(gauss);
    w->import(landau);
    w->import(cb);
    w->import(landau_X_gauss);
    w->import(cb_X_gauss);
    w->import(cb_gauss);

    // Pick out signal function based on preference
    std::string sig_pdf_name_unique = Form("%s%s",sig_pdf_name.c_str(),ws_unique_id.c_str());
    RooAbsPdf *sig = w->pdf(sig_pdf_name_unique.c_str());

    // Consruct background parameters and function
    RooRealVar b1("b1","b_{1}",  0.72,-10.0,10.0);
    RooRealVar b2("b2","b_{2}", -0.17,-10.0,10.0);
    RooRealVar b3("b3","b_{3}",  0.05,-10.0,10.0);
    RooRealVar b4("b4","b_{4}", -0.01,-10.0,10.0);
    std::string bg_pdf_name_unique = Form("bg%s",ws_unique_id.c_str());
    RooChebychev bg(bg_pdf_name_unique.c_str(),bg_pdf_name_unique.c_str(),*m,(use_poly4_bg==1 ? RooArgList(b1,b2,b3,b4) : RooArgList(b1,b2)));
    
    // Combine signal and background functions
    double sgfrac = 0.1;
    double sgYield_init = sgfrac * count;
    double bgYield_init = (1.0-sgfrac) * count;
    RooRealVar sgYield(Form("%s%s",sgYield_name.c_str(),ws_unique_id.c_str()), "fitted yield for signal", sgYield_init, 0., 2.0*count);
    RooRealVar bgYield(Form("%s%s",bgYield_name.c_str(),ws_unique_id.c_str()), "fitted yield for background", bgYield_init, 0., 2.0*count);
    RooAddPdf model(model_name.c_str(), Form("%s+%s",sig_pdf_name_unique.c_str(),bg_pdf_name_unique.c_str()), RooArgList(*sig,bg), RooArgList(sgYield,bgYield)); //NOTE: N-1 Coefficients!  Unless you want extended ML Fit

    // Fit invariant mass spectrum
    std::unique_ptr<RooFitResult> fit_result_data{model.fitTo(*rooDataSetResult, Save(), PrintLevel(-1))};

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
    double i_ds = rooDataSetResult->sumEntries(signal_cut.c_str());
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
    rooDataSetResult->plotOn(mframe_1d);
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
    c_massfit->SaveAs(Form("%s_%s_%s_%s.pdf",c_massfit->GetName(),sig_pdf_name.c_str(),ws_unique_id.c_str(),bin_id.c_str()));

    // Add yield variables to workspace
    w->import(sgYield);
    w->import(bgYield);
    w->import(cbFrac);

    // Add model to workspace
    w->import(model);

    // Return background fraction and error
    std::vector<double> result;
    result.push_back(eps_sg_th1_int);
    result.push_back(eps_sg_th1_int_err);
    return result;
}

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
    auto& data1 = static_cast<RooDataSet&>(*rooDataSetResult);
    RooDataSet data_sg_sw{dataset_sg_name.c_str(), data1.GetTitle(), &data1, *data1.get(), nullptr, "sgYield_sw"};
    RooDataSet data_bg_sw{dataset_bg_name.c_str(), data1.GetTitle(), &data1, *data1.get(), nullptr, "bgYield_sw"};

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
* @param outdir Name of output directory
* @param outroot Name of output ROOT file
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param pol Luminosity averaged beam polarization
* @param helicity Name of helicity variable
* @param binvars List of kinematic binning variables
* @param bincut Kinematic variable cut for bin
* @param binid Bin unique id
* @param depolvars List of depolarization variables
* @param depolvarbins List of number of bins in each depolarization variable for binned fit
* @param fitvars List of names for each fit variable
* @param fitvarbins List of number of bins in each fit variable for plotting asymmetry and binned fit
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
        std::string                      outdir,
        TFile                           *outroot,
        RooWorkspace                    *w,
        std::string                      dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
        double                           pol,
        std::string                      helicity,
        std::vector<std::string>         binvars,
        std::string                      bincut,
        int                              binid,
        std::vector<std::string>         depolvars,
        std::vector<int>                 depolvarbins,
        std::vector<std::string>         fitvars,
        std::vector<int>                 fitvarbins,
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

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

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
    std::string model_name = Form("model_%s_%d",method_name.c_str(),binid); //TODO: Make model names more specific above to avoid naming conflicts...
    if (use_extended_nll) { model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf", {{"plus", &model_pos}, {"minus", &model_neg}}, *h); } //TODO: Set these from helicity states...
    else { model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf", {{"plus", &_model_pos}, {"minus", &_model_neg}}, *h); }

    // Fit the pdf to data
    std::unique_ptr<RooFitResult> r;
    if (use_binned_fit) {

        // Set bin numbers for each fit variable
        for (int idx=0; idx<fitvars.size(); idx++) {
            f[idx]->setBins(fitvarbins[idx]);
        }
        if (!use_average_depol) {
            for (int idx=0; idx<depolvars.size(); idx++) {
                d[idx]->setBins(depolvarbins[idx]);
            }
        }

        // Create binned data
        std::unique_ptr<RooDataHist> dh = (std::unique_ptr<RooDataHist>)bin_ds->binnedClone();

        // Fit pdf
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*dh, RooFit::Save(), RooFit::SumW2Error(use_sumw2error), RooFit::PrintLevel(-1));

    } else {

        // Fit pdf
        r = (std::unique_ptr<RooFitResult>)model->fitTo(*bin_ds, RooFit::Save(), RooFit::SumW2Error(use_sumw2error), RooFit::PrintLevel(-1));

        // Set bin numbers for each fit variable
        for (int idx=0; idx<fitvars.size(); idx++) {
            f[idx]->setBins(fitvarbins[idx]);
        }
        if (!use_average_depol) {
            for (int idx=0; idx<depolvars.size(); idx++) {
                d[idx]->setBins(depolvarbins[idx]);
            }
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

    // Define the asymmetry as a function
    RooFormulaVar f_asym("f_asym","Asymmetry function",Form("%.3f*%s",pol,fitformula.c_str()), *argset);//NOTE: NEED TO CORRECT FOR POLARIZATION FACTOR.

    // Loop fit variables and plot fit projections
    for (int idx=0; idx<fitvars.size(); idx++) {

        // Plot projection of fitted distribution in fit variable
        RooPlot *xframe = f[idx]->frame(RooFit::Bins(fitvarbins[idx]), RooFit::Title(Form("%s Projection, Bin: %s",f[idx]->GetTitle(),bincut.c_str())));
        bin_ds->plotOn(xframe, RooFit::Asymmetry(*h));
        f_asym.plotOn(xframe, RooFit::LineColor(kRed));

        // Draw the frame on the canvas
        std::string c1_x_name = Form("c1_%s__fitvar_%s__binid_%d",outdir.c_str(),fitvars[idx].c_str(),binid);
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

    // Go back to parent directory
    outroot->cd("..");

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
* @brief Compute an asymmetry using an unbinned maximum likelihood fit with 1 fit variable.
* 
* Compute the bin count, bin variable values and errors, depolarization variable values and errors,
* and fit parameter values and errors using an unbinned maximum likelihood fit to the asymmetry.
* Note that the given asymmetry formula \f$ A(x, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) \f$ will be converted internally to a PDF of the form
*
* @f[
* \begin{aligned}
* PDF(h, x, &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + h \cdot P \cdot A(x, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...),
* \end{aligned}
* @f]
* and a simultaneous fit will be applied over the data subsets distinguished by the helicity states.  The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the corresponding depolarization factors.
*
* The variable names should be replaced in the fit formula by `x`\f$\rightarrow\f$`x[0]`, `a_0`\f$\rightarrow\f$`x[1]`, etc.
*
* @param outdir Name of output directory
* @param outroot Name of output ROOT file
* @param frame ROOT RDataframe from which to get bin count
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param binvars List of kinematic binning variables
* @param binvarlims List of minimum and maximum bounds for each kinematic binning variable
* @param bincut Kinematic variable cut for bin
* @param bintitle Bin title
* @param pol Luminosity averaged beam polarization
* @param depolvars List of depolarization variables
* @param helicity Name of helicity variable
* @param helicity_states Map of state names to helicity values
* @param fitformula The asymmetry formula
* @param nparams Number of parameters in the asymmetry formula (up to 5)
* @param params List of initial values for asymmetry parameters
* @param fitvarx Name of fit variable
* @param fitvarxtitle Title of fit variable
* @param xbins Number of bins in fit variable for plotting asymmetry
* @param xmin Minimum bound for fit variable
* @param xmax Maximum bound for fit variable
* @param use_sumW2Error Option to use RooFit::SumW2Error(true) option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param out Output stream
*
* @return List of bin count, bin variable means and errors, depolarization variable means and errors, fit parameters and errors
*/
std::vector<double> getKinBinAsymUBML1D(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
    RooWorkspace *w,
    std::string dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
    std::vector<std::string> binvars,
    std::vector<std::vector<double>> binvarlims,
    std::string bincut,
    std::string bintitle,
    double      pol,
    std::vector<std::string>   depolvars,
    std::string  helicity = "heli",
    std::map<std::string,int> helicity_states = {},
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(5),
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    int xbins                  = 16,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    bool use_sumW2Error        = true,
    bool use_average_depol     = false,
    bool use_extended_nll      = false,
    std::ostream &out          = std::cout
    ) {

    // Set plotting title for bin
    std::string title    = Form("%s %s",fitvarxtitle.c_str(),bincut.c_str());

    // Define the helicity variable
    RooCategory h(helicity.c_str(), helicity.c_str());
    for (auto it = helicity_states.begin(); it != helicity_states.end(); it++) {
        h.defineType(it->first.c_str(), it->second);
    }

    // TODO: Load fit avariables from workspace
    RooRealVar *x = w->var(fitvarx.c_str());

    //TODO Load dataset from workspace
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Apply bin cuts
    auto binframe = frame.Filter(bincut.c_str());
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(bincut.c_str());

    // Get count
    auto count = (int)bin_ds->sumEntries();
   
    // Get bin variable means
    std::vector<double> binvar_means;
    std::vector<double> binvar_errs;
    for (int i=0; i<binvars.size(); i++) {
        // auto mean   = (double)*binframe.Mean(binvars[i].c_str());
        // auto stddev = (double)*binframe.StdDev(binvar[i].c_str());
        RooRealVar binvar(binvars[i].c_str(), binvars[i].c_str(), binvarlims[i][0], binvarlims[i][1]);
        double mean   = bin_ds->mean(binvar);
        double stddev = TMath::Sqrt(bin_ds->moment(binvar,2.0));
        binvar_means.push_back(mean);
        binvar_errs.push_back(stddev);
    }

    // Get depolarization factors
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

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create fit parameters
    std::vector<std::string> aNames;
    std::vector<std::vector<double>> parlims; //TODO: Set entries from function argument.
    double parLimit = 0.5;
    RooRealVar *a[nparams];
    for (int aa=0; aa<nparams; aa++) {
        parlims.push_back({-parLimit,parLimit});
        std::string aName = Form("a%d",aa);
        aNames.push_back(aName);
        a[aa] = new RooRealVar(aNames[aa].c_str(),aNames[aa].c_str(),params[aa],parlims[aa][0],parlims[aa][1]);
    }

    // Add parameters to argument list in order
    RooArgSet *arglist = new RooArgSet();
    arglist->add(*x); // Fit variables
    for (int aa=0; aa<nparams; aa++) { // Fit asymmetry parameters
        arglist->add(*a[aa]);
    }
    if (!use_average_depol) {
        for (int dd=0; dd<depolvars.size(); dd++) { // Fit depolarization variables
            arglist->add(*d[dd]);
        }
    }

    // Create 1D PDF positive helicity
    std::string fitformula_pos = Form("1.0+%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf _gen_pos("_gen_pos", fitformula_pos.c_str(), *arglist);

    // Create extended pdf positive helicity
    RooRealVar nsig_pos("nsig_pos", "number of signal events", count/2, 0.0, count);
    RooExtendPdf gen_pos("gen_pos", "extended signal pdf", _gen_pos, nsig_pos);

    // Create 1D PDF negative helicity
    std::string fitformula_neg = Form("1.0-%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf _gen_neg("_gen_neg", fitformula_neg.c_str(), *arglist);

    // Create extended pdf negative helicity
    RooRealVar nsig_neg("nsig_neg", "number of signal events", count/2, 0.0, count);
    RooExtendPdf gen_neg("gen_neg", "extended signal pdf", _gen_neg, nsig_neg);

    // Create simultaneous pdf
    RooSimultaneous * gen;
    if (use_extended_nll) { gen = new RooSimultaneous("gen", "simultaneous pdf", {{"plus", &gen_pos}, {"minus", &gen_neg}}, h); } //TODO: Set these from helicity states...
    else { gen = new RooSimultaneous("gen", "simultaneous pdf", {{"plus", &_gen_pos}, {"minus", &_gen_neg}}, h); }

    // Fit the pdf to data
    std::unique_ptr<RooFitResult> r{gen->fitTo(*bin_ds, RooFit::Save(), RooFit::SumW2Error(use_sumW2Error), RooFit::PrintLevel(-1))};

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
    RooFormulaVar f_asym("f_asym","Asymmetry function",Form("%.3f*%s",pol,fitformula.c_str()), *arglist);//NOTE: NEED TO CORRECT FOR POLARIZATION FACTOR.

    // Plot projection of fitted distribution in x.
    RooPlot *xframe = x->frame(RooFit::Bins(xbins), RooFit::Title(Form("%s Projection, Bin: %s",fitvarxtitle.c_str(),bincut.c_str())));
    bin_ds->plotOn(xframe, RooFit::Asymmetry(h));
    f_asym.plotOn(xframe, RooFit::LineColor(kRed));

    // Draw the frame on the canvas
    std::string c1_x_name = Form("c1_%s__fitvarx_%s",outdir.c_str(),fitvarx.c_str());
    TCanvas *c1_x = new TCanvas(c1_x_name.c_str(), c1_x_name.c_str());
    gPad->SetLeftMargin(0.15);
    xframe->GetYaxis()->SetTitleOffset(1.4);
    xframe->Draw();
    c1_x->Print(Form("%s.pdf",c1_x_name.c_str()));

    // Get fit parameters
    std::vector<double> pars;
    std::vector<double> Epars;
    for (int aa=0; aa<nparams; aa++) {
        pars.push_back((double)a[aa]->getVal());
        Epars.push_back((double)a[aa]->getError());
    }

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinAsymUBML1D():" << std::endl;
    out << " bincut     = " << bincut.c_str() << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvar_means[idx] << "±" << binvar_errs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " pol        = " << pol << std::endl;
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
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " initial params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " nsig_pos = " << (double)nsig_pos.getVal() << "±" << (double)nsig_pos.getError() << std::endl;
    out << " nsig_neg = " << (double)nsig_neg.getVal() << "±" << (double)nsig_neg.getError() << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    std::vector<double> arr; //NOTE: Dimension = 1+2*binvars.size()+2*depolvars.size()+2*nparams
    arr.push_back(count);
    for (int idx=0; idx<binvars.size(); idx++) {
        arr.push_back(binvar_means[idx]);
        arr.push_back(binvar_errs[idx]);
    }
    for (int idx=0; idx<depolvars.size(); idx++) {
        arr.push_back(depols[idx]);
        arr.push_back(depolerrs[idx]);
    }
    for (int idx=0; idx<nparams; idx++) {
        arr.push_back(pars[idx]);
        arr.push_back(Epars[idx]);
    }

    return arr;

} // std::vector<double> getKinBinAsymUBML1D()

/**
* @brief Loop kinematic bins and fit a 1D asymmetry.
*
* Loop bins and compute an asymmetry using an unbinned maximum likelihood fit.
* Optionally apply an invariant mass fit and background correction using the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* Note that the given asymmetry formula \f$ A(x, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) \f$ will be converted internally to a PDF of the form
*
* @f[
* \begin{aligned}
* PDF(h, x, &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + h \cdot P \cdot A(x, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...),
* \end{aligned}
* @f]
* and a simultaneous fit will be applied over the data subsets distinguished by the helicity states.  The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the corresponding depolarization factors.
*
* The variable names should be replaced in the fit formula by `x`\f$\rightarrow\f$`x[0]`, `a_0`\f$\rightarrow\f$`x[1]`, etc.
*
*
* @param outdir Name of output directory
* @param outroot Name of output ROOT file
* @param frame ROOT RDataFrame
* @param method Asymmetry fit method
* @param binvar Kinematic binning variable
* @param nbins Number of kinematic variable bins
* @param bins List of bin limits, which must have dimension nbins+1
* @param pol Luminosity averaged beam polarization
* @param depolvars List of depolarization variables (up to 5)
* @param workspace_name Name of workspace in which to work
* @param workspace_title Title of workspace in which to work
* @param dataset_name Dataset name
* @param dataset_title Dataset title
* @param helicity Name of helicity variable
* @param fitvarx Name of fit variable
* @param xmin Minimum bound for fit variable
* @param xmax Maximum bound for fit variable
* @param massvar Invariant mass variable name
* @param mmin Minimum bound for invariant mass variable
* @param mmax Maximum bound for invariant mass variable
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param model_name Full PDF name
* @param fitformula The asymmetry formula
* @param nparams Number of parameters in the asymmetry formula (up to 5)
* @param params List of initial values for asymmetry parameters
* @param fitvarxtitle Title of fit variable
* @param xbins Number of bins in fit variable for plotting asymmetry
* @param use_sumW2Error Option to use RooFit::SumW2Error(true) option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_splot Option to use sPlot method and perform fit with sWeighted dataset
* @param graph_title Title of kinematically binned graph of asymmetry parameters
* @param marker_color Graph ROOT marker color
* @param marker_style Graph ROOT marker style
* @param out Output stream
*/
void getKinBinnedAsymUBML1D(
        std::string outdir,
        TFile      *outroot,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        std::string     method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
        std::string binvar, // Variable name to bin in
        int         nbins, // Number of bins
        double      *bins, // Bin limits (length=nbins+1)
        double      pol,
        std::vector<std::string> depolvars,
        std::string workspace_name  = "w",
        std::string workspace_title = "workspace",
        std::string dataset_name    = "dataset",
        std::string dataset_title   = "dataset",
        std::string helicity   = "heli",
        std::string fitvarx         = "x",
        double      xmin            = 0.0,
        double      xmax            = 1.0,
        std::string massvar         = "mass_ppim",
        double      mmin            = 1.08,
        double      mmax            = 1.24,
        std::string sgYield_name    = "sgYield",
        std::string bgYield_name    = "bgYield",
        std::string model_name      = "model",
        std::string fitformula      = "[0]*sin(x)+[1]*sin(2*x)",
        int         nparams         = 2,
        std::vector<double> params  = std::vector<double>(5),
        std::string fitvarxtitle    = "#phi_{h p#pi^{-}}",
        int         xbins           = 16,
        bool        use_sumW2Error  = true,
        bool        use_average_depol = false,
        bool        use_extended_nll = false,
        bool        use_splot       = true,
        std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
        int         marker_color    = 4, // 4 is blue
        int         marker_style    = 20, // 20 is circle
        std::ostream &out           = std::cout
    ) {

    // Check arguments
    if (method != "BSA1D") {std::cerr<<" *** ERROR *** Method must be BSA1D.  Exiting...\n"; return;}
    if (nbins<1) {std::cerr<<" *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}
    if (depolvars.size()!=nparams) {std::cerr<<" *** WARNING *** depolvars.size() does not match the number of parameters injected."<<std::endl;}

    // Starting message
    out << "----------------------- getKinBinnedAsymUBML1D ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double depols[nparams][nbins];
    double edepols[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];

    // Create workspace
    RooWorkspace *w = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());

    // Create bin var vector and outer bin lims vector
    std::vector binvars = {binvar};
    std::vector<std::vector<double>> binvarlims_outer = {{bins[0],bins[nbins]}};
    std::vector<std::vector<double>> depolvarlims;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvarlims.push_back({0.0,1.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.
    }
    std::map<std::string,int> helicity_states;
    helicity_states["plus"]  =  1;
    helicity_states["zero"]  =  0;
    helicity_states["minus"] = -1;

    // Create dataset
    data::createDataset(
        frame,
        w,
        dataset_name,
        dataset_title,
        helicity,
        helicity_states,
        {fitvarx,massvar},
        {{xmin,xmax},{mmin,mmax}},
        binvars,
        binvarlims_outer,
        depolvars,
        depolvarlims
    );

    // Apply Lambda mass fit
    if (use_splot) {
        applyLambdaMassFit(
            w,
            massvar,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name
        );
    }

    // Apply sPlot
    std::string fit_dataset_name = dataset_name;
    if (use_splot) {
        std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
        std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
        applySPlot(
            w,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name,
            dataset_sg_name,
            dataset_bg_name
        );
        fit_dataset_name = dataset_sg_name;
    }

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {

        // Get bin limits
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bincut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto binframe = frame.Filter(bincut.c_str());

        // Compute bin results
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        std::vector<double> bin_data = getKinBinAsymUBML1D(
                                binoutdir,
                                outroot,
                                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                w,
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                helicity,
                                helicity_states,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xbins,
                                xmin,
                                xmax,
                                use_sumW2Error,
                                use_average_depol,
                                use_extended_nll,
                                out
                            );

        // Get bin data
        int k = 0;
        counts[binidx] = (int)bin_data[k++];
        xs[binidx]     = bin_data[k++]; //NOTE: ASSUME THIS IS 1D BINNING FOR NOW
        exs[binidx]    = bin_data[k++];
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx][binidx] = bin_data[k++];
            edepols[idx][binidx] = bin_data[k++];
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = bin_data[k++];
            eys[idx][binidx] = bin_data[k++];
        }

        // Divide out depolarization factors
        if (use_average_depol) {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx] / depols[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx] / depols[idx][binidx];
            }
        } else {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        }

        // Output message
        out << "--- Accpetance, background, and depolarization corrected signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            if (use_average_depol) {
                out << " ys["<< idx <<"]["<<binidx<<"]             = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]            = " << eys[idx][binidx] << "\n";
                out << " depols["<< idx <<"]["<<binidx<<"]         = " << depols[idx][binidx] << "\n";
                out << " edepols["<< idx <<"]["<<binidx<<"]        = " << edepols[idx][binidx] << "\n";
            }
            out << " ys_corrected["<< idx <<"]["<<binidx<<"]   = " << ys_corrected[idx][binidx] << "\n";
            out << " eys_corrected["<< idx <<"]["<<binidx<<"]  = " << eys_corrected[idx][binidx] << "\n";
        }
        out << "---------------------------\n";
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

        // Plot results graph
        TCanvas *c1 = new TCanvas();
        c1->SetBottomMargin(0.125);
        c1->cd(0);

        // Stylistic choices that aren't really necessary
        gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
        gr->SetMarkerSize(1.25);
        gr->GetXaxis()->SetTitleSize(0.05);
        gr->GetXaxis()->SetTitleOffset(0.9);
        gr->GetYaxis()->SetTitleSize(0.05);
        gr->GetYaxis()->SetTitleOffset(0.9);

        // More necessary stylistic choices
        gr->SetTitle(graph_title.c_str());
        gr->SetMarkerColor(marker_color); // 4  blue
        gr->SetMarkerStyle(marker_style); // 20 circle
        gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
        gr->GetXaxis()->SetTitle(binvar.c_str());
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
        gr->Draw("PA");

        // Add CLAS12 Preliminary watermark
        TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary"); //TODO: Add option for watermark
        lt->SetTextAngle(45);
        lt->SetTextColor(18);
        lt->SetTextSize(0.1);
        lt->SetNDC();
        lt->Draw();

        // Add zero line
        TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
        f2->SetLineColor(1); // 1 black
        f2->SetLineWidth(1); // 1 thinnest
        f2->Draw("SAME");

        // Set outname and save
        TString fname;
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",method.c_str(),fitvarx.c_str(),binvar.c_str(),bins[0],bins[nbins],idx); //TODO: Replace all bin naming schemes with bin id after introducing arbitrary dimensional binning
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedAsymUBML1D -------------------\n";

} // getKinBinnedAsymUBML1D()

/**
* @brief Loop kinematic bins and fit a 1D asymmetry, correcting for background with sideband subtraction or <a href="http://arxiv.org/abs/physics/0402083">sPlots</a>.
*
* Loop bins and compute an asymmetry using an unbinned maximum likelihood fit.
* Optionally apply an invariant mass fit and background correction using the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* or the sideband subtraction method.  Note that the asymmetry fit formula \f$ A(x, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) \f$ will be converted internally to a PDF of the form
*
* @f[
* \begin{aligned}
* PDF(h, x, &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + h \cdot P \cdot A(x, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...),
* \end{aligned}
* @f]
* and a simultaneous fit will be applied over the data subsets distinguished by the helicity states.  The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the corresponding depolarization factors.
*
* The variable names should be replaced in the fit formula by `x`\f$\rightarrow\f$`x[0]`, `a_0`\f$\rightarrow\f$`x[1]`, etc.
*
* @param outdir Name of output directory
* @param outroot Name of output ROOT file
* @param frame ROOT RDataframe from which to create RooFit datasets
* @param method Asymmetry fit method
* @param binvar Kinematic binning variable
* @param nbins Number of kinematic variable bins
* @param bins List of bin limits, which must have dimension nbins+1
* @param pol Luminosity averaged beam polarization
* @param depolvars List of depolarization variables (up to 5)
* @param workspace_name Name of workspace in which to work
* @param workspace_title Title of workspace in which to work
* @param dataset_name Dataset name
* @param dataset_title Dataset title
* @param helicity Name of helicity variable
* @param fitvarx Name of fit variable
* @param xmin Minimum bound for fit variable
* @param xmax Maximum bound for fit variable
* @param massvar Invariant mass variable name
* @param mass_min Minimum bound for invariant mass variable
* @param mass_max Maximum bound for invariant mass variable
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param model_name Full PDF name
* @param fitformula The asymmetry formula
* @param nparams Number of parameters in the asymmetry formula (up to 5)
* @param params List of initial values for asymmetry parameters
* @param fitvarxtitle Title of fit variable
* @param xbins Number of bins in fit variable for plotting asymmetry
* @param use_sumW2Error Option to use RooFit::SumW2Error() option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_splot Option to use sPlot method and perform fit with sWeighted dataset
* @param graph_title Title of kinematically binned graph of asymmetry parameters
* @param marker_color Graph ROOT marker color
* @param marker_style Graph ROOT style color
* @param out Output stream
* @param sgcut Signal region cut for sideband subtraction background correction
* @param bgcut Signal region cut for sideband subtraction background correction
* @param mass_nbins_hist Number of bins for the invariant mass histogram
* @param mass_nbins_conv Number of invariant mass bins for convolution of PDFs
* @param sig_pdf_name Signal PDF name
* @param sg_region_min Invariant mass signal region lower bound
* @param sg_region_max Invariant mass signal region upper bound
* @param use_sb_subtraction Option to use sideband subtraction for background correction
*/
void getKinBinnedAsym1D(
        std::string outdir,
        TFile      *outroot,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        std::string     method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
        std::string binvar, // Variable name to bin in
        int         nbins, // Number of bins
        double      *bins, // Bin limits (length=nbins+1)
        double      pol,
        std::vector<std::string> depolvars,
        std::string workspace_name  = "w",
        std::string workspace_title = "workspace",
        std::string dataset_name    = "dataset",
        std::string dataset_title   = "dataset",
        std::string helicity   = "heli",
        std::string fitvarx         = "x",
        double      xmin            = 0.0,
        double      xmax            = 1.0,
        std::string massvar         = "mass_ppim",
        double      mass_min        = 1.08,
        double      mass_max        = 1.24,
        std::string sgYield_name    = "sgYield",
        std::string bgYield_name    = "bgYield",
        std::string model_name      = "model",
        std::string fitformula      = "[0]*sin(x)+[1]*sin(2*x)",
        int         nparams         = 2,
        std::vector<double> params  = std::vector<double>(5),
        std::string fitvarxtitle    = "#phi_{h p#pi^{-}}",
        int         xbins           = 16,
        bool        use_sumW2Error  = true,
        bool use_average_depol      = false,
        bool use_extended_nll       = false,
        bool        use_splot       = true,
        std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
        int         marker_color    = 4, // 4 is blue
        int         marker_style    = 20, // 20 is circle
        std::ostream &out           = std::cout,
        std::string sgcut           = "Q2>1",
        std::string bgcut           = "Q2>1",
        int mass_nbins_hist         = 100,
        int mass_nbins_conv         = 1000,
        std::string sig_pdf_name    = "cb", //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")
        double sg_region_min        = 1.11,
        double sg_region_max        = 1.13,
        bool   use_sb_subtraction   = false
    ) {

    // Check arguments
    if (method != "BSA1D") {std::cerr<<" *** ERROR *** Method must be BSA1D.  Exiting...\n"; return;}
    if (nbins<1) {std::cerr<<" *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}
    if (depolvars.size()!=nparams) {std::cerr<<" *** WARNING *** depolvars.size() does not match the number of parameters injected."<<std::endl;}
    if (use_sb_subtraction && use_splot) {std::cerr<<" *** ERROR *** Cannot simultaneously use sideband subtraction and sPlot.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedAsym1D ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_sb[nparams][nbins];
    double eys_sb[nparams][nbins];
    double depols[nparams][nbins];
    double edepols[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];

    // Create bin var vector and outer bin lims vector
    std::vector binvars = {binvar};
    std::vector<std::vector<double>> binvarlims_outer = {{bins[0],bins[nbins]}};
    std::vector<std::vector<double>> depolvarlims;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvarlims.push_back({0.0,1.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.
    }

    // Filter frames for signal and sideband
    auto frame_sg = frame.Filter(sgcut.c_str());
    auto frame_sb = frame.Filter(bgcut.c_str());

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {

        // Get bin limits
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Create workspace
        RooWorkspace *ws    = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());
        RooWorkspace *ws_sg = new RooWorkspace(Form("%s_sg",workspace_name.c_str()),Form("%s_signal",workspace_title.c_str()));
        RooWorkspace *ws_sb = new RooWorkspace(Form("%s_sb",workspace_name.c_str()),Form("%s_sideband",workspace_title.c_str())); //NOTE: Use separate signal and sideband workspaces for dataset, variable, and pdf name uniqueness.

        // Make bin cut on frame
        std::string  bincut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto binframe = frame.Filter(bincut.c_str());
        auto binframe_sg = frame_sg.Filter(bincut.c_str());

        std::map<std::string,int> helicity_states;
        helicity_states["plus"]  =  1;
        helicity_states["zero"]  =  0;
        helicity_states["minus"] = -1;

        // Create bin dataset
        data::createDataset(
            binframe,
            ws,
            dataset_name,
            dataset_title,
            helicity,
            helicity_states,
            {fitvarx,massvar},
            {{xmin,xmax},{mass_min,mass_max}},
            binvars,
            binvarlims_outer,
            depolvars,
            depolvarlims
        );

        // Apply Lambda mass fit to FULL bin frame
        std::string bin_id = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
        std::vector<double> epss = applyLambdaMassFit(
                ws,
                massvar,
                dataset_name,
                sgYield_name,
                bgYield_name,
                binframe,
                mass_nbins_hist,
                mass_min,
                mass_max,
                mass_nbins_conv,
                model_name,
                sig_pdf_name,
                sg_region_min,
                sg_region_max,
                "",//ws_unique_id->This changes pdf,yieldvar names, but NOT (bin,depol,mass,fit)vars,pdf parameters which are saved internally.
                1,//use_poly4_bg
                bin_id
            );

        // Apply sPlot
        std::string fit_dataset_name = dataset_name; // -> Use this for sPlot
        if (use_splot) {
            std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
            std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
            applySPlot(
                ws,
                dataset_name,
                sgYield_name,
                bgYield_name,
                model_name,
                dataset_sg_name,
                dataset_bg_name
            );
            fit_dataset_name = dataset_sg_name;
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
                {fitvarx,massvar},
                {{xmin,xmax},{mass_min,mass_max}},
                binvars,
                binvarlims_outer,
                depolvars,
                depolvarlims
            );
        }

        // Compute signal region bin results
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        std::vector<double> bin_data = getKinBinAsymUBML1D(
                                binoutdir,
                                outroot,
                                (use_sb_subtraction ? binframe_sg : binframe), //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                (use_sb_subtraction ? ws_sg : ws),
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                helicity,
                                helicity_states,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xbins,
                                xmin,
                                xmax,
                                use_sumW2Error,
                                use_average_depol,
                                use_extended_nll,
                                out
                            );

        // Compute sideband region bin results
        std::string  sbbinoutdir = Form("method_%s_sbbin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        std::vector<double> bin_data_sb;
        if (use_sb_subtraction) {

            // Make bin cut on sideband frame
            auto binframe_sb = frame_sb.Filter(bincut.c_str());

            // Create sideband dataset
            data::createDataset(
                binframe_sb,
                ws_sb, //NOTE: Use separate sideband workspace for dataset, variable, and pdf name uniqueness.
                dataset_name,
                dataset_title,
                helicity,
                helicity_states,
                {fitvarx,massvar},
                {{xmin,xmax},{mass_min,mass_max}},
                binvars,
                binvarlims_outer,
                depolvars,
                depolvarlims
            );

            // Compute sideband bin results
            bin_data_sb = getKinBinAsymUBML1D(
                                binoutdir,
                                outroot,
                                binframe_sb, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                ws_sb,
                                dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                helicity,
                                helicity_states,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xbins,
                                xmin,
                                xmax,
                                use_sumW2Error,
                                use_average_depol,
                                use_extended_nll,
                                out
                            );
        } 

        // Get bin data
        int k = 0;
        counts[binidx] = (int)bin_data[k++];
        xs[binidx]     = bin_data[k++]; //NOTE: ASSUME THIS IS 1D BINNING FOR NOW
        exs[binidx]    = bin_data[k++];
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx][binidx] = bin_data[k++];
            edepols[idx][binidx] = bin_data[k++];
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = bin_data[k++];
            eys[idx][binidx] = bin_data[k++];
        }

        // Apply sideband subtraction to asymmetries assuming that depolarization factors do not vary much
        double epsilon, epsilon_err;
        if (use_sb_subtraction) {
            int k2 = 3 + depolvars.size();
            epsilon = epss[0];
            epsilon_err = epss[1];
            for (int idx=0; idx<nparams; idx++) {
                ys_sb[idx][binidx] = bin_data_sb[k2++];
                eys_sb[idx][binidx] = bin_data_sb[k2++];
                ys[idx][binidx]  = (ys[idx][binidx] - epsilon * ys_sb[idx][binidx]) / (1.0 - epsilon);
                eys[idx][binidx] = TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx] + epsilon * epsilon * eys_sb[idx][binidx]*eys_sb[idx][binidx]) / (1.0 - epsilon) ;
            }
        }

        // Divide out depolarization factors
        if (use_average_depol) {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx] / depols[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx] / depols[idx][binidx];
            }
        } else {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        }

        // Output message
        out << "--- Accpetance, background, and depolarization corrected signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            if (use_average_depol) {
                out << " ys["<< idx <<"]["<<binidx<<"]             = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]            = " << eys[idx][binidx] << "\n";
                out << " depols["<< idx <<"]["<<binidx<<"]         = " << depols[idx][binidx] << "\n";
                out << " edepols["<< idx <<"]["<<binidx<<"]        = " << edepols[idx][binidx] << "\n";
            } if (use_sb_subtraction) {
                out << " ys_sb["<< idx <<"]["<<binidx<<"]       = " << ys_sb[idx][binidx] << "\n";
                out << " eys_sb["<< idx <<"]["<<binidx<<"]      = " << eys_sb[idx][binidx] << "\n";
            }
            out << " ys_corrected["<< idx <<"]["<<binidx<<"]   = " << ys_corrected[idx][binidx] << "\n";
            out << " eys_corrected["<< idx <<"]["<<binidx<<"]  = " << eys_corrected[idx][binidx] << "\n";
        }
        out << "---------------------------\n";
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

        // Plot results graph
        TCanvas *c1 = new TCanvas();
        c1->SetBottomMargin(0.125);
        c1->cd(0);

        // Stylistic choices that aren't really necessary
        gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
        gr->SetMarkerSize(1.25);
        gr->GetXaxis()->SetTitleSize(0.05);
        gr->GetXaxis()->SetTitleOffset(0.9);
        gr->GetYaxis()->SetTitleSize(0.05);
        gr->GetYaxis()->SetTitleOffset(0.9);

        // More necessary stylistic choices
        gr->SetTitle(graph_title.c_str());
        gr->SetMarkerColor(marker_color); // 4  blue
        gr->SetMarkerStyle(marker_style); // 20 circle
        gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
        gr->GetXaxis()->SetTitle(binvar.c_str());
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
        gr->Draw("PA");

        // Add CLAS12 Preliminary watermark
        TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
        lt->SetTextAngle(45);
        lt->SetTextColor(18);
        lt->SetTextSize(0.1);
        lt->SetNDC();
        lt->Draw();

        // Add zero line
        TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
        f2->SetLineColor(1); // 1 black
        f2->SetLineWidth(1); // 1 thinnest
        f2->Draw("SAME");

        // Set outname and save
        TString fname;
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",method.c_str(),fitvarx.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedAsym1D -------------------\n";

} // getKinBinnedAsym1D()

/**
* @brief Compute an asymmetry using an unbinned maximum likelihood fit with 2 fit variables.
*
* Compute the bin count, bin variable values and errors, depolarization variable values and errors,
* and fit parameter values and errors using an unbinned maximum likelihood fit and an asymmetry fit.
* Note that the given asymmetry formula \f$ A(x, y, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) \f$ will be converted internally to a PDF of the form
*
* @f[
* \begin{aligned}
* PDF(h, x, y, &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + h \cdot P \cdot A(x, y, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...),
* \end{aligned}
* @f]
* and a simultaneous fit will be applied over the data subsets distinguished by the helicity states.  The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the corresponding depolarization factors.
*
* The variable names should be replaced in the fit formula by `x`\f$\rightarrow\f$`x[0]`, `y`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[2]`, etc.
*
* @param outdir Name of output directory
* @param outroot Name of output ROOT file
* @param frame ROOT RDataframe from which to create a RooDataSet
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param binvars List of kinematic binning variables
* @param binvarlims List of minimum and maximum bounds for each kinematic binning variable
* @param bincut Kinematic variable cut for bin
* @param bintitle Bin title
* @param pol Luminosity averaged beam polarization
* @param depolvars List of depolarization variables (up to 5)
* @param helicity Name of helicity variable
* @param helicity_states Map of state names to helicity values
* @param fitformula The asymmetry formula
* @param nparams Number of parameters in the asymmetry formula (up to 5)
* @param params List of initial values for asymmetry parameters
* @param fitvarx Name of fit variable 1
* @param fitvarxtitle Title of fit variable 1
* @param xbins Number of bins in fit variable 1 for plotting asymmetry
* @param xmin Minimum bound for fit variable 1
* @param xmax Maximum bound for fit variable 1
* @param fitvary Name of fit variable 2
* @param fitvarytitle Title of fit variable 2
* @param ybins Number of bins in fit variable 2 for plotting asymmetry
* @param ymin Minimum bound for fit variable 2
* @param ymax Maximum bound for fit variable 2
* @param use_sumW2Error Option to use RooFit::SumW2Error(true) option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param out Output stream
*
* @return List of bin count, bin variable means and errors, depolarization variable means and errors, fit parameters and errors
*/
std::vector<double> getKinBinAsymUBML2D(
    std::string  outdir,
    TFile      * outroot,
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
    RooWorkspace *w,
    std::string dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
    std::vector<std::string> binvars,
    std::vector<std::vector<double>> binvarlims,
    std::string bincut,
    std::string bintitle,
    double      pol,
    std::vector<std::string>   depolvars,
    std::string  helicity = "heli",
    std::map<std::string,int> helicity_states = {},
    std::string  fitformula    = "[0]*sin(x)+[1]*sin(2*x)",
    int          nparams       = 2,
    std::vector<double> params = std::vector<double>(5),
    std::string  fitvarx       = "phi_h",
    std::string  fitvarxtitle  = "#phi_{h p#pi^{-}}",
    int xbins                  = 16,
    double xmin                = 0.0,
    double xmax                = 2*TMath::Pi(),
    std::string  fitvary       = "phi_h",
    std::string  fitvarytitle  = "#phi_{h p#pi^{-}}",
    int ybins                  = 16,
    double ymin                = 0.0,
    double ymax                = 2*TMath::Pi(),
    bool use_sumW2Error        = true,
    bool use_average_depol     = false,
    bool use_extended_nll      = false,
    std::ostream &out          = std::cout
    ) {

    // Set plotting title for bin
    std::string title    = Form("%s %s",fitvarxtitle.c_str(),bincut.c_str());

    // Define the helicity variable
    RooCategory h(helicity.c_str(), helicity.c_str());
    for (auto it = helicity_states.begin(); it != helicity_states.end(); it++) {
        h.defineType(it->first.c_str(), it->second);
    }

    // TODO: Load fit avariables from workspace
    RooRealVar *x = w->var(fitvarx.c_str());
    RooRealVar *y = w->var(fitvary.c_str());

    //TODO Load dataset from workspace
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Apply bin cuts
    auto binframe = frame.Filter(bincut.c_str());
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(bincut.c_str());

    // Get count
    auto count = (int)bin_ds->sumEntries();

    // Get bin variable means
    std::vector<double> binvar_means;
    std::vector<double> binvar_errs;
    for (int i=0; i<binvars.size(); i++) {
        // auto mean   = (double)*binframe.Mean(binvars[i].c_str());
        // auto stddev = (double)*binframe.StdDev(binvar[i].c_str());
        RooRealVar binvar(binvars[i].c_str(), binvars[i].c_str(), binvarlims[i][0], binvarlims[i][1]);
        double mean   = bin_ds->mean(binvar);
        double stddev = TMath::Sqrt(bin_ds->moment(binvar,2.0));
        binvar_means.push_back(mean);
        binvar_errs.push_back(stddev);
    }

    // Get depolarization factors
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

    // Make subdirectory
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Switch off histogram stats
    gStyle->SetOptStat(0);

    // Create fit parameters
    std::vector<std::string> aNames;
    std::vector<std::vector<double>> parlims; //TODO: Set entries from function argument.
    double parLimit = 0.5;
    RooRealVar *a[nparams];
    for (int aa=0; aa<nparams; aa++) {
        parlims.push_back({-parLimit,parLimit});
        std::string aName = Form("a%d",aa);
        aNames.push_back(aName);
        a[aa] = new RooRealVar(aNames[aa].c_str(),aNames[aa].c_str(),params[aa],parlims[aa][0],parlims[aa][1]);
    }

    // Add parameters to argument list in order
    RooArgSet *arglist = new RooArgSet();
    arglist->add(*x); // Fit variables
    arglist->add(*y);
    for (int aa=0; aa<nparams; aa++) { // Fit asymmetry parameters
        arglist->add(*a[aa]);
    }
    if (!use_average_depol) {
        for (int dd=0; dd<depolvars.size(); dd++) { // Fit depolarization variables
            arglist->add(*d[dd]);
        }
    }

    // Create 1D PDF positive helicity
    std::string fitformula_pos = Form("1.0+%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf _gen_pos("_gen_pos", fitformula_pos.c_str(), *arglist);

    // Create extended pdf positive helicity
    RooRealVar nsig_pos("nsig_pos", "number of signal events", count/2, 0.0, count);
    RooExtendPdf gen_pos("gen_pos", "extended signal pdf", _gen_pos, nsig_pos);

    // Create 1D PDF negative helicity
    std::string fitformula_neg = Form("1.0-%.3f*%s",pol,fitformula.c_str());
    RooGenericPdf _gen_neg("_gen_neg", fitformula_neg.c_str(), *arglist);

    // Create extended pdf negative helicity
    RooRealVar nsig_neg("nsig_neg", "number of signal events", count/2, 0.0, count);
    RooExtendPdf gen_neg("gen_neg", "extended signal pdf", _gen_neg, nsig_neg);

    // Create simultaneous pdf
    RooSimultaneous * gen;
    if (use_extended_nll) { gen = new RooSimultaneous("gen", "simultaneous pdf", {{"plus", &gen_pos}, {"minus", &gen_neg}}, h); } //TODO: Set these from helicity states...
    else { gen = new RooSimultaneous("gen", "simultaneous pdf", {{"plus", &_gen_pos}, {"minus", &_gen_neg}}, h); }

    // Fit the pdf to data
    std::unique_ptr<RooFitResult> r{gen->fitTo(*bin_ds, RooFit::Save(), RooFit::SumW2Error(use_sumW2Error), RooFit::PrintLevel(-1))};

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
    RooFormulaVar f_asym("f_asym","Asymmetry function",Form("%.3f*%s",pol,fitformula.c_str()), *arglist);//NOTE: NEED TO CORRECT FOR POLARIZATION FACTOR.

    // Plot projection of fitted distribution in x.
    RooPlot *xframe = x->frame(RooFit::Bins(xbins), RooFit::Title(Form("%s Projection, Bin: %s",fitvarxtitle.c_str(),bincut.c_str())));
    bin_ds->plotOn(xframe, RooFit::Asymmetry(h));
    f_asym.plotOn(xframe, RooFit::LineColor(kRed));

    // Draw the frame on the canvas
    std::string c1_x_name = Form("c1_%s__fitvarx_%s",outdir.c_str(),fitvarx.c_str());
    TCanvas *c1_x = new TCanvas(c1_x_name.c_str(), c1_x_name.c_str());
    gPad->SetLeftMargin(0.15);
    xframe->GetYaxis()->SetTitleOffset(1.4);
    xframe->Draw();
    c1_x->Print(Form("%s.pdf",c1_x_name.c_str()));

    // Plot projection of fitted distribution in y.
    RooPlot *yframe = y->frame(RooFit::Bins(ybins), RooFit::Title(Form("%s Projection, Bin: %s",fitvarytitle.c_str(),bincut.c_str())));
    bin_ds->plotOn(yframe, RooFit::Asymmetry(h));
    f_asym.plotOn(xframe, RooFit::LineColor(kRed));

    // Draw the frame on the canvas
    std::string c1_y_name = Form("c1_%s__fitvary_%s",outdir.c_str(),fitvary.c_str());
    TCanvas *c1_y = new TCanvas(c1_y_name.c_str(), c1_y_name.c_str());
    gPad->SetLeftMargin(0.15);
    yframe->GetYaxis()->SetTitleOffset(1.4);
    yframe->Draw();
    c1_y->Print(Form("%s.pdf",c1_y_name.c_str()));

    // Get fit parameters
    std::vector<double> pars;
    std::vector<double> Epars;
    for (int aa=0; aa<nparams; aa++) {
        pars.push_back((double)a[aa]->getVal());
        Epars.push_back((double)a[aa]->getError());
    }

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " getKinBinAsymUBML2D():" << std::endl;
    out << " bincut     = " << bincut.c_str() << std::endl;
    out << " bincount   = " << count << std::endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvar_means[idx] << "±" << binvar_errs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " pol        = " << pol << std::endl;
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
    out << " fitformula = " << fitformula.c_str() << std::endl;
    out << " nparams    = " << nparams <<std::endl;
    out << " initial params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << pars[idx] << "±" << Epars[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << " nsig_pos = " << (double)nsig_pos.getVal() << "±" << (double)nsig_pos.getError() << std::endl;
    out << " nsig_neg = " << (double)nsig_neg.getVal() << "±" << (double)nsig_neg.getError() << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Go back to parent directory
    outroot->cd("..");

    // Fill return array
    std::vector<double> arr; //NOTE: Dimension = 1+2*binvars.size()+2*depolvars.size()+2*nparams
    arr.push_back(count);
    for (int idx=0; idx<binvars.size(); idx++) {
        arr.push_back(binvar_means[idx]);
        arr.push_back(binvar_errs[idx]);
    }
    for (int idx=0; idx<depolvars.size(); idx++) {
        arr.push_back(depols[idx]);
        arr.push_back(depolerrs[idx]);
    }
    for (int idx=0; idx<nparams; idx++) {
        arr.push_back(pars[idx]);
        arr.push_back(Epars[idx]);
    }

    return arr;

} // std::vector<double> getKinBinAsymUBML2D()

/**
* @brief Loop kinematic bins and fit a 2D asymmetry.
*
* Loop bins and compute an asymmetry using an unbinned maximum likelihood fit.
* Optionally apply an invariant mass fit and background correction using the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* Note that the given asymmetry formula \f$ A(x, y, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) \f$ will be converted internally to a PDF of the form
*
* @f[
* \begin{aligned}
* PDF(h, x, y, &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* 1 + h \cdot P \cdot A(x, y, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...).
* \end{aligned}
* @f]
* and a simultaneous fit will be applied over the data subsets distinguished by the helicity states.  The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the corresponding depolarization factors.
*
* The variable names should be replaced in the fit formula by `x`\f$\rightarrow\f$`x[0]`, `y`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[2]`, etc.
*
* @param outdir Name of output directory
* @param outroot Name of output ROOT file
* @param frame ROOT RDataFrame
* @param method Asymmetry fit method
* @param binvar Kinematic binning variable
* @param nbins Number of kinematic variable bins
* @param bins List of bin limits, which must have dimension nbins+1
* @param pol Luminosity averaged beam polarization
* @param depolvars List of depolarization variables (up to 5)
* @param workspace_name Name of workspace in which to work
* @param workspace_title Title of workspace in which to work
* @param dataset_name Dataset name
* @param dataset_title Dataset title
* @param helicity Name of helicity variable
* @param fitvarx Name of fit variable 1
* @param xmin Minimum bound for fit variable 1
* @param xmax Maximum bound for fit variable 1
* @param fitvary Name of fit variable 2
* @param ymin Minimum bound for fit variable 2
* @param ymax Maximum bound for fit variable 2
* @param massvar Invariant mass variable name
* @param mmin Minimum bound for invariant mass variable
* @param mmax Maximum bound for invariant mass variable
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param model_name Full PDF name
* @param fitformula The asymmetry formula
* @param nparams Number of parameters in the asymmetry formula (up to 5)
* @param params List of initial values for asymmetry parameters
* @param fitvarxtitle Title of fit variable 1
* @param fitvarytitle Title of fit variable 2
* @param xbins Number of bins in fit variable 1 for plotting asymmetry
* @param ybins Number of bins in fit variable 2 for plotting asymmetry
* @param use_sumW2Error Option to use RooFit::SumW2Error(true) option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_splot Option to use sPlot method and perform fit with sWeighted dataset
* @param graph_title Title of kinematically binned graph of asymmetry parameters
* @param marker_color Graph ROOT marker color
* @param marker_style Graph ROOT marker style
* @param out Output stream
*/
void getKinBinnedAsymUBML2D(
        std::string outdir,
        TFile      *outroot,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        std::string     method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
        std::string binvar, // Variable name to bin in
        int         nbins, // Number of bins
        double      *bins, // Bin limits (length=nbins+1)
        double      pol,
        std::vector<std::string> depolvars,
        std::string workspace_name  = "w",
        std::string workspace_title = "workspace",
        std::string dataset_name    = "dataset",
        std::string dataset_title   = "dataset",
        std::string helicity   = "heli",
        std::string fitvarx         = "x",
        double      xmin            = 0.0,
        double      xmax            = 1.0,
        std::string fitvary         = "y",
        double      ymin            = 0.0,
        double      ymax            = 1.0,
        std::string massvar         = "mass_ppim",
        double      mmin            = 1.08,
        double      mmax            = 1.24,
        std::string sgYield_name    = "sgYield",
        std::string bgYield_name    = "bgYield",
        std::string model_name      = "model",
        std::string fitformula      = "[0]*sin(x)+[1]*sin(2*x)",
        int         nparams         = 2,
        std::vector<double> params  = std::vector<double>(5),
        std::string fitvarxtitle    = "x",
        std::string fitvarytitle    = "y",
        int         xbins           = 16,
        int         ybins           = 16,
        bool        use_sumW2Error  = true,
        bool        use_average_depol = false,
        bool        use_extended_nll = false,
        bool        use_splot       = true,
        std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
        int         marker_color    = 4, // 4 is blue
        int         marker_style    = 20, // 20 is circle
        std::ostream &out           = std::cout
    ) {

    // Check arguments
    if (method != "BSA2D") {std::cerr<<" *** ERROR *** Method must be BSA2D.  Exiting...\n"; return;}
    if (nbins<1) {std::cerr<<" *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}
    if (depolvars.size()!=nparams) {std::cerr<<" *** WARNING *** depolvars.size() does not match the number of parameters injected."<<std::endl;}

    // Starting message
    out << "----------------------- getKinBinnedAsymUBML2D ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double depols[nparams][nbins];
    double edepols[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];

    // Create workspace
    RooWorkspace *w = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());

    // Create bin var vector and outer bin lims vector
    std::vector binvars = {binvar};
    std::vector<std::vector<double>> binvarlims_outer = {{bins[0],bins[nbins]}};
    std::vector<std::vector<double>> depolvarlims;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvarlims.push_back({0.0,1.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.
    }
    std::map<std::string,int> helicity_states;
    helicity_states["plus"]  =  1;
    helicity_states["zero"]  =  0;
    helicity_states["minus"] = -1;

    // Create dataset
    data::createDataset(
        frame,
        w,
        dataset_name,
        dataset_title,
        helicity,
        helicity_states,
        {fitvarx,fitvary,massvar},
        {{xmin,xmax},{ymin,ymax},{mmin,mmax}},
        binvars,
        binvarlims_outer,
        depolvars,
        depolvarlims
    );

    // Apply Lambda mass fit
    if (use_splot) {
        applyLambdaMassFit(
            w,
            massvar,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name
        );
    }

    // Apply sPlot
    std::string fit_dataset_name = dataset_name;
    if (use_splot) {
        std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
        std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
        applySPlot(
            w,
            dataset_name,
            sgYield_name,
            bgYield_name,
            model_name,
            dataset_sg_name,
            dataset_bg_name
        );
        fit_dataset_name = dataset_sg_name;
    }

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {

        // Get bin limits
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Make bin cut on frame
        std::string  bincut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto binframe = frame.Filter(bincut.c_str());

        // Compute bin results
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        std::vector<double> bin_data = getKinBinAsymUBML2D(
                                binoutdir,
                                outroot,
                                frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                w,
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                helicity,
                                helicity_states,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xbins,
                                xmin,
                                xmax,
                                fitvary,
                                fitvarytitle,
                                ybins,
                                ymin,
                                ymax,
                                use_sumW2Error,
                                use_average_depol,
                                use_extended_nll,
                                out
                            );

        // Get bin data
        int k = 0;
        counts[binidx] = (int)bin_data[k++];
        xs[binidx]     = bin_data[k++]; //NOTE: ASSUME THIS IS 1D BINNING FOR NOW
        exs[binidx]    = bin_data[k++];
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx][binidx] = bin_data[k++];
            edepols[idx][binidx] = bin_data[k++];
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = bin_data[k++];
            eys[idx][binidx] = bin_data[k++];
        }

        // Divide out depolarization factors
        if (use_average_depol) {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx] / depols[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx] / depols[idx][binidx];
            }
        } else {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        }

        // Output message
        out << "--- Accpetance, background, and depolarization corrected signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            if (use_average_depol) {
                out << " ys["<< idx <<"]["<<binidx<<"]             = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]            = " << eys[idx][binidx] << "\n";
                out << " depols["<< idx <<"]["<<binidx<<"]         = " << depols[idx][binidx] << "\n";
                out << " edepols["<< idx <<"]["<<binidx<<"]        = " << edepols[idx][binidx] << "\n";
            }
            out << " ys_corrected["<< idx <<"]["<<binidx<<"]   = " << ys_corrected[idx][binidx] << "\n";
            out << " eys_corrected["<< idx <<"]["<<binidx<<"]  = " << eys_corrected[idx][binidx] << "\n";
        }
        out << "---------------------------\n";
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

        // Plot results graph
        TCanvas *c1 = new TCanvas();
        c1->SetBottomMargin(0.125);
        c1->cd(0);

        // Stylistic choices that aren't really necessary
        gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
        gr->SetMarkerSize(1.25);
        gr->GetXaxis()->SetTitleSize(0.05);
        gr->GetXaxis()->SetTitleOffset(0.9);
        gr->GetYaxis()->SetTitleSize(0.05);
        gr->GetYaxis()->SetTitleOffset(0.9);

        // More necessary stylistic choices
        gr->SetTitle(graph_title.c_str());
        gr->SetMarkerColor(marker_color); // 4  blue
        gr->SetMarkerStyle(marker_style); // 20 circle
        gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
        gr->GetXaxis()->SetTitle(binvar.c_str());
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
        gr->Draw("PA");

        // Add CLAS12 Preliminary watermark
        TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
        lt->SetTextAngle(45);
        lt->SetTextColor(18);
        lt->SetTextSize(0.1);
        lt->SetNDC();
        lt->Draw();

        // Add zero line
        TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
        f2->SetLineColor(1); // 1 black
        f2->SetLineWidth(1); // 1 thinnest
        f2->Draw("SAME");

        // Set outname and save
        TString fname;
        fname.Form("%s_%s_%s_%s_%.3f_%.3f_A%d",method.c_str(),fitvarx.c_str(),fitvary.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedAsymUBML2D -------------------\n";

} // getKinBinnedAsymUBML2D()

/**
* @brief Loop kinematic bins and fit a 1D asymmetry, correcting for background with sideband subtraction or <a href="http://arxiv.org/abs/physics/0402083">sPlots</a>.
*
* Loop bins and compute an asymmetry using an unbinned maximum likelihood fit.
* Optionally apply an invariant mass fit and background correction using the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* or the sideband subtraction method.  Note that the asymmetry fit formula \f$ A(x, y, a_0, a_1, a_2, ..., d_0, d_1, d_2, ... ) \f$ will be converted internally to a PDF of the form
*
* @f[
* \begin{aligned}
* PDF(h, x, y, &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + h \cdot P \cdot A(x, y, a_0, a_1, a_2, ..., d_0, d_1, d_2, ...),
* \end{aligned}
* @f]
*
* The variable names should be replaced in the fit formula by `x`\f$\rightarrow\f$`x[0]`, `y`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[2]`, etc.
*
* @param outdir Name of output directory
* @param outroot Name of output ROOT file
* @param frame ROOT RDataframe from which to create RooFit datasets
* @param method Asymmetry fit method
* @param binvar Kinematic binning variable
* @param nbins Number of kinematic variable bins
* @param bins List of bin limits, which must have dimension nbins+1
* @param pol Luminosity averaged beam polarization
* @param depolvars List of depolarization variables (up to 5)
* @param workspace_name Name of workspace in which to work
* @param workspace_title Title of workspace in which to work
* @param dataset_name Dataset name
* @param dataset_title Dataset title
* @param helicity Name of helicity variable
* @param fitvarx Name of first fit variable
* @param xmin Minimum bound for first fit variable
* @param xmax Maximum bound for first fit variable
* @param fitvary Name of second fit variable
* @param ymin Minimum bound for second fit variable
* @param ymax Maximum bound for second fit variable
* @param massvar Invariant mass variable name
* @param mass_min Minimum bound for invariant mass variable
* @param mass_max Maximum bound for invariant mass variable
* @param sgYield_name Signal yield variable name
* @param bgYield_name Background yield variable name
* @param model_name Full PDF name
* @param fitformula The asymmetry formula
* @param nparams Number of parameters in the asymmetry formula (up to 5)
* @param params List of initial values for asymmetry parameters
* @param fitvarxtitle Title of fit variable
* @param fitvarytitle Title of fit variable
* @param xbins Number of bins in fit variable 1 for plotting asymmetry
* @param ybins Number of bins in fit variable 2 for plotting asymmetry
* @param use_sumW2Error Option to use RooFit::SumW2Error() option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_splot Option to use sPlot method and perform fit with sWeighted dataset
* @param graph_title Title of kinematically binned graph of asymmetry parameters
* @param marker_color Graph ROOT marker color
* @param marker_style Graph ROOT style color
* @param out Output stream
* @param sgcut Signal region cut for sideband subtraction background correction
* @param bgcut Signal region cut for sideband subtraction background correction
* @param mass_nbins_hist Number of bins for the invariant mass histogram
* @param mass_nbins_conv Number of invariant mass bins for convolution of PDFs
* @param sig_pdf_name Signal PDF name
* @param sg_region_min Invariant mass signal region lower bound
* @param sg_region_max Invariant mass signal region upper bound
* @param use_sb_subtraction Option to use sideband subtraction for background correction
*/
void getKinBinnedAsym2D(
        std::string outdir,
        TFile      *outroot,
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        std::string     method, // ONLY getKinBinBSAGeneric ('BSA') is allowed at the moment
        std::string binvar, // Variable name to bin in
        int         nbins, // Number of bins
        double      *bins, // Bin limits (length=nbins+1)
        double      pol,
        std::vector<std::string> depolvars,
        std::string workspace_name  = "w",
        std::string workspace_title = "workspace",
        std::string dataset_name    = "dataset",
        std::string dataset_title   = "dataset",
        std::string helicity   = "heli",
        std::string fitvarx         = "x",
        double      xmin            = 0.0,
        double      xmax            = 1.0,
        std::string fitvary         = "y",
        double      ymin            = 0.0,
        double      ymax            = 1.0,
        std::string massvar         = "mass_ppim",
        double      mass_min        = 1.08,
        double      mass_max        = 1.24,
        std::string sgYield_name    = "sgYield",
        std::string bgYield_name    = "bgYield",
        std::string model_name      = "model",
        std::string fitformula      = "x[0]*(x[2]*x[4]*cos(x[1])+x[3]*x[5]*cos(2.0*x[1]))",
        int         nparams         = 2,
        std::vector<double> params  = std::vector<double>(2),
        std::string fitvarxtitle    = "x",
        std::string fitvarytitle    = "y",
        int         xbins           = 16,
        int         ybins           = 16,
        bool        use_sumW2Error  = true,
        bool use_average_depol      = false,
        bool use_extended_nll       = false,
        bool        use_splot       = true,
        std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
        int         marker_color    = 4, // 4 is blue
        int         marker_style    = 20, // 20 is circle
        std::ostream &out           = std::cout,
        std::string sgcut           = "Q2>1",
        std::string bgcut           = "Q2>1",
        int mass_nbins_hist         = 100,
        int mass_nbins_conv         = 1000,
        std::string sig_pdf_name    = "cb", //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")
        double sg_region_min        = 1.11,
        double sg_region_max        = 1.13,
        bool   use_sb_subtraction   = false
    ) {

    // Check arguments
    if (method != "BSA2D") {std::cerr<<" *** ERROR *** Method must be BSA2D.  Exiting...\n"; return;}
    if (nbins<1) {std::cerr<<" *** ERROR *** Number of " << binvar << " bins is too small.  Exiting...\n"; return;}
    if (depolvars.size()!=nparams) {std::cerr<<" *** WARNING *** depolvars.size() does not match the number of parameters injected."<<std::endl;}
    if (use_sb_subtraction && use_splot) {std::cerr<<" *** ERROR *** Cannot simultaneously use sideband subtraction and sPlot.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedAsym2D ----------------------\n";
    out << "Getting " << binvar << " binned hist...\n";
    out << "bins = { "; for (int i=0; i<nbins; i++) { out << bins[i] << " , ";} out << bins[nbins] << " }\n";

    // Make output directory in ROOT file and cd
    outroot->mkdir(outdir.c_str());
    outroot->cd(outdir.c_str());

    // Initialize data arrays
    double xs[nbins];
    double exs[nbins];
    int    counts[nbins];

    double ys[nparams][nbins];
    double eys[nparams][nbins];
    double ys_sb[nparams][nbins];
    double eys_sb[nparams][nbins];
    double depols[nparams][nbins];
    double edepols[nparams][nbins];
    double ys_corrected[nparams][nbins];
    double eys_corrected[nparams][nbins];

    // Create bin var vector and outer bin lims vector
    std::vector binvars = {binvar};
    std::vector<std::vector<double>> binvarlims_outer = {{bins[0],bins[nbins]}};
    std::vector<std::vector<double>> depolvarlims;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvarlims.push_back({0.0,1.0});//DEBUGGING: FOR NOW ASSUME ALL DEPOLARIZATION VARIABLES ARE IN THIS RANGE.
    }

    // Filter frames for signal and sideband
    auto frame_sg = frame.Filter(sgcut.c_str());
    auto frame_sb = frame.Filter(bgcut.c_str());

    // Loop bins and get data
    for (int binidx=0; binidx<nbins; binidx++) {

        // Get bin limits
        double bin_min = bins[binidx];
        double bin_max = bins[binidx+1];

        // Create workspace
        RooWorkspace *ws    = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());
        RooWorkspace *ws_sg = new RooWorkspace(Form("%s_sg",workspace_name.c_str()),Form("%s_signal",workspace_title.c_str()));
        RooWorkspace *ws_sb = new RooWorkspace(Form("%s_sb",workspace_name.c_str()),Form("%s_sideband",workspace_title.c_str())); //NOTE: Use separate signal and sideband workspaces for dataset, variable, and pdf name uniqueness.

        // Make bin cut on frame
        std::string  bincut = Form("(%s>=%.16f && %s<%.16f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
        auto binframe = frame.Filter(bincut.c_str());
        auto binframe_sg = frame_sg.Filter(bincut.c_str());

        std::map<std::string,int> helicity_states;
        helicity_states["plus"]  =  1;
        helicity_states["zero"]  =  0;
        helicity_states["minus"] = -1;

        // Create bin dataset
        data::createDataset(
            binframe,
            ws,
            dataset_name,
            dataset_title,
            helicity,
            helicity_states,
            {fitvarx,fitvary,massvar},
            {{xmin,xmax},{ymin,ymax},{mass_min,mass_max}},
            binvars,
            binvarlims_outer,
            depolvars,
            depolvarlims
        );

        // Apply Lambda mass fit to FULL bin frame
        std::string bin_id = Form("%s_%.3f_%.3f",binvar.c_str(),bin_min,bin_max);
        std::vector<double> epss = applyLambdaMassFit(
                ws,
                massvar,
                dataset_name,
                sgYield_name,
                bgYield_name,
                binframe,
                mass_nbins_hist,
                mass_min,
                mass_max,
                mass_nbins_conv,
                model_name,
                sig_pdf_name,
                sg_region_min,
                sg_region_max,
                "",//ws_unique_id->This changes pdf,yieldvar names, but NOT (bin,depol,mass,fit)vars,pdf parameters which are saved internally.
                1,//use_poly4_bg
                bin_id
            );

        // Apply sPlot
        std::string fit_dataset_name = dataset_name; // -> Use this for sPlot
        if (use_splot) {
            std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
            std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
            applySPlot(
                ws,
                dataset_name,
                sgYield_name,
                bgYield_name,
                model_name,
                dataset_sg_name,
                dataset_bg_name
            );
            fit_dataset_name = dataset_sg_name;
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
                {fitvarx,fitvary,massvar},
                {{xmin,xmax},{ymin,ymax},{mass_min,mass_max}},
                binvars,
                binvarlims_outer,
                depolvars,
                depolvarlims
            );
        }

        // Compute signal region bin results
        std::string  binoutdir = Form("method_%s_bin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        std::vector<double> bin_data = getKinBinAsymUBML2D(
                                binoutdir,
                                outroot,
                                (use_sb_subtraction ? binframe_sg : binframe), //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                (use_sb_subtraction ? ws_sg : ws),
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                helicity,
                                helicity_states,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xbins,
                                xmin,
                                xmax,
                                fitvary,
                                fitvarytitle,
                                ybins,
                                ymin,
                                ymax,
                                use_sumW2Error,
                                use_average_depol,
                                use_extended_nll,
                                out
                            );

        // Compute sideband region bin results
        std::string  sbbinoutdir = Form("method_%s_sbbin_%s_%.3f_%.3f",method.c_str(),binvar.c_str(),bin_min,bin_max);
        std::vector<double> bin_data_sb;
        if (use_sb_subtraction) {

            // Make bin cut on sideband frame
            auto binframe_sb = frame_sb.Filter(bincut.c_str());

            // Create sideband dataset
            data::createDataset(
                binframe_sb,
                ws_sb, //NOTE: Use separate sideband workspace for dataset, variable, and pdf name uniqueness.
                dataset_name,
                dataset_title,
                helicity,
                helicity_states,
                {fitvarx,fitvary,massvar},
                {{xmin,xmax},{ymin,ymax},{mass_min,mass_max}},
                binvars,
                binvarlims_outer,
                depolvars,
                depolvarlims
            );

            // Compute sideband bin results
            bin_data_sb = getKinBinAsymUBML2D(
                                binoutdir,
                                outroot,
                                binframe_sb, //NOTE: FRAME SHOULD ALREADY BE FILTERED WITH OVERALL CUTS
                                ws_sb,
                                dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                binvars,
                                binvarlims_outer,
                                bincut,
                                (std::string)Form("bin_%d",binidx),//Bintitle
                                pol,
                                depolvars,
                                helicity,
                                helicity_states,
                                fitformula,
                                nparams,
                                params,
                                fitvarx,
                                fitvarxtitle,
                                xbins,
                                xmin,
                                xmax,
                                fitvary,
                                fitvarytitle,
                                ybins,
                                ymin,
                                ymax,
                                use_sumW2Error,
                                use_average_depol,
                                use_extended_nll,
                                out
                            );
        } 

        // Get bin data
        int k = 0;
        counts[binidx] = (int)bin_data[k++];
        xs[binidx]     = bin_data[k++]; //NOTE: ASSUME THIS IS 1D BINNING FOR NOW
        exs[binidx]    = bin_data[k++];
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx][binidx] = bin_data[k++];
            edepols[idx][binidx] = bin_data[k++];
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx][binidx] = bin_data[k++];
            eys[idx][binidx] = bin_data[k++];
        }

        // Apply sideband subtraction to asymmetries assuming that depolarization factors do not vary much
        double epsilon, epsilon_err;
        if (use_sb_subtraction) {
            int k2 = 3 + depolvars.size();
            epsilon = epss[0];
            epsilon_err = epss[1];
            for (int idx=0; idx<nparams; idx++) {
                ys_sb[idx][binidx] = bin_data_sb[k2++];
                eys_sb[idx][binidx] = bin_data_sb[k2++];
                ys[idx][binidx]  = (ys[idx][binidx] - epsilon * ys_sb[idx][binidx]) / (1.0 - epsilon);
                eys[idx][binidx] = TMath::Sqrt(eys[idx][binidx]*eys[idx][binidx] + epsilon * epsilon * eys_sb[idx][binidx]*eys_sb[idx][binidx]) / (1.0 - epsilon) ;
            }
        }

        // Divide out depolarization factors
        if (use_average_depol) {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx] / depols[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx] / depols[idx][binidx];
            }
        } else {
            for (int idx=0; idx<nparams; idx++) {
                ys_corrected[idx][binidx] = ys[idx][binidx];
                eys_corrected[idx][binidx] = eys[idx][binidx];
            }
        }

        // Output message
        out << "--- Accpetance, background, and depolarization corrected signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            if (use_average_depol) {
                out << " ys["<< idx <<"]["<<binidx<<"]             = " << ys[idx][binidx] << "\n";
                out << " eys["<< idx <<"]["<<binidx<<"]            = " << eys[idx][binidx] << "\n";
                out << " depols["<< idx <<"]["<<binidx<<"]         = " << depols[idx][binidx] << "\n";
                out << " edepols["<< idx <<"]["<<binidx<<"]        = " << edepols[idx][binidx] << "\n";
            } if (use_sb_subtraction) {
                out << " ys_sb["<< idx <<"]["<<binidx<<"]       = " << ys_sb[idx][binidx] << "\n";
                out << " eys_sb["<< idx <<"]["<<binidx<<"]      = " << eys_sb[idx][binidx] << "\n";
            }
            out << " ys_corrected["<< idx <<"]["<<binidx<<"]   = " << ys_corrected[idx][binidx] << "\n";
            out << " eys_corrected["<< idx <<"]["<<binidx<<"]  = " << eys_corrected[idx][binidx] << "\n";
        }
        out << "---------------------------\n";
    }

    // Loop results and plot
    for (int idx=0; idx<nparams; idx++) {

        // Create graph of results binned in binvar
        TGraphErrors *gr = new TGraphErrors(nbins,xs,ys_corrected[idx],exs,eys_corrected[idx]);
        gr->Write("gr");

        // Plot results graph
        TCanvas *c1 = new TCanvas();
        c1->SetBottomMargin(0.125);
        c1->cd(0);

        // Stylistic choices that aren't really necessary
        gStyle->SetEndErrorSize(5); gStyle->SetTitleFontSize(0.05);
        gr->SetMarkerSize(1.25);
        gr->GetXaxis()->SetTitleSize(0.05);
        gr->GetXaxis()->SetTitleOffset(0.9);
        gr->GetYaxis()->SetTitleSize(0.05);
        gr->GetYaxis()->SetTitleOffset(0.9);

        // More necessary stylistic choices
        gr->SetTitle(graph_title.c_str());
        gr->SetMarkerColor(marker_color); // 4  blue
        gr->SetMarkerStyle(marker_style); // 20 circle
        gr->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);                                                       
        gr->GetXaxis()->SetTitle(binvar.c_str());
        std::string ytitle = Form("BSA A%d",idx);
        gr->GetYaxis()->SetTitle(ytitle.c_str());
        gr->Draw("PA");

        // Add CLAS12 Preliminary watermark
        TLatex *lt = new TLatex(0.3,0.2,"CLAS12 Preliminary");
        lt->SetTextAngle(45);
        lt->SetTextColor(18);
        lt->SetTextSize(0.1);
        lt->SetNDC();
        lt->Draw();

        // Add zero line
        TF1 *f2 = new TF1("f2","0",bins[0],bins[nbins]);
        f2->SetLineColor(1); // 1 black
        f2->SetLineWidth(1); // 1 thinnest
        f2->Draw("SAME");

        // Set outname and save
        TString fname;
        fname.Form("%s_%s_%s_%.3f_%.3f_A%d",method.c_str(),fitvarx.c_str(),binvar.c_str(),bins[0],bins[nbins],idx);
        c1->Print(fname+".pdf");
        gr->SaveAs(fname+".root","recreate");

        // Output message
        out << " Saved graph to " << fname << ".root\n";
    }

    // Cd out of outdir
    outroot->cd("..");

    // Ending message
    out << "------------------- END of getKinBinnedAsym2D -------------------\n";

} // getKinBinnedAsym2D()

} // namespace analysis {