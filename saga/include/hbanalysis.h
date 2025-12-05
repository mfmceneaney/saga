#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <map>

// ROOT Includes
#include <TStyle.h>
#include <TCanvas.h>
// #include <TAxis.h>
// #include <TLegend.h>
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
#include <RooFormulaVar.h>
// #include <RooProduct.h>
#include <RooDataSet.h>
// #include <RooPlot.h>
// #include <RooAbsDataHelper.h>
// #include <RooDataHist.h>
#include <RooArgList.h>
// #include <RooAddPdf.h>
// #include <RooGenericPdf.h>
// #include <RooExtendPdf.h>
// #include <RooSimultaneous.h>
// #include <RooFFTConvPdf.h>
// #include <RooCrystalBall.h>
// #include <RooLandau.h>
// #include <RooGaussian.h>
// #include <RooChebychev.h>
// #include <RooFitResult.h>
#include <RooWorkspace.h>

// // RooStats Includes
// #include <RooStats/SPlot.h>

// Local Includes
#include <data.h>
#include <bins.h>
#include <util.h>
#include <signal.h>

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 5/Dec./2025
* @version 0.0.0
* @brief Fit asymmetries using the Helicity Balance (HB) method and sideband subtraction 
* or the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a> for background correction.
*/

namespace saga {

namespace hbanalysis {

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
* @brief Fit an asymmetry with the Helicity Balance (HB) method.
*
* Compute the bin count, bin variable mean values and variances, depolarization variable values and errors,
* and fit the asymmetry with the Helicity Balance (HB) method.  The asymmetry parameter will be computed with:
* @f[
* D^{\Lambda}_{LL'} = \frac{1}{\alpha_{\Lambda} \overline{\lambda_{\ell}^2}}\frac{\sum^{N_{\Lambda}}_{i=1}\lambda_{\ell,i}
\cos{\theta_{LL'}^i}}{\sum^{N_{\Lambda}}_{i=1}D(y_i) \cos^2{\theta_{LL'}^i}} \,,
* @f]
* where `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the depolarization factors.
* Here, \f$\lambda_{\ell,i}\f$ indicates the beam helicity for a given event \f$i\f$i,
* and \f$\overline{\lambda^{2}_{\ell}}\f$ is the luminosity averaged beam polarization.
* The method relies on the assumption that the luminosity averaged helicity \f$\overline{\lambda_{\ell}}=0\f$
* to allow the acceptance method to cancel out.
* See <a href="https://cds.cern.ch/record/732977">Gunar Schnell's thesis</a> from New Mexico State University, 1999 for a full derivation.
* \f$N_{\Lambda}\f$ is the number of \f$\Lambda\f$ events in the bin,
* and \f$D(y_i)\f$ and \f$\cos{\theta_{LL'}^i}\f$ are the depolarization factor and the decay angle
* in the \f$\Lambda\f$ CM frame respectively for the given event.
* Similarly, the error will be computed as follows.  Letting
* @f[
* \begin{aligned}
*     A& = \lambda_{\ell,i} \cos{\theta_{LL'}^i}, \\
*     B& =D(y_i) \cos^2{\theta_{LL'}^i} \,,
* \end{aligned}
* @f]
* the statistical scale uncertainty may be expressed as
* @f[
* \begin{aligned}
*     \bigg{(}\frac{\delta D^{\Lambda}_{LL'}}{D^{\Lambda}_{LL'}}\bigg{)}^2& = \bigg{[}\delta\bigg{(}\frac{\text{Sum}[A]}{\text{Sum}[B]}\bigg{)}\bigg{/}\bigg{(}\frac{\text{Sum}[A]}{\text{Sum}[B]}\bigg{)}\bigg{)}\bigg{]}^2  \\
*     & = \bigg{(} \frac{\text{Var}[A]^2}{\text{Sum}[A]^2} + \frac{\text{Var}[B]^2}{\text{Sum}[B]^2} - 2 \frac{\text{Var}[AB]}{\text{Sum}[A]\text{Sum}[B]}\bigg{)} \,.
* \end{aligned}
* @f]
* Here, \f$\text{Var}[X_i] = \sum_i (X_i - \text{Mean}[X_i])^2\f$ denotes the variance of
* a quantity \f$X_i\f$, and Sum and Mean are exactly the operations named.
* Since both \f$A\f$ and \f$B\f$ are polynomial functions of \f$\cos{\theta_{LL'}}\f$, \f$A\f$ and \f$B\f$ are correlated.
* Thus, we include the covariance term \f$\text{Var[AB]}\f$ in the uncertainty calculation.
*
* The returned vector will have the following entries:
*
* - Bin count
*
* - For each bin variable:
*
*   - Bin variable mean value
*
*   - Bin variable standard deviation
*
* - For each depolarization variable:
*
*   - Depolarization variable mean value
*
*   - Depolarization variable standard deviation
*
* - The raw asymmetries and errors using the actual counts
*   **or**, in the case of an extended fit, using the fitted counts, for each of
*
*   - Beam helicity \f$\lambda_{\ell}\f$
*
*   - Target spin \f$S\f$
*
*   - Beam helicity times target spin \f$\lambda_{\ell}\cdot S\f$
*
* - For the (only) Helicity Balance parameter:
*
*   - Helicity Balance parameter mean value
*
*   - Helicity Balance parameter error
*
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param bpol Luminosity averaged beam polarization \f$\overline{\lambda_{\ell}^2}\f$
* @param tpol Luminosity averaged target polarization \f$\overline{S^2}\f$
* @param alpha Lambda decay asymmetry parameter \f$\alpha_{\Lambda}\f$
* @param helicity Name of the helicity variable
* @param tspin Name of the target spin variable
* @param htspin Name of the beam helicity times target spin variable
* @param combined_spin_state Name of the combined spin state variable
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param depolvars List of depolarization variables
* @param fitvars List of asymmetry fit variables
* @param out Output stream
*
* @return List of bin count, bin variable means and errors, depolarization variable means and errors, fit parameters and errors
*
* @throws Runtime error
*/
vector<double> fitHB(
        RooWorkspace               *w,
        string                      dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
        double                      bpol,
        double                      tpol,
        double                      alpha,
        string                      helicity,
        string                      tspin,
        string                      htspin,
        string                      combined_spin_state,
        string                      binid,
        string                      bincut,
        vector<string>              binvars,
        vector<string>              depolvars,
        vector<string>              fitvars,
        ostream &out                = cout
    ) {

    // Set method name
    string method_name = "fitHB";
    int nparams = 1;

    // Check arguments
    if (depolvars.size()!=1) {
        string msg = Form("[%s]: Number of depolarization variables must be 1", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }
    if (fitvars.size()!=1) {
        string msg = Form("[%s]: Number of fit variables must be 1", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }

    // Load helicity variable from workspace
    LOG_DEBUG(Form("[%s] Loading RooCategory variables %s, %s, %s, and %s from workspace...", method_name.c_str(), helicity.c_str(), tspin.c_str(), htspin.c_str(), combined_spin_state.c_str()));
    RooCategory * h  = w->cat(helicity.c_str());
    RooCategory * t  = w->cat(tspin.c_str());
    RooCategory * ht = w->cat(htspin.c_str());
    RooCategory * ss = w->cat(combined_spin_state.c_str());

    // Load fit variables from workspace
    RooRealVar * f[(const int)fitvars.size()];
    for (int i=0; i<fitvars.size(); i++) {
        LOG_DEBUG(Form("[%s] Loading RooRealVar fit variable %s from workspace...", method_name.c_str(), fitvars[i].c_str()));
        f[i] = w->var(fitvars[i].c_str());
    }

    // Load dataset from workspace
    LOG_DEBUG(Form("[%s] Loading RooDataSet %s from workspace...", method_name.c_str(), dataset_name.c_str()));
    RooDataSet *ds = (RooDataSet*)w->data(dataset_name.c_str());

    // Apply bin cuts
    LOG_DEBUG(Form("[%s] Applying bin cut: %s", method_name.c_str(), bincut.c_str()));
    RooDataSet *bin_ds = (RooDataSet*)ds->reduce(bincut.c_str());

    // Get count
    LOG_DEBUG(Form("[%s] Getting bin count...", method_name.c_str()));
    auto count = (int)bin_ds->sumEntries();

    // Get bin variable means and errors
    vector<double> binvarmeans;
    vector<double> binvarerrs;
    RooRealVar * b[(const int)binvars.size()];
    for (int i=0; i<binvars.size(); i++) {
        LOG_DEBUG(Form("[%s] Getting mean and error for bin variable %s...", method_name.c_str(), binvars[i].c_str()));
        b[i] = w->var(binvars[i].c_str());
        double mean   = bin_ds->mean(*b[i]);
        double stddev = sqrt(bin_ds->moment(*b[i],2.0));
        binvarmeans.push_back(mean);
        binvarerrs.push_back(stddev);
    }

    // Get depolarization factor means and errors
    vector<double> depols;
    vector<double> depolerrs;
    RooRealVar * d[(const int)depolvars.size()];
    for (int i=0; i<depolvars.size(); i++) {
        LOG_DEBUG(Form("[%s] Getting mean and error for depolarization variable %s...", method_name.c_str(), depolvars[i].c_str()));
        d[i] = w->var(depolvars[i].c_str());
        double mean   = bin_ds->mean(*d[i]);
        double stddev = sqrt(bin_ds->moment(*d[i],2.0));
        depols.push_back(mean);
        depolerrs.push_back(stddev);
    }

    // Add helicity*depolvar*fitvar (hdf)
    string hdf_var_name = Form("%s%s_hdf", method_name.c_str(), binid.c_str());
    string hdf_var_formula = Form("%s*%s*%s", h->GetName(), d[0]->GetName(), f[0]->GetName());
    RooFormulaVar hdf_fvar(
        hdf_var_name.c_str(),
        hdf_var_formula.c_str(),
        RooArgSet(*h, *d[0], *f[0])
    );
    RooRealVar* hdf_var = (RooRealVar*)bin_ds->addColumn(hdf_fvar);
    double hdf_var_mean   = bin_ds->mean(*hdf_var);
    double hdf_var_stddev = sqrt(bin_ds->moment(*hdf_var,2.0));

    // Add depolvar^2*fitvar^2 (d2f2)
    string d2f2_var_name = Form("%s%s_d2f2", method_name.c_str(), binid.c_str());
    string d2f2_var_formula = Form("%s*%s*%s*%s", d[0]->GetName(), f[0]->GetName(), d[0]->GetName(), f[0]->GetName());
    RooFormulaVar d2f2_fvar(
        d2f2_var_name.c_str(),
        d2f2_var_formula.c_str(),
        RooArgSet(*d[0], *f[0])
    );
    RooRealVar* d2f2_var = (RooRealVar*)bin_ds->addColumn(d2f2_fvar);
    double d2f2_var_mean   = bin_ds->mean(*d2f2_var);
    double d2f2_var_stddev = sqrt(bin_ds->moment(*d2f2_var,2.0));

    // Add covariance (covar_hdf_d2f2)
    string covar_hdf_d2f2_var_name = Form("%s%s_covar_hdf_d2f2", method_name.c_str(), binid.c_str());
    string covar_hdf_d2f2_var_formula = Form("(%s-%.8f)*(%s-%.8f)", hdf_var->GetName(), hdf_var_mean, d2f2_var->GetName(), d2f2_var_mean);
    RooFormulaVar covar_hdf_d2f2_fvar(
        covar_hdf_d2f2_var_name.c_str(),
        covar_hdf_d2f2_var_formula.c_str(),
        RooArgSet(*hdf_var, *d2f2_var)
    );
    RooRealVar* covar_hdf_d2f2_var = (RooRealVar*)bin_ds->addColumn(covar_hdf_d2f2_fvar);
    double covar_hdf_d2f2_var_mean   = bin_ds->mean(*covar_hdf_d2f2_var);
    // double covar_hdf_d2f2_var_stddev = sqrt(bin_ds->moment(*covar_hdf_d2f2_var,2.0));

    // Compute asymmetry with helicity balance method
    double a = 0.0;
    if (d2f2_var_mean==0) {
        LOG_WARN(Form("[%s]: d2f2_var_mean=0.  Setting asymmetry to zero to avoid divison errors.",method_name.c_str()));
    } else {
        a = hdf_var_mean / (d2f2_var_mean * alpha * bpol);
    }

    // Compute asymmetry error with helicity balance method
    // NOTE: #sigma^2 = variance / count but #mu^2 = (sum / count)^2
    // so need extra factor of 1/count if dividing
    double a_err = 0.0;
    if (d2f2_var_mean==0) {
        LOG_WARN(Form("[%s]: d2f2_var_mean=0 || hdf_var_mean==0.  Setting asymmetry error to zero to avoid divison errors.",method_name.c_str()));
    } else {
        a_err = abs(a) * sqrt(
            hdf_var_stddev*hdf_var_stddev / (count * hdf_var_mean * hdf_var_mean)
            + d2f2_var_stddev*d2f2_var_stddev / (count * d2f2_var_mean * d2f2_var_mean)
            - 2 * covar_hdf_d2f2_var_mean / (count * hdf_var_mean * d2f2_var_mean)
        );
    }

    // Get fit parameter values and errors
    vector<double> params = { a };
    vector<double> paramerrs = { a_err };

    // Get the raw counts and poissonian errors
    LOG_DEBUG(Form("[%s] Getting raw counts and poissonian errors for helicity and spin states...", method_name.c_str()));
    vector<double> counts;
    vector<double> counterrs;
    double count_h_pos  = (double)bin_ds->reduce(Form("%s>0",h->GetName()))->sumEntries();
    double count_h_neg  = (double)bin_ds->reduce(Form("%s<0",h->GetName()))->sumEntries();
    double count_t_pos  = (double)bin_ds->reduce(Form("%s>0",t->GetName()))->sumEntries();
    double count_t_neg  = (double)bin_ds->reduce(Form("%s<0",t->GetName()))->sumEntries();
    double count_ht_pos = (double)bin_ds->reduce(Form("%s>0",ht->GetName()))->sumEntries();
    double count_ht_neg = (double)bin_ds->reduce(Form("%s<0",ht->GetName()))->sumEntries();
    double counterr_h_pos   = (double)sqrt(count_h_pos);
    double counterr_h_neg   = (double)sqrt(count_h_neg);
    double counterr_t_pos   = (double)sqrt(count_t_pos);
    double counterr_t_neg   = (double)sqrt(count_t_neg);
    double counterr_ht_pos  = (double)sqrt(count_ht_pos);
    double counterr_ht_neg  = (double)sqrt(count_ht_neg);

    // Set the raw asymmetries
    LOG_DEBUG(Form("[%s] Calculating raw asymmetries and errors...", method_name.c_str()));
    vector<double> rawasyms;
    vector<double> rawasymerrs;
    double asym_h  = (count_h_pos-count_h_neg)/(count_h_pos+count_h_neg);
    double asym_t  = (count_t_pos-count_t_neg)/(count_t_pos+count_t_neg);
    double asym_ht = (count_ht_pos-count_ht_neg)/(count_ht_pos+count_ht_neg);
    rawasyms.push_back(asym_h);
    rawasyms.push_back(asym_t);
    rawasyms.push_back(asym_ht);
    double asymerr_h  = (double)sqrt(4.0*count_h_pos*count_h_neg/TMath::Power(count,3)); //NOTE: Use binomial error assuming correlated counts from: http://blast.lns.mit.edu/BlastTalk/archive/att-5707/01-asymmetry_calculations.pdf
    double asymerr_t  = (double)sqrt(4.0*count_t_pos*count_t_neg/TMath::Power(count,3)); //NOTE: Assume q=N^+/N is your random variable following a binomial distribution
    double asymerr_ht = (double)sqrt(4.0*count_ht_pos*count_ht_neg/TMath::Power(count,3)); //NOTE: A = q - (q-1) = 2q-1 ---> (dN^{+})^2 = N*q*(1-q) ---> dA^2 = 4*N^{+}*N^{-}/N^3
    rawasymerrs.push_back(asymerr_h);
    rawasymerrs.push_back(asymerr_t);
    rawasymerrs.push_back(asymerr_ht);

    // Print out fit info
    LOG_DEBUG(Form("[%s] Printing fit information...", method_name.c_str()));
    out << "--------------------------------------------------" << endl;
    out << " "<<method_name.c_str()<<"():" << endl;
    out << " bpol        = " << bpol << endl;
    out << " tpol        = " << tpol << endl;
    out << " alpha       = " << alpha << endl;
    out << " bincut     = " << bincut.c_str() << endl;
    out << " bincount   = " << count << endl;
    out << " binvar means = [" ;
    for (int idx=0; idx<binvars.size(); idx++) {
        out << binvarmeans[idx] << "±" << binvarerrs[idx];
        if (idx<binvars.size()-1) { out << " , "; }
    }
    out << "]" << endl;
    out << " depolvars  = [" ;
    for (int idx=0; idx<depolvars.size(); idx++) {
        out << depolvars[idx];
        if (idx<depolvars.size()-1) { out << " , "; }
    }
    out << "]" << endl;
    out << " depols  = [" ;
    for (int idx=0; idx<depols.size(); idx++) {
        out << depols[idx];
        if (idx<depols.size()-1) { out << " , "; }
    }
    out << "]" << endl;
    out << " fitvars  = [" ;
    for (int idx=0; idx<fitvars.size(); idx++) {
        out << fitvars[idx];
        if (idx<fitvars.size()-1) { out << " , "; }
    }
    out << "]" << endl;
    out << " params = [" ;
    for (int idx=0; idx<nparams; idx++) {
        out << params[idx] << "±" << paramerrs[idx];
        if (idx<nparams-1) { out << " , "; }
    }
    out << "]" << endl;
    out << " rawasyms = [" ;
    for (int idx=0; idx<rawasyms.size(); idx++) {
        out << rawasyms[idx] << "±" << rawasymerrs[idx];
        if (idx<rawasyms.size()-1) { out << " , "; }
    }
    out << "]" << endl;
    out << "--------------------------------------------------" << endl;

    // Fill return array
    LOG_DEBUG(Form("[%s] Filling return array...", method_name.c_str()));
    vector<double> arr; //NOTE: Dimension = 1+2*binvars.size()+2*depolvars.size()+2*rawasyms.size()+2*nparams
    arr.push_back(count);
    for (int idx=0; idx<binvars.size(); idx++) {
        arr.push_back(binvarmeans[idx]);
        arr.push_back(binvarerrs[idx]);
    }
    for (int idx=0; idx<depolvars.size(); idx++) {
        arr.push_back(depols[idx]);
        arr.push_back(depolerrs[idx]);
    }
    for (int idx=0; idx<rawasyms.size(); idx++) {
        arr.push_back(rawasyms[idx]);
        arr.push_back(rawasymerrs[idx]);
    }
    for (int idx=0; idx<nparams; idx++) {
        arr.push_back(params[idx]);
        arr.push_back(paramerrs[idx]);
    }

    return arr;

} // vector<double> fitHB()

} // namespace hbanalysis {

} // namespace saga {
