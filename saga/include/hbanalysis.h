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

/**
* @brief Loop kinematic bins and fit an asymmetry, correcting for background with sideband subtraction or <a href="http://arxiv.org/abs/physics/0402083">sPlots</a>.
*
* Loop bins cuts and fit an asymmetry with the `saga::hbanalysis::fitHB()` method.  Optionally, apply an invariant mass fit and background correction using the
* sideband subtraction method or the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a>.
* The mass fit will be applied with `saga::signal::fitMass()` and the sPlot method will use `saga::signal::applySPlot()`.
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
* - For each depolarization variable `depolvar`
*
*   - `<depolvar>`: Mean value
*
*   - `<depolvar>_err`: Standard deviation
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
* - For each asymmetry fit parameter `asymfitpar` (only one allowed for the HB method)
*
*   - `<asymfitpar>`: Final parameter value
*
*   - `<asymfitpar>_err`: Final parameter error
*
* The following columns will be added in the case of a single mass fit for applied to the entire bin:
*
* - `int_sg_pdf_val`: Signal PDF integral \f$N_{SG}^{PDF}\f$ value in the signal region
*
* - `int_sg_pdf_err`: Signal PDF integral error  \f$\delta N_{SG}^{PDF}\f$ in the signal region
*
* - `int_bg_pdf_val`: Background PDF integral \f$N_{BG}^{PDF}\f$ in the signal region
*
* - `int_bg_pdf_err`: Background PDF integral error \f$\delta N_{BG}^{PDF}\f$ in the signal region
*
* - `int_model_pdf_val`: Full PDF integral \f$N^{PDF}\f$ in the signal region
*
* - `int_model_pdf_err`: Full PDF integral error \f$\delta N^{PDF}\f$ in the signal region
*
* - `int_ds_val` Full dataset sum \f$N^{DS}\f$ in the signal region
*
* - `int_ds_err`: Poissonian error \f$\sqrt{N^{DS}}\f$ of the full dataset sum in the signal region
*
* - `eps_bg_pdf`: Background fraction \f$\varepsilon_{1} = \frac{N_{BG}^{PDF}}{N^{DS}}\f$
*
* - `eps_bg_pdf_err`: Background fraction error \f$\delta\varepsilon_{1}\f$
*
* - `eps_sg_pdf`: Background fraction \f$\varepsilon_{2} = 1 - \frac{N_{SG}^{PDF}}{N^{DS}}\f$
*
* - `eps_sg_pdf_err`: Background fraction error \f$\delta\varepsilon_{2}\f$
*
* - `eps_pdf`: Background fraction \f$\varepsilon_{3} = 1 - \frac{N_{SG}^{PDF}}{N^{PDF}}\f$
*
* - `eps_pdf_err`: Background fraction error \f$\delta\varepsilon_{3}\f$
*
* - For each mass fit variable:
*
*   - `<chi2>`: \f$\chi^2\f$ value of the 1D projection of the full PDF in that variable
*
* - For each mass fit signal PDF parameter `massfitpar_sg`
*
*   - `<massfitpar_sg>`: Final parameter value
*
*   - `<massfitpar_sg>_err`: Final parameter error
*
* - For each mass fit background PDF parameter `massfitpar_bg`
*
*   - `<massfitpar_bg>`: Final parameter value
*
*   - `<massfitpar_bg>_err`: Final parameter error
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
* @param asymfitvars List of asymmetry fit variables names
* @param asymfitvar_titles List of asymmetry fit variables titles
* @param asymfitvar_lims List asymmetry fit variable minimum and maximum bounds 
* @param asymfitvar_bins List of asymmetry fit variables bins
* @param massfitvars List of invariant mass fit variables names
* @param massfitvar_titles List of invariant mass fit variables titles
* @param massfitvar_lims List invariant mass fit variable minimum and maximum bounds 
* @param massfitvar_bins List of invariant mass fit variables bins

* @param bpol Luminosity averaged beam polarization \f$\overline{\lambda_{\ell}^2}\f$
* @param tpol Luminosity averaged target polarization \f$\overline{S^2}\f$
* @param alpha Lambda decay asymmetry parameter \f$\alpha_{\Lambda}\f$

* @param massfit_yamlfile_map Map of bin ids to the paths of yaml files specifying the remaining mass fit arguments.  Note that the values specified here will function as the defaults.
* @param massfit_pdf_name Base name of the combined signal and background pdf.  Note, the actual PDF name will be: `<pdf_name>_<binid>`.
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

* @param use_splot Option to use sPlot method and perform fit with sWeighted dataset

* @param massfit_sgcut Signal region cut for sideband subtraction background correction.  Note, this will automatically be formed from `massfit_sgregion_lims` if not specified.
* @param massfit_bgcut Background region cut for sideband subtraction background correction
* @param use_sb_subtraction Option to use sideband subtraction for background correction
* @param use_binned_sb_bgfracs Option to use background fractions from invariant mass fits binned in the asymmetry fit variable for background correction
* @param asymfitvar_bincuts Map of unique bin id ints to bin variable cuts for asymmetry fit variable bins
* @param bgfracvar Name of binned background fraction variable
* @param bgfracvar_lims List of binned background fraction variable minimum and maximum bounds
* @param bgfrac_idx Index to select which formulation to use for the background fraction in `saga::signal::setBinnedBGFractions()`

* @param massfit_plot_bg_pars Option to plot background pdf parameters on TLegend for the signal and background mass fit
* @param massfit_lg_text_size Size of TLegend text for the signal and background mass fit
* @param massfit_lg_margin Margin of TLegend for the signal and background mass fit
* @param massfit_lg_ncols Number of columns in TLegend for the signal and background mass fit
* @param massfit_use_sumw2error Option to use `RooFit::SumW2Error(true)` option for the signal and background mass fit which is necessary if using a weighted dataset 
* @param massfit_use_extended_nll Option to use an extended Negative Log Likelihood function for minimization for the signal and background mass fit
* @param massfit_use_binned_fit Option to use a binned fit to the data for the signal and background mass fit

* @param out Output stream

* @throws runtime_error if invalid arguments are provided
*/
void getKinBinnedHB(
        string                      scheme_name,
        RNode                            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        string                      workspace_name,
        string                      workspace_title,

        // parameters passed to data::createDataset()
        string                      dataset_name,
        string                      dataset_title,
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
        vector<string>         asymfitvars,
        vector<string>         asymfitvar_titles,
        vector<vector<double>> asymfitvar_lims,
        vector<int>                 asymfitvar_bins,
        vector<string>         massfitvars,
        vector<string>         massfitvar_titles,
        vector<vector<double>> massfitvar_lims,
        vector<int>                 massfitvar_bins,

        // parameterss passed to analysis::fitHB()
        double                           bpol,
        double                           tpol,
        double                           alpha,

        // parameters passed to saga::signal::fitMass()
        map<string,string> massfit_yamlfile_map,
        string                       massfit_pdf_name,
        string                       massfit_formula_sg,
        string                       massfit_formula_bg,
        string                       massfit_sgYield_name,
        string                       massfit_bgYield_name,
        double                            massfit_initsgfrac,
        vector<double>               massfit_parinits_sg,
        vector<string>          massfit_parnames_sg,
        vector<string>          massfit_partitles_sg,
        vector<string>          massfit_parunits_sg,
        vector<vector<double>>  massfit_parlims_sg,
        vector<double>               massfit_parinits_bg,
        vector<string>          massfit_parnames_bg,
        vector<string>          massfit_partitles_bg,
        vector<string>          massfit_parunits_bg,
        vector<vector<double>>  massfit_parlims_bg,
        vector<vector<double>>  massfit_sgregion_lims,

        // Parameters passed to analysis::applySPlots()
        bool                             use_splot,

        // Parameters used for sb subtraction
        string                      massfit_sgcut,
        string                      massfit_bgcut,
        bool                             use_sb_subtraction,
        bool                             use_binned_sb_bgfracs,
        map<int,string>        asymfitvar_bincuts,
        string                      bgfracvar,
        vector<double>              bgfracvar_lims,
        int                              bgfrac_idx               = 0,

        // Parameters passed to signal::fitMass()
        double                           massfit_lg_text_size     = 0.04,
        double                           massfit_lg_margin        = 0.1,
        int                              massfit_lg_ncols         = 1,
        bool                             massfit_plot_bg_pars     = false,
        bool                             massfit_use_sumw2error   = false,
        bool                             massfit_use_extended_nll = true,
        bool                             massfit_use_binned_fit   = false,

        // Ouput stream
        ostream                    &out                      = cout
    ) {

    string method_name = "getKinBinnedHB";
    int nparams = 1;

    // Check arguments
    if (binvars.size()<1) {
        string msg = Form("[%s]: Number of bin variables is <1", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }
    if (depolvars.size()!=1) {
        string msg = Form("[%s]: Number of depolarization variables must be 1", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }
    if (asymfitvars.size()!=1) {
        string msg = Form("[%s]: Number of asymmetry fit variables must be 1", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }
    if ((use_sb_subtraction && use_binned_sb_bgfracs) || (use_sb_subtraction && use_splot) || (use_binned_sb_bgfracs && use_splot)) {
        string msg = Form("[%s]: Sideband subtraction, sideband subtraction with binned background fractions, and the sPlot method are all mutually exclusive.", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }
    if (use_binned_sb_bgfracs) {
        string msg = Form("[%s]: Sideband subtraction with binned background fractions is not implemented with Helicity Balance method.", method_name.c_str());
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

    // Filter frames for signal and sideband
    LOG_DEBUG(Form("[%s]: Filtering frames for signal and sideband regions...", method_name.c_str()));
    massfit_sgcut = (massfit_sgcut.size()>0) ? massfit_sgcut : saga::util::addLimitCuts("",massfitvars,massfit_sgregion_lims);
    auto frame_sg = (massfit_sgcut.size()>0) ? frame.Filter(massfit_sgcut.c_str()) : frame;
    auto frame_sb = (massfit_bgcut.size()>0) ? frame.Filter(massfit_bgcut.c_str()) : frame;

    // Set condition for single mass fit
    bool single_massfit = (massfit_pdf_name!="" && !use_binned_sb_bgfracs && (use_splot || use_sb_subtraction));

    // Open output CSV
    LOG_DEBUG(Form("[%s]: Opening output CSV file...", method_name.c_str()));
    string csvpath = Form("%s.csv",scheme_name.c_str());
    ofstream csvoutf; csvoutf.open(csvpath.c_str());
    ostream &csvout = csvoutf;
    string csv_separator = ",";
    vector<string> rawasymvars = { "bsa", "tsa", "dsa"};

    // Set CSV column headers
    // COLS: bin_id,count,{binvarmean,binvarerr},{depolvarmean,depolvarerr},{rawasym,rawasymerr},{asymfitvar,asymfitvarerr},{fitvar_info if requested}
    LOG_DEBUG(Form("[%s]: Setting CSV headers...", method_name.c_str()));
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
    for (int rr=0; rr<rawasymvars.size(); rr++) {
        csvout << rawasymvars[rr].c_str() << csv_separator.c_str();
        csvout << rawasymvars[rr].c_str() << "_err" << csv_separator.c_str();
    }
    for (int aa=0; aa<nparams; aa++) {
        csvout << Form("a%d",aa) << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
        csvout << Form("a%d",aa) << "_err";
        if (aa<nparams-1 || single_massfit || use_binned_sb_bgfracs) csvout << csv_separator.c_str();
        else csvout << endl;//NOTE: IMPORTANT!
    }

    // Optionally add background asymmetries
    if (use_binned_sb_bgfracs || use_sb_subtraction) {
        for (int aa=0; aa<nparams; aa++) {
            csvout << Form("a_bg%d",aa) << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
            csvout << Form("a_bg%d",aa) << "_err";
            if (aa<nparams-1 || single_massfit) csvout << csv_separator.c_str();
            else csvout << endl;//NOTE: IMPORTANT!
        }
    }

    // Optionally add mass fit outputs
    if (single_massfit) {

        // Add signal region integration values
        csvout << "int_sg_pdf_val" << csv_separator.c_str();
        csvout << "int_sg_pdf_err" << csv_separator.c_str();
        csvout << "int_bg_pdf_val" << csv_separator.c_str();
        csvout << "int_bg_pdf_err" << csv_separator.c_str();
        csvout << "int_model_pdf_val" << csv_separator.c_str();
        csvout << "int_model_pdf_err" << csv_separator.c_str();
        csvout << "int_ds_val" << csv_separator.c_str();
        csvout << "int_ds_err" << csv_separator.c_str();
        csvout << "eps_bg_pdf" << csv_separator.c_str();
        csvout << "eps_bg_pdf_err" << csv_separator.c_str();
        csvout << "eps_sg_pdf" << csv_separator.c_str();
        csvout << "eps_sg_pdf_err" << csv_separator.c_str();
        csvout << "eps_pdf" << csv_separator.c_str();
        csvout << "eps_pdf_err" << csv_separator.c_str();

        // Add chi2 / ndf fit values
        for (int idx=0; idx<massfitvars.size(); idx++) {
            csvout << Form("chi2ndf_1d_%s", massfitvars[idx].c_str()) << csv_separator.c_str();
        }

        // Add mass fit signal PDF parameters and errors
        for (int aa=0; aa<massfit_parinits_sg.size(); aa++) {
            csvout << massfit_parnames_sg[aa].c_str() << csv_separator.c_str();
            csvout << massfit_parnames_sg[aa].c_str() << "_err" << csv_separator.c_str();
        }

        // Add mass fit background PDF parameters and errors
        for (int aa=0; aa<massfit_parinits_bg.size(); aa++) {
            csvout << massfit_parnames_bg[aa].c_str() << csv_separator.c_str();
            csvout << massfit_parnames_bg[aa].c_str() << "_err";
            if (aa<massfit_parinits_bg.size()-1) csvout << csv_separator.c_str();
            else csvout << endl;//NOTE: IMPORTANT!
        }

    }

    // Loop bins and get data
    for (auto it = bincuts.begin(); it != bincuts.end(); it++) {

        // Get bin id and cut
        int         bin_id  = it->first;
        string bin_cut = it->second;

        // Set bin id string
        string scheme_binid = Form("scheme_%s_bin_%d",scheme_name.c_str(),bin_id);
        LOG_DEBUG(Form("[%s]: Processing bin id %d with cut: %s", method_name.c_str(), bin_id, bin_cut.c_str()));

        // Create workspace
        LOG_DEBUG(Form("[%s]: Creating workspaces...", method_name.c_str()));
        RooWorkspace *ws    = new RooWorkspace(workspace_name.c_str(),workspace_title.c_str());
        RooWorkspace *ws_sg = new RooWorkspace(Form("%s_sg",workspace_name.c_str()),Form("%s_signal",workspace_title.c_str()));
        RooWorkspace *ws_sb = new RooWorkspace(Form("%s_sb",workspace_name.c_str()),Form("%s_sideband",workspace_title.c_str())); //NOTE: Use separate signal and sideband workspaces for dataset, variable, and pdf name uniqueness.

        // Make bin cut on frame
        LOG_DEBUG(Form("[%s]: Filtering frame for bin...", method_name.c_str()));
        auto binframe = frame.Filter(bin_cut.c_str());
        auto binframe_sg = frame_sg.Filter(bin_cut.c_str());

        // Create bin dataset
        LOG_DEBUG(Form("[%s]: Creating dataset...", method_name.c_str()));
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
            asymfitvars,
            asymfitvar_titles,
            asymfitvar_lims,
            asymfitvar_bins,
            massfitvars,
            massfitvar_titles,
            massfitvar_lims,
            massfitvar_bins
        );

        // Apply a generic mass fit to the FULL bin dataset
        vector<double> massfit_result;
        if (single_massfit) {  //NOTE: A mass fit in each bin is needed for basic sideband subtraction and splots.

            // Set yaml path for mass fit parameters
            string yamlfile = massfit_yamlfile_map[scheme_binid];

            // Fit the mass spectrum
            LOG_DEBUG(Form("[%s]: Fitting mass spectrum...", method_name.c_str()));
            vector<double> massfit_result = saga::signal::fitMass(
                    ws, // RooWorkspace                    *w,
                    dataset_name, // string                      dataset_name,
                    scheme_binid, // string                      binid,
                    bin_cut, // string                      bincut,
                    binvars, // vector<string>         binvars,
                    massfitvars, // vector<string>         fitvars,
                    yamlfile, // string                      yamlfile,
                    massfit_pdf_name, // string                      massfit_pdf_name,
                    massfit_formula_sg, // string                      massfit_formula_sg,
                    massfit_formula_bg, // string                      massfit_formula_bg,
                    massfit_sgYield_name, // string                      massfit_sgYield_name,
                    massfit_bgYield_name, // string                      massfit_bgYield_name,
                    massfit_initsgfrac, // double                           massfit_initsgfrac,
                    massfit_parinits_sg, // vector<double>              massfit_parinits_sg,
                    massfit_parnames_sg, // vector<string>         massfit_parnames_sg,
                    massfit_partitles_sg, // vector<string>         massfit_partitles_sg,
                    massfit_parunits_sg, // vector<string>         massfit_parunits_sg,
                    massfit_parlims_sg, // vector<vector<double>> massfit_parlims_sg,
                    massfit_parinits_bg, // vector<double>              massfit_parinits_bg,
                    massfit_parnames_bg, // vector<string>         massfit_parnames_bg,
                    massfit_partitles_bg, // vector<string>         massfit_partitles_bg,
                    massfit_parunits_bg, // vector<string>         massfit_parunits_bg,
                    massfit_parlims_bg, // vector<vector<double>> massfit_parlims_bg,
                    massfit_sgregion_lims, // vector<vector<double>> massfit_sgregion_lims,
                    massfit_lg_text_size, // double                           massfit_lg_text_size     = 0.04,
                    massfit_lg_margin, // double                           massfit_lg_margin        = 0.1,
                    massfit_lg_ncols, // int                              massfit_lg_ncols         = 1,
                    massfit_plot_bg_pars, // bool                             massfit_plot_bg_pars     = false,
                    massfit_use_sumw2error, // bool                             massfit_use_sumw2error   = false,
                    massfit_use_extended_nll, // bool                             massfit_use_extended_nll = true,
                    massfit_use_binned_fit, // bool                             massfit_use_binned_fit   = false,
                    out // ostream                    &out              = cout
            );
        }

        // Apply sPlot
        string fit_dataset_name = dataset_name; // -> Use this for sPlot
        if (use_splot) {
            LOG_DEBUG(Form("[%s]: Applying sPlot to dataset...", method_name.c_str()));
            string dataset_sg_name = (string)Form("%s_sg_sw",dataset_name.c_str());
            string dataset_bg_name = (string)Form("%s_bg_sw",dataset_name.c_str());
            saga::signal::applySPlot(
                ws,
                dataset_name,
                Form("%s_%s",massfit_sgYield_name.c_str(),scheme_binid.c_str()),//NOTE: getGenAsymPdf() renames these variables to ensure workspace uniqueness
                Form("%s_%s",massfit_bgYield_name.c_str(),scheme_binid.c_str()),
                Form("%s_%s",massfit_pdf_name.c_str(),scheme_binid.c_str()),
                dataset_sg_name,
                dataset_bg_name
            );
            fit_dataset_name = dataset_sg_name;
        }

        // Weight dataset from binned mass fits
        string fit_sb_dataset_name = ""; // -> Use this for binned sideband backgrounds
        if (use_binned_sb_bgfracs) {
            string rds_out_name = (string)Form("%s_sg",dataset_name.c_str());
            string sb_rds_out_name = (string)Form("%s_sb",dataset_name.c_str());
            LOG_DEBUG(Form("[%s]: Setting binned background fractions for dataset %s...", method_name.c_str(), dataset_name.c_str()));
            saga::signal::setBinnedBGFractions(
                ws, // RooWorkspace                    *w,
                dataset_name, // string                      dataset_name,
                scheme_binid, // string                      binid,
                bin_cut, // string                      bincut,
                binvars, // vector<string>         binvars,
                massfitvars, // vector<string>         fitvars,
                massfit_yamlfile_map, // map<string,string> yamlfile_map
                massfit_pdf_name, // string                      massfit_pdf_name,
                massfit_formula_sg, // string                      massfit_formula_sg,
                massfit_formula_bg, // string                      massfit_formula_bg,
                massfit_sgYield_name, // string                      massfit_sgYield_name,
                massfit_bgYield_name, // string                      massfit_bgYield_name,
                massfit_initsgfrac, // double                           massfit_initsgfrac,
                massfit_parinits_sg, // vector<double>              massfit_parinits_sg,
                massfit_parnames_sg, // vector<string>         massfit_parnames_sg,
                massfit_partitles_sg, // vector<string>         massfit_partitles_sg,
                massfit_parunits_sg, // vector<string>         massfit_parunits_sg,
                massfit_parlims_sg, // vector<vector<double>> massfit_parlims_sg,
                massfit_parinits_bg, // vector<double>              massfit_parinits_bg,
                massfit_parnames_bg, // vector<string>         massfit_parnames_bg,
                massfit_partitles_bg, // vector<string>         massfit_partitles_bg,
                massfit_parunits_bg, // vector<string>         massfit_parunits_bg,
                massfit_parlims_bg, // vector<vector<double>> massfit_parlims_bg,
                massfit_sgregion_lims, // vector<vector<double>> massfit_sgregion_lims,

                binframe, // RNode                            frame, // arguments for this method
                massfit_bgcut, // string                      bgcut, 
                asymfitvars, // vector<string>         asymfitvars,
                asymfitvar_bincuts, // map<int,string>        asymfitvar_bincuts,
                rds_out_name, // string                      rds_out_name,
                sb_rds_out_name, // string                   sb_rds_out_name,
                bgfracvar, // string                      bgfracvar,
                bgfracvar_lims, // vector<double>              bgfracvar_lims,

                massfit_lg_text_size, // double                           massfit_lg_text_size     = 0.04,
                massfit_lg_margin, // double                           massfit_lg_margin        = 0.1,
                massfit_lg_ncols, // int                              massfit_lg_ncols         = 1,
                massfit_plot_bg_pars, // bool                             massfit_plot_bg_pars     = false,
                massfit_use_sumw2error, // bool                             massfit_use_sumw2error   = false,
                massfit_use_extended_nll, // bool                             massfit_use_extended_nll = true,
                massfit_use_binned_fit, // bool                             massfit_use_binned_fit   = false,

                bgfrac_idx, // int                               bgfrac_idx               = 0,
                0.0, // double                           bgfracs_default  = 0.0 // arguments for this method
                out // ostream                    &out              = cout
            );
            fit_dataset_name = rds_out_name;
            fit_sb_dataset_name = sb_rds_out_name;
        }

        // Create signal region dataset for sideband subtraction
        if (use_sb_subtraction) {
            LOG_DEBUG(Form("[%s]: Creating sideband dataset...", method_name.c_str()));
            data::createDataset(
                binframe_sg,
                ws_sg, //NOTE: Use separate sideband workspace for dataset, variable, and pdf name uniqueness.
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
        LOG_DEBUG(Form("[%s]: Fitting asymmetry...", method_name.c_str()));
        vector<double> asymfit_result = fitHB(
                                (use_sb_subtraction ? ws_sg : ws),
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                bpol,
                                tpol,
                                alpha,
                                helicity,
                                tspin,
                                htspin,
                                combined_spin_state,
                                scheme_binid,
                                bin_cut,
                                binvars,
                                depolvars,
                                asymfitvars,
                                out
                            );

        // Compute sideband region bin results
        string sb_scheme_binid = Form("sb_%s",scheme_binid.c_str());
        vector<double> asymfit_result_sb;
        if (use_sb_subtraction) {

            // Make bin cut on sideband frame
            auto binframe_sb = frame_sb.Filter(bin_cut.c_str());

            // Create sideband dataset
            LOG_DEBUG(Form("[%s]: Creating sideband dataset...", method_name.c_str()));
            data::createDataset(
                binframe_sb,
                ws_sb, //NOTE: Use separate sideband workspace for dataset, variable, and pdf name uniqueness.
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
            LOG_DEBUG(Form("[%s]: Fitting sideband asymmetry...", method_name.c_str()));
            asymfit_result_sb = fitHB(
                                ws_sb,
                                dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                bpol,
                                tpol,
                                alpha,
                                helicity,
                                tspin,
                                htspin,
                                combined_spin_state,
                                sb_scheme_binid,
                                bin_cut,
                                binvars,
                                depolvars,
                                asymfitvars,
                                out
                            );
        }

        // Initialize data
        LOG_DEBUG(Form("[%s]: Extracting results...", method_name.c_str()));
        int nbinvars = binvars.size();
        int nparams  = nparams;
        double xs[nbinvars];
        double exs[nbinvars];
        int    count;

        double ys[nparams];
        double eys[nparams];
        double ys_sb[nparams];
        double eys_sb[nparams];
        double depols[nparams];
        double edepols[nparams];
        double rawasyms[(const int)rawasymvars.size()];
        double rawasymerrs[(const int)rawasymvars.size()];

        // Get asymmetry fit bin data
        int k = 0;
        count = (int)asymfit_result[k++];
        for (int idx=0; idx<binvars.size(); idx++) {
            xs[idx]     = asymfit_result[k++];
            exs[idx]    = asymfit_result[k++];
        }
        for (int idx=0; idx<depolvars.size(); idx++) {
            depols[idx] = asymfit_result[k++];
            edepols[idx] = asymfit_result[k++];
        }
        for (int idx=0; idx<rawasymvars.size(); idx++) {
            rawasyms[idx] = asymfit_result[k++];
            rawasymerrs[idx] = asymfit_result[k++];
        }
        for (int idx=0; idx<nparams; idx++) {
            ys[idx] = asymfit_result[k++];
            eys[idx] = asymfit_result[k++];
        }
        if (use_binned_sb_bgfracs) {
            for (int idx=0; idx<nparams; idx++) {
                ys_sb[idx] = asymfit_result[k++];
                eys_sb[idx] = asymfit_result[k++];
            }
        }


        // Get mass fit bin data
        double int_sg_pdf_val;
        double int_sg_pdf_err;
        double int_bg_pdf_val;
        double int_bg_pdf_err;
        double int_model_pdf_val;
        double int_model_pdf_err;
        double int_ds_val;
        double int_ds_err;
        double eps_bg_pdf;
        double eps_bg_pdf_err;
        double eps_sg_pdf;
        double eps_sg_pdf_err;
        double eps_pdf;
        double eps_pdf_err;
        vector<double> chi2ndfs;
        vector<double> massfit_pars_sg;
        vector<double> massfit_parerrs_sg;
        vector<double> massfit_pars_bg;
        vector<double> massfit_parerrs_bg;
        if (massfit_result.size()>0) {

            // Start counter
            int m = 1; //NOTE: Ignore count which should first entry

            // Add signal region integration values
            int_sg_pdf_val    = massfit_result[m++];
            int_sg_pdf_err    = massfit_result[m++];
            int_bg_pdf_val    = massfit_result[m++];
            int_bg_pdf_err    = massfit_result[m++];
            int_model_pdf_val = massfit_result[m++];
            int_model_pdf_err = massfit_result[m++];
            int_ds_val        = massfit_result[m++];
            int_ds_err        = massfit_result[m++];

            // Add background fractions
            eps_bg_pdf     = massfit_result[m++];
            eps_bg_pdf_err = massfit_result[m++];
            eps_sg_pdf     = massfit_result[m++];
            eps_sg_pdf_err = massfit_result[m++];
            eps_pdf        = massfit_result[m++];
            eps_pdf_err    = massfit_result[m++];

            // Add chi2/ndfs
            for (int idx=0; idx<massfitvars.size(); idx++) {
                chi2ndfs.push_back(massfit_result[m++]);
            }

            // Skip bin variables
            m += binvars.size();

            // Add signal PDF parameters and errors
            for (int idx=0; idx<massfit_parinits_sg.size(); idx++) {
                massfit_pars_sg.push_back(massfit_result[m++]);
                massfit_parerrs_sg.push_back(massfit_result[m++]);
            }

            // Add background PDF parameters and errors
            for (int idx=0; idx<massfit_parinits_bg.size(); idx++) {
                massfit_pars_bg.push_back(massfit_result[m++]);
                massfit_parerrs_bg.push_back(massfit_result[m++]);
            }
        }

        // Apply sideband subtraction to asymmetries
        double epsilon, epsilon_err;
        if (use_sb_subtraction) {
            LOG_DEBUG(Form("[%s]: Applying sideband subtraction...", method_name.c_str()));
            int k2 = 1 + binvars.size() + depolvars.size();
            epsilon = eps_bg_pdf;
            epsilon_err = eps_bg_pdf_err;
            for (int idx=0; idx<nparams; idx++) {
                ys_sb[idx] = asymfit_result_sb[k2++];
                eys_sb[idx] = asymfit_result_sb[k2++];
                ys[idx]  = (ys[idx] - epsilon * ys_sb[idx]) / (1.0 - epsilon);
                eys[idx] = TMath::Sqrt(eys[idx]*eys[idx] + epsilon * epsilon * eys_sb[idx]*eys_sb[idx]) / (1.0 - epsilon);
            }
        }

        // Output message
        out << "--- Acceptance, background, and depolarization corrected signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            if (use_sb_subtraction || use_binned_sb_bgfracs) {
                out << " ys_sb["<< idx <<"]       = " << ys_sb[idx] << "\n";
                out << " eys_sb["<< idx <<"]      = " << eys_sb[idx] << "\n";
            }
            out << " ys["<< idx <<"]   = " << ys[idx] << "\n";
            out << " eys["<< idx <<"]  = " << eys[idx] << "\n";
        }
        out << "---------------------------\n";

        // Write out a row of data to csv
        // COLS: bin_id,count,{binvarmean,binvarerr},{depolvarmean,depolvarerr},{rawasym,rawasymerr},{asymfitvar,asymfitvarerr}(,{bg_asymfitvar,bg_asymfitvarerr})
        LOG_DEBUG(Form("[%s]: Writing results to CSV...", method_name.c_str()));
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
        for (int rr=0; rr<rawasymvars.size(); rr++) {
            csvout << rawasyms[rr] << csv_separator.c_str();
            csvout << rawasymerrs[rr] << csv_separator.c_str();
        }
        for (int aa=0; aa<nparams; aa++) {
            csvout << ys[aa] << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
            csvout << eys[aa];
            if (aa<nparams-1 || single_massfit || use_binned_sb_bgfracs) csvout << csv_separator.c_str();
            else csvout << endl;//NOTE: IMPORTANT!
        }

        // Optionally add background asymmetries
        if (use_binned_sb_bgfracs || use_sb_subtraction) {
            for (int aa=0; aa<nparams; aa++) {
                csvout << ys_sb[aa] << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
                csvout << eys_sb[aa];
                if (aa<nparams-1 || single_massfit) csvout << csv_separator.c_str();
                else csvout << endl;//NOTE: IMPORTANT!
            }
        }

        // Optionally add mass fit outputs
        // COLS: {integration values and errors in signal region},{background fraction values and errors},{chi2/ndf},{signal PDF parameters and errors},{background PDF parameters and errors}
        if (single_massfit) {

            // Add signal region integration values
            csvout << int_sg_pdf_val << csv_separator.c_str();
            csvout << int_sg_pdf_err << csv_separator.c_str();
            csvout << int_bg_pdf_val << csv_separator.c_str();
            csvout << int_bg_pdf_err << csv_separator.c_str();
            csvout << int_model_pdf_val << csv_separator.c_str();
            csvout << int_model_pdf_err << csv_separator.c_str();
            csvout << int_ds_val << csv_separator.c_str();
            csvout << int_ds_err << csv_separator.c_str();
            csvout << eps_bg_pdf << csv_separator.c_str();
            csvout << eps_bg_pdf_err << csv_separator.c_str();
            csvout << eps_sg_pdf << csv_separator.c_str();
            csvout << eps_sg_pdf_err << csv_separator.c_str();
            csvout << eps_pdf << csv_separator.c_str();
            csvout << eps_pdf_err << csv_separator.c_str();

            // Add chi2 / ndf fit values
            for (int idx=0; idx<chi2ndfs.size(); idx++) {
                csvout << chi2ndfs[idx] << csv_separator.c_str();
            }

            // Add mass fit signal PDF parameters and errors
            for (int aa=0; aa<massfit_pars_sg.size(); aa++) {
                csvout << massfit_pars_sg[aa] << csv_separator.c_str();
                csvout << massfit_parerrs_sg[aa] << csv_separator.c_str();
            }

            // Add mass fit background PDF parameters and errors
            for (int aa=0; aa<massfit_pars_bg.size(); aa++) {
                csvout << massfit_pars_sg[aa] << csv_separator.c_str();
                csvout << massfit_parerrs_bg[aa];
                if (aa<massfit_pars_bg.size()-1) csvout << csv_separator.c_str();
                else csvout << endl;//NOTE: IMPORTANT!
            }

        }
    }

    csvoutf.close();
    out << " Saved asymmetry fit results to " << csvpath.c_str() << endl;

    // Ending message
    out << "------------------- END of " << method_name.c_str() << " -------------------\n";

} // getKinBinnedHB()

} // namespace hbanalysis {

} // namespace saga {
