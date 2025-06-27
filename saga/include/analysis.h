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
#include <RooPlot.h>
// #include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
// #include <RooAddPdf.h>
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

// // RooStats includes
// #include <RooStats/SPlot.h>

// Local includes
#include <data.h>
#include <bins.h>
#include <util.h>
#include <signal.h>

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 12/Dec./2024
* @version 0.0.0
* @brief Fit asymmetries using RooFit unbinned Maximum Likelihood methods and sideband subtraction 
* or the sPlot method from <a href="http://arxiv.org/abs/physics/0402083">arXiv:physics/0402083</a> for background correction.
*/

namespace saga {

namespace analysis {

using namespace RooFit;
using RNode = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;

/**
* @brief Create a subset of a RooArgSet for a given pdf formula
*
* @param argset RooArgSet containing RooRealVar* variables
* @param fitformula PDF formula passed to RooGenericPdf
* @param varformulas List of variable formulas passed to RooGenericPdf for all arguments in the RooArgSet
* @param varnames List of variable names in RooArgSet directly corresponding to the list of variable formulas
*
* @return RooArgSet*
*/
RooArgSet* getSubRooArgSet(
    RooArgSet* argset,
    std::string fitformula,
    std::vector<std::string> varformulas,
    std::vector<std::string> varnames
) {

    // Isolate the argset for the target spin dependent terms
    RooArgSet *subargset = new RooArgSet();

    // Loop variable formulas and check if they are in the fit formula
    for (int idx = 0; idx<varformulas.size(); idx++) {

        // Check if fit formula contains variable formula
        if (fitformula.find(varformulas[idx]) != std::string::npos) {

            // Find the RooRealVar in the argset
            RooRealVar *var = (RooRealVar*)argset->find(varnames[idx].c_str());
            if (var==nullptr) {
                std::cerr << "ERROR: RooRealVar \"" << varnames[idx].c_str() << "\" with formula \"" << varformulas[idx].c_str() << "\" not found in argset" << std::endl;
                continue;
            }
            subargset->add(*var);
        }
    }

    return subargset;

} // RooArgSet* getSubRooArgSet()


/**
* @brief Adjust the parameter indexing of a pdf formula for a given subset of parameters
*
* Note that the output formula will have the variable notation `x[<idx>]` where `idx`
* indicates the integer index of the variable in the appropriate RooArgSet.
*
* @param fitformula PDF formula passed to RooGenericPdf
* @param varformulas List of variable formulas passed to RooGenericPdf
* @param max_idx If this parameter is \f$>0\f$, then the formulas will be substituted in descending order starting at `idx==max_idx`.
*
* @return std::string
*/
std::string getSubFormula(
    std::string fitformula,
    std::vector<std::string> varformulas,
    int max_idx = 0
) {

    // Loop variable formulas and map old formulas to new formulas
    std::string subfitformula = fitformula;

    // If you are doing forward ordering
    if (max_idx==0) {
        int add_idx = max_idx;
        for (int idx = 0; idx<varformulas.size(); idx++) {

            // Check if fit formula contains variable formula
            if (fitformula.find(varformulas[idx]) != std::string::npos) {

                // Replace variable and increment index of variables added
                saga::util::replaceAll(subfitformula, varformulas[idx], Form("x[%d]", add_idx));
                add_idx++;
            }
        }
    } else { // Or backward ordering
        int add_idx = max_idx;
        for (int idx = varformulas.size()-1; idx>=0; idx--) {

            // Check if fit formula contains variable formula
            if (fitformula.find(varformulas[idx]) != std::string::npos) {

                // Replace variable and increment index of variables added
                saga::util::replaceAll(subfitformula, varformulas[idx], Form("x[%d]", add_idx));
                add_idx--;
            }
        }
    }

    return subfitformula;

} // RooArgSet* getSubRooArgSet()

/**
* @brief Create a PDF for fitting a generic asymmetry with a maximum likelihood fit.
* 
* Create a PDF given the formulas for the asymmetries coupling to each combination of beam helicity
* and target spin states.  The PDF will be constructed internally using <a href="https://root.cern.ch/doc/master/classRooGenericPdf.html">RooGenericPdf</a>
* in the form:
*
* @f[
* \begin{aligned}
* PDF(\lambda_{\ell}, S_{||}, x_0, x_1, ..., &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + A_{UU,UT}(\vec{x}, \vec{a}, \vec{d}) \\
* &   + \lambda_{\ell} \cdot \overline{\lambda_{\ell}^2} \cdot A_{LU,LT}(\vec{x}, \vec{a}, \vec{d}) \\
* &   + S_{||} \cdot \overline{S^2} \cdot A_{UL}(\vec{x}, \vec{a}, \vec{d}) \\
* &   + \lambda_{\ell} \cdot \overline{\lambda_{\ell}^2} \cdot S_{||} \cdot \overline{S^2} \cdot A_{LL}(\vec{x}, \vec{a}, \vec{d}), \\
* \end{aligned}
* @f]
* where the appropriate terms will be dropped if there is no dependence on beam helicity or target spin.
* The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the depolarization factors.
*
* If `categories_as_float` contains the beam helicity or target spin variable names,
* the PDF will use these as independent variables.  Otherwise, a simultaneous PDF will be
* formed over the various helicity and spin states.
*
* Note that in the case of an \f$A_{UT}\f$ or \f$A_{LT}\f$ asymmetry, the relevant formula should be included in the argument for the \f$A_{UU}\f$ or \f$A_{LU}\f$ formula
* respectively since \f$A_{UT}\f$ and \f$A_{LT}\f$ should only have kinematic dependence on \f$\phi_{S}\f$ rather than categorical dependence on \f$S_{\perp}\f$.
*
* The variable names in the fit formulas should follow the <a href="https://root.cern.ch/doc/master/classTFormula.html">TFormula</a> notation, e.g.,
* `x_0`\f$\rightarrow\f$`x[0]`, `x_1`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[N_x]`, `a_1`\f$\rightarrow\f$`x[N_x+1]`, etc.
* 
* Note that in the case that all three fit formula terms are used, the formulas and corresponding
* argument sets for the PDFs that depend only on either beam helicity or target spin will be reduced
* to the appropriate subset of variables used in the corresponding fit formula.  This ensures
* that the PDFs will still compile correctly even when uploaded to the RooWorkspace.
*
* When using a single fit formula the indexing of parameters in the fit formula is straightforward.
* However, when using all three fit formulas, the indexing is global across all three formulas.
*
* The returned list contains:
*
* - The name of the full model in the workspace
*
* - The name of each yield variable in the workspace in the case of an extended NLL fit
* 
* @param w RooWorkspace in which to work
* @param categories_as_float List of category variables to include use as asymmetry fit variables and automatically add to PDF formula
* @param h Beam helicity \f$\lambda_{\ell}\in(-1,0,1)\f$
* @param t Target spin \f$S\in(-1,0,1)\f$
* @param ht Beam helicity times target spin \f$\lambda_{\ell} \cdot S\in(-1,0,1)\f$
* @param ss Combined beam helicity and target spin state \f$ss = (\lambda_{\ell}+1)\cdot10 + (S+1)\f$
* @param argset Argument set for PDF
* @param argnames Argument names for PDF
* @param method_name Method name, used to name PDF
* @param binid Unique bin id, used to name PDF
* @param fitformula_uu Fit formula for the asymmetry terms \f$A_{UU,UT}\f$
* @param fitformula_pu Fit formula for the beam helicity dependent asymmetry terms \f$A_{LU,LT}\f$
* @param fitformula_up Fit formula for the target spin dependent asymmetry terms \f$A_{UL}\f$
* @param fitformula_pp Fit formula for the beam helicity and target spin dependent asymmetry terms \f$A_{LL}\f$
* @param bpol Luminosity averaged beam polarization \f$\overline{\lambda_{\ell}^2}\f$
* @param tpol Luminosity averaged target polarization \f$\overline{S^2}\f$
* @param count Bin count
* @param use_extended_nll Option to use an extended likelihood term
* 
* @return List of model name and all yield variable names
*/
std::vector<std::string> getGenAsymPdf(
    RooWorkspace *w,
    std::vector<std::string> categories_as_float,
    RooCategory *h,
    RooCategory *t,
    RooCategory *ht,
    RooCategory *ss,
    RooArgSet *argset,
    std::vector<std::string> argnames,
    std::string method_name,
    std::string binid,
    std::string fitformula_uu,
    std::string fitformula_pu,
    std::string fitformula_up,
    std::string fitformula_pp,
    double bpol,
    double tpol,
    int count,
    bool use_extended_nll
) {

    // Get the total number of states and set the starting count for the extended case
    int nstates = 0;
    if (fitformula_pu!="") nstates += 1;
    if (fitformula_up!="") nstates += 1;
    if (fitformula_pp!="") nstates += 1;
    nstates *= 3;
    if (fitformula_uu!="") nstates += 1; //NOTE: There is only one unpolarized PDF so add AFTER multiplying by three for the spin dependent PDFs.
    double ninit = count / nstates;

    // Set variable formulas list
    std::vector<std::string> varformulas;
    for (int idx=0; idx<argset->size(); idx++) {
        varformulas.push_back(Form("x[%d]", (int)(idx-categories_as_float.size())));
    }
    //NOTE: The user does not need to specify the helicity and target spin variables in the fit formula,
    // but they are prepended to the asymmetry fit variables (helicity, then tspin).

    // Set model and yield names
    std::string model_name = Form("model_%s_%s",method_name.c_str(),binid.c_str()); //TODO: Make model names more specific above to avoid naming conflicts...
    std::vector<std::string> model_and_yield_names;

    // Create simple pdf here if not using simultaneous PDF
    if (categories_as_float.size()>0) {

        // Check whether you have helicity and/or tspin and set formulas
        std::string helicity_formula = (categories_as_float.size()==2) ? "x[-2]" : "x[-1]";
        std::string tspin_formula = "x[-1]";

        // Create the PDF formula
        std::string fitformula_full = "";

        // Create the PDF formula
        if (fitformula_uu!="") {
            fitformula_full = fitformula_uu;
        }
        if (fitformula_pu!="") {
            std::string fitformula_new = Form("%s*%.3f*(%s)",helicity_formula.c_str(),bpol,fitformula_pu.c_str());
            if (fitformula_full=="") fitformula_full = fitformula_new;
            else fitformula_full = Form("%s + %s",fitformula_full.c_str(),fitformula_new.c_str());
        }
        if (fitformula_up!="") {
            std::string fitformula_new = Form("%s*%.3f*(%s)",tspin_formula.c_str(),tpol,fitformula_up.c_str());
            if (fitformula_full=="") fitformula_full = fitformula_new;
            else fitformula_full = Form("%s + %s",fitformula_full.c_str(),fitformula_new.c_str());
        }
        if (fitformula_pp!="") {
            std::string fitformula_new = Form("%s*%s*%.3f*%.3f*(%s)",helicity_formula.c_str(),tspin_formula.c_str(),bpol,tpol,fitformula_pp.c_str());
            if (fitformula_full=="") fitformula_full = fitformula_new;
            else fitformula_full = Form("%s + %s",fitformula_full.c_str(),fitformula_new.c_str());
        }

        // Isolate the argset for the target spin dependent terms
        RooArgSet *argset_full = getSubRooArgSet(argset, fitformula_full.c_str(), varformulas, argnames);
        std::string subfitformula_full = getSubFormula(fitformula_full.c_str(), varformulas, argset->size()-1);

        // Create PDF
        fitformula_full = fitformula_full!="" ? Form("1.0+(%s)",subfitformula_full.c_str()): "1.0";
        RooGenericPdf _model_full(Form("_%s_full",model_name.c_str()), fitformula_full.c_str(), *argset_full);

        // Create extended PDF
        RooRealVar nsig_full(Form("nsig_%s_full",model_name.c_str()), "number of signal events", count, 0.0, 2.0*count);
        RooExtendPdf model_full(Form("%s_full",model_name.c_str()), "extended signal pdf", _model_full, nsig_full);

        // Import the PDF and return the model and yield names
        if (use_extended_nll) {
            w->import(model_full);
            model_and_yield_names.push_back(model_full.GetName());
            model_and_yield_names.push_back(nsig_full.GetName());
        }
        else {
            model_and_yield_names.push_back(_model_full.GetName());
            w->import(_model_full);
        }
        return model_and_yield_names;
    }

    // Create simultaneous pdf and name
    RooSimultaneous * model;
    model_and_yield_names.push_back(model_name);

    //NOTE: `_<int><int>` on the variables below correspond to beam helicity and target spin states (+1) respectively, i.e., (-1,0,1) -> (0,1,2).

    // Dummy formula so that unused generic pdfs still compile
    std::string fitformula_unused = "0.0";

    //----- Unpolarized and transverse target spin dependent terms -----//

    // Isolate the argset for the target spin dependent terms
    RooArgSet *argset_uu = getSubRooArgSet(argset, fitformula_uu!="" ? fitformula_uu.c_str() : fitformula_unused.c_str(), varformulas, argnames);
    std::string subfitformula_uu = getSubFormula(fitformula_uu!="" ? fitformula_uu.c_str() : fitformula_unused.c_str(), varformulas);

    // Create pdf helicity==0
    std::string fitformula_11 = fitformula_uu!="" ? Form("1.0+(%s)",subfitformula_uu.c_str()): "1.0";
    RooGenericPdf _model_11(Form("_%s_11",model_name.c_str()), fitformula_11.c_str(), *argset_uu);

    // Create extended pdf helicity==0
    RooRealVar nsig_11(Form("nsig_%s_11",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_11(Form("%s_11",model_name.c_str()), "extended signal pdf", _model_11, nsig_11);

    //----- Create UU/UT starting fit formula -----//
    std::string fitformula_uu_ut = fitformula_uu!="" ? Form("1.0+(%s)",fitformula_uu.c_str()): "1.0";

    //----- Beam helicity and target spin dependent terms for fit of these terms ONLY -----//

    // Set the fit formulas
    std::string fitformula_22_00 = Form("%s+%.3f*(%s)",fitformula_uu_ut.c_str(),bpol*tpol,fitformula_pp!="" ? fitformula_pp.c_str() : fitformula_unused.c_str());
    std::string fitformula_20_02 = Form("%s-%.3f*(%s)",fitformula_uu_ut.c_str(),bpol*tpol,fitformula_pp!="" ? fitformula_pp.c_str() : fitformula_unused.c_str());

    // Isolate the argset for the target spin dependent terms
    RooArgSet *argset_pp = getSubRooArgSet(argset, fitformula_22_00, varformulas, argnames);
    std::string subfitformula_22_00 = getSubFormula(fitformula_22_00, varformulas);
    std::string subfitformula_20_02 = getSubFormula(fitformula_20_02, varformulas);

    // Create pdf htspin==+1
    RooGenericPdf _model_22_00(Form("_%s_22_00",model_name.c_str()), subfitformula_22_00.c_str(), *argset_pp);

    // Create extended pdf htspin==+1
    RooRealVar nsig_22_00(Form("nsig_%s_22_00",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_22_00(Form("%s_22_00",model_name.c_str()), "extended signal pdf", _model_22_00, nsig_22_00);

    // Create pdf htspin==-1
    RooGenericPdf _model_20_02(Form("_%s_20_02",model_name.c_str()), subfitformula_20_02.c_str(), *argset_pp);

    // Create extended pdf htspin==-1
    RooRealVar nsig_20_02(Form("nsig_%s_20_02",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_20_02(Form("%s_20_02",model_name.c_str()), "extended signal pdf", _model_20_02, nsig_20_02);

    //----- Beam helicity dependent terms -----//

    // Set the fit formulas
    std::string fitformula_21 = Form("%s+%.3f*(%s)",fitformula_uu_ut.c_str(),bpol,fitformula_pu!="" ? fitformula_pu.c_str() : fitformula_unused.c_str());
    std::string fitformula_01 = Form("%s-%.3f*(%s)",fitformula_uu_ut.c_str(),bpol,fitformula_pu!="" ? fitformula_pu.c_str() : fitformula_unused.c_str());

    // Isolate the argset for the target spin dependent terms
    RooArgSet *argset_pu = getSubRooArgSet(argset, fitformula_21, varformulas, argnames);
    std::string subfitformula_21 = getSubFormula(fitformula_21, varformulas);
    std::string subfitformula_01 = getSubFormula(fitformula_01, varformulas);

    // Create pdf (h,t,ht) -> ( 1, 0, 0)
    RooGenericPdf _model_21(Form("_%s_21",model_name.c_str()), subfitformula_21.c_str(), *argset_pu);

    // Create extended pdf (h,t,ht) -> ( 1, 0, 0)
    RooRealVar nsig_21(Form("nsig_%s_21",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_21(Form("%s_21",model_name.c_str()), "extended signal pdf", _model_21, nsig_21);

    // Create pdf (h,t,ht) -> (-1, 0, 0)
    RooGenericPdf _model_01(Form("_%s_01",model_name.c_str()), subfitformula_01.c_str(), *argset_pu);

    // Create extended pdf (h,t,ht) -> (-1, 0, 0)
    RooRealVar nsig_01(Form("nsig_%s_01",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_01(Form("%s_01",model_name.c_str()), "extended signal pdf", _model_01, nsig_01);

    //----- Target spin dependent terms -----//

    // Set the fit formulas
    std::string fitformula_12 = Form("%s+%.3f*(%s)",fitformula_uu_ut.c_str(),tpol,fitformula_up!="" ? fitformula_up.c_str() : fitformula_unused.c_str());
    std::string fitformula_10 = Form("%s-%.3f*(%s)",fitformula_uu_ut.c_str(),tpol,fitformula_up!="" ? fitformula_up.c_str() : fitformula_unused.c_str());

    // Isolate the argset for the target spin dependent terms
    RooArgSet *argset_up = getSubRooArgSet(argset, fitformula_12, varformulas, argnames);
    std::string subfitformula_12 = getSubFormula(fitformula_12, varformulas);
    std::string subfitformula_10 = getSubFormula(fitformula_10, varformulas);

    // Create pdf (h,t,ht) -> ( 0, 1, 0)
    RooGenericPdf _model_12(Form("_%s_12",model_name.c_str()), subfitformula_12.c_str(), *argset_up);

    // Create extended pdf (h,t,ht) -> ( 0, 1, 0)
    RooRealVar nsig_12(Form("nsig_%s_12",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_12(Form("%s_12",model_name.c_str()), "extended signal pdf", _model_12, nsig_12);

    // Create pdf (h,t,ht) -> ( 0,-1, 0)
    RooGenericPdf _model_10(Form("_%s_10",model_name.c_str()), subfitformula_10.c_str(), *argset_up);

    // Create extended pdf (h,t,ht) -> ( 0,-1, 0)
    RooRealVar nsig_10(Form("nsig_%s_10",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_10(Form("%s_10",model_name.c_str()), "extended signal pdf", _model_10, nsig_10);

    //----- Beam helicity and target spin dependent terms -----//
    // Create pdf (h,t,ht) -> ( 1, 1, 1)
    std::string fitformula_22 = Form("%s+%.3f*(%s)+%.3f*(%s)+%.3f*(%s)",fitformula_uu_ut.c_str(),bpol,fitformula_pu!="" ? fitformula_pu.c_str() : fitformula_unused.c_str(),tpol,fitformula_up!="" ? fitformula_up.c_str() : fitformula_unused.c_str(),bpol*tpol,fitformula_pp!="" ? fitformula_pp.c_str() : fitformula_unused.c_str());
    RooGenericPdf _model_22(Form("_%s_22",model_name.c_str()), fitformula_22.c_str(), *argset);

    // Create extended pdf (h,t,ht) -> ( 1, 1, 1)
    RooRealVar nsig_22(Form("nsig_%s_22",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_22(Form("%s_22",model_name.c_str()), "extended signal pdf", _model_22, nsig_22);

    // Create pdf (h,t,ht) -> (-1,-1, 1)
    std::string fitformula_00 = Form("%s-%.3f*(%s)-%.3f*(%s)+%.3f*(%s)",fitformula_uu_ut.c_str(),bpol,fitformula_pu!="" ? fitformula_pu.c_str() : fitformula_unused.c_str(),tpol,fitformula_up!="" ? fitformula_up.c_str() : fitformula_unused.c_str(),bpol*tpol,fitformula_pp!="" ? fitformula_pp.c_str() : fitformula_unused.c_str());
    RooGenericPdf _model_00(Form("_%s_00",model_name.c_str()), fitformula_00.c_str(), *argset);

    // Create extended pdf (h,t,ht) -> (-1,-1, 1)
    RooRealVar nsig_00(Form("nsig_%s_00",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_00(Form("%s_00",model_name.c_str()), "extended signal pdf", _model_00, nsig_00);

    // Create pdf (h,t,ht) -> (-1, 1,-1)
    std::string fitformula_02 = Form("%s-%.3f*(%s)+%.3f*(%s)-%.3f*(%s)",fitformula_uu_ut.c_str(),bpol,fitformula_pu!="" ? fitformula_pu.c_str() : fitformula_unused.c_str(),tpol,fitformula_up!="" ? fitformula_up.c_str() : fitformula_unused.c_str(),bpol*tpol,fitformula_pp!="" ? fitformula_pp.c_str() : fitformula_unused.c_str());
    RooGenericPdf _model_02(Form("_%s_02",model_name.c_str()), fitformula_02.c_str(), *argset);

    // Create extended pdf (h,t,ht) -> (-1, 1,-1)
    RooRealVar nsig_02(Form("nsig_%s_02",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_02(Form("%s_02",model_name.c_str()), "extended signal pdf", _model_02, nsig_02);

    // Create pdf (h,t,ht) -> ( 1,-1,-1)
    std::string fitformula_20 = Form("%s+%.3f*(%s)-%.3f*(%s)-%.3f*(%s)",fitformula_uu_ut.c_str(),bpol,fitformula_pu!="" ? fitformula_pu.c_str() : fitformula_unused.c_str(),tpol,fitformula_up!="" ? fitformula_up.c_str() : fitformula_unused.c_str(),bpol*tpol,fitformula_pp!="" ? fitformula_pp.c_str() : fitformula_unused.c_str());
    RooGenericPdf _model_20(Form("_%s_20",model_name.c_str()), fitformula_20.c_str(), *argset);

    // Create extended pdf (h,t,ht) -> ( 1,-1,-1)
    RooRealVar nsig_20(Form("nsig_%s_20",model_name.c_str()), "number of signal events", ninit, 0.0, count);
    RooExtendPdf model_20(Form("%s_20",model_name.c_str()), "extended signal pdf", _model_20, nsig_20);

    // Construct helicity dependent pdf
    if (fitformula_pu!="" && fitformula_up=="" && fitformula_pp=="") {

        // Create the simultaneous pdf
        if (use_extended_nll) {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {h->lookupName(1), &model_21}, {h->lookupName(0), &model_11}, {h->lookupName(-1), &model_01}
            },
            *h);
            model_and_yield_names.push_back(nsig_21.GetName());
            model_and_yield_names.push_back(nsig_11.GetName());
            model_and_yield_names.push_back(nsig_01.GetName());
        }
        else {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {h->lookupName(1), &_model_21}, {h->lookupName(0), &_model_11}, {h->lookupName(-1), &_model_01}
            },
            *h);
        }
    }

    // Construct target spin dependent pdf
    if (fitformula_pu=="" && fitformula_up!="" && fitformula_pp=="") {

        // Create the simultaneous pdf
        if (use_extended_nll) {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {t->lookupName(1), &model_12}, {t->lookupName(0), &model_11}, {t->lookupName(-1), &model_10}
            },
            *t);
            model_and_yield_names.push_back(nsig_12.GetName());
            model_and_yield_names.push_back(nsig_11.GetName());
            model_and_yield_names.push_back(nsig_10.GetName());
        }
        else {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {t->lookupName(1), &_model_12}, {t->lookupName(0), &_model_11}, {t->lookupName(-1), &_model_10}
            },
            *t);
        }
    }

    // Construct beam helicity and target spin dependent pdf
    if (fitformula_pu=="" && fitformula_up=="" && fitformula_pp!="") {

        // Create the simultaneous pdf
        if (use_extended_nll) {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {ht->lookupName(1), &model_22_00}, {ht->lookupName(0), &model_11}, {ht->lookupName(-1), &model_20_02}
            },
            *ht);
            model_and_yield_names.push_back(nsig_22_00.GetName());
            model_and_yield_names.push_back(nsig_11.GetName());
            model_and_yield_names.push_back(nsig_20_02.GetName());
        }
        else {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {ht->lookupName(1), &_model_22_00}, {ht->lookupName(0), &_model_11}, {ht->lookupName(-1), &_model_20_02}
            },
            *ht);
        }
    }

    // Construct FULL beam helicity and target spin dependent pdf **WITHOUT** PU asymmetries
    if (fitformula_pu=="" && fitformula_up!="" && fitformula_pp!="") {

        // Create the simultaneous pdf
        if (use_extended_nll) {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {ss->lookupName(11), &model_11}, // Polarization states: UU
                {ss->lookupName(12), &model_12}, {ss->lookupName(10), &model_10}, // Polarization states: UP
                {ss->lookupName(22), &model_22}, {ss->lookupName(0),  &model_00}, // Polarization states: PP
                {ss->lookupName(2),  &model_02}, {ss->lookupName(20), &model_20}  // Polarization states: PP
            },
            *ss);
            model_and_yield_names.push_back(nsig_11.GetName());
            model_and_yield_names.push_back(nsig_12.GetName()); model_and_yield_names.push_back(nsig_10.GetName());
            model_and_yield_names.push_back(nsig_22.GetName()); model_and_yield_names.push_back(nsig_00.GetName());
            model_and_yield_names.push_back(nsig_02.GetName()); model_and_yield_names.push_back(nsig_20.GetName());            
        }
        else {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {ss->lookupName(11), &_model_11}, // Polarization states: UU
                {ss->lookupName(12), &_model_12}, {ss->lookupName(10), &_model_10}, // Polarization states: UP
                {ss->lookupName(22), &_model_22}, {ss->lookupName(0),  &_model_00}, // Polarization states: PP
                {ss->lookupName(2),  &_model_02}, {ss->lookupName(20), &_model_20}  // Polarization states: PP
            },
            *ss);
        }
    }

    // Construct FULL beam helicity and target spin dependent pdf
    if (fitformula_pu!="" && fitformula_up!="" && fitformula_pp!="") {

        // Create the simultaneous pdf
        if (use_extended_nll) {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {ss->lookupName(11), &model_11}, // Polarization states: UU
                {ss->lookupName(21), &model_21}, {ss->lookupName(1),  &model_01}, // Polarization states: PU
                {ss->lookupName(12), &model_12}, {ss->lookupName(10), &model_10}, // Polarization states: UP
                {ss->lookupName(22), &model_22}, {ss->lookupName(0),  &model_00}, // Polarization states: PP
                {ss->lookupName(2),  &model_02}, {ss->lookupName(20), &model_20}  // Polarization states: PP
            },
            *ss);
            model_and_yield_names.push_back(nsig_11.GetName());
            model_and_yield_names.push_back(nsig_21.GetName()); model_and_yield_names.push_back(nsig_01.GetName());
            model_and_yield_names.push_back(nsig_12.GetName()); model_and_yield_names.push_back(nsig_10.GetName());
            model_and_yield_names.push_back(nsig_22.GetName()); model_and_yield_names.push_back(nsig_00.GetName());
            model_and_yield_names.push_back(nsig_02.GetName()); model_and_yield_names.push_back(nsig_20.GetName());            
        }
        else {
            model = new RooSimultaneous(model_name.c_str(), "simultaneous pdf",
            {
                {ss->lookupName(11), &_model_11}, // Polarization states: UU
                {ss->lookupName(21), &_model_21}, {ss->lookupName(1),  &_model_01}, // Polarization states: PU
                {ss->lookupName(12), &_model_12}, {ss->lookupName(10), &_model_10}, // Polarization states: UP
                {ss->lookupName(22), &_model_22}, {ss->lookupName(0),  &_model_00}, // Polarization states: PP
                {ss->lookupName(2),  &_model_02}, {ss->lookupName(20), &_model_20}  // Polarization states: PP
            },
            *ss);
        }
    }

    w->import(*model);
    return model_and_yield_names;

} // std::vector<std::string> getGenAsymPdf()

/**
* @brief Fit an asymmetry.
*
* Compute the bin count, bin variable mean values and variances, depolarization variable values and errors,
* and fit the asymmetry with a binned or unbinned dataset using a maximum likelihood fit method with an optional extended likelihood term.
* Note that for the maximum likelihood fit, the given asymmetry formulas \f$ A_{(UU,PU,UP,PP)}(x_0, x_1, ..., a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) \f$
* will be used internally by `getGenAsymPdf()` to construct a simultaneous PDF of the form:
* @f[
* \begin{aligned}
* PDF(\lambda_{\ell}, S_{||}, x_0, x_1, ..., &a_0, a_1, a_2, ..., d_0, d_1, d_2, ...) = \\
* & 1 + A_{UU=(UU,UT)}(\vec{x}, \vec{a}, \vec{d}) \\
* &   + \lambda_{\ell} \cdot \overline{\lambda_{\ell}^2} \cdot A_{PU=(LU,LT)}(\vec{x}, \vec{a}, \vec{d}) \\
* &   + S_{||} \cdot \overline{S^2} \cdot A_{UP=UL}(\vec{x}, \vec{a}, \vec{d}) \\
* &   + \lambda_{\ell} \cdot \overline{\lambda_{\ell}^2} \cdot S_{||} \cdot \overline{S^2} \cdot A_{PP=LL}(\vec{x}, \vec{a}, \vec{d}), \\
* \end{aligned}
* @f]
* where the appropriate terms will be dropped if there is no dependence on beam helicity or target spin.
* A simultaneous fit will be applied over the data subsets distinguished by the beam helicity, target spin, and product of beam helicity and target spin.
* The `a_<int>` denote the asymmetry amplitudes and the `d_<int>` denote the depolarization factors.
*
* The variable names in the fit formulas should follow the <a href="https://root.cern.ch/doc/master/classTFormula.html">TFormula</a> notation, e.g.,
* `x_0`\f$\rightarrow\f$`x[0]`, `x_1`\f$\rightarrow\f$`x[1]`, `a_0`\f$\rightarrow\f$`x[N_x]`, `a_1`\f$\rightarrow\f$`x[N_x+1]`, etc.
*
* In the case that a sideband dataset is supplied via `sb_dataset_name`, the initial dataset (`dataset_name`) is taken to be the signal region dataset.
* Then, a background PDF \f$A_{BG}(\vec{x}, \vec{a}_{BG}, \vec{d})\f$ is created identically to the signal PDF \f$A_{SG}(\vec{x}, \vec{a}_{SG}, \vec{d})\f$.
* The datasets should be created from `saga::signal::setBinnedBGFractions()` so that
* the background fraction variable \f$\varepsilon\f$ is already be loaded in the workspace and present in the datasets.
* A simultaneous PDF will be constructed over the combined signal region (\f$SG\f$) and sideband region (\f$SB\f$) datasets with the form:
* @f[
* PDF(\vec{x}, \vec{a}, \vec{d}) = 
* 1 + \bigg{\{}
* \begin{array}
* e \varepsilon(\vec{x}) \cdot A_{BG}(\vec{x}, \vec{a}_{BG}, \vec{d}) + (1 - \varepsilon(\vec{x})) \cdot A_{SG}(\vec{x}, \vec{a}_{SG}, \vec{d}), & \text{ } \vec{x} \in SG \\
* A_{BG}(\vec{x}, \vec{a}_{BG}, \vec{d}), & \text{ } \vec{x} \in SB \\
* \end{array}.
@f]
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
* - For each asymmetry fit parameter:
*
*   - Asymmetry fit parameter mean value
*
*   - Asymmetry fit parameter error
*
* The following entries will be appended if using sideband subtraction with binned background fractions:
*
* - For each background asymmetry fit parameter:
*
*   - Asymmetry fit parameter mean value
*
*   - Asymmetry fit parameter error
*
* @param w RooWorkspace in which to work
* @param dataset_name Dataset name
* @param bpol Luminosity averaged beam polarization \f$\overline{\lambda_{\ell}^2}\f$
* @param tpol Luminosity averaged target polarization \f$\overline{S^2}\f$
* @param categories_as_float List of category variables to treat as asymmetry fit variables and automatically add to PDF formulas
* @param helicity Name of the helicity variable
* @param tspin Name of the target spin variable
* @param htspin Name of the beam helicity times target spin variable
* @param combined_spin_state Name of the combined spin state variable
* @param binid Bin unique id
* @param bincut Kinematic variable cut for bin
* @param binvars List of kinematic binning variables
* @param depolvars List of depolarization variables
* @param fitvars List of asymmetry fit variables
* @param fitformula_uu The asymmetry formula in ROOT TFormula format for unpolarized and transverse target spin (\f$\phi_{S}\f$) dependent terms
* @param fitformula_pu The asymmetry formula in ROOT TFormula format for beam helicity dependent terms
* @param fitformula_up The asymmetry formula in ROOT TFormula format for target spin dependent terms
* @param fitformula_pp The asymmetry formula in ROOT TFormula format for beam helicity and target spin dependent terms
* @param initparams List of initial values for asymmetry parameters
* @param initparamlims List of initial asymmetry parameter minimum and maximum bounds
* @param use_sumw2error Option to use `RooFit::SumW2Error(true)` option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_binned_fit Option to use a binned fit to the data
* @param sb_dataset_name Name of the sideband dataset to use for the simultaneous fit of the signal + background and background PDFs
* @param bgfracvar Name of binned background fraction variable passed to `saga::signal::setBinnedBGFractions()`
* @param out Output stream
*
* @return List of bin count, bin variable means and errors, depolarization variable means and errors, fit parameters and errors
*/
std::vector<double> fitAsym(
        RooWorkspace                    *w,
        std::string                      dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
        double                           bpol,
        double                           tpol,
        std::vector<std::string>         categories_as_float,
        std::string                      helicity,
        std::string                      tspin,
        std::string                      htspin,
        std::string                      combined_spin_state,
        std::string                      binid,
        std::string                      bincut,
        std::vector<std::string>         binvars,
        std::vector<std::string>         depolvars,
        std::vector<std::string>         fitvars,
        std::string                      fitformula_uu,
        std::string                      fitformula_pu,
        std::string                      fitformula_up,
        std::string                      fitformula_pp,
        std::vector<double>              initparams,
        std::vector<std::vector<double>> initparamlims,
        bool use_sumw2error              = true,
        bool use_average_depol           = false,
        bool use_extended_nll            = false,
        bool use_binned_fit              = false,
        std::string sb_dataset_name      = "", //NOTE: If this is non-empty, a simultaneous fit of sg+bg and bg PDFs will be applied to the signal and sideband datasets respectively.
        std::string bgfracvar            = "",
        std::ostream &out                = std::cout
    ) {

    // Set method name
    std::string method_name = "fitAsym";

    // Load helicity variable from workspace
    RooCategory * h  = w->cat(helicity.c_str());
    RooCategory * t  = w->cat(tspin.c_str());
    RooCategory * ht = w->cat(htspin.c_str());
    RooCategory * ss = w->cat(combined_spin_state.c_str());

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
    std::vector<std::string> argnames;
    for (int ff=0; ff<fitvars.size(); ff++) { // Fit independent variables
        argset->add(*f[ff]);
        argnames.push_back(f[ff]->GetName());
    }
    for (int aa=0; aa<nparams; aa++) { // Fit asymmetry amplitude parameters
        argset->add(*a[aa]);
        argnames.push_back(a[aa]->GetName());
    }
    if (!use_average_depol) {
        for (int dd=0; dd<depolvars.size(); dd++) { // Fit depolarization factor variables
            argset->add(*d[dd]);
            argnames.push_back(d[dd]->GetName());
        }
    }

    // Optionally load sideband dataset from workspace and apply bin cuts //NOTE: DO THIS AFTER GETING THE COUNT, BINVARS, AND DEPOLVARS VALUES AND ERRORS.
    RooDataSet *ds_sb;
    RooDataSet *bin_ds_sb;
    RooCategory *region;
    RooRealVar *bgf;
    if (sb_dataset_name!="") {

        // Load bgfrac variable from workspace
        bgf = w->var(bgfracvar.c_str());

        // Initialize region category
        region = new RooCategory("region", "region");
        region->defineType("signal");
        region->defineType("sideband");

        // Load sideband dataset and apply bin cut
        ds_sb = (RooDataSet*)w->data(sb_dataset_name.c_str());
        bin_ds_sb = (RooDataSet*)ds_sb->reduce(bincut.c_str());

        // Construct combined dataset indexed on signal and background
        bin_ds = new RooDataSet(bin_ds->GetName(), bin_ds->GetTitle(), *bin_ds->get(), Index(*region),
                            Import({{"signal", bin_ds}, {"sideband", bin_ds_sb}}));
    }

    // Create asymmetry amplitude parameters for the background pdf
    std::vector<std::string> anames_sb;
    RooRealVar *a_sb[nparams];
    if (sb_dataset_name!="") {
        for (int aa=0; aa<nparams; aa++) {
            std::string aname = Form("a_bg%d",aa);
            anames_sb.push_back(aname);
            a_sb[aa] = new RooRealVar(anames_sb[aa].c_str(),anames_sb[aa].c_str(),initparams[aa],initparamlims[aa][0],initparamlims[aa][1]);
        }
    }

    // Add parameters to argument list in order for the background pdf
    RooArgSet *argset_sb = new RooArgSet();
    std::vector<std::string> argnames_sb;
    if (sb_dataset_name!="") {
        for (int ff=0; ff<fitvars.size(); ff++) { // Fit independent variables
            argset_sb->add(*f[ff]);
            argnames_sb.push_back(f[ff]->GetName());
        }
        for (int aa=0; aa<nparams; aa++) { // Fit asymmetry amplitude parameters
            argset_sb->add(*a_sb[aa]);
            argnames_sb.push_back(a_sb[aa]->GetName());
        }
        if (!use_average_depol) {
            for (int dd=0; dd<depolvars.size(); dd++) { // Fit depolarization factor variables
                argset_sb->add(*d[dd]);
                argnames_sb.push_back(d[dd]->GetName());
            }
        }
    }

    // Create and load asymmetry PDF
    std::vector<std::string> model_and_yield_names = getGenAsymPdf(
        w,
        categories_as_float,
        h,
        t,
        ht,
        ss,
        argset,
        argnames,
        method_name,
        binid,
        fitformula_uu,
        fitformula_pu,
        fitformula_up,
        fitformula_pp,
        bpol,
        tpol,
        count,
        use_extended_nll
    );
    std::string model_name = model_and_yield_names[0];

    // Create the asymmetry PDF for the background data in the case of sideband subtraction
    std::vector<std::string> model_and_yield_names_bg;
    std::string model_name_bg;
    std::string binid_sb;
    if (sb_dataset_name!="") {
        binid_sb = Form("%s_sb",binid.c_str());
        model_and_yield_names_bg = getGenAsymPdf(
            w,
            categories_as_float,
            h,
            t,
            ht,
            ss,
            argset_sb,
            argnames_sb,
            method_name,
            binid_sb,
            fitformula_uu,
            fitformula_pu,
            fitformula_up,
            fitformula_pp,
            bpol,
            tpol,
            count,
            use_extended_nll
        );
        model_name_bg = model_and_yield_names_bg[0];
    }

    // Load PDFs
    RooAbsPdf *model;
    RooAbsPdf *model_sg;
    RooAbsPdf *model_bg;
    RooAbsPdf *model_bg_plus_sg;

    // Load the signal PDF
    if (sb_dataset_name=="") {

        model = w->pdf(model_name.c_str());

    } else {

        // Load signal and background PDFs
        model_sg = w->pdf(model_name.c_str());
        model_bg = w->pdf(model_name_bg.c_str());

        // Add signal and background PDFs using the bgfrac variable as the coefficient of the background
        std::string model_bg_plus_sg_name = Form("%s_plus_bg", model_name.c_str());
        model_bg_plus_sg = new RooAddPdf(model_bg_plus_sg_name.c_str(), model_bg_plus_sg_name.c_str(), RooArgList(*model_bg, *model_sg), RooArgList(*bgf)); //NOTE: ORDER IS IMPORTANT HERE!

        // Construct a simultaneous PDF for the signal and background regions
        model = new RooSimultaneous(Form("%s_sb_bgfracs",model_name.c_str()), "simultaneous pdf",
                {
                    {"signal", model_bg_plus_sg},
                    {"sideband", model_bg}
                },
                *region
        );
    }

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

    // // Define the asymmetry as a function
    // RooFormulaVar f_asym("f_asym","Asymmetry function",Form("%.3f*%s",bpol,fitformula.c_str()), *argset);//NOTE: NEED TO CORRECT FOR POLARIZATION FACTOR.

    // // Switch off histogram stats
    // gStyle->SetOptStat(0);

    // // Loop fit variables and plot fit projections
    // for (int idx=0; idx<fitvars.size(); idx++) {

    //     // Plot projection of fitted distribution in fit variable
    //     RooPlot *xframe = f[idx]->frame(RooFit::Title(Form("%s Projection, Bin: %s",f[idx]->GetTitle(),bincut.c_str())));
    //     bin_ds->plotOn(xframe, RooFit::Asymmetry(*h));
    //     f_asym.plotOn(xframe, RooFit::LineColor(kRed));

    //     // Draw the frame on the canvas
    //     std::string c1_x_name = Form("c1_%s__fitvar_%s",binid.c_str(),fitvars[idx].c_str());
    //     TCanvas *c1_x = new TCanvas(c1_x_name.c_str(), c1_x_name.c_str());
    //     gPad->SetLeftMargin(0.15);
    //     xframe->GetYaxis()->SetTitleOffset(1.4);
    //     xframe->Draw();
    //     c1_x->Print(Form("%s.pdf",c1_x_name.c_str()));
    // }

    // Get fit parameter values and errors
    std::vector<double> params;
    std::vector<double> paramerrs;
    for (int aa=0; aa<nparams; aa++) {
        RooRealVar *avar = (RooRealVar*)w->var(a[aa]->GetName()); //NOTE: Load from workspace since parameters are copied to work space when you import the PDF.
        params.push_back((double)avar->getVal());
        paramerrs.push_back((double)avar->getError());
    }

    // Get fit parameter values and errors for background values
    std::vector<double> params_bg;
    std::vector<double> paramerrs_bg;
    if (sb_dataset_name!="") {
        for (int aa=0; aa<nparams; aa++) {
            RooRealVar *avar_sb = (RooRealVar*)w->var(a_sb[aa]->GetName()); //NOTE: Load from workspace since parameters are copied to work space when you import the PDF.
            params_bg.push_back((double)avar_sb->getVal());
            paramerrs_bg.push_back((double)avar_sb->getError());
        }
    }

    // Get the raw counts and poissonian errors
    std::vector<double> counts;
    std::vector<double> counterrs;
    double count_h_pos  = (double)bin_ds->reduce(Form("%s>0",h->GetName()))->sumEntries();
    double count_h_neg  = (double)bin_ds->reduce(Form("%s<0",h->GetName()))->sumEntries();
    double count_t_pos  = (double)bin_ds->reduce(Form("%s>0",t->GetName()))->sumEntries();
    double count_t_neg  = (double)bin_ds->reduce(Form("%s<0",t->GetName()))->sumEntries();
    double count_ht_pos = (double)bin_ds->reduce(Form("%s>0",ht->GetName()))->sumEntries();
    double count_ht_neg = (double)bin_ds->reduce(Form("%s<0",ht->GetName()))->sumEntries();
    double counterr_h_pos   = (double)TMath::Sqrt(count_h_pos);
    double counterr_h_neg   = (double)TMath::Sqrt(count_h_neg);
    double counterr_t_pos   = (double)TMath::Sqrt(count_t_pos);
    double counterr_t_neg   = (double)TMath::Sqrt(count_t_neg);
    double counterr_ht_pos  = (double)TMath::Sqrt(count_ht_pos);
    double counterr_ht_neg  = (double)TMath::Sqrt(count_ht_neg);

    // Get the fitted yields and errors if using an extended fit
    if (use_extended_nll && categories_as_float.size()==0) {

        // Loop fitted yields and get values and errors
        for (int nn=1; nn<model_and_yield_names.size(); nn++) {
            RooRealVar *avar = (RooRealVar*)w->var(model_and_yield_names[nn].c_str()); //NOTE: Load from workspace since parameters are copied to work space when you import the PDF.
            counts.push_back((double)avar->getVal());
            counterrs.push_back((double)avar->getError());
        }

        // Get yields for helicity dependent pdf
        if (fitformula_pu!="" && fitformula_up=="" && fitformula_pp=="") {

            int j = 0;
            double nsig_21 = counts[j++];
            double nsig_11 = counts[j++];
            double nsig_01 = counts[j++];
            int k = 0;
            double nsigerr_21 = counterrs[k++];
            double nsigerr_11 = counterrs[k++];
            double nsigerr_01 = counterrs[k++];

            count_h_pos    = nsig_21;
            count_h_neg    = nsig_01;
            counterr_h_pos = nsigerr_21;
            counterr_h_neg = nsigerr_01;
        }

        // Get yields for target spin dependent pdf
        if (fitformula_pu=="" && fitformula_up!="" && fitformula_pp=="") {

            int j = 0;
            double nsig_12 = counts[j++];
            double nsig_11 = counts[j++];
            double nsig_10 = counts[j++];
            int k = 0;
            double nsigerr_12 = counterrs[k++];
            double nsigerr_11 = counterrs[k++];
            double nsigerr_10 = counterrs[k++];

            count_t_pos    = nsig_12;
            count_t_neg    = nsig_10;
            counterr_t_pos = nsigerr_12;
            counterr_t_neg = nsigerr_10;
        }

        // Get yields for beam helicity and target spin dependent pdf
        if (fitformula_pu=="" && fitformula_up=="" && fitformula_pp!="") {

            int j = 0;
            double nsig_22_00 = counts[j++];
            double nsig_11    = counts[j++];
            double nsig_20_02 = counts[j++];
            int k = 0;
            double nsigerr_22_00 = counterrs[k++];
            double nsigerr_11    = counterrs[k++];
            double nsigerr_20_02 = counterrs[k++];

            count_ht_pos    = nsig_22_00;
            count_ht_neg    = nsig_20_02;
            counterr_ht_pos = nsigerr_22_00;
            counterr_ht_neg = nsigerr_20_02;
        }

        // Get yields for FULL beam helicity and target spin dependent pdf **WITHOUT** PU asymmetries
        if (fitformula_pu=="" && fitformula_up!="" && fitformula_pp!="") {

            int j = 0;
            double nsig_11 = counts[j++];
            double nsig_12 = counts[j++]; double nsig_10 = counts[j++];
            double nsig_22 = counts[j++]; double nsig_00 = counts[j++];
            double nsig_02 = counts[j++]; double nsig_20 = counts[j++];
            int k = 0;
            double nsigerr_11 = counterrs[k++];
            double nsigerr_12 = counterrs[k++]; double nsigerr_10 = counterrs[k++];
            double nsigerr_22 = counterrs[k++]; double nsigerr_00 = counterrs[k++];
            double nsigerr_02 = counterrs[k++]; double nsigerr_20 = counterrs[k++];
        }

        // Get yields for FULL beam helicity and target spin dependent pdf
        if (fitformula_pu!="" && fitformula_up!="" && fitformula_pp!="") {

            int j = 0;
            double nsig_11 = counts[j++];
            double nsig_21 = counts[j++]; double nsig_01 = counts[j++];
            double nsig_12 = counts[j++]; double nsig_10 = counts[j++];
            double nsig_22 = counts[j++]; double nsig_00 = counts[j++];
            double nsig_02 = counts[j++]; double nsig_20 = counts[j++];
            int k = 0;
            double nsigerr_11 = counterrs[k++];
            double nsigerr_21 = counterrs[k++]; double nsigerr_01 = counterrs[k++];
            double nsigerr_12 = counterrs[k++]; double nsigerr_10 = counterrs[k++];
            double nsigerr_22 = counterrs[k++]; double nsigerr_00 = counterrs[k++];
            double nsigerr_02 = counterrs[k++]; double nsigerr_20 = counterrs[k++];

            count_h_pos     = nsig_21;
            count_h_neg     = nsig_01;
            counterr_h_pos  = nsigerr_21;
            counterr_h_neg  = nsigerr_01;
            count_t_pos     = nsig_12;
            count_t_neg     = nsig_10;
            counterr_t_pos  = nsigerr_12;
            counterr_t_neg  = nsigerr_10;
            count_ht_pos    = nsig_22 + nsig_00;
            count_ht_neg    = nsig_20 + nsig_02;
            counterr_ht_pos = (double)TMath::Sqrt(nsigerr_22*nsigerr_22 + nsigerr_00*nsigerr_00);
            counterr_ht_neg = (double)TMath::Sqrt(nsigerr_20*nsigerr_20 + nsigerr_02*nsigerr_02);
        }

        // Set counts knowing how they are added in saga::analysis::getGenAsymPdf()
    }

    // Set the raw asymmetries
    std::vector<double> rawasyms;
    std::vector<double> rawasymerrs;
    double asym_h  = (count_h_pos-count_h_neg)/(count_h_pos+count_h_neg);
    double asym_t  = (count_t_pos-count_t_neg)/(count_t_pos+count_t_neg);
    double asym_ht = (count_ht_pos-count_ht_neg)/(count_ht_pos+count_ht_neg);
    rawasyms.push_back(asym_h);
    rawasyms.push_back(asym_t);
    rawasyms.push_back(asym_ht);
    double asymerr_h  = (1.0 + asym_h )*(counterr_h_pos*counterr_h_pos   + counterr_h_neg*counterr_h_neg  )/(count_h_pos+count_h_neg);
    double asymerr_t  = (1.0 + asym_t )*(counterr_t_pos*counterr_t_pos   + counterr_t_neg*counterr_t_neg  )/(count_t_pos+count_t_neg);
    double asymerr_ht = (1.0 + asym_ht)*(counterr_ht_pos*counterr_ht_pos + counterr_ht_neg*counterr_ht_neg)/(count_ht_pos+count_ht_neg);
    rawasymerrs.push_back(asymerr_h);
    rawasymerrs.push_back(asymerr_t);
    rawasymerrs.push_back(asymerr_ht);

    // Print out fit info
    out << "--------------------------------------------------" << std::endl;
    out << " "<<method_name.c_str()<<"():" << std::endl;
    out << " bpol        = " << bpol << std::endl;
    out << " tpol        = " << tpol << std::endl;
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
    out << " fitformula_uu  = " << fitformula_uu.c_str() << std::endl;
    out << " fitformula_pu  = " << fitformula_pu.c_str() << std::endl;
    out << " fitformula_up  = " << fitformula_up.c_str() << std::endl;
    out << " fitformula_pp  = " << fitformula_pp.c_str() << std::endl;
    out << " nparams        = " << nparams <<std::endl;
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
    if (sb_dataset_name!="") {
        out << " params_bg = [" ;
        for (int idx=0; idx<nparams; idx++) {
            out << params_bg[idx] << "±" << paramerrs_bg[idx];
            if (idx<nparams-1) { out << " , "; }
        }
        out << "]" << std::endl;
    }
    if (use_extended_nll) {
        out << " yields = [" ;
        for (int idx=1; idx<model_and_yield_names_bg.size(); idx++) {
            RooRealVar *yield = (RooRealVar*)w->var(model_and_yield_names_bg[idx].c_str());
            out << yield->getVal() << "±" << yield->getError();
            if (idx<nparams-1) { out << " , "; }
        }
        out << "]" << std::endl;
    }
    out << " rawasyms = [" ;
    for (int idx=0; idx<rawasyms.size(); idx++) {
        out << rawasyms[idx] << "±" << rawasymerrs[idx];
        if (idx<rawasyms.size()-1) { out << " , "; }
    }
    out << "]" << std::endl;
    out << "--------------------------------------------------" << std::endl;

    // Fill return array
    std::vector<double> arr; //NOTE: Dimension = 1+2*binvars.size()+2*depolvars.size()+2*rawasyms.size()+2*nparams(+2*nparams)
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
    if (sb_dataset_name!="") {
        for (int idx=0; idx<nparams; idx++) {
            arr.push_back(params_bg[idx]);
            arr.push_back(paramerrs_bg[idx]);
        }
    }

    return arr;

} // std::vector<double> fitAsym()

/**
* @brief Loop kinematic bins and fit an asymmetry, correcting for background with sideband subtraction or <a href="http://arxiv.org/abs/physics/0402083">sPlots</a>.
*
* Loop bins cuts and fit an asymmetry with the `saga::analysis::fitAsym()` method.  Optionally, apply an invariant mass fit and background correction using the
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
* - For each asymmetry fit parameter `asymfitpar`
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
* @param asymfit_formula_uu The asymmetry formula in ROOT TFormula format for the unpolarized and transverse target spin (\f$\phi_{S}\f$) dependent asymmetries
* @param asymfit_formula_pu The asymmetry formula in ROOT TFormula format for the beam helicity dependent asymmetries
* @param asymfit_formula_up The asymmetry formula in ROOT TFormula format for the target spin dependent asymmetries
* @param asymfit_formula_pp The asymmetry formula in ROOT TFormula format for the beam helicity and target spin dependent asymmetries
* @param asymfitpar_inits List of initial values for asymmetry fit variables
* @param asymfitpar_initlims List of initial asymmetry fit variables minimum and maximum bounds
* @param use_sumw2error Option to use `RooFit::SumW2Error(true)` option when fitting to dataset which is necessary if using a weighted dataset
* @param use_average_depol Option to divide out average depolarization in bin instead of including depolarization as an independent variable in the fit
* @param use_extended_nll Option to use an extended Negative Log Likelihood function for minimization
* @param use_binned_fit Option to use a binned fit to the data

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
*/
void getKinBinnedAsym(
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
        std::vector<std::string>         asymfitvars,
        std::vector<std::string>         asymfitvar_titles,
        std::vector<std::vector<double>> asymfitvar_lims,
        std::vector<int>                 asymfitvar_bins,
        std::vector<std::string>         massfitvars,
        std::vector<std::string>         massfitvar_titles,
        std::vector<std::vector<double>> massfitvar_lims,
        std::vector<int>                 massfitvar_bins,

        // parameterss passed to analysis::fitAsym()
        double                           bpol,
        double                           tpol,
        std::string                      asymfit_formula_uu,
        std::string                      asymfit_formula_pu,
        std::string                      asymfit_formula_up,
        std::string                      asymfit_formula_pp,
        std::vector<double>              asymfitpar_inits,
        std::vector<std::vector<double>> asymfitpar_initlims,
        bool                             use_sumw2error,
        bool                             use_average_depol,
        bool                             use_extended_nll,
        bool                             use_binned_fit,

        // parameters passed to saga::signal::fitMass() //TODO: Add init fit parameter value and limits arguments here...assuming you always want a chebychev polynomial background...
        std::map<std::string,std::string> massfit_yamlfile_map,
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

        // Parameters passed to analysis::applySPlots()
        bool                             use_splot,

        // Parameters used for sb subtraction
        std::string                      massfit_sgcut,
        std::string                      massfit_bgcut,
        bool                             use_sb_subtraction,
        bool                             use_binned_sb_bgfracs,
        std::map<int,std::string>        asymfitvar_bincuts,
        std::string                      bgfracvar,
        std::vector<double>              bgfracvar_lims,
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
        std::ostream                    &out                      = std::cout
    ) {

    // Check arguments
    if (binvars.size()<1) {std::cerr<<"ERROR: Number of bin variables is <1.  Exiting...\n"; return;}
    if (depolvars.size()!=asymfitvars.size()) {std::cerr<<"WARNING: depolvars.size() does not match the number of parameters injected."<<std::endl;}
    if ((use_sb_subtraction && use_binned_sb_bgfracs) || (use_sb_subtraction && use_splot) || (use_binned_sb_bgfracs && use_splot)) {std::cerr<<"ERROR: Sideband subtraction, sideband subtraction with binned background fractions, and the sPlot method are all mutually exclusive.  Exiting...\n"; return;}

    // Starting message
    out << "----------------------- getKinBinnedAsym ----------------------\n";
    out << "bincuts = { ";
    for (auto it = bincuts.begin(); it != bincuts.end(); it++) {
        out << it->first << " : " << it->second.c_str() << " , ";
    }
    out << " }\n";

    // Filter frames for signal and sideband
    massfit_sgcut = (massfit_sgcut.size()>0) ? massfit_sgcut : saga::util::addLimitCuts("",massfitvars,massfit_sgregion_lims);
    auto frame_sg = (massfit_sgcut.size()>0) ? frame.Filter(massfit_sgcut.c_str()) : frame;
    auto frame_sb = (massfit_bgcut.size()>0) ? frame.Filter(massfit_bgcut.c_str()) : frame;

    // Set condition for single mass fit
    bool single_massfit = (massfit_pdf_name!="" && !use_binned_sb_bgfracs && (use_splot || use_sb_subtraction));

    // Open output CSV
    std::string csvpath = Form("%s.csv",scheme_name.c_str());
    std::ofstream csvoutf; csvoutf.open(csvpath.c_str());
    std::ostream &csvout = csvoutf;
    std::string csv_separator = ",";
    std::vector<std::string> rawasymvars = { "asym_h", "asym_t", "asym_ht"}; //TODO: OPTIONALLY CORRECT THESE BY A GIVEN DEPOLVAR

    // Set CSV column headers
    // COLS: bin_id,count,{binvarmean,binvarerr},{depolvarmean,depolvarerr},{rawasym,rawasymerr},{asymfitvar,asymfitvarerr},{fitvar_info if requested}
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
    for (int aa=0; aa<asymfitpar_inits.size(); aa++) {
        csvout << Form("a%d",aa) << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
        csvout << Form("a%d",aa) << "_err";
        if (aa<asymfitpar_inits.size()-1 || single_massfit || use_binned_sb_bgfracs) csvout << csv_separator.c_str();
        else csvout << std::endl;//NOTE: IMPORTANT!
    }

    // Optionally add background asymmetries
    if (use_binned_sb_bgfracs || use_sb_subtraction) {
        for (int aa=0; aa<asymfitpar_inits.size(); aa++) {
            csvout << Form("a_bg%d",aa) << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
            csvout << Form("a_bg%d",aa) << "_err";
            if (aa<asymfitpar_inits.size()-1 || single_massfit) csvout << csv_separator.c_str();
            else csvout << std::endl;//NOTE: IMPORTANT!
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
            else csvout << std::endl;//NOTE: IMPORTANT!
        }

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
        std::vector<double> massfit_result;
        if (single_massfit) {  //NOTE: A mass fit in each bin is needed for basic sideband subtraction and splots.

            // Set yaml path for mass fit parameters
            std::string yamlfile = massfit_yamlfile_map[scheme_binid];

            // Fit the mass spectrum
            std::vector<double> massfit_result = saga::signal::fitMass(
                    ws, // RooWorkspace                    *w,
                    dataset_name, // std::string                      dataset_name,
                    scheme_binid, // std::string                      binid,
                    bin_cut, // std::string                      bincut,
                    binvars, // std::vector<std::string>         binvars,
                    massfitvars, // std::vector<std::string>         fitvars,
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
                    massfit_use_sumw2error, // bool                             massfit_use_sumw2error   = false,
                    massfit_use_extended_nll, // bool                             massfit_use_extended_nll = true,
                    massfit_use_binned_fit, // bool                             massfit_use_binned_fit   = false,
                    out // std::ostream                    &out              = std::cout
            );
        }

        // Apply sPlot
        std::string fit_dataset_name = dataset_name; // -> Use this for sPlot
        if (use_splot) {
            std::string dataset_sg_name = (std::string)Form("%s_sg_sw",dataset_name.c_str());
            std::string dataset_bg_name = (std::string)Form("%s_bg_sw",dataset_name.c_str());
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
        std::string fit_sb_dataset_name = ""; // -> Use this for binned sideband backgrounds
        if (use_binned_sb_bgfracs) {
            std::string rds_out_name = (std::string)Form("%s_sg",dataset_name.c_str());
            std::string sb_rds_out_name = (std::string)Form("%s_sb",dataset_name.c_str());
            saga::signal::setBinnedBGFractions(
                ws, // RooWorkspace                    *w,
                dataset_name, // std::string                      dataset_name,
                scheme_binid, // std::string                      binid,
                bin_cut, // std::string                      bincut,
                binvars, // std::vector<std::string>         binvars,
                massfitvars, // std::vector<std::string>         fitvars,
                massfit_yamlfile_map, // std::map<std::string,std::string> yamlfile_map
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

                binframe, // RNode                            frame, // arguments for this method
                massfit_bgcut, // std::string                      bgcut, 
                asymfitvars, // std::vector<std::string>         asymfitvars,
                asymfitvar_bincuts, // std::map<int,std::string>        asymfitvar_bincuts,
                rds_out_name, // std::string                      rds_out_name,
                sb_rds_out_name, // std::string                   sb_rds_out_name,
                bgfracvar, // std::string                      bgfracvar,
                bgfracvar_lims, // std::vector<double>              bgfracvar_lims,

                massfit_lg_text_size, // double                           massfit_lg_text_size     = 0.04,
                massfit_lg_margin, // double                           massfit_lg_margin        = 0.1,
                massfit_lg_ncols, // int                              massfit_lg_ncols         = 1,
                massfit_plot_bg_pars, // bool                             massfit_plot_bg_pars     = false,
                massfit_use_sumw2error, // bool                             massfit_use_sumw2error   = false,
                massfit_use_extended_nll, // bool                             massfit_use_extended_nll = true,
                massfit_use_binned_fit, // bool                             massfit_use_binned_fit   = false,

                bgfrac_idx, // int                               bgfrac_idx               = 0,
                0.0, // double                           bgfracs_default  = 0.0 // arguments for this method
                out // std::ostream                    &out              = std::cout
            );
            fit_dataset_name = rds_out_name;
            fit_sb_dataset_name = sb_rds_out_name;
        }

        // Create signal region dataset for sideband subtraction
        if (use_sb_subtraction) {
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
        std::vector<double> asymfit_result = fitAsym(
                                (use_sb_subtraction ? ws_sg : ws),
                                fit_dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                bpol,
                                tpol,
                                categories_as_float,
                                helicity,
                                tspin,
                                htspin,
                                combined_spin_state,
                                scheme_binid,
                                bin_cut,
                                binvars,
                                depolvars,
                                asymfitvars,
                                asymfit_formula_uu,
                                asymfit_formula_pu,
                                asymfit_formula_up,
                                asymfit_formula_pp,
                                asymfitpar_inits,
                                asymfitpar_initlims,
                                use_sumw2error,
                                use_average_depol,
                                use_extended_nll,
                                use_binned_fit,
                                fit_sb_dataset_name,
                                bgfracvar,
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
            asymfit_result_sb = fitAsym(
                                ws_sb,
                                dataset_name, //NOTE: DATASET SHOULD ALREADY BE FILTERED WITH OVERALL CUTS AND CONTAIN WEIGHT VARIABLE IF NEEDED
                                bpol,
                                tpol,
                                categories_as_float,
                                helicity,
                                tspin,
                                htspin,
                                combined_spin_state,
                                sb_scheme_binid,
                                bin_cut,
                                binvars,
                                depolvars,
                                asymfitvars,
                                asymfit_formula_uu,
                                asymfit_formula_pu,
                                asymfit_formula_up,
                                asymfit_formula_pp,
                                asymfitpar_inits,
                                asymfitpar_initlims,
                                use_sumw2error,
                                use_average_depol,
                                use_extended_nll,
                                use_binned_fit,
                                fit_sb_dataset_name,
                                bgfracvar,
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
        std::vector<double> chi2ndfs;
        std::vector<double> massfit_pars_sg;
        std::vector<double> massfit_parerrs_sg;
        std::vector<double> massfit_pars_bg;
        std::vector<double> massfit_parerrs_bg;
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
        out << "--- Acceptance, background, and depolarization corrected signal ---\n";
        for (int idx=0; idx<nparams; idx++) {
            if (use_average_depol) {
                out << " ys["<< idx <<"]             = " << ys[idx] << "\n";
                out << " eys["<< idx <<"]            = " << eys[idx] << "\n";
                out << " depols["<< idx <<"]         = " << depols[idx] << "\n";
                out << " edepols["<< idx <<"]        = " << edepols[idx] << "\n";
            } if (use_sb_subtraction || use_binned_sb_bgfracs) {
                out << " ys_sb["<< idx <<"]       = " << ys_sb[idx] << "\n";
                out << " eys_sb["<< idx <<"]      = " << eys_sb[idx] << "\n";
            }
            out << " ys_corrected["<< idx <<"]   = " << ys_corrected[idx] << "\n";
            out << " eys_corrected["<< idx <<"]  = " << eys_corrected[idx] << "\n";
        }
        out << "---------------------------\n";

        // Write out a row of data to csv
        // COLS: bin_id,count,{binvarmean,binvarerr},{depolvarmean,depolvarerr},{rawasym,rawasymerr},{asymfitvar,asymfitvarerr}(,{bg_asymfitvar,bg_asymfitvarerr})
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
            csvout << rawasyms[rr] << csv_separator.c_str();
        }
        for (int aa=0; aa<nparams; aa++) {
            csvout << ys_corrected[aa] << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
            csvout << eys_corrected[aa];
            if (aa<nparams-1 || single_massfit || use_binned_sb_bgfracs) csvout << csv_separator.c_str();
            else csvout << std::endl;//NOTE: IMPORTANT!
        }

        // Optionally add background asymmetries
        if (use_binned_sb_bgfracs || use_sb_subtraction) {
            for (int aa=0; aa<asymfitpar_inits.size(); aa++) {
                csvout << ys_sb[aa] << csv_separator.c_str();//NOTE: This is the default naming from analysis::fitAsym()
                csvout << eys_sb[aa];
                if (aa<nparams-1 || single_massfit) csvout << csv_separator.c_str();
                else csvout << std::endl;//NOTE: IMPORTANT!
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
                else csvout << std::endl;//NOTE: IMPORTANT!
            }

        }
    }

    csvoutf.close();
    out << " Saved asymmetry fit results to " << csvpath.c_str() << std::endl;

    // Ending message
    out << "------------------- END of getKinBinnedAsym -------------------\n";

} // getKinBinnedAsym()

} // namespace analysis {

} // namespace saga {
