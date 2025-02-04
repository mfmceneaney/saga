#include <iostream>
#include <memory>
#include <fstream>
#include <string>

// YAML Includes
#include <yaml-cpp/yaml.h>

// // ARGPARSE Includes
// #include "argparse.h"

// ROOT Includes
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <TRandom.h>
#include <TMath.h>
#include <TF3.h>

// Project Includes
#include <analysis.h>

// Define find and replace function following solution here: https://stackoverflow.com/questions/5878775/how-to-find-and-replace-string
void replace_all(
    std::string& s,
    std::string const& toReplace,
    std::string const& replaceWith
) {
    std::string buf;
    std::size_t pos = 0;
    std::size_t prevPos;

    // Reserves rough estimate of final size of string.
    buf.reserve(s.size());

    while (true) {
        prevPos = pos;
        pos = s.find(toReplace, pos);
        if (pos == std::string::npos)
            break;
        buf.append(s, prevPos, pos - prevPos);
        buf += replaceWith;
        pos += toReplace.size();
    }

    buf.append(s, prevPos, s.size() - prevPos);
    s.swap(buf);
}

void execute(const YAML::Node& node) {

    // Process arguments
    std::string outdir = "";
    if (node["outdir"]) {
        outdir = node["outdir"].as<std::string>();
    }
    std::cout << "INFO: outdir: " << outdir << std::endl;

    std::string outpath = "out.root";
    if (node["outpath"]) {
        outpath = node["outpath"].as<std::string>();
    }
    std::cout << "INFO: outpath: " << outpath << std::endl;

    std::string inpath = "";
    if (node["inpath"]) {
        inpath = node["inpath"].as<std::string>();
    }
    std::cout << "INFO: inpath: " << inpath << std::endl;

    std::string tree = "";
    if (node["tree"]) {
        tree = node["tree"].as<std::string>();
    }
    std::cout << "INFO: tree: " << tree << std::endl;

    int nthreads = 1;
    if (node["nthreads"]) {
        nthreads = node["nthreads"].as<int>();
    }
    std::cout << "INFO: nthreads: " << nthreads << std::endl;

    std::string cuts = "";
    if (node["cuts"]) {
        cuts = node["cuts"].as<std::string>();
    }
    std::cout << "INFO: cuts: " << cuts << std::endl;

    std::string method = "BSA1D"; //NOTE: This can stay a std::string since it's a TString later...
    if (node["method"]) {
        method = node["method"].as<std::string>();
    }
    std::cout << "INFO: method: " << method << std::endl;

    std::string fitvar1 = "";
    if (node["fitvar1"]) {
        fitvar1 = node["fitvar1"].as<std::string>();
    }
    std::cout << "INFO: fitvar1: " << fitvar1 << std::endl;

    std::string fitvar1formula = "";
    if (node["fitvar1formula"]) {
        fitvar1formula = node["fitvar1formula"].as<std::string>();
    }
    std::cout << "INFO: fitvar1formula: " << fitvar1formula << std::endl;

    std::string fitvar1formulamc = "";
    if (node["fitvar1formulamc"]) {
        fitvar1formulamc = node["fitvar1formulamc"].as<std::string>();
    }
    std::cout << "INFO: fitvar1formulamc: " << fitvar1formulamc << std::endl;

    std::string fitvar1title = "dphi";
    if (node["fitvar1title"]) {
        fitvar1title = node["fitvar1title"].as<std::string>();
    }
    std::cout << "INFO: fitvar1title: " << fitvar1title << std::endl;

    int fitvar1bins = 16;
    if (node["fitvar1bins"]) {
        fitvar1bins = node["fitvar1bins"].as<int>();
    }
    std::cout << "INFO: fitvar1bins: " << fitvar1bins << std::endl;

    // Define fit / injection function formulas
    std::string fitformula = "";
    if (node["fitformula"]) {
        fitformula = node["fitformula"].as<std::string>();
    }
    std::cout << "INFO: fitformula: " << fitformula << std::endl;

    int nparams = 2; //NOTE: THIS IS ONLY NECESSARILY FOR  fitformula.  fsgasyms and fbgasyms should have nparameters equal in number to the number of provided asymmetries to inject in sgasyms and bgasyms.
    if (node["nparams"]) {
        nparams = node["nparams"].as<int>();
    }
    std::cout << "INFO: nparams: " << nparams << std::endl;

    std::vector<double> params = std::vector<double>(nparams);
    if (node["params"]) {
        params = node["params"].as<std::vector<double>>();
    }
    std::cout << "INFO: params: [ ";
    for (int idx=0; idx<params.size(); idx++) {
        if (idx!=params.size()-1) { std::cout << params[idx]<<", "; }
        else { std::cout << params[idx]; }
    }
    std::cout << " ]" << std::endl;

    std::string fitopt = "S";
    if (node["fitopt"]) {
        fitopt = node["fitopt"].as<std::string>();
    }
    std::cout << "INFO: fitopt: " << fitopt << std::endl;

    std::string gammavar = "gamma";
    if (node["gammavar"]) {
        gammavar = node["gammavar"].as<std::string>();
    }
    std::cout << "INFO: gammavar: " << gammavar << std::endl;

    std::string gammavarformula = "";
    if (node["gammavarformula"]) {
        gammavarformula = node["gammavarformula"].as<std::string>();
    }
    std::cout << "INFO: gammavarformula: " << gammavarformula << std::endl;

    std::string gammavarformulamc = "";
    if (node["gammavarformulamc"]) {
        gammavarformulamc = node["gammavarformulamc"].as<std::string>();
    }
    std::cout << "INFO: gammavarformulamc: " << gammavarformulamc << std::endl;

    std::string epsilonvar = "epsilon";
    if (node["epsilonvar"]) {
        epsilonvar = node["epsilonvar"].as<std::string>();
    }
    std::cout << "INFO: epsilonvar: " << epsilonvar << std::endl;

    std::string epsilonvarformula = "";
    if (node["epsilonvarformula"]) {
        epsilonvarformula = node["epsilonvarformula"].as<std::string>();
    }
    std::cout << "INFO: epsilonvarformula: " << epsilonvarformula << std::endl;

    std::string epsilonvarformulamc = "";
    if (node["epsilonvarformulamc"]) {
        epsilonvarformulamc = node["epsilonvarformulamc"].as<std::string>();
    }
    std::cout << "INFO: epsilonvarformulamc: " << epsilonvarformulamc << std::endl;

    std::vector<std::string> depolvars;
    if (node["depolvars"]) {
        depolvars = node["depolvars"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: depolvars: [ ";
    for (int idx=0; idx<depolvars.size(); idx++) {
        if (idx!=depolvars.size()-1) { std::cout << depolvars[idx]<<", "; }
        else { std::cout << depolvars[idx]; }
    }
    std::cout << " ]" << std::endl;

    std::vector<std::string> depolvarformulas;
    if (node["depolvarformulas"]) {
        depolvarformulas = node["depolvarformulas"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: depolvarformulas: [ ";
    for (int idx=0; idx<depolvarformulas.size(); idx++) {
        if (idx!=depolvarformulas.size()-1) { std::cout << depolvarformulas[idx]<<", "; }
        else { std::cout << depolvarformulas[idx]; }
    }
    std::cout << " ]" << std::endl;

    std::vector<std::string> depolvarformulasmc;
    if (node["depolvarformulasmc"]) {
        depolvarformulasmc = node["depolvarformulasmc"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: depolvarformulasmc: [ ";
    for (int idx=0; idx<depolvarformulasmc.size(); idx++) {
        if (idx!=depolvarformulasmc.size()-1) { std::cout << depolvarformulasmc[idx]<<", "; }
        else { std::cout << depolvarformulasmc[idx]; }
    }
    std::cout << " ]" << std::endl;

    //----------------------------------------------------------------------------------------------------//
    // Get bin variables and poly4 bg maps
    std::map<std::string,std::vector<double>> binvars;
    std::map<std::string,std::vector<int>> poly4map;
    if (node["binvars"]) {
        if (node["binvars"].IsMap()) {
            std::cout << "binvars:" << std::endl;//DEBUGGING
            for (auto it = node["binvars"].begin(); it != node["binvars"].end(); ++it) { //TODO: How to check if too many binning variables...

                // Get bin variable name
                std::string name = it->first.as<std::string>();

                // Get nbins and bins from yaml
                int nbins = 0;
                if (node["binvars"][name]["nbins"]) {
                    nbins = node["binvars"][name]["nbins"].as<int>();
                }
                std::vector<double> bins;
                if (node["binvars"][name]["bins"]) {
                    bins = node["binvars"][name]["bins"].as<std::vector<double>>();
                }
                std::vector<int> poly4bins;
                if (node["binvars"][name]["poly4bins"]) {
                    poly4bins = node["binvars"][name]["poly4bins"].as<std::vector<int>>();
                }
                
                // Set bin limits if just given nbins and outer limits
                std::vector<double> vec = bins;
                if (nbins>0 && bins.size()==2) {
                    vec = {}; //NOTE: IMPORTANT!  RESET VEC IF INFERRING BINWIDTH.
                    double binwidth = (bins[1] - bins[0])/nbins;
                    for (int bin=0; bin<nbins+1; bin++) {
                        double binval = bins[0] + binwidth * bin;
                        vec.push_back(binval);
                    }
                } else if (nbins==0) { vec = bins; }
                else { std::cerr<<"*** ERROR *** COULD NOT READ BINS" << std::endl; }

                // Add to bin variables map
                binvars.insert(std::pair<std::string, std::vector<double>>(name, vec));
                std::cout << "\tINFO: "<<name<<": [ ";//DEBUGGING
                for (int bin=0; bin<vec.size(); bin++) {
                    if (bin!=vec.size()-1) { std::cout << vec[bin]<<", "; }
                    else { std::cout << vec[bin]; }
                }
                std::cout << " ]" << std::endl;

                // Add to poly4 bin variables map
                poly4map.insert(std::pair<std::string, std::vector<int>>(name, poly4bins));
                std::cout << "\tINFO: "<<name<<": [ ";//DEBUGGING
                for (int bin=0; bin<poly4bins.size(); bin++) {
                    if (bin!=poly4bins.size()-1) { std::cout << poly4bins[bin]<<", "; }
                    else { std::cout << poly4bins[bin]; }
                }
                std::cout << " ]" << std::endl;
            }
        }
    }
    double bgfraction = 1.0;
    if (node["bgfraction"]) {
        bgfraction = node["bgfraction"].as<double>();
    }
    std::cout << "INFO: bgfraction: " << bgfraction << std::endl;

    bool use_bgfraction = false;
    if (node["use_bgfraction"]) {
        use_bgfraction = node["use_bgfraction"].as<bool>();
    }
    std::cout << "INFO: use_bgfraction: " << use_bgfraction << std::endl;

    int seed = 2;
    if (node["inject_seed"]) {
        seed = node["inject_seed"].as<int>();
    }
    std::cout << "INFO: inject_seed: " << seed << std::endl;

    bool inject_asym = false;
    if (node["inject_asym"]) {
        inject_asym = node["inject_asym"].as<bool>();
    }
    std::cout << "INFO: inject_asym: " << inject_asym << std::endl;

    std::vector<double> sgasyms;
    if (node["sgasyms"]) {
        sgasyms = node["sgasyms"].as<std::vector<double>>();
    }
    std::cout << "INFO: sgasyms: [ ";
    for (int idx=0; idx<sgasyms.size(); idx++) {
        if (idx!=sgasyms.size()-1) { std::cout << sgasyms[idx]<<", "; }
        else { std::cout << sgasyms[idx]; }
    }
    std::cout << " ]" << std::endl;

    std::vector<double> bgasyms;
    if (node["bgasyms"]) {
        bgasyms = node["bgasyms"].as<std::vector<double>>();
    }
    std::cout << "INFO: bgasyms: [ ";
    for (int idx=0; idx<bgasyms.size(); idx++) {
        if (idx!=bgasyms.size()-1) { std::cout << bgasyms[idx]<<", "; }
        else { std::cout << bgasyms[idx]; }
    }
    std::cout << " ]" << std::endl;

    double beam_polarization   = 0.8922; // Average Polarization for Fall 2018 Outbending data runs >= 5331
    if (node["beam_polarization"]) {
        beam_polarization = node["beam_polarization"].as<double>();
    }
    std::cout << "INFO: beam_polarization: " << beam_polarization << std::endl;

    std::string mass_name = "mass_ppim";
    if (node["mass_name"]) {
        mass_name = node["mass_name"].as<std::string>();   
    }
    std::cout << "INFO: mass_name: " << mass_name << std::endl;

    int n_mass_bins = 100;
    if (node["n_mass_bins"]) {
        n_mass_bins = node["n_mass_bins"].as<int>();
    }
    std::cout << "INFO: n_mass_bins: " << n_mass_bins << std::endl;

    double mass_min = 1.08;
    if (node["mass_min"]) {
        mass_min = node["mass_min"].as<double>();
    }
    std::cout << "INFO: mass_min: " << mass_min << std::endl;

    double mass_max = 1.24;
    if (node["mass_max"]) {
        mass_max = node["mass_max"].as<double>();
    }
    std::cout << "INFO: mass_max: " << mass_max << std::endl;

    std::string mass_draw_opt = "";
    if (node["mass_draw_opt"]) {
        mass_draw_opt = node["mass_draw_opt"].as<std::string>();
    }
    std::cout << "INFO: mass_draw_opt: " << mass_draw_opt << std::endl;

    std::string graph_title = "graph_title";
    if (node["graph_title"]) {
        graph_title = node["graph_title"].as<std::string>();
    }
    std::cout << "INFO: graph_title: " << graph_title << std::endl;

    int marker_color = 4;
    if (node["marker_color"]) {
        marker_color = node["marker_color"].as<int>();
    }
    std::cout << "INFO: marker_color: " << marker_color << std::endl;

    int marker_style = 20;
    if (node["marker_style"]) {
        marker_style = node["marker_style"].as<int>(); 
    }
    std::cout << "INFO: marker_style: " << marker_style << std::endl;

    std::string logpath = "out.txt";
    if (node["logpath"]) {
        logpath = node["logpath"].as<std::string>();
    }
    std::cout << "INFO: logpath: " << logpath << std::endl;

    int n_fitvar1_bins = 10;
    if (node["n_fitvar1_bins"]) {
        n_fitvar1_bins = node["n_fitvar1_bins"].as<int>();
    }
    std::cout << "INFO: n_fitvar1_bins: " << n_fitvar1_bins << std::endl;

    double fitvar1_min = 0.0;
    if (node["fitvar1_min"]) {
        fitvar1_min = node["fitvar1_min"].as<double>();
    }
    std::cout << "INFO: fitvar1_min: " << fitvar1_min << std::endl;

    double fitvar1_max = 2*TMath::Pi();
    if (node["fitvar1_max"]) {
        fitvar1_max = node["fitvar1_max"].as<double>();
    }
    std::cout << "INFO: fitvar1_max: " << fitvar1_max << std::endl;

    bool use_sumW2Error = true;
    if (node["use_sumW2Error"]) {
        use_sumW2Error = node["use_sumW2Error"].as<bool>();
    }
    std::cout << "INFO: use_sumW2Error: " << use_sumW2Error << std::endl;

    bool use_splot = true;
    if (node["use_splot"]) {
        use_splot = node["use_splot"].as<bool>();
    }
    std::cout << "INFO: use_splot: " << use_splot << std::endl;

    std::string helicity_name = "heli";
    if (node["helicity_name"]) {
        helicity_name = node["helicity_name"].as<std::string>();
    }
    std::cout << "INFO: helicity_name: " << helicity_name << std::endl;

    std::string helicity_formula = "-helicity"; //NOTE: Make sure to flip helicity for RGA fall 2018 data and check if needed for other datasets.
    if (node["helicity_formula"]) {
        helicity_formula = node["helicity_formula"].as<std::string>();
    }
    std::cout << "INFO: helicity_formula: " << helicity_formula << std::endl;

    std::string mc_cuts = "Q2>1"; //NOTE: This may not be empty!
    if (node["mc_cuts"]) {
        mc_cuts = node["mc_cuts"].as<std::string>();
    }
    std::cout << "INFO: mc_cuts: " << mc_cuts << std::endl;

    // Get particle suffixes for MC matching branches
    std::vector<std::string> particle_suffixes;
    if (node["particle_suffixes"]) {
        particle_suffixes = node["particle_suffixes"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: particle_suffixes: [ ";
    for (int idx=0; idx<particle_suffixes.size(); idx++) {
        if (idx!=particle_suffixes.size()-1) { std::cout << particle_suffixes[idx]<<", "; }
        else { std::cout << particle_suffixes[idx]; }
    }
    std::cout << " ]" << std::endl;

    std::string mc_sg_match_name = "mc_sg_match"; //NOTE: This may not be empty!
    if (node["mc_sg_match_name"]) {
        mc_sg_match_name = node["mc_sg_match_name"].as<std::string>();
    }
    std::cout << "INFO: mc_sg_match_name: " << mc_sg_match_name << std::endl;

    std::string mc_sg_match_formula = "ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc && p"; //NOTE: This may not be empty!
    if (node["mc_sg_match_formula"]) {
        mc_sg_match_formula = node["mc_sg_match_formula"].as<std::string>();
    }
    std::cout << "INFO: mc_sg_match_formula: " << mc_sg_match_formula << std::endl;

    std::string fsgasyms_xs_name = "fsgasyms_xs"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_name"]) {
        fsgasyms_xs_name = node["fsgasyms_xs_name"].as<std::string>();
    }
    std::cout << "INFO: fsgasyms_xs_name: " << fsgasyms_xs_name << std::endl;

    std::string fsgasyms_xs_formula = "0.747*depol0mc*sgasym0*fitvar1_mc"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_formula"]) {
        fsgasyms_xs_formula = node["fsgasyms_xs_formula"].as<std::string>();
    }
    std::cout << "INFO: fsgasyms_xs_formula: " << fsgasyms_xs_formula << std::endl;

    std::string fbgasyms_xs_name = "fbgasyms_xs"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_name"]) {
        fbgasyms_xs_name = node["fbgasyms_xs_name"].as<std::string>();
    }
    std::cout << "INFO: fbgasyms_xs_name: " << fbgasyms_xs_name << std::endl;

    std::string fbgasyms_xs_formula = "(float)1.0"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_formula"]) {
        fbgasyms_xs_formula = node["fbgasyms_xs_formula"].as<std::string>();
    }
    std::cout << "INFO: fbgasyms_xs_formula: " << fbgasyms_xs_formula << std::endl;

    std::string randvar_name = "randvar"; //NOTE: This may not be empty!
    if (node["randvar_name"]) {
        randvar_name = node["randvar_name"].as<std::string>();
    }
    std::cout << "INFO: randvar_name: " << randvar_name << std::endl;

    std::string xs_name = "XS"; //NOTE: This may not be empty!
    if (node["xs_name"]) {
        xs_name = node["xs_name"].as<std::string>();
    }
    std::cout << "INFO: xs_name: " << xs_name << std::endl;

    // Additional parameters for getKinBinnedAsym1D below
    double sg_region_min = 1.11;
    if (node["sg_region_min"]) {
        sg_region_min = node["sg_region_min"].as<double>();
    }
    std::cout << "INFO: sg_region_min: " << sg_region_min << std::endl;

    double sg_region_max = 1.13;
    if (node["sg_region_max"]) {
        sg_region_max = node["sg_region_max"].as<double>();
    }
    std::cout << "INFO: sg_region_max: " << sg_region_max << std::endl;
    
    std::string sgcut = Form("%s>%.8f && %s<%.8f",mass_name.c_str(),sg_region_min,mass_name.c_str(),sg_region_max); //NOTE: THIS DEFAULT NEEDS TO OCCUR AFTER mass_name IS SET!
    if (node["sgcut"]) {
        sgcut = node["sgcut"].as<std::string>();
    }
    std::cout << "INFO: sgcut: " << sgcut << std::endl;

    std::string bgcut = "(mass_ppim>1.08 && mass_ppim<1.11)  || (mass_ppim>1.15 && mass_ppim<1.18)";
    if (node["bgcut"]) {
        bgcut = node["bgcut"].as<std::string>();
    }
    std::cout << "INFO: bgcut: " << bgcut << std::endl;

    std::string sig_pdf_name = "cb";
    if (node["sig_pdf_name"]) {
        sig_pdf_name = node["sig_pdf_name"].as<std::string>();
    }
    std::cout << "INFO: sig_pdf_name: " << sig_pdf_name << std::endl; //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")

    int mass_nbins_hist = 100;
    if (node["mass_nbins_hist"]) {
        mass_nbins_hist = node["mass_nbins_hist"].as<int>();
    }
    std::cout << "INFO: mass_nbins_hist: " << mass_nbins_hist << std::endl;

    int mass_nbins_conv = 1000;
    if (node["mass_nbins_conv"]) {
        mass_nbins_conv = node["mass_nbins_conv"].as<int>();
    }
    std::cout << "INFO: mass_nbins_conv: " << mass_nbins_conv << std::endl;

    bool use_sb_subtraction = false;
    if (node["use_sb_subtraction"]) {
        use_sb_subtraction = node["use_sb_subtraction"].as<bool>();
    }
    std::cout << "INFO: use_sb_subtraction: " << use_sb_subtraction << std::endl;

    bool use_average_depol = false;
    if (node["use_average_depol"]) {
        use_average_depol = node["use_average_depol"].as<bool>();
    }
    std::cout << "INFO: use_average_depol: " << use_average_depol << std::endl;

    bool use_extended_nll = false;
    if (node["use_extended_nll"]) {
        use_extended_nll = node["use_extended_nll"].as<bool>();
    }
    std::cout << "INFO: use_extended_nll: " << use_extended_nll << std::endl;

    //----------------------------------------------------------------------------------------------------//
    // ANALYSIS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Create random number generator for MC asymmetry injection
    TRandom *gRandom = new TRandom(seed); //NOTE: IMPORTANT: Need `new` here to get a pointer.

    // Add all absolute bin limits to overall cuts
    std::string binlims_cuts = "";
    for (auto it = binvars.begin(); it != binvars.end(); ++it) { //TODO: How to check if too many binning variables...

        // Get bin variable name and bin limits
        std::string binvar = it->first;
        std::vector<double> bins_ = it->second;
        double binmin = bins_.at(0);
        double binmax = bins_.at(bins_.size()-1);

        // Add to bin limits cuts
        if (binlims_cuts.size()>0) {
            binlims_cuts = Form("%s && %s>=%.16f && %s<%.16f",binlims_cuts.c_str(),binvar.c_str(),binmin,binvar.c_str(),binmax);
        } else {
            binlims_cuts = Form("%s>=%.16f && %s<%.16f",binvar.c_str(),binmin,binvar.c_str(),binmax);
        }

    } // for (auto it = binvars.begin(); it != binvars.end(); ++it) {
    std::cout << "INFO: binlims_cuts = "<<binlims_cuts.c_str() << std::endl;

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Create branch names for mc variables assuming they all append `_mc` to the corresponding data branch name
    std::string fitvar1_mc = Form("%s_mc",fitvar1.c_str());
    std::string gammavar_mc = Form("%s_mc",gammavar.c_str());
    std::string epsilonvar_mc = Form("%s_mc",epsilonvar.c_str());
    std::vector<std::string> depolvars_mc;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvars_mc.push_back(Form("%s_mc",depolvars[idx].c_str()));
    }
    std::cout << "INFO: Defined MC variable : " << fitvar1_mc.c_str() << std::endl;
    std::cout << "INFO: Defined MC variable : " << gammavar_mc.c_str() << std::endl;
    std::cout << "INFO: Defined MC variable : " << epsilonvar_mc.c_str() << std::endl;
    for (int idx=0; idx<depolvars.size(); idx++) {
        std::cout << "INFO: Defined MC variable : " << depolvars_mc[idx].c_str() << std::endl;
    }

    if (inject_asym) {
        // Find and replace asymmetry names with injected values in fsgasyms_xs : example string fsgasyms_xs="0.747*depolvars_mc0*sgasym0*fitvar1_mc"
        for (int idx=0; idx<depolvars.size(); idx++) {
            replace_all(fsgasyms_xs_formula, Form("depolvars_mc%d",idx), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
        }
        replace_all(fsgasyms_xs_formula, "fitvar1_mc", fitvar1_mc.c_str()); // Replace fitvar1_mc
        for (int idx=0; idx<sgasyms.size(); idx++) {
            replace_all(fsgasyms_xs_formula, Form("sgasyms%d",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
        }
        std::cout << "INFO: Updated " << fsgasyms_xs_name.c_str() << " = " << fsgasyms_xs_formula.c_str() << std::endl;

        // Find and replace placeholder variable names with actual values in fbgasyms_xs : example string fbgasyms_xs="0.747*depolvars_mc0*bgasym0*fitvar1_mc"
        for (int idx=0; idx<depolvars.size(); idx++) {
            replace_all(fbgasyms_xs_formula, Form("depolvars_mc%d",idx), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
        }
        replace_all(fbgasyms_xs_formula, "fitvar1_mc", fitvar1_mc.c_str()); // Replace fitvar1_mc
        for (int idx=0; idx<bgasyms.size(); idx++) {
            replace_all(fbgasyms_xs_formula, Form("bgasyms%d",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
        }
        std::cout << "INFO: Updated " << fbgasyms_xs_name.c_str() << " = " << fbgasyms_xs_formula.c_str() << std::endl;
    }

    // Pre-define depolarization variables.
    auto d2 = (!inject_asym) ? d
                .Define(gammavar.c_str(),gammavarformula.c_str())
                .Define(epsilonvar.c_str(),epsilonvarformula.c_str())
                .Define(depolvars[0].c_str(),depolvarformulas[0].c_str()) :
                d
                .Define(gammavar.c_str(),gammavarformula.c_str())
                .Define(gammavar_mc.c_str(),gammavarformulamc.c_str())
                .Define(epsilonvar.c_str(),epsilonvarformula.c_str())
                .Define(epsilonvar_mc.c_str(),epsilonvarformulamc.c_str())
                .Define(depolvars[0].c_str(),depolvarformulas[0].c_str())
                .Define(depolvars_mc[0].c_str(),depolvarformulasmc[0].c_str()); 
    for (int idx=1; idx<depolvars.size(); idx++) { //NOTE: START AT 1 HERE BECAUSE FIRST DEPOLARIZATION VARIABLE IS DEFINED ABOVE.
        d2 = (!inject_asym) ? d2.
                                Define(depolvars[idx].c_str(),depolvarformulas[idx].c_str()) :
                                d2
                                .Define(depolvars[idx].c_str(),depolvarformulas[idx].c_str())
                                .Define(depolvars_mc[idx].c_str(),depolvarformulasmc[idx].c_str());
    }
    //TODO: Add output message about defined branches

    // Define fit variables if fit formulas are not empty
    bool define_fitvars = (fitvar1formula.size()!=0 && (inject_asym==(fitvar1formulamc.size()!=0)));
    if (define_fitvars) {
        d2 = (!inject_asym) ? d2
                                .Define(fitvar1.c_str(),fitvar1formula.c_str()) :
                                d2
                                .Define(fitvar1.c_str(),fitvar1formula.c_str())
                                .Define(fitvar1_mc.c_str(),fitvar1formulamc.c_str());
    } // if (define_fitvars) {
    //TODO: Add output message about defined branches

    // Apply overall cuts AFTER defining depolarization and fit variables
    auto d2_filtered = (!inject_asym) ? d2.Filter(cuts.c_str()) :
                    d2.Filter(Form("(%s) && (%s)",cuts.c_str(),mc_cuts.c_str()));

    // Define MC matching variable names
    std::vector<std::string> theta_vars;
    std::vector<std::string> phi_vars;
    std::vector<std::string> theta_mc_vars;
    std::vector<std::string> phi_mc_vars;
    std::vector<std::string> dtheta_vars;
    std::vector<std::string> dphi_vars;
    if (inject_asym) {
        for (int idx=0; idx<particle_suffixes.size(); idx++) {
            theta_vars.push_back(Form("theta%s",particle_suffixes[idx].c_str()));
            phi_vars.push_back(Form("phi%s",particle_suffixes[idx].c_str()));
            theta_mc_vars.push_back(Form("theta%s_mc",particle_suffixes[idx].c_str()));
            phi_mc_vars.push_back(Form("phi%s_mc",particle_suffixes[idx].c_str()));
            dtheta_vars.push_back(Form("dtheta%s",particle_suffixes[idx].c_str()));
            dphi_vars.push_back(Form("dphi%s",particle_suffixes[idx].c_str()));
        }
    }

    // Define MC matching angular difference variable branches
    if (inject_asym) {
        for (int idx=0; idx<particle_suffixes.size(); idx++) {
            d2_filtered = d2_filtered.Define(dtheta_vars[idx].c_str(),[](float theta, float theta_mc){ return TMath::Abs(theta-theta_mc); },{theta_vars[idx].c_str(),theta_mc_vars[idx].c_str()})
                .Define(dphi_vars[idx].c_str(),[](float phi, float phi_mc){
                    return (float) (TMath::Abs(phi-phi_mc)<TMath::Pi()
                    ? TMath::Abs(phi-phi_mc) : 2*TMath::Pi() - TMath::Abs(phi-phi_mc));
                    },{phi_vars[idx].c_str(),phi_mc_vars[idx].c_str()});
        }
    }
    //TODO: Add output message about defined branches

    // Define signal matching condition, and XS values branches
    if (inject_asym) {
        d2_filtered = d2_filtered.Define(mc_sg_match_name.c_str(),mc_sg_match_formula.c_str()); //TODO: Throw error if formulas are empty
        d2_filtered = d2_filtered.Define(fsgasyms_xs_name.c_str(),fsgasyms_xs_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_name.c_str(),fbgasyms_xs_formula.c_str());
    }
    //TODO: Add output message about defined branches

    // Define helicity variable, injecting and applying MC matching cuts if requested
    auto frame = (!inject_asym) ? d2_filtered.Define(helicity_name.c_str(), helicity_formula.c_str()) :
                    d2_filtered.Define(randvar_name.c_str(),[&gRandom](){ return (float)gRandom->Rndm(); },{})
                    .Define(xs_name.c_str(), [&beam_polarization]
                        (bool mc_sg_match, float fsgasyms_xs, float fbgasyms_xs) {
                            return (float)((mc_sg_match) ?
                            0.5*(1.0 + beam_polarization*fsgasyms_xs) :
                            0.5*(1.0 + beam_polarization*fbgasyms_xs));
                        },
                        {mc_sg_match_name.c_str(),fsgasyms_xs_name.c_str(),fbgasyms_xs_name.c_str()})
                    .Define(helicity_name.c_str(), [](float my_rand_var, float XS) {
                        return (float)(my_rand_var<XS ? 1.0 : -1.0);
                    },
                    {randvar_name.c_str(),xs_name.c_str()});
    //TODO: Add output message about defined branches

    // Make sure injection values are all computed before running analysis
    if (inject_asym) {
        double my_testvar  = (double)*frame.Mean(randvar_name.c_str());
        double my_testvar1 = (double)*frame.Mean(xs_name.c_str());
        double my_testvar2 = (double)*frame.Mean(helicity_name.c_str());
    }

    // Create output log
    std::ofstream outf; outf.open(logpath.c_str());
    std::ostream &out = outf; //std::cout;

    // Create output ROOT file
    TFile * outroot = TFile::Open(outpath.c_str(),"RECREATE");

    // Loop variables to bin in
    for (auto it = binvars.begin(); it != binvars.end(); ++it) { //TODO: How to check if too many binning variables...

        // Get bin variable name and bin limits
        std::string binvar = it->first;
        std::vector<double> bins_ = it->second;
        const int nbins = bins_.size()-1; //NOTE: IMPORTANT: -1 is because you give bin limits!
        double bins[nbins];

        // Set poly4 mask
        int poly4bins[nbins];
        for (int entry=0; entry<nbins; entry++) { poly4bins[entry] = 0; }
        for (int entry=0; entry<poly4map[binvar.c_str()].size(); entry++) {
            int bin = poly4map[binvar.c_str()][entry]-1;
            if (bin<nbins) poly4bins[bin] = 1;
        }
        for (int bin=0; bin<bins_.size(); bin++) { bins[bin] = bins_[bin]; }

        // Set binvar outdir name
        std::string binvar_outdir = Form("binvar_%s",binvar.c_str());

        // Produce graphs of asymmetry fit parameters corrected for depolarization and background binned in given kinematic variable
        analysis::getKinBinnedAsym1D(
            outdir, //std::string outdir,
            outroot, //TFile      *outroot,
            frame, //ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            method, // std::string     method, // ONLY getKinBinAsymUBML1D ('BSA1D') is allowed at the moment
            binvar, // std::string binvar, // Variable name to bin in
            nbins, // int         nbins, // Number of bins
            bins, // double      *bins, // Bin limits (length=nbins+1)
            beam_polarization, //double      pol,
            depolvars, //std::vector<std::string> depolvars,
            "w", //std::string workspace_name  = "w",
            "workspace", //std::string workspace_title = "workspace",
            "dataset", //std::string dataset_name    = "dataset",
            "dataset", //std::string dataset_title   = "dataset",
            helicity_name, //std::string helicity_name   = "heli",
            fitvar1, //std::string fitvarx         = "x",
            fitvar1_min, //double      xmin            = 0.0,
            fitvar1_max, //double      xmax            = 1.0,
            mass_name, //std::string massvar         = "mass_ppim",
            mass_min, //double      mmin            = 1.08,
            mass_max, //double      mmax            = 1.24,
            "sgYield", //std::string sgYield_name    = "sgYield",
            "bgYield", //std::string bgYield_name    = "bgYield",
            "model", //std::string model_name      = "model",
            fitformula, //std::string fitformula      = "[0]*sin(x)+[1]*sin(2*x)",
            nparams, //int         nparams         = 2,
            params, //std::vector<double> params  = std::vector<double>(5),
            fitvar1title, //std::string fitvarxtitle    = "#phi_{h p#pi^{-}}",
            fitvar1bins, //int         xbins           = 16,
            use_sumW2Error, //bool        use_sumW2Error  = true,
            use_average_depol, //bool use_average_depol      = false,
            use_extended_nll, // bool use_extended_nll       = false,
            use_splot, //bool        use_splot       = true,
            graph_title, //std::string graph_title     = "BSA A_{LU} vs. #Delta#phi", // Histogram title
            marker_color, // int         marker_color    = 4, // 4 is blue
            marker_style, // int         marker_style    = 20, // 20 is circle
            out, //std::ostream &out           = std::cout
            sgcut, //std::string sgcut           = "Q2>1",
            bgcut, //std::string bgcut           = "Q2>1",
            mass_nbins_hist, //int mass_nbins_hist         = 100,
            mass_nbins_conv, //int mass_nbins_conv         = 1000,
            sig_pdf_name, //std::string sig_pdf_name    = "cb", //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")
            sg_region_min, //double sg_region_min        = 1.11,
            sg_region_max, //double sg_region_max        = 1.13,
            use_sb_subtraction //bool   use_sb_subtraction   = false
        );
   
    }// loop over bin variables

} // void execute()

int main(int argc, char** argv) {

    // Check # of arguments
    if (argc==1) {
        std::cout << "Usage: " << argv[0] << " /path/to/config.yaml" << std::endl;
        return 0;
    }

    // Load YAML file
    const char * yamlfile = argv[1];
    YAML::Node config = YAML::LoadFile(yamlfile);

    // Process arguments
    execute(config);

    return 0;

} // int main()
