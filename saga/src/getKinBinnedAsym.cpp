#include <iostream>
#include <memory>
#include <fstream>
#include <string>
#include <map>

// YAML Includes
#include <yaml-cpp/yaml.h>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>
#include <TRandom.h>
#include <TMath.h>

// Project Includes
#include <analysis.h>
#include <bins.h>
#include <data.h>
#include <util.h>

void execute(const YAML::Node& node) {

    // Process arguments

    // OUTPATH
    std::string baseoutpath = "";//NOTE: This will be prepended to the default output path like so: `<baseoutpath><binscheme_name>.csv`.
    if (node["baseoutpath"]) {
        baseoutpath = node["baseoutpath"].as<std::string>();
    }
    std::cout << "INFO: baseoutpath: " << baseoutpath << std::endl;

    // INPATH
    std::string inpath = "";
    if (node["inpath"]) {
        inpath = node["inpath"].as<std::string>();
    }
    std::cout << "INFO: inpath: " << inpath << std::endl;

    // TREE
    std::string tree = "";
    if (node["tree"]) {
        tree = node["tree"].as<std::string>();
    }
    std::cout << "INFO: tree: " << tree << std::endl;

    // NTHREADS
    int nthreads = 1;
    if (node["nthreads"]) {
        nthreads = node["nthreads"].as<int>();
    }
    std::cout << "INFO: nthreads: " << nthreads << std::endl;

    // CUTS
    std::string cuts = "";
    if (node["cuts"]) {
        cuts = node["cuts"].as<std::string>();
    }
    std::cout << "INFO: cuts: " << cuts << std::endl;

    //----------------------------------------------------------------------//
    // BEGIN ASYMMETRY INJECTION ARGUMENTS

    // INJECT_ASYM
    bool inject_asym = false;
    if (node["inject_asym"]) {
        inject_asym = node["inject_asym"].as<bool>();
    }
    std::cout << "INFO: inject_asym: " << inject_asym << std::endl;

    // SEED
    int seed = 2;
    if (node["inject_seed"]) {
        seed = node["inject_seed"].as<int>();
    }
    std::cout << "INFO: inject_seed: " << seed << std::endl;

    // MC_CUTS
    std::string mc_cuts = "Q2>1"; //NOTE: This may not be empty!
    if (node["mc_cuts"]) {
        mc_cuts = node["mc_cuts"].as<std::string>();
    }
    std::cout << "INFO: mc_cuts: " << mc_cuts << std::endl;

    // SGASYMS
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

    // BGASYMS
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

    // PARTICLE_SUFFIXES -> For MC matching with mc_sg_match cut
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

    // MC_SG_MATCH
    std::string mc_sg_match_name = "mc_sg_match"; //NOTE: This may not be empty!
    if (node["mc_sg_match_name"]) {
        mc_sg_match_name = node["mc_sg_match_name"].as<std::string>();
    }
    std::cout << "INFO: mc_sg_match_name: " << mc_sg_match_name << std::endl;

    // MC_SG_MATCH_FORMULA
    std::string mc_sg_match_formula = "ppid_p_mc==3122 && pidx_p_mc==pidx_pim_mc"; //NOTE: This may not be empty!
    if (node["mc_sg_match_formula"]) {
        mc_sg_match_formula = node["mc_sg_match_formula"].as<std::string>();
    }
    std::cout << "INFO: mc_sg_match_formula: " << mc_sg_match_formula << std::endl;

    // FSGASYMS_XS_PU_NAME
    std::string fsgasyms_xs_pu_name = "fsgasyms_xs_pu"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_pu_name"]) {
        fsgasyms_xs_pu_name = node["fsgasyms_xs_pu_name"].as<std::string>();
    }
    std::cout << "INFO: fsgasyms_xs_pu_name: " << fsgasyms_xs_pu_name << std::endl;

    // FSGASYMS_XS_PU_FORMULA
    std::string fsgasyms_xs_pu_formula = "(float)0.0"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_pu_formula"]) {
        fsgasyms_xs_pu_formula = node["fsgasyms_xs_pu_formula"].as<std::string>();
        if (fsgasyms_xs_pu_formula=="") { fsgasyms_xs_pu_formula = "(float)0.0"; }
    }
    std::cout << "INFO: fsgasyms_xs_pu_formula: " << fsgasyms_xs_pu_formula << std::endl;

    // FSGASYMS_XS_UP_NAME
    std::string fsgasyms_xs_up_name = "fsgasyms_xs_up"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_up_name"]) {
        fsgasyms_xs_up_name = node["fsgasyms_xs_up_name"].as<std::string>();
    }
    std::cout << "INFO: fsgasyms_xs_up_name: " << fsgasyms_xs_up_name << std::endl;

    // FSGASYMS_XS_UP_FORMULA
    std::string fsgasyms_xs_up_formula = "(float)0.0"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_up_formula"]) {
        fsgasyms_xs_up_formula = node["fsgasyms_xs_up_formula"].as<std::string>();
        if (fsgasyms_xs_up_formula=="") { fsgasyms_xs_up_formula = "(float)0.0"; }
    }
    std::cout << "INFO: fsgasyms_xs_up_formula: " << fsgasyms_xs_up_formula << std::endl;

    // FSGASYMS_XS_PP_NAME
    std::string fsgasyms_xs_pp_name = "fsgasyms_xs_pp"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_pp_name"]) {
        fsgasyms_xs_pp_name = node["fsgasyms_xs_pp_name"].as<std::string>();
    }
    std::cout << "INFO: fsgasyms_xs_pp_name: " << fsgasyms_xs_pp_name << std::endl;

    // FSGASYMS_XS_PP_FORMULA
    std::string fsgasyms_xs_pp_formula = "(float)0.0"; //NOTE: This may not be empty!
    if (node["fsgasyms_xs_pp_formula"]) {
        fsgasyms_xs_pp_formula = node["fsgasyms_xs_pp_formula"].as<std::string>();
        if (fsgasyms_xs_pp_formula=="") { fsgasyms_xs_pp_formula = "(float)0.0"; }
    }
    std::cout << "INFO: fsgasyms_xs_pp_formula: " << fsgasyms_xs_pp_formula << std::endl;

    // FBGASYMS_XS_PU_NAME
    std::string fbgasyms_xs_pu_name = "fbgasyms_xs_pu"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_pu_name"]) {
        fbgasyms_xs_pu_name = node["fbgasyms_xs_pu_name"].as<std::string>();
    }
    std::cout << "INFO: fbgasyms_xs_pu_name: " << fbgasyms_xs_pu_name << std::endl;

    // FBGASYMS_XS_PU_FORMULA
    std::string fbgasyms_xs_pu_formula = "(float)0.0"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_pu_formula"]) {
        fbgasyms_xs_pu_formula = node["fbgasyms_xs_pu_formula"].as<std::string>();
        if (fbgasyms_xs_pu_formula=="") { fbgasyms_xs_pu_formula = "(float)0.0"; }
    }
    std::cout << "INFO: fbgasyms_xs_pu_formula: " << fbgasyms_xs_pu_formula << std::endl;

    // FBGASYMS_XS_UP_NAME
    std::string fbgasyms_xs_up_name = "fbgasyms_xs_up"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_up_name"]) {
        fbgasyms_xs_up_name = node["fbgasyms_xs_up_name"].as<std::string>();
    }
    std::cout << "INFO: fbgasyms_xs_up_name: " << fbgasyms_xs_up_name << std::endl;

    // FBGASYMS_XS_UP_FORMULA
    std::string fbgasyms_xs_up_formula = "(float)0.0"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_up_formula"]) {
        fbgasyms_xs_up_formula = node["fbgasyms_xs_up_formula"].as<std::string>();
        if (fbgasyms_xs_up_formula=="") { fbgasyms_xs_up_formula = "(float)0.0"; }
    }
    std::cout << "INFO: fbgasyms_xs_up_formula: " << fbgasyms_xs_up_formula << std::endl;

    // FBGASYMS_XS_PP_NAME
    std::string fbgasyms_xs_pp_name = "fbgasyms_xs_pp"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_pp_name"]) {
        fbgasyms_xs_pp_name = node["fbgasyms_xs_pp_name"].as<std::string>();
    }
    std::cout << "INFO: fbgasyms_xs_pp_name: " << fbgasyms_xs_pp_name << std::endl;

    // FBGASYMS_XS_PP_FORMULA
    std::string fbgasyms_xs_pp_formula = "(float)0.0"; //NOTE: This may not be empty!
    if (node["fbgasyms_xs_pp_formula"]) {
        fbgasyms_xs_pp_formula = node["fbgasyms_xs_pp_formula"].as<std::string>();
        if (fbgasyms_xs_pp_formula=="") { fbgasyms_xs_pp_formula = "(float)0.0"; }
    }
    std::cout << "INFO: fbgasyms_xs_pp_formula: " << fbgasyms_xs_pp_formula << std::endl;

    // COMBINED_SPIN_STATE
    std::string combined_spin_state = "ss"; //NOTE: This may not be empty!
    if (node["combined_spin_state"]) {
        combined_spin_state = node["combined_spin_state"].as<std::string>();
    }
    std::cout << "INFO: combined_spin_state: " << combined_spin_state << std::endl;

    // END MC ASYMMETRY INJECTION ARGUMENTS
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN RUN-DEPENDENT CSV VARIABLE ARGUMENTS

    // RDF_KEY_COLS
    std::vector<std::string> rdf_key_cols;
    if (node["rdf_key_cols"]) {
        rdf_key_cols = node["rdf_key_cols"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: rdf_key_cols: [ ";
    for (int idx=0; idx<rdf_key_cols.size(); idx++) {
        if (idx!=rdf_key_cols.size()-1) { std::cout << rdf_key_cols[idx]<<", "; }
        else { std::cout << rdf_key_cols[idx]; }
    }
    std::cout << " ]" << std::endl;

    // CSV_PATHS
    std::vector<std::string> csv_paths;
    if (node["csv_paths"]) {
        csv_paths = node["csv_paths"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: csv_paths: [ ";
    for (int idx=0; idx<csv_paths.size(); idx++) {
        if (idx!=csv_paths.size()-1) { std::cout << csv_paths[idx]<<", "; }
        else { std::cout << csv_paths[idx]; }
    }
    std::cout << " ]" << std::endl;

    // CSV_KEY_COLS
    std::vector<std::string> csv_key_cols;
    if (node["csv_key_cols"]) {
        csv_key_cols = node["csv_key_cols"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: csv_key_cols: [ ";
    for (int idx=0; idx<csv_key_cols.size(); idx++) {
        if (idx!=csv_key_cols.size()-1) { std::cout << csv_key_cols[idx]<<", "; }
        else { std::cout << csv_key_cols[idx]; }
    }
    std::cout << " ]" << std::endl;

    // COL_NAMES
    std::vector<std::vector<std::string>> col_names;
    if (node["col_names"]) {
        col_names = node["col_names"].as<std::vector<std::vector<std::string>>>();
    }
    std::cout << "INFO: col_names: [ \n";
    for (int idx=0; idx<col_names.size(); idx++) {
        std::cout << "\t[ ";
        for (int j=0; j<col_names[idx].size()-1; j++) { std::cout << col_names[idx][j].c_str() << " , "; }
        std::cout << col_names[idx][col_names[idx].size()-1].c_str() << " ],\n";
    }
    std::cout << " ]" << std::endl;

    // COL_ALIASES
    std::map<std::string,std::string> col_aliases;
    if (node["col_aliases"]) {
        col_aliases = node["col_aliases"].as<std::map<std::string,std::string>>();
    }
    std::cout << "INFO: col_aliases: { ";
    for (auto it = col_aliases.begin(); it != col_aliases.end(); ++it) {
        std::cout << it->first.c_str()<<" : "<<it->second.c_str()<<", ";
    }
    std::cout << " }" << std::endl;

    // END RUN-DEPENDENT CSV VARIABLE ARGUMENTS
    //----------------------------------------------------------------------//

    // VAR_FORMULAS
    std::vector<std::vector<std::string>> var_formulas;
    if (node["var_formulas"]) {
        var_formulas = node["var_formulas"].as<std::vector<std::vector<std::string>>>();
    }
    std::cout << "INFO: var_formulas: { \n";
    for (int idx=0; idx<var_formulas.size(); idx++) {
        std::cout << "\t[ " << var_formulas[idx][0].c_str() << " : " << var_formulas[idx][1].c_str() << " ],\n";
    }
    std::cout << " }" << std::endl;

    //----------------------------------------------------------------------//
    // BEGIN HELICITY AND SPIN VARIABLES
    
    // HELICITY_NAME
    std::string helicity_name = "heli";
    if (node["helicity_name"]) {
        helicity_name = node["helicity_name"].as<std::string>();
    }
    std::cout << "INFO: helicity_name: " << helicity_name << std::endl;

    // HELICITY_FORMULA
    std::string helicity_formula = "-helicity"; //NOTE: Make sure to flip helicity for RGA fall 2018 data and check if needed for other datasets.
    if (node["helicity_formula"]) {
        helicity_formula = node["helicity_formula"].as<std::string>();
    }
    std::cout << "INFO: helicity_formula: " << helicity_formula << std::endl;

    // HELICITY_STATES
    std::map<std::string,int> helicity_states = {{"plus",1}, {"zero",0}, {"minus",-1}};
    if (node["helicity_states"]) {
        helicity_states = node["helicity_states"].as<std::map<std::string,int>>();
    }
    std::cout << "INFO: helicity_states: { ";
    for (auto it = helicity_states.begin(); it != helicity_states.end(); ++it) {
        std::cout << it->first<<" : "<<it->second<<", ";
    }
    std::cout << " }" << std::endl;

    // TSPIN_NAME
    std::string tspin_name = "heli";
    if (node["tspin_name"]) {
        tspin_name = node["tspin_name"].as<std::string>();
    }
    std::cout << "INFO: tspin_name: " << tspin_name << std::endl;

    // TSPIN_FORMULA
    std::string tspin_formula = "-tspin"; //NOTE: Make sure to flip tspin for RGA fall 2018 data and check if needed for other datasets.
    if (node["tspin_formula"]) {
        tspin_formula = node["tspin_formula"].as<std::string>();
    }
    std::cout << "INFO: tspin_formula: " << tspin_formula << std::endl;

    // TSPIN_STATES
    std::map<std::string,int> tspin_states = {{"plus",1}, {"zero",0}, {"minus",-1}};
    if (node["tspin_states"]) {
        tspin_states = node["tspin_states"].as<std::map<std::string,int>>();
    }
    std::cout << "INFO: tspin_states: { ";
    for (auto it = tspin_states.begin(); it != tspin_states.end(); ++it) {
        std::cout << it->first<<" : "<<it->second<<", ";
    }
    std::cout << " }" << std::endl;

    // HTSPIN_NAME
    std::string htspin_name = "heli";
    if (node["htspin_name"]) {
        htspin_name = node["htspin_name"].as<std::string>();
    }
    std::cout << "INFO: htspin_name: " << htspin_name << std::endl;

    // HTSPIN_STATES
    std::map<std::string,int> htspin_states = {{"plus",1}, {"zero",0}, {"minus",-1}};
    if (node["htspin_states"]) {
        htspin_states = node["htspin_states"].as<std::map<std::string,int>>();
    }
    std::cout << "INFO: htspin_states: { ";
    for (auto it = htspin_states.begin(); it != htspin_states.end(); ++it) {
        std::cout << it->first<<" : "<<it->second<<", ";
    }
    std::cout << " }" << std::endl;


    // END HELICITY AND SPIN VARIABLES
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN BINNING SCHEME ARGUMENTS

    // BINSCHEMES
    std::map<std::string,std::map<int,std::string>> bincuts_map;
    std::map<std::string,std::vector<std::string>> binschemes_vars;
    if (node["binschemes"] && node["binschemes"].IsMap()) {

        // Get bin scheme node and get bin cuts maps
        auto node_binschemes = node["binschemes"];
        bincuts_map = saga::bins::getBinCutsMap(node["binschemes"]);

        // Get list of bin variables used in scheme
        binschemes_vars = saga::bins::getBinSchemesVars(node["binschemes"]);

    } else if (node["binschemes_paths"]) {
        
        // Get list of paths to yamls containing bin scheme definitions 
        std::vector<std::string> binschemes_paths = node["binschemes_paths"].as<std::vector<std::string>>();

        // Loop paths and add bin schemes
        for (int idx=0; idx<binschemes_paths.size(); idx++) {

            // Load YAML file and get bin cuts maps
            std::cout<<"INFO: Loading bin scheme from : "<<binschemes_paths[idx].c_str()<<std::endl;
            YAML::Node bincut_config = YAML::LoadFile(binschemes_paths[idx].c_str());
            std::map<std::string,std::map<int,std::string>> new_bincuts_map = saga::bins::getBinCutsMap(bincut_config);
            bincuts_map.insert(new_bincuts_map.begin(), new_bincuts_map.end());
            std::map<std::string,std::vector<std::string>> new_binschemes_vars = saga::bins::getBinSchemesVars(bincut_config);
            binschemes_vars.insert(new_binschemes_vars.begin(), new_binschemes_vars.end());
        }
    }

    // NBATCHES AND IBATCH
    if (node["nbatches"] && node["ibatch"]) {
        int nbatches = node["nbatches"].as<int>();
        std::cout << "INFO: nbatches: " << nbatches << std::endl;
        int ibatch   = node["ibatch"].as<int>();
        std::cout << "INFO: ibatch: " << ibatch << std::endl;

        // Reduce bin cuts map into a single batch for parallelization
        if (nbatches>1 && ibatch>=0 && ibatch<nbatches) bincuts_map = saga::bins::getBinCutsMapBatch(bincuts_map, nbatches, ibatch);
    }

    // Show bin cuts map
    for (auto it = bincuts_map.begin(); it != bincuts_map.end(); ++it) {
        std::cout << "INFO: bincuts_map["<<it->first<<"]: { ";
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            std::cout <<"\t\t"<< it2->first<<": "<<it2->second<<", \n";
        }
        std::cout << "}" << std::endl;
    }
    // END BINNING SCHEME ARGUMENTS
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN BIN VARIABLES

    // BINVARS
    std::vector<std::string> binvars;
    if (node["binvars"]) {
        binvars = node["binvars"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: binvars: [ ";
    for (int idx=0; idx<binvars.size(); idx++) {
        if (idx!=binvars.size()-1) { std::cout << binvars[idx]<<", "; }
        else { std::cout << binvars[idx]; }
    }
    std::cout << " ]" << std::endl;

    // BINVAR_TITLES
    std::vector<std::string> binvar_titles = binvars; //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    if (node["binvar_titles"]) {
        binvar_titles = node["binvar_titles"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: binvar_titles: [ ";
    for (int idx=0; idx<binvar_titles.size(); idx++) {
        if (idx!=binvar_titles.size()-1) { std::cout << binvar_titles[idx]<<", "; }
        else { std::cout << binvar_titles[idx]; }
    }
    std::cout << " ]" << std::endl;

    // BINVAR_LIMS
    std::vector<std::vector<double>> binvar_lims;
    if (node["binvar_lims"]) {
        binvar_lims = node["binvar_lims"].as<std::vector<std::vector<double>>>();
    }
    std::cout << "INFO: binvar_lims: [ ";
    for (int idx=0; idx<binvar_lims.size(); idx++) {
        if (idx!=binvar_lims.size()-1) { std::cout << "[ " << binvar_lims[idx][0] << ", " << binvar_lims[idx][1] << " ], "; }
        else { std::cout << "[ " << binvar_lims[idx][0] << ", " << binvar_lims[idx][1] << " ] "; }
    }
    std::cout << " ]" << std::endl;

    // BINVAR_BINS
    std::vector<int> binvar_bins;
    if (node["binvar_bins"]) {
        binvar_bins = node["binvar_bins"].as<std::vector<int>>();
    }
    std::cout << "INFO: binvar_bins: [ ";
    for (int idx=0; idx<binvar_bins.size(); idx++) {
        if (idx!=binvar_bins.size()-1) { std::cout << binvar_bins[idx]<<", "; }
        else { std::cout << binvar_bins[idx]; }
    }
    std::cout << " ]" << std::endl;

    // END BIN VARIABLES
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN DEPOLARIZATION VARIABLES

    // DEPOLVARS
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

    // DEPOLVAR_TITLES
    std::vector<std::string> depolvar_titles = depolvars; //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    if (node["depolvar_titles"]) {
        depolvar_titles = node["depolvar_titles"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: depolvar_titles: [ ";
    for (int idx=0; idx<depolvar_titles.size(); idx++) {
        if (idx!=depolvar_titles.size()-1) { std::cout << depolvar_titles[idx]<<", "; }
        else { std::cout << depolvar_titles[idx]; }
    }
    std::cout << " ]" << std::endl;

    // DEPOLVAR_LIMS
    std::vector<std::vector<double>> depolvar_lims;
    if (node["depolvar_lims"]) {
        depolvar_lims = node["depolvar_lims"].as<std::vector<std::vector<double>>>();
    }
    std::cout << "INFO: depolvar_lims: [ ";
    for (int idx=0; idx<depolvar_lims.size(); idx++) {
        if (idx!=depolvar_lims.size()-1) { std::cout << "[ " << depolvar_lims[idx][0] << ", " << depolvar_lims[idx][1] << " ], "; }
        else { std::cout << "[ " << depolvar_lims[idx][0] << ", " << depolvar_lims[idx][1] << " ] "; }
    }
    std::cout << " ]" << std::endl;

    // DEPOLVAR_BINS
    std::vector<int> depolvar_bins;
    if (node["depolvar_bins"]) {
        depolvar_bins = node["depolvar_bins"].as<std::vector<int>>();
    }
    std::cout << "INFO: depolvar_bins: [ ";
    for (int idx=0; idx<depolvar_bins.size(); idx++) {
        if (idx!=depolvar_bins.size()-1) { std::cout << depolvar_bins[idx]<<", "; }
        else { std::cout << depolvar_bins[idx]; }
    }
    std::cout << " ]" << std::endl;

    // END DEPOLARIZATION VARIABLES
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN ASYMMETRY FIT VARIABLES

    // ASYMFITVARS
    std::vector<std::string> asymfitvars;
    if (node["asymfitvars"]) {
        asymfitvars = node["asymfitvars"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: asymfitvars: [ ";
    for (int idx=0; idx<asymfitvars.size(); idx++) {
        if (idx!=asymfitvars.size()-1) { std::cout << asymfitvars[idx]<<", "; }
        else { std::cout << asymfitvars[idx]; }
    }
    std::cout << " ]" << std::endl;

    // ASYMFITVAR_TITLES
    std::vector<std::string> asymfitvar_titles = asymfitvars; //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    if (node["asymfitvar_titles"]) {
        asymfitvar_titles = node["asymfitvar_titles"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: asymfitvar_titles: [ ";
    for (int idx=0; idx<asymfitvar_titles.size(); idx++) {
        if (idx!=asymfitvar_titles.size()-1) { std::cout << asymfitvar_titles[idx]<<", "; }
        else { std::cout << asymfitvar_titles[idx]; }
    }
    std::cout << " ]" << std::endl;

    // ASYMFITVAR_LIMS
    std::vector<std::vector<double>> asymfitvar_lims;
    if (node["asymfitvar_lims"]) {
        asymfitvar_lims = node["asymfitvar_lims"].as<std::vector<std::vector<double>>>();
    }
    std::cout << "INFO: asymfitvar_lims: [ ";
    for (int idx=0; idx<asymfitvar_lims.size(); idx++) {
        if (idx!=asymfitvar_lims.size()-1) { std::cout << "[ " << asymfitvar_lims[idx][0] << ", " << asymfitvar_lims[idx][1] << " ], "; }
        else { std::cout << "[ " << asymfitvar_lims[idx][0] << ", " << asymfitvar_lims[idx][1] << " ] "; }
    }
    std::cout << " ]" << std::endl;

    // ASYMFITVAR_BINS
    std::vector<int> asymfitvar_bins;
    if (node["asymfitvar_bins"]) {
        asymfitvar_bins = node["asymfitvar_bins"].as<std::vector<int>>();
    }
    std::cout << "INFO: asymfitvar_bins: [ ";
    for (int idx=0; idx<asymfitvar_bins.size(); idx++) {
        if (idx!=asymfitvar_bins.size()-1) { std::cout << asymfitvar_bins[idx]<<", "; }
        else { std::cout << asymfitvar_bins[idx]; }
    }
    std::cout << " ]" << std::endl;

    // END ASYMMETRY FIT VARIABLES
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN MASS FIT VARIABLES

    // MASSFITVARS
    std::vector<std::string> massfitvars;
    if (node["massfitvars"]) {
        massfitvars = node["massfitvars"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: massfitvars: [ ";
    for (int idx=0; idx<massfitvars.size(); idx++) {
        if (idx!=massfitvars.size()-1) { std::cout << massfitvars[idx]<<", "; }
        else { std::cout << massfitvars[idx]; }
    }
    std::cout << " ]" << std::endl;

    // MASSFITVAR_TITLES
    std::vector<std::string> massfitvar_titles = massfitvars; //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    if (node["massfitvar_titles"]) {
        massfitvar_titles = node["massfitvar_titles"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: massfitvar_titles: [ ";
    for (int idx=0; idx<massfitvar_titles.size(); idx++) {
        if (idx!=massfitvar_titles.size()-1) { std::cout << massfitvar_titles[idx]<<", "; }
        else { std::cout << massfitvar_titles[idx]; }
    }
    std::cout << " ]" << std::endl;

    // MASSFITVAR_LIMS
    std::vector<std::vector<double>> massfitvar_lims;
    if (node["massfitvar_lims"]) {
        massfitvar_lims = node["massfitvar_lims"].as<std::vector<std::vector<double>>>();
    }
    std::cout << "INFO: massfitvar_lims: [ ";
    for (int idx=0; idx<massfitvar_lims.size(); idx++) {
        if (idx!=massfitvar_lims.size()-1) { std::cout << "[ " << massfitvar_lims[idx][0] << ", " << massfitvar_lims[idx][1] << " ], "; }
        else { std::cout << "[ " << massfitvar_lims[idx][0] << ", " << massfitvar_lims[idx][1] << " ] "; }
    }
    std::cout << " ]" << std::endl;

    // MASSFITVAR_BINS
    std::vector<int> massfitvar_bins;
    if (node["massfitvar_bins"]) {
        massfitvar_bins = node["massfitvar_bins"].as<std::vector<int>>();
    }
    std::cout << "INFO: massfitvar_bins: [ ";
    for (int idx=0; idx<massfitvar_bins.size(); idx++) {
        if (idx!=massfitvar_bins.size()-1) { std::cout << massfitvar_bins[idx]<<", "; }
        else { std::cout << massfitvar_bins[idx]; }
    }
    std::cout << " ]" << std::endl;

    // END MASS FIT VARIABLES
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN FIT PARAMETERS

    // BPOL
    double bpol = 0.8922; // Average Beam Polarization \overline{P_B^2} for Fall 2018 Outbending data runs >= 5331 is 0.8922
    if (node["bpol"]) {
        bpol = node["bpol"].as<double>();
    }
    std::cout << "INFO: bpol: " << bpol << std::endl;

    // TPOL
    double tpol = 1.0; // Target polarization
    if (node["tpol"]) {
        tpol = node["tpol"].as<double>();
    }
    std::cout << "INFO: tpol: " << tpol << std::endl;

    // ASYMFIT_FORMULA_PU
    std::string asymfit_formula_pu = "";
    if (node["asymfit_formula_pu"]) {
        asymfit_formula_pu = node["asymfit_formula_pu"].as<std::string>();
    }
    std::cout << "INFO: asymfit_formula_pu: " << asymfit_formula_pu << std::endl;

    // ASYMFIT_FORMULA_UP
    std::string asymfit_formula_up = "";
    if (node["asymfit_formula_up"]) {
        asymfit_formula_up = node["asymfit_formula_up"].as<std::string>();
    }
    std::cout << "INFO: asymfit_formula_up: " << asymfit_formula_up << std::endl;

    // ASYMFIT_FORMULA_PP
    std::string asymfit_formula_pp = "";
    if (node["asymfit_formula_pp"]) {
        asymfit_formula_pp = node["asymfit_formula_pp"].as<std::string>();
    }
    std::cout << "INFO: asymfit_formula_pp: " << asymfit_formula_pp << std::endl;

    // ASYMFITPAR_INITS
    std::vector<double> asymfitpar_inits;
    if (node["asymfitpar_inits"]) {
        asymfitpar_inits = node["asymfitpar_inits"].as<std::vector<double>>();
    }
    std::cout << "INFO: asymfitpar_inits: [ ";
    for (int idx=0; idx<asymfitpar_inits.size(); idx++) {
        if (idx!=asymfitpar_inits.size()-1) { std::cout << asymfitpar_inits[idx]<<", "; }
        else { std::cout << asymfitpar_inits[idx]; }
    }
    std::cout << " ]" << std::endl;

    // ASYMFITPAR_INITS
    std::vector<std::vector<double>> asymfitpar_initlims;
    if (node["asymfitpar_initlims"]) {
        asymfitpar_initlims = node["asymfitpar_initlims"].as<std::vector<std::vector<double>>>();
    }
    std::cout << "INFO: asymfitpar_initlims: [ ";
    for (int idx=0; idx<asymfitpar_initlims.size(); idx++) {
        if (idx!=asymfitpar_initlims.size()-1) { std::cout << "[ " << asymfitpar_initlims[idx][0] << ", " << asymfitpar_initlims[idx][1] << " ], "; }
        else { std::cout << "[ " << asymfitpar_initlims[idx][0] << ", " << asymfitpar_initlims[idx][1] << " ] "; }
    }
    std::cout << " ]" << std::endl;

    // USE_SUMW2ERROR
    bool use_sumw2error = true;
    if (node["use_sumw2error"]) {
        use_sumw2error = node["use_sumw2error"].as<bool>();
    }
    std::cout << "INFO: use_sumw2error: " << use_sumw2error << std::endl;

    // USE_AVERAGE_DEPOL
    bool use_average_depol = false;
    if (node["use_average_depol"]) {
        use_average_depol = node["use_average_depol"].as<bool>();
    }
    std::cout << "INFO: use_average_depol: " << use_average_depol << std::endl;

    // USE_EXTENDED_NLL
    bool use_extended_nll = false;
    if (node["use_extended_nll"]) {
        use_extended_nll = node["use_extended_nll"].as<bool>();
    }
    std::cout << "INFO: use_extended_nll: " << use_extended_nll << std::endl;

    // USE_BINNED_FIT
    bool use_binned_fit = false;
    if (node["use_binned_fit"]) {
        use_binned_fit = node["use_binned_fit"].as<bool>();
    }
    std::cout << "INFO: use_binned_fit: " << use_binned_fit << std::endl;
    // END FIT PARAMETERS
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN FIT ARGUMENTS

    // MASSFIT_MODEL_NAME
    std::string massfit_model_name = "model";//NOTE: JUST FIX FOR NOW

    // MASSFIT_NBINS_CONV
    int massfit_nbins_conv = 1000;
    if (node["massfit_nbins_conv"]) {
        massfit_nbins_conv = node["massfit_nbins_conv"].as<int>();
    }
    std::cout << "INFO: massfit_nbins_conv: " << massfit_nbins_conv << std::endl;

    // MASSFIT_SIG_PDF_NAME
    std::string massfit_sig_pdf_name = "cb";
    if (node["massfit_sig_pdf_name"]) {
        massfit_sig_pdf_name = node["massfit_sig_pdf_name"].as<std::string>();
    }
    std::cout << "INFO: massfit_sig_pdf_name: " << massfit_sig_pdf_name << std::endl; //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")

    // MASSFIT_SG_REGION_MIN //TODO: CONVERT TO STRING CUT SINCE ALLOWING MULTIPLE INVARIANT MASS FIT VARIABLES
    double massfit_sg_region_min = 1.11;
    if (node["massfit_sg_region_min"]) {
        massfit_sg_region_min = node["massfit_sg_region_min"].as<double>();
    }
    std::cout << "INFO: massfit_sg_region_min: " << massfit_sg_region_min << std::endl;

    // MASSFIT_SG_REGION_MAX 
    double massfit_sg_region_max = 1.13;
    if (node["massfit_sg_region_max"]) {
        massfit_sg_region_max = node["massfit_sg_region_max"].as<double>();
    }
    std::cout << "INFO: massfit_sg_region_max: " << massfit_sg_region_max << std::endl;

    // END FIT ARGUMENTS
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN SPLOT ARGUMENTS

    // SGYIELD_NAME
    std::string sgyield_name = "sgyield";
    if (node["sgyield_name"]) {
        sgyield_name = node["sgyield_name"].as<std::string>();
    }
    std::cout << "INFO: sgyield_name: " << sgyield_name << std::endl;

    // BGYIELD_NAME
    std::string bgyield_name = "bgyield";
    if (node["bgyield_name"]) {
        bgyield_name = node["bgyield_name"].as<std::string>();
    }
    std::cout << "INFO: bgyield_name: " << bgyield_name << std::endl;

    // USE_SPLOT
    bool use_splot = true;
    if (node["use_splot"]) {
        use_splot = node["use_splot"].as<bool>();
    }
    std::cout << "INFO: use_splot: " << use_splot << std::endl;

    // END SPLOT ARGUMENTS
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN SIDEBAND SUBTRACTION ARGUMENTS

    // MASSFIT_SGCUT //TODO: Rename to SB_SGCUT
    std::string massfit_sgcut = "";
    if (node["massfit_sgcut"]) {
        massfit_sgcut = node["massfit_sgcut"].as<std::string>();
    }
    std::cout << "INFO: massfit_sgcut: " << massfit_sgcut << std::endl;

    // MASSFIT_BGCUT //TODO: Rename to SB_BGCUT
    std::string massfit_bgcut = "";
    if (node["massfit_bgcut"]) {
        massfit_bgcut = node["massfit_bgcut"].as<std::string>();
    }
    std::cout << "INFO: massfit_bgcut: " << massfit_bgcut << std::endl;

    // USE_SB_SUBTRACTION
    bool use_sb_subtraction = false;
    if (node["use_sb_subtraction"]) {
        use_sb_subtraction = node["use_sb_subtraction"].as<bool>();
    }
    std::cout << "INFO: use_sb_subtraction: " << use_sb_subtraction << std::endl;

    // USE_BINNED_SB_WEIGHTS
    bool use_binned_sb_weights = false;
    if (node["use_binned_sb_weights"]) {
        use_binned_sb_weights = node["use_binned_sb_weights"].as<bool>();
    }
    std::cout << "INFO: use_binned_sb_weights: " << use_binned_sb_weights << std::endl;

    // ASYMFITVAR_BINSCHEME
    std::map<std::string,std::map<int,std::string>> asymfitvar_bincuts_map;
    if (node["asymfitvar_binschemes"] && node["asymfitvar_binschemes"].IsMap()) {

        // Get bin scheme node and get bin cuts maps
        asymfitvar_bincuts_map = saga::bins::getBinCutsMap(node["asymfitvar_binschemes"]);
    }

    // END SIDEBAND SUBTRACTION ARGUMENTS
    //----------------------------------------------------------------------//

    // LOGPATH
    std::string logpath = "out.txt";
    if (node["logpath"]) {
        logpath = node["logpath"].as<std::string>();
    }
    std::cout << "INFO: logpath: " << logpath << std::endl;

    // DUMP_VARS
    std::vector<std::string> dump_vars; //NOTE: If empty ALL variables from dataset will be dumped if dump_dataset==true.
    if (node["dump_vars"]) {
        dump_vars = node["dump_vars"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: dump_vars: [ ";
    for (int idx=0; idx<dump_vars.size(); idx++) {
        if (idx!=dump_vars.size()-1) { std::cout << dump_vars[idx].c_str()<<", "; }
        else { std::cout << dump_vars[idx].c_str(); }
    }
    std::cout << " ]" << std::endl;

    // DUMP_DATASET
    bool dump_dataset = false;
    if (node["dump_dataset"]) {
        dump_dataset = node["dump_dataset"].as<bool>();
    }
    std::cout << "INFO: dump_dataset: " << dump_dataset << std::endl;

    //----------------------------------------------------------------------------------------------------//
    // ANALYSIS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Create random number generator for MC asymmetry injection
    TRandom *gRandom = new TRandom(seed); //NOTE: IMPORTANT: Need `new` here to get a pointer.

    // Add all absolute variable limits to overall cuts
    cuts = saga::util::addLimitCuts(cuts,binvars,binvar_lims);
    cuts = saga::util::addLimitCuts(cuts,depolvars,depolvar_lims);
    cuts = saga::util::addLimitCuts(cuts,asymfitvars,asymfitvar_lims);
    cuts = saga::util::addLimitCuts(cuts,massfitvars,massfitvar_lims);
    std::cout << "INFO: cuts: "<<cuts.c_str() << std::endl;

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Create branch names for mc variables assuming they all append `_mc` to the corresponding data branch name
    std::vector<std::string> asymfitvars_mc;
    for (int idx=0; idx<asymfitvars.size(); idx++) {
        asymfitvars_mc.push_back(Form("%s_mc",asymfitvars[idx].c_str()));
        std::cout << "INFO: Defined MC variable : " << asymfitvars_mc[idx].c_str() << std::endl;
    }
    std::vector<std::string> depolvars_mc;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvars_mc.push_back(Form("%s_mc",depolvars[idx].c_str()));
        std::cout << "INFO: Defined MC variable : " << depolvars_mc[idx].c_str() << std::endl;
    }

    // Replace variable names in XS formulas with their MC counter parts and asymmetries with their injected values
    if (inject_asym) {
        // Find and replace asymmetry names with injected values in fsgasyms_xs : example string fsgasyms_xs="0.747*depolvars_mc0*sgasym0*fitvar1_mc"
        for (int idx=0; idx<asymfitvars.size(); idx++) {
            saga::util::replaceAll(fsgasyms_xs_pu_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_up_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_pp_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
        }
        for (int idx=0; idx<sgasyms.size(); idx++) {
            saga::util::replaceAll(fsgasyms_xs_pu_formula, Form("sgasyms[%d]",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fsgasyms_xs_up_formula, Form("sgasyms[%d]",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fsgasyms_xs_pp_formula, Form("sgasyms[%d]",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
        }
        for (int idx=0; idx<depolvars.size(); idx++) {
            saga::util::replaceAll(fsgasyms_xs_pu_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_up_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_pp_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
        }
        std::cout << "INFO: Updated " << fsgasyms_xs_pu_name.c_str() << " = " << fsgasyms_xs_pu_formula.c_str() << std::endl;
        std::cout << "INFO: Updated " << fsgasyms_xs_up_name.c_str() << " = " << fsgasyms_xs_up_formula.c_str() << std::endl;
        std::cout << "INFO: Updated " << fsgasyms_xs_pp_name.c_str() << " = " << fsgasyms_xs_pp_formula.c_str() << std::endl;

        // Find and replace placeholder variable names with actual values in fbgasyms_xs : example string fbgasyms_xs="0.747*depolvars_mc0*bgasym0*fitvar1_mc"
        for (int idx=0; idx<asymfitvars.size(); idx++) {
            saga::util::replaceAll(fbgasyms_xs_pu_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_up_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_pp_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
        }
        for (int idx=0; idx<bgasyms.size(); idx++) {
            saga::util::replaceAll(fbgasyms_xs_pu_formula, Form("bgasyms[%d]",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fbgasyms_xs_up_formula, Form("bgasyms[%d]",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fbgasyms_xs_pp_formula, Form("bgasyms[%d]",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
        }
        for (int idx=0; idx<depolvars.size(); idx++) {
            saga::util::replaceAll(fbgasyms_xs_pu_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_up_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_pp_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
        }
        std::cout << "INFO: Updated " << fbgasyms_xs_pu_name.c_str() << " = " << fbgasyms_xs_pu_formula.c_str() << std::endl;
        std::cout << "INFO: Updated " << fbgasyms_xs_up_name.c_str() << " = " << fbgasyms_xs_up_formula.c_str() << std::endl;
        std::cout << "INFO: Updated " << fbgasyms_xs_pp_name.c_str() << " = " << fbgasyms_xs_pp_formula.c_str() << std::endl;
    }

    // Define variables from formulas
    auto d2 = d.Define("__dummyvar__","(float)0.0"); //NOTE: Define a dummy variable to declare the data frame in this scope.
    for (int idx=0; idx<var_formulas.size(); idx++) {
        d2 = d2.Define(var_formulas[idx][0].c_str(),var_formulas[idx][1].c_str());
        std::cout<<"INFO: Defined branch "<<var_formulas[idx][0].c_str()<<std::endl;
    }

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
        d2_filtered = d2_filtered.Define(fsgasyms_xs_pu_name.c_str(),fsgasyms_xs_pu_formula.c_str());
        d2_filtered = d2_filtered.Define(fsgasyms_xs_up_name.c_str(),fsgasyms_xs_up_formula.c_str());
        d2_filtered = d2_filtered.Define(fsgasyms_xs_pp_name.c_str(),fsgasyms_xs_pp_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_pu_name.c_str(),fbgasyms_xs_pu_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_up_name.c_str(),fbgasyms_xs_up_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_pp_name.c_str(),fbgasyms_xs_pp_formula.c_str());
    }
    //TODO: Add output message about defined branches

    // Define run-dependent columns from CSV
    for (int idx=0; idx<csv_paths.size(); idx++) {
        d2_filtered = saga::data::mapDataFromCSV<Long64_t,double>(
            d2_filtered,
            rdf_key_cols[idx],
            csv_paths[idx],
            csv_key_cols[idx],
            col_names[idx],
            col_aliases,
            true,
            ','
        );
    }
    //TODO: Add output messages about defined branches

    // Define helicity variable, injecting and applying MC matching cuts if requested
    auto frame = (!inject_asym) ?
                    d2_filtered.Define(helicity_name.c_str(), helicity_formula.c_str())
                                .Define(tspin_name.c_str(), tspin_formula.c_str()) :
                    saga::data::injectAsym(
                        d2_filtered,
                        seed,
                        bpol,
                        tpol,
                        mc_sg_match_name,
                        fsgasyms_xs_pu_name,
                        fsgasyms_xs_up_name,
                        fsgasyms_xs_pp_name,
                        fbgasyms_xs_pu_name,
                        fbgasyms_xs_up_name,
                        fbgasyms_xs_pp_name,
                        combined_spin_state,
                        helicity_name,
                        tspin_name
                    );
    //TODO: Add output message about defined branches

    // Make sure injection values are all computed before running analysis
    if (inject_asym) {
        double my_testvar  = (double)*frame.Mean(combined_spin_state.c_str());
        double my_testvar1 = (double)*frame.Mean(helicity_name.c_str());
        double my_testvar2 = (double)*frame.Mean(tspin_name.c_str());
    }

    // Dump dataset to ROOT file and exit
    if (dump_dataset) {
        std::string out_ds_path = Form("%sdataset.root", baseoutpath.c_str());
        std::cout<<"INFO: Dumping dataset to: "<<out_ds_path.c_str()<<std::endl;
        if (dump_vars.size()==0) frame.Snapshot(tree.c_str(), out_ds_path.c_str());
        else frame.Snapshot(tree.c_str(), out_ds_path.c_str(), dump_vars);

        return;
    }

    // Create output log
    std::ofstream outf; outf.open(logpath.c_str());
    std::ostream &out = outf; //std::cout;

    // Loop bin schemes
    for (auto it = bincuts_map.begin(); it != bincuts_map.end(); ++it) {

        // Get bin scheme name and bin cut map
        std::string binscheme_name = it->first;
        std::map<int,std::string> bincuts = it->second;

        // Get bin scheme bin variables
        std::vector<std::string> scheme_binvars = binschemes_vars[binscheme_name];

        // Get asymmetry fit variables bin scheme
        std::map<int,std::string> asymfitvar_bincuts;
        if (use_binned_sb_weights) {
            asymfitvar_bincuts = asymfitvar_bincuts_map[binscheme_name];
        }

        // Get additional bin variable attributes needed for this binning scheme
        std::vector<std::string> scheme_binvar_titles;
        std::vector<std::vector<double>> scheme_binvar_lims;
        std::vector<int> scheme_binvar_bins;
        for (int idx=0; idx<binvars.size(); idx++) {
            int count = std::count(scheme_binvars.begin(), scheme_binvars.end(), binvars[idx]);
            if (count>0) {
                scheme_binvar_titles.push_back(binvar_titles[idx]);
                scheme_binvar_lims.push_back(binvar_lims[idx]);
                scheme_binvar_bins.push_back(binvar_bins[idx]);
            }
        }

        // Produce graphs of asymmetry fit parameters corrected for depolarization and background binned in given kinematic variable
        std::string scheme_name = Form("%s%s",baseoutpath.c_str(),binscheme_name.c_str());
        saga::analysis::getKinBinnedAsym(
            scheme_name, //std::string                      scheme_name,
            frame, //ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            "w", //std::string                      workspace_name,
            "workspace", //std::string                      workspace_title,

            // parameters passed to saga::data::createDataset()
            "dataset", //std::string                      dataset_name,
            "dataset", //std::string                      dataset_title,
            helicity_name, // std::string                      helicity,
            helicity_states, //std::map<std::string,int>        helicity_states,
            tspin_name, // std::string                         tspin,
            tspin_states, //std::map<std::string,int>           tspin_states,
            htspin_name, // std::string                         htspin,
            htspin_states, //std::map<std::string,int>           htspin_states,
            combined_spin_state, //std::string                      combined_spin_state,
            bincuts, //std::map<int,std::string>        bincuts,
            scheme_binvars, //std::vector<std::string>         binvars,
            scheme_binvar_titles, //std::vector<std::string>         binvar_titles,
            scheme_binvar_lims, //std::vector<std::vector<double>> binvar_lims,
            scheme_binvar_bins, //std::vector<int>                 binvar_bins,
            depolvars, //std::vector<std::string>         depolvars,
            depolvar_titles, //std::vector<std::string>         depolvar_titles,
            depolvar_lims, //std::vector<std::vector<double>> depolvar_lims,
            depolvar_bins, //std::vector<int>                 depolvar_bins,
            asymfitvars, //std::vector<std::string>         asymfitvars,
            asymfitvar_titles, //std::vector<std::string>         asymfitvar_titles,
            asymfitvar_lims, //std::vector<std::vector<double>> asymfitvar_lims,
            asymfitvar_bins, //std::vector<int>                 asymfitvar_bins,
            massfitvars, //std::vector<std::string>         massfitvars,
            massfitvar_titles, //std::vector<std::string>         massfitvar_titles,
            massfitvar_lims, //std::vector<std::vector<double>> massfitvar_lims,
            massfitvar_bins, //std::vector<int>                 massfitvar_bins,

            // parameterss passed to analysis::fitAsym()
            bpol, //double                           bpol,
            tpol, //double                           tpol,
            asymfit_formula_pu, //std::string                      asymfit_formula_pu,
            asymfit_formula_up, //std::string                      asymfit_formula_up,
            asymfit_formula_pp, //std::string                      asymfit_formula_pp,
            asymfitpar_inits, //std::vector<double>              asymfitpar_inits,
            asymfitpar_initlims, //std::vector<std::vector<double>> asymfitpar_initlims,
            use_sumw2error, //bool                             use_sumw2error,
            use_average_depol, //bool                             use_average_depol,
            use_extended_nll, //bool                             use_extended_nll,
            use_binned_fit, //bool                             use_binned_fit,

            // parameters passed to saga::data::createDataset() and analysis::applyLambdaMassFit() //TODO: Add init fit parameter value and limits arguments here...assuming you always want a chebychev polynomial background...
            massfit_model_name, //std::string                      massfit_model_name,
            massfit_nbins_conv, //int                              massfit_nbins_conv,
            massfit_sig_pdf_name, //std::string                      massfit_sig_pdf_name, //NOTE: This must be one of ("gauss","landau","cb","landau_X_gauss","cb_X_gauss")
            massfit_sg_region_min, //double                           massfit_sg_region_min,
            massfit_sg_region_max, //double                           massfit_sg_region_max,

            // Parameters passed to analysis::applySPlots()
            sgyield_name, //std::string                      sgYield_name,
            bgyield_name, //std::string                      bgYield_name,
            use_splot, //bool                             use_splot,

            // Parameters used for sb subtraction
            massfit_sgcut, //std::string                      massfit_sgcut,
            massfit_bgcut, //std::string                      massfit_bgcut,
            use_sb_subtraction, //bool                             use_sb_subtraction,
            use_binned_sb_weights, //bool                             use_binned_sb_weights,
            asymfitvar_bincuts, //std::map<int,std::string>        asymfitvar_bincuts,

            // Output stream
            out //std::ostream &out                = std::cout
        );
    } // for (auto it = bincuts_map.begin(); it != bincuts_map.end(); ++it) {

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
