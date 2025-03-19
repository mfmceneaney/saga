#include <iostream>
#include <memory>
#include <fstream>
#include <string>
#include <map>

// YAML Includes
#include <yaml-cpp/yaml.h>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>
#include <TMath.h>

// Project Includes
#include <bins.h>
#include <util.h>

void execute(const YAML::Node& node) {

    // Process arguments

    // OUTPATH
    std::string baseoutpath = "";//NOTE: This will be prepended to the default output path like so: `<baseoutpath><binscheme_name>_kinematics.csv`.
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

    // MC_CUTS
    std::string mc_cuts = "Q2>1"; //NOTE: This may not be empty!
    if (node["mc_cuts"]) {
        mc_cuts = node["mc_cuts"].as<std::string>();
    }
    std::cout << "INFO: mc_cuts: " << mc_cuts << std::endl;

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

    // END ASYMMETRY INJECTION ARGUMENTS
    //----------------------------------------------------------------------//

    // VAR_FORMULAS
    std::map<std::string,std::string> var_formulas;
    if (node["var_formulas"]) {
        var_formulas = node["var_formulas"].as<std::map<std::string,std::string>>();
    }
    std::cout << "INFO: var_formulas: { \n";
    for (auto it = var_formulas.begin(); it != var_formulas.end(); ++it) {
        std::cout << "\t[ " << it->first.c_str() << " : " << it->second.c_str() << " ],\n";
    }
    std::cout << " }" << std::endl;

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

    // END MASS FIT VARIABLES
    //----------------------------------------------------------------------//


    //----------------------------------------------------------------------//
    // BEGIN KINEMATIC VARIABLES

    // KINVARS
    std::vector<std::vector<std::string>> kinvars;
    if (node["kinvars"]) {
        kinvars = node["kinvars"].as<std::vector<std::vector<std::string>>>();
    }
    std::cout << "INFO: kinvars: [ ";
    for (int idx=0; idx<kinvars.size(); idx++) {
        if (idx!=kinvars.size()-1) { std::cout << "\n\t[" << kinvars[idx][0].c_str() << " , " << kinvars[idx][1].c_str() <<" ],"; }
        else { std::cout << "\n\t[" << kinvars[idx][0].c_str() << " , " << kinvars[idx][1].c_str() <<" ]\n"; }
    }
    std::cout << " ]" << std::endl;

    // KINVAR_LIMS
    std::vector<std::vector<std::vector<double>>> kinvar_lims;
    if (node["kinvar_lims"]) {
        kinvar_lims = node["kinvar_lims"].as<std::vector<std::vector<std::vector<double>>>>();
    }
    std::cout << "INFO: kinvar_lims: [ ";
    for (int idx=0; idx<kinvar_lims.size(); idx++) {
        if (idx!=kinvar_lims.size()-1) { std::cout << "\n\t[" << "\n\t\t[ " << kinvar_lims[idx][0][0] << ", " << kinvar_lims[idx][0][1] << " ]," << "\n\t\t[ " << kinvar_lims[idx][1][0] << ", " << kinvar_lims[idx][1][1] << " ], " << "\n\t],"; }
        else { std::cout << "\n\t[" << "\n\t\t[ " << kinvar_lims[idx][0][0] << ", " << kinvar_lims[idx][0][1] << " ], " << "\n\t\t[ " << kinvar_lims[idx][1][0] << ", " << kinvar_lims[idx][1][1] << " ], " << "\n"; }
    }
    std::cout << " ]" << std::endl;

    // KINVAR_BINS
    std::vector<std::vector<int>> kinvar_bins;
    if (node["kinvar_bins"]) {
        kinvar_bins = node["kinvar_bins"].as<std::vector<std::vector<int>>>();
    }
    std::cout << "INFO: kinvar_bins: [ ";
    for (int idx=0; idx<kinvar_bins.size(); idx++) {
        if (idx!=kinvar_bins.size()-1) { std::cout << "\n\t[" << kinvar_bins[idx][0] << " , " << kinvar_bins[idx][1] <<" ],"; }
        else { std::cout << "\n\t[" << kinvar_bins[idx][0] << " , " << kinvar_bins[idx][1] <<" ]\n"; }
    }
    std::cout << " ]" << std::endl;

    // END KINEMATIC VARIABLES
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------------------------------------//
    // ANALYSIS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Add all absolute variable limits to overall cuts
    cuts = saga::util::addLimitCuts(cuts,binvars,binvar_lims);
    cuts = saga::util::addLimitCuts(cuts,depolvars,depolvar_lims);
    cuts = saga::util::addLimitCuts(cuts,asymfitvars,asymfitvar_lims);
    cuts = saga::util::addLimitCuts(cuts,massfitvars,massfitvar_lims);
    std::cout << "INFO: cuts: "<<cuts.c_str() << std::endl;

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Define variables from formulas
    auto d2 = d.Define("__dummyvar__","(float)0.0"); //NOTE: Define a dummy variable to declare the data frame in this scope.
    for (auto it = var_formulas.begin(); it != var_formulas.end(); ++it) {
        d2 = d2.Define(it->first.c_str(),it->second.c_str());
        std::cout<<"INFO: Defined branch "<<it->first.c_str()<<std::endl;
    }

    // Apply overall cuts AFTER defining depolarization and fit variables
    auto d2_filtered = d2.Filter(Form("(%s) && (%s)",cuts.c_str(),mc_cuts.c_str()));

    // Define MC matching variable names
    std::vector<std::string> theta_vars;
    std::vector<std::string> phi_vars;
    std::vector<std::string> theta_mc_vars;
    std::vector<std::string> phi_mc_vars;
    std::vector<std::string> dtheta_vars;
    std::vector<std::string> dphi_vars;
    for (int idx=0; idx<particle_suffixes.size(); idx++) {
        theta_vars.push_back(Form("theta%s",particle_suffixes[idx].c_str()));
        phi_vars.push_back(Form("phi%s",particle_suffixes[idx].c_str()));
        theta_mc_vars.push_back(Form("theta%s_mc",particle_suffixes[idx].c_str()));
        phi_mc_vars.push_back(Form("phi%s_mc",particle_suffixes[idx].c_str()));
        dtheta_vars.push_back(Form("dtheta%s",particle_suffixes[idx].c_str()));
        dphi_vars.push_back(Form("dphi%s",particle_suffixes[idx].c_str()));
    }

    // Define MC matching angular difference variable branches
    for (int idx=0; idx<particle_suffixes.size(); idx++) {
        d2_filtered = d2_filtered.Define(dtheta_vars[idx].c_str(),[](float theta, float theta_mc){ return TMath::Abs(theta-theta_mc); },{theta_vars[idx].c_str(),theta_mc_vars[idx].c_str()})
            .Define(dphi_vars[idx].c_str(),[](float phi, float phi_mc){
                return (float) (TMath::Abs(phi-phi_mc)<TMath::Pi()
                ? TMath::Abs(phi-phi_mc) : 2*TMath::Pi() - TMath::Abs(phi-phi_mc));
                },{phi_vars[idx].c_str(),phi_mc_vars[idx].c_str()});
    }
    //TODO: Add output message about defined branches
 
    // Loop bin schemes
    for (auto it = bincuts_map.begin(); it != bincuts_map.end(); ++it) {

        // Get bin scheme name and bin cut map
        std::string binscheme_name = it->first;
        std::map<int,std::string> bincuts = it->second;

        // Get statistics and average values and errors of kinematic variables in each bin
        std::string scheme_name = Form("%s%s",baseoutpath.c_str(),binscheme_name.c_str());
        saga::bins::getBinKinematicsTH2Ds(
            d2_filtered, // ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
            scheme_name, // std::string                                                   scheme_name,
            bincuts, // std::map<int,std::string>                                     bincuts,
            kinvars, // std::vector<std::vector<std::string>>                         kinvars,
            kinvar_lims, // std::vector<std::vector<std::vector<double>>>                 kinvar_lims,
            kinvar_bins // std::vector<std::vector<int>>                                 kinvar_bins
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
