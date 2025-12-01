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
#include <log.h>
#include <bins.h>
#include <data.h>
#include <util.h>

void execute(const YAML::Node& node) {

    std::string   message_prefix  = "[getBinKinematics]: ";
    bool          verbose         = true;
    std::ostream &yamlargout      = std::cout;

    //----------------------------------------------------------------------//
    // BEGIN ARGUMENTS
    saga::log::Logger::instance().setOutputStream(yamlargout);
    std::string log_level = saga::util::getYamlArg<std::string>(node, "log_level", "INFO", message_prefix, verbose);
    saga::log::Logger::instance().setLogLevelFromString(log_level);
    std::string baseoutpath = saga::util::getYamlArg<std::string>(node,"baseoutpath","",message_prefix,verbose); //NOTE: This will be prepended to the default output path like so: `<baseoutpath><binscheme_name>.csv`.
    std::string inpath = saga::util::getYamlArg<std::string>(node,"inpath","",message_prefix,verbose);
    std::string tree = saga::util::getYamlArg<std::string>(node,"tree","t",message_prefix,verbose);
    int nthreads = saga::util::getYamlArg<int>(node,"nthreads",1,message_prefix,verbose);
    std::string cuts = saga::util::getYamlArg<std::string>(node,"cuts","",message_prefix,verbose);

    //----------------------------------------------------------------------//
    // BEGIN MC MATCHING ARGUMENTS
    std::string mc_cuts = saga::util::getYamlArg<std::string>(node,"mc_cuts","Q2>1",message_prefix,verbose); //NOTE: This may not be empty!
    std::vector<std::string> particle_suffixes = saga::util::getYamlArg<std::vector<std::string>>(node, "particle_suffixes", {}, message_prefix, verbose); // -> For MC matching with mc_sg_match cut

    //----------------------------------------------------------------------//
    // VAR_FORMULAS
    std::vector<std::vector<std::string>> var_formulas = saga::util::getYamlArg<std::vector<std::vector<std::string>>>(node, "var_formulas", {}, message_prefix, verbose);

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
        std::vector<std::string> binschemes_paths = saga::util::getYamlArg<std::vector<std::string>>(node, "binschemes_paths", {}, message_prefix, verbose);

        // Loop paths and add bin schemes
        for (int idx=0; idx<binschemes_paths.size(); idx++) {

            // Load YAML file and get bin cuts maps
            if (verbose) yamlargout << message_prefix.c_str() << "Loading bin scheme from : " <<binschemes_paths[idx].c_str() << std::endl;
            YAML::Node bincut_config = YAML::LoadFile(binschemes_paths[idx].c_str());
            std::map<std::string,std::map<int,std::string>> new_bincuts_map = saga::bins::getBinCutsMap(bincut_config);
            bincuts_map.insert(new_bincuts_map.begin(), new_bincuts_map.end());
            std::map<std::string,std::vector<std::string>> new_binschemes_vars = saga::bins::getBinSchemesVars(bincut_config);
            binschemes_vars.insert(new_binschemes_vars.begin(), new_binschemes_vars.end());
        }
    }

    // Reduce bin cuts map into a single batch for parallelization
    if (node["nbatches"] && node["ibatch"]) {
        int nbatches = saga::util::getYamlArg<int>(node, "nbatches", 1, message_prefix, verbose);
        int ibatch   = saga::util::getYamlArg<int>(node, "ibatch", 0, message_prefix, verbose);
        if (nbatches>1 && ibatch>=0 && ibatch<nbatches) bincuts_map = saga::bins::getBinCutsMapBatch(bincuts_map, nbatches, ibatch);
    }

    // Show bin cuts map
    if (verbose) {
        for (auto it = bincuts_map.begin(); it != bincuts_map.end(); ++it) {
            yamlargout << message_prefix.c_str() << "bincuts_map["<<it->first<<"]: { \n";
            for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                yamlargout <<"\t\t"<< it2->first<<": "<<it2->second<<", \n";
            }
            yamlargout << "}" << std::endl;
        }
    }

    // END BINNING SCHEME ARGUMENTS
    //----------------------------------------------------------------------//

    //----------------------------------------------------------------------//
    // BEGIN BIN VARIABLES
    std::vector<std::string> binvars = saga::util::getYamlArg<std::vector<std::string>>(node, "binvars", {}, message_prefix, verbose);
    std::vector<std::vector<double>> binvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "binvar_lims", {}, message_prefix, verbose);

    //----------------------------------------------------------------------//
    // BEGIN DEPOLARIZATION VARIABLES
    std::vector<std::string> depolvars = saga::util::getYamlArg<std::vector<std::string>>(node, "depolvars", {}, message_prefix, verbose);
    std::vector<std::vector<double>> depolvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "depolvar_lims", {}, message_prefix, verbose);

    //----------------------------------------------------------------------//
    // BEGIN ASYMMETRY FIT VARIABLES
    std::vector<std::string> asymfitvars = saga::util::getYamlArg<std::vector<std::string>>(node, "asymfitvars", {}, message_prefix, verbose);
    std::vector<std::vector<double>> asymfitvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "asymfitvar_lims", {}, message_prefix, verbose);

    //----------------------------------------------------------------------//
    // BEGIN MASS FIT VARIABLES
    std::vector<std::string> massfitvars = saga::util::getYamlArg<std::vector<std::string>>(node, "massfitvars", {}, message_prefix, verbose);
    std::vector<std::vector<double>> massfitvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfitvar_lims", {}, message_prefix, verbose);

    //----------------------------------------------------------------------//
    // BEGIN KINEMATIC VARIABLES
    std::vector<std::string> kinvars = saga::util::getYamlArg<std::vector<std::string>>(node, "kinvars", {}, message_prefix, verbose);

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
    for (int idx=0; idx<var_formulas.size(); idx++) {
        d2 = d2.Define(var_formulas[idx][0].c_str(),var_formulas[idx][1].c_str());
        yamlargout << message_prefix.c_str() << "Defined branch "<<var_formulas[idx][0].c_str()<<std::endl;
    }

    // Apply overall cuts AFTER defining depolarization and fit variables
    std::string all_cuts = (mc_cuts.size()>0) ? Form("(%s) && (%s)",cuts.c_str(),mc_cuts.c_str()) : cuts;
    auto d2_filtered = d2.Filter(all_cuts.c_str());

    // Define angular difference variables
    if (mc_cuts.size()>0) d2_filtered = saga::data::defineAngularDiffVars(d2_filtered, particle_suffixes, "theta", "phi", "_mc");
    //TODO: Add output message about defined branches
 
    // Loop bin schemes
    for (auto it = bincuts_map.begin(); it != bincuts_map.end(); ++it) {

        // Get bin scheme name and bin cut map
        std::string binscheme_name = it->first;
        std::map<int,std::string> bincuts = it->second;

        // Get statistics and average values and errors of kinematic variables in each bin
        std::string scheme_name = Form("%s%s",baseoutpath.c_str(),binscheme_name.c_str());
        saga::bins::getBinKinematics(
            d2_filtered, //ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
            scheme_name, //std::string                                                   scheme_name,
            bincuts, //std::map<int,std::string>                                     bincuts,
            kinvars //std::vector<std::string>                                      kinvarrs,
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
