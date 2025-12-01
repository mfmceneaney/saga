#include <iostream>
#include <memory>
#include <fstream>
#include <string>

// YAML Includes
#include <yaml-cpp/yaml.h>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>

// Project Includes
#include <log.h>
#include <bins.h>

void execute(const YAML::Node& node) {

    // Process arguments
    std::string   message_prefix  = "[findBinLims]: ";
    bool          verbose         = true;
    std::ostream &yamlargout      = std::cout;

    //----------------------------------------------------------------------//
    // BEGIN ARGUMENTS
    saga::log::Logger::instance().setOutputStream(yamlargout);
    std::string log_level = saga::util::getYamlArg<std::string>(node, "log_level", "INFO", message_prefix, verbose);
    saga::log::Logger::instance().setLogLevelFromString(log_level);
    std::string outpath = saga::util::getYamlArg<std::string>(node,"outpath","out.yaml",message_prefix,verbose);
    std::string inpath = saga::util::getYamlArg<std::string>(node,"inpath","",message_prefix,verbose);
    std::string tree = saga::util::getYamlArg<std::string>(node,"tree","t",message_prefix,verbose);
    int nthreads = saga::util::getYamlArg<int>(node,"nthreads",1,message_prefix,verbose);
    std::string cuts = saga::util::getYamlArg<std::string>(node,"cuts","",message_prefix,verbose);

    //----------------------------------------------------------------------//
    // BEGIN BIN VARIABLES
    std::vector<std::string> binvars = saga::util::getYamlArg<std::vector<std::string>>(node, "binvars", {}, message_prefix, verbose);

    // NBINS_LIST
    std::vector<int> nbins_list = saga::util::getYamlArg<std::vector<int>>(node, "nbins_list", {}, message_prefix, verbose);
    if (binvars.size()!=nbins_list.size()) {
        std::cerr << "ERROR: binvars.size() must match nbins_list.size()" << std::endl;
    }

    //----------------------------------------------------------------------//
    // VAR_FORMULAS
    std::vector<std::vector<std::string>> var_formulas = saga::util::getYamlArg<std::vector<std::vector<std::string>>>(node, "var_formulas", {}, message_prefix, verbose);

    //----------------------------------------------------------------------------------------------------//
    // FINDBINLIMS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Define variables from formulas
    auto d2 = d.Define("__dummyvar__","(float)0.0"); //NOTE: Define a dummy variable to declare the data frame in this scope.
    for (int idx=0; idx<var_formulas.size(); idx++) {
        d2 = d2.Define(var_formulas[idx][0].c_str(),var_formulas[idx][1].c_str());
        yamlargout << message_prefix.c_str() << "Defined branch "<<var_formulas[idx][0].c_str()<<std::endl;
    }
    
    // Apply overall cuts AFTER defining variables
    auto d2_filtered = d2.Filter(cuts.c_str());

    // Check for nested bin scheme and then default to grid
    std::string node_binscheme_name = "binscheme";
    auto node_binscheme = node[node_binscheme_name];
    if (node_binscheme && node_binscheme.IsMap()) {
        saga::bins::findNestedBinLims(
            d2_filtered, // ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
            node_binscheme, // YAML::Node node,
            "", // std::string node_name = "",
            "nbins", // std::string nbins_key = "nbins",
            "lims", // std::string lims_key = "lims",
            "nested", // std::string nested_key = "nested",
            {} // std::vector<std::string> bin_cuts = {}
        );

        // Save bin limits to a yaml file
        std::ofstream outf; outf.open(outpath.c_str());
        std::ostream &out = outf;
        YAML::Emitter yout(out);
        YAML::Node node_out;
        node_out[node_binscheme_name] = node_binscheme;
        yout << node_out;

    } else {

        // Initialize map that will be converted to yaml
        std::map<std::string,std::vector<double>> binlimits;

        // Loop variables to bin in
        for (int binvar_idx = 0; binvar_idx<binvars.size(); binvar_idx++) {

            // Get desired number of bins
            int nbins = nbins_list[binvar_idx];

            // Find bin limits
            std::vector<double> _binlimits = saga::bins::findBinLims(
                d2_filtered,
                binvars[binvar_idx],
                nbins
            );

            // Add bin limits to map
            binlimits[binvars[binvar_idx]] = _binlimits;
    
        }// loop over bin variables

        // Save bin limits to a yaml file
        std::ofstream outf; outf.open(outpath.c_str());
        std::ostream &out = outf;
        YAML::Emitter yout(out);
        yout << YAML::Flow << binlimits;
    }

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
