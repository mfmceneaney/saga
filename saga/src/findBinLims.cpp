#include <iostream>
#include <memory>
#include <fstream>
#include <string>

// YAML Includes
#include <yaml-cpp/yaml.h>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>

// Project Includes
#include <bins.h>

void execute(const YAML::Node& node) {

    // Process arguments

    // INPATH
    std::string inpath = "";
    if (node["inpath"]) {
        inpath = node["inpath"].as<std::string>();
    }
    std::cout << "INFO: inpath: " << inpath << std::endl;

    // TREREGet tree name
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

    // OUTPATH
    std::string outpath = "out.yaml";
    if (node["outpath"]) {
        outpath = node["outpath"].as<std::string>();
    }
    std::cout << "INFO: outpath: " << outpath << std::endl;
    
    // BINVVARS
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

    // NBINS_LIST
    std::vector<int> nbins_list;
    if (node["nbins_list"]) {
        nbins_list = node["nbins_list"].as<std::vector<int>>();
    }
    std::cout << "INFO: nbins_list: [ ";
    for (int idx=0; idx<nbins_list.size(); idx++) {
        if (idx!=nbins_list.size()-1) { std::cout << nbins_list[idx]<<", "; }
        else { std::cout << nbins_list[idx]; }
    }
    std::cout << " ]" << std::endl;
    if (binvars.size()!=nbins_list.size()) {
        std::cerr << "ERROR: binvars.size() must match nbins_list.size()" << std::endl;
    }

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

    //----------------------------------------------------------------------------------------------------//
    // FINDBINLIMS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Define variables from formulas
    auto d2 = d.Define("__dummyvar__","(float)0.0"); //NOTE: Define a dummy variable to declare the data frame in this scope.
    for (auto it = var_formulas.begin(); it != var_formulas.end(); ++it) {
        d2 = d2.Define(it->first.c_str(),it->second.c_str());
        std::cout<<"INFO: Defined branch "<<it->first.c_str()<<std::endl;
    }
    
    // Apply overall cuts AFTER defining variables
    auto d2_filtered = d2.Filter(cuts.c_str());

    // Check for nested bin scheme and then default to grid
    auto node_binscheme = node["binscheme"];
    if (node_binscheme && node_binscheme.IsMap()) {
        std::map<std::vector<int>,std::vector<double>> nested_bin_lims; //TODO: See if can drop these arguments
        std::vector<std::string> new_binvars = {};
        saga::bins::findNestedBinLims(
            d2_filtered, // ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
            nested_bin_lims, // std::map<std::vector<int>,std::vector<double>> nested_bin_lims,
            {}, // std::vector<std::string> binvars,
            node_binscheme, // YAML::Node node_nested,
            "", // std::string node_nested_name,
            {}, // std::vector<int> coordinates = {},
            "nbins" // std::string nbins_key = "nbins"
        );

        // Save bin limits to a yaml file
        std::ofstream outf; outf.open(outpath.c_str());
        std::ostream &out = outf;
        YAML::Emitter yout(out);
        yout << node_binscheme;

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
