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

    // Get input ROOT tree path
    std::string inpath = "";
    if (node["inpath"]) {
        inpath = node["inpath"].as<std::string>();
    }
    std::cout << "INFO: inpath: " << inpath << std::endl;

    // Get tree name
    std::string tree = "";
    if (node["tree"]) {
        tree = node["tree"].as<std::string>();
    }
    std::cout << "INFO: tree: " << tree << std::endl;

    // Get number of threads for frame
    int nthreads = 1;
    if (node["nthreads"]) {
        nthreads = node["nthreads"].as<int>();
    }
    std::cout << "INFO: nthreads: " << nthreads << std::endl;

    // Get frame cuts
    std::string cuts = "";
    if (node["cuts"]) {
        cuts = node["cuts"].as<std::string>();
    }
    std::cout << "INFO: cuts: " << cuts << std::endl;

    // Get output path
    std::string outpath = "out.yaml";
    if (node["outpath"]) {
        outpath = node["outpath"].as<std::string>();
    }
    std::cout << "INFO: outpath: " << outpath << std::endl;
    
    // Get bin variable names
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

    // Get defined variable names
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

    // Get defined variable names
    std::vector<std::string> defvars;
    if (node["defvars"]) {
        defvars = node["defvars"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: defvars: [ ";
    for (int idx=0; idx<defvars.size(); idx++) {
        if (idx!=defvars.size()-1) { std::cout << defvars[idx]<<", "; }
        else { std::cout << defvars[idx]; }
    }
    std::cout << " ]" << std::endl;

    // Get defined variable names
    std::vector<std::string> defvar_formulas;
    if (node["defvar_formulas"]) {
        defvar_formulas = node["defvar_formulas"].as<std::vector<std::string>>();
    }
    std::cout << "INFO: defvar_formulas: [ ";
    for (int idx=0; idx<defvar_formulas.size(); idx++) {
        if (idx!=defvar_formulas.size()-1) { std::cout << defvar_formulas[idx]<<", "; }
        else { std::cout << defvar_formulas[idx]; }
    }
    std::cout << " ]" << std::endl;
    if (defvars.size()!=defvar_formulas.size()) {
        std::cerr << "ERROR: defvars.size() must match defvar_formulas.size()" << std::endl;
    }

    //----------------------------------------------------------------------------------------------------//
    // FINDBINLIMITS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Pre-define additional variables.
    auto d2 = d.Define("_randvar_","(float)1.0");//NOTE: Need the type of d2 declared outside the loop to match the type assigned inside the loop below.
    for (int idx=0; idx<defvars.size(); idx++) {
        d2 = d2.Define(defvars[idx].c_str(),defvar_formulas[idx].c_str());
    }
    
    // Apply overall cuts AFTER defining variables
    auto d2_filtered = d2.Filter(cuts.c_str());

    // Initialize map that will be converted to yaml
    std::map<std::string,std::vector<double>> binlimits;

    // Loop variables to bin in
    for (int binvar_idx = 0; binvar_idx<binvars.size(); binvar_idx++) {

        // Get desired number of bins
        int nbins = nbins_list[binvar_idx];

        // Find bin limits
        std::vector<double> _binlimits = saga::bins::findBinLimits(
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
