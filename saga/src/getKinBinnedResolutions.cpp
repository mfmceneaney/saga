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
#include <resolution.h>

void execute(const YAML::Node& node) {

    // Process arguments
    std::string   message_prefix  = "INFO: ";
    bool          verbose         = true;
    std::ostream &yamlargout      = std::cout;

    //----------------------------------------------------------------------//
    // BEGIN ARGUMENTS
    saga::log::Logger::instance().setOutputStream(yamlargout);
    std::string log_level = saga::util::getYamlArg<std::string>(node, "log_level", "INFO", message_prefix, verbose, yamlargout);
    saga::log::Logger::instance().setLogLevelFromString(log_level);
    std::string baseoutpath = saga::util::getYamlArg<std::string>(node,"baseoutpath","",message_prefix,verbose,yamlargout); //NOTE: This will be prepended to the default output path like so: `<baseoutpath><binscheme_name>.csv`.
    std::string inpath = saga::util::getYamlArg<std::string>(node,"inpath","",message_prefix,verbose,yamlargout);
    std::string tree = saga::util::getYamlArg<std::string>(node,"tree","t",message_prefix,verbose,yamlargout);
    int nthreads = saga::util::getYamlArg<int>(node,"nthreads",1,message_prefix,verbose,yamlargout);
    std::string cuts = saga::util::getYamlArg<std::string>(node,"cuts","",message_prefix,verbose,yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN MC MATCHING ARGUMENTS
    std::string mc_cuts = saga::util::getYamlArg<std::string>(node,"mc_cuts","Q2>1",message_prefix,verbose,yamlargout); //NOTE: This may not be empty!
    std::vector<std::string> particle_suffixes = saga::util::getYamlArg<std::vector<std::string>>(node, "particle_suffixes", {}, message_prefix, verbose, yamlargout); // -> For MC matching with mc_match cut
    std::string combined_spin_state = saga::util::getYamlArg<std::string>(node, "combined_spin_state", "ss", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!

    //----------------------------------------------------------------------//
    // VAR_FORMULAS
    std::vector<std::vector<std::string>> var_formulas = saga::util::getYamlArg<std::vector<std::vector<std::string>>>(node, "var_formulas", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN HELICITY AND SPIN VARIABLES
    std::vector<std::string> categories_as_float = saga::util::getYamlArg<std::vector<std::string>>(node, "categories_as_float", {}, message_prefix, verbose, yamlargout);
    std::string helicity_name = saga::util::getYamlArg<std::string>(node, "helicity_name", "heli", message_prefix, verbose, yamlargout);
    std::string helicity_formula = saga::util::getYamlArg<std::string>(node, "helicity_formula", "-helicity", message_prefix, verbose, yamlargout); //NOTE: Make sure to flip helicity for RGA fall 2018 data and check if needed for other datasets.
    std::map<std::string,int> helicity_states = saga::util::getYamlArg<std::map<std::string,int>>(node, "helicity_states", {{"plus",1}, {"zero",0}, {"minus",-1}}, message_prefix, verbose, yamlargout);
    std::string tspin_name = saga::util::getYamlArg<std::string>(node, "tspin_name", "heli", message_prefix, verbose, yamlargout);
    std::string tspin_formula = saga::util::getYamlArg<std::string>(node, "tspin_formula", "-tspin", message_prefix, verbose, yamlargout); //NOTE: You will probably need to define this from a run-dependent branch loaded CSV.
    std::map<std::string,int> tspin_states = saga::util::getYamlArg<std::map<std::string,int>>(node, "tspin_states", {{"plus", 1}, {"zero", 0}, {"minus", -1}}, message_prefix, verbose, yamlargout);
    std::string htspin_name = saga::util::getYamlArg<std::string>(node, "htspin_name", "htspin", message_prefix, verbose, yamlargout);
    std::map<std::string,int> htspin_states = saga::util::getYamlArg<std::map<std::string,int>>(node, "htspin_states", {{"plus", 1}, {"zero", 0}, {"minus", -1}}, message_prefix, verbose, yamlargout);

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
        std::vector<std::string> binschemes_paths = saga::util::getYamlArg<std::vector<std::string>>(node, "binschemes_paths", {}, message_prefix, verbose, yamlargout);

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
        int nbatches = saga::util::getYamlArg<int>(node, "nbatches", 1, message_prefix, verbose, yamlargout);
        int ibatch   = saga::util::getYamlArg<int>(node, "ibatch", 0, message_prefix, verbose, yamlargout);
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
    std::vector<std::string> binvars = saga::util::getYamlArg<std::vector<std::string>>(node, "binvars", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> binvar_titles = saga::util::getYamlArg<std::vector<std::string>>(node, "binvar_titles", {}, message_prefix, verbose, yamlargout); //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    std::vector<std::vector<double>> binvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "binvar_lims", {}, message_prefix, verbose, yamlargout);
    std::vector<int> binvar_bins = saga::util::getYamlArg<std::vector<int>>(node, "binvar_bins", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN DEPOLARIZATION VARIABLES
    std::vector<std::string> depolvars = saga::util::getYamlArg<std::vector<std::string>>(node, "depolvars", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> depolvar_titles = saga::util::getYamlArg<std::vector<std::string>>(node, "depolvar_titles", {}, message_prefix, verbose, yamlargout); //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    std::vector<std::vector<double>> depolvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "depolvar_lims", {}, message_prefix, verbose, yamlargout);
    std::vector<int> depolvar_bins = saga::util::getYamlArg<std::vector<int>>(node, "depolvar_bins", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN ASYMMETRY FIT VARIABLES
    std::vector<std::string> resfitvars = saga::util::getYamlArg<std::vector<std::string>>(node, "resfitvars", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> resfitvar_titles = saga::util::getYamlArg<std::vector<std::string>>(node, "resfitvar_titles", {}, message_prefix, verbose, yamlargout); //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    std::vector<std::vector<double>> resfitvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "resfitvar_lims", {}, message_prefix, verbose, yamlargout);
    std::vector<int> resfitvar_bins = saga::util::getYamlArg<std::vector<int>>(node, "resfitvar_bins", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN MASS FIT VARIABLES
    std::vector<std::string> massfitvars = saga::util::getYamlArg<std::vector<std::string>>(node, "massfitvars", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfitvar_titles = saga::util::getYamlArg<std::vector<std::string>>(node, "massfitvar_titles", {}, message_prefix, verbose, yamlargout); //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    std::vector<std::vector<double>> massfitvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfitvar_lims", {}, message_prefix, verbose, yamlargout);
    std::vector<int> massfitvar_bins = saga::util::getYamlArg<std::vector<int>>(node, "massfitvar_bins", {}, message_prefix, verbose, yamlargout);

    // MC_SUFFIX
    std::string mc_suffix = saga::util::getYamlArg<std::string>(node, "mc_suffix", "_mc", message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN RESOLUTION FIT ARGUMENTS
    std::map<std::string,std::string> yamlfile_map = saga::util::getYamlArg<std::map<std::string,std::string>>(node, "yamlfile_map", {}, message_prefix, verbose, yamlargout);
    std::string pdf_name = saga::util::getYamlArg<std::string>(node, "pdf_name", "gauss", message_prefix, verbose, yamlargout); //NOTE: A mass fit and background correction will only be run if this is non-empty!
    std::string fitformula = saga::util::getYamlArg<std::string>(node, "fitformula", "gaus(x[0],x[1],x[2],x[3])", message_prefix, verbose, yamlargout); //NOTE: This is parsed by RooGenericPdf using TFormula
    std::vector<std::string> parnames = saga::util::getYamlArg<std::vector<std::string>>(node, "parnames", {"constant","mu","sigma"}, message_prefix, verbose, yamlargout);
    std::vector<std::string> partitles = saga::util::getYamlArg<std::vector<std::string>>(node, "partitles", {"C","#mu","#sigma"}, message_prefix, verbose, yamlargout);
    std::vector<std::string> parunits = saga::util::getYamlArg<std::vector<std::string>>(node, "parunits", {"","",""}, message_prefix, verbose, yamlargout);
    std::vector<double> parinits = saga::util::getYamlArg<std::vector<double>>(node, "parinits", {1.0, 0.0, 0.1}, message_prefix, verbose, yamlargout);
    std::vector<std::vector<double>> parlims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "parlims", {{1.0, 1.0}, {-1.0, 1.0}, {0.0, 1.0}}, message_prefix, verbose, yamlargout);
    double lg_text_size = saga::util::getYamlArg<double>(node, "lg_text_size", 0.04, message_prefix, verbose, yamlargout);
    double lg_margin = saga::util::getYamlArg<double>(node, "lg_margin", 0.1, message_prefix, verbose, yamlargout);
    double lg_ncols = saga::util::getYamlArg<double>(node, "lg_ncols", 1, message_prefix, verbose, yamlargout);
    bool use_sumw2error = saga::util::getYamlArg<bool>(node, "use_sumw2error", false, message_prefix, verbose, yamlargout);
    bool use_extended_nll = saga::util::getYamlArg<bool>(node, "use_extended_nll", true, message_prefix, verbose, yamlargout);
    bool use_binned_fit = saga::util::getYamlArg<bool>(node, "use_binned_fit", false, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------------------------------------//
    // Additional arguments
    std::string logpath = saga::util::getYamlArg<std::string>(node, "logpath", "out.txt", message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------------------------------------//
    // ANALYSIS
    //----------------------------------------------------------------------------------------------------//
    
    // Allow multithreading
    ROOT::EnableImplicitMT(nthreads);

    // Add all absolute variable limits to overall cuts
    cuts = saga::util::addLimitCuts(cuts,binvars,binvar_lims);
    cuts = saga::util::addLimitCuts(cuts,depolvars,depolvar_lims);
    cuts = saga::util::addLimitCuts(cuts,massfitvars,massfitvar_lims);
    yamlargout << message_prefix.c_str() << "cuts: "<<cuts.c_str() << std::endl;

    // Add resolution fit variable cuts to MC cuts since they depend on MC variables.
    mc_cuts = saga::util::addLimitCuts(mc_cuts,resfitvars,resfitvar_lims);
    yamlargout << message_prefix.c_str() << "mc_cuts: "<<mc_cuts.c_str() << std::endl;

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Define variables from formulas
    auto d2 = d.Define("__dummyvar__","(float)0.0"); //NOTE: Define a dummy variable to declare the data frame in this scope.
    for (int idx=0; idx<var_formulas.size(); idx++) {
        d2 = d2.Define(var_formulas[idx][0].c_str(),var_formulas[idx][1].c_str());
        yamlargout << message_prefix.c_str() << "Defined branch "<<var_formulas[idx][0].c_str()<<std::endl;
    }

    // Apply overall cuts AFTER defining depolarization and fit variables
    auto d2_filtered = d2.Filter(cuts.c_str());

    // Define angular difference variables
    if (mc_cuts.size()>0) {
        d2_filtered = saga::data::defineAngularDiffVars(d2_filtered, particle_suffixes, "theta", "phi", "_mc");
        d2_filtered = d2_filtered.Filter(mc_cuts.c_str());//TODO: Add output message about defined branches
    }

    // Define beam helicity and target spin variables
    std::string combined_spin_state_formula  = Form("(int)(10*(%s+1)+%s+1)",helicity_name.c_str(),tspin_name.c_str());
    d2_filtered = d2_filtered.Define(helicity_name.c_str(), helicity_formula.c_str())
                                .Define(tspin_name.c_str(), tspin_formula.c_str())
                                .Define(combined_spin_state.c_str(), combined_spin_state_formula.c_str());

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

        // Produce graphs of resolution fit parameters binned in given kinematic variable
        std::string scheme_name = Form("%s%s",baseoutpath.c_str(),binscheme_name.c_str());
        saga::resolution::getKinBinnedResolutions(
            scheme_name, // std::string                      scheme_name,
            d2_filtered, // RNode                            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
            "w", // std::string                      workspace_name,
            "workspace", // std::string                      workspace_title,

            // // parameters passed to data::createDataset()
            "dataset", // std::string                      dataset_name,
            "dataset", // std::string                      dataset_title,
            categories_as_float, // std::vector<std::string>         categories_as_float,
            helicity_name, // std::string                      helicity,
            helicity_states, // std::map<std::string,int>        helicity_states,
            tspin_name, // std::string                      tspin,
            tspin_states, // std::map<std::string,int>        tspin_states,
            htspin_name, // std::string                      htspin,
            htspin_states, // std::map<std::string,int>        htspin_states,
            combined_spin_state, // std::string                      combined_spin_state,
            bincuts, // std::map<int,std::string>        bincuts,
            scheme_binvars, // std::vector<std::string>         binvars,
            scheme_binvar_titles, // std::vector<std::string>         binvar_titles,
            scheme_binvar_lims, // std::vector<std::vector<double>> binvar_lims,
            scheme_binvar_bins, // std::vector<int>                 binvar_bins,
            depolvars, // std::vector<std::string>         depolvars,
            depolvar_titles, // std::vector<std::string>         depolvar_titles,
            depolvar_lims, // std::vector<std::vector<double>> depolvar_lims,
            depolvar_bins, // std::vector<int>                 depolvar_bins,
            resfitvars, // std::vector<std::string>         resfitvars,
            resfitvar_titles, // std::vector<std::string>         resfitvar_titles,
            resfitvar_lims, // std::vector<std::vector<double>> resfitvar_lims,
            resfitvar_bins, // std::vector<int>                 resfitvar_bins,
            massfitvars, // std::vector<std::string>         massfitvars,
            massfitvar_titles, // std::vector<std::string>         massfitvar_titles,
            massfitvar_lims, // std::vector<std::vector<double>> massfitvar_lims,
            massfitvar_bins, // std::vector<int>                 massfitvar_bins,

            // // parameters passed to saga::signal::fitResolution() //TODO: Add init fit parameter value and limits arguments here...assuming you always want a chebychev polynomial background...
            yamlfile_map, // std::map<std::string,std::string> yamlfile_map,
            pdf_name, // std::string                       pdf_name,
            fitformula, // std::string                       fitformula,
            parnames, // std::vector<std::string>          parnames,
            partitles, // std::vector<std::string>          partitles,
            parunits, // std::vector<std::string>          parunits,
            parinits, // std::vector<double>               parinits,
            parlims, // std::vector<std::vector<double>>  parlims,

            // // Parameters passed to signal::fitResolution()
            lg_text_size, // double                           lg_text_size     = 0.04,
            lg_margin, // double                           lg_margin        = 0.1,
            lg_ncols, // int                              lg_ncols         = 1,
            use_sumw2error, // bool                             use_sumw2error   = false,
            use_extended_nll, // bool                             use_extended_nll = true,
            use_binned_fit, // bool                             use_binned_fit   = false,

            // // Ouput stream
            out // std::ostream                    &out              = std::cout
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
