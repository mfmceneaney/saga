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
#include <analysis.h>
#include <bins.h>
#include <data.h>
#include <util.h>

void execute(const YAML::Node& node) {

    // Process arguments
    std::string   message_prefix  = "INFO: ";
    bool          verbose         = true;
    std::ostream &yamlargout      = std::cout;

    //----------------------------------------------------------------------//
    // BEGIN ARGUMENTS
    std::string baseoutpath = saga::util::getYamlArg<std::string>(node,"baseoutpath","",message_prefix,verbose,yamlargout); //NOTE: This will be prepended to the default output path like so: `<baseoutpath><binscheme_name>.csv`.
    std::string inpath = saga::util::getYamlArg<std::string>(node,"inpath","",message_prefix,verbose,yamlargout);
    std::string tree = saga::util::getYamlArg<std::string>(node,"tree","t",message_prefix,verbose,yamlargout);
    int nthreads = saga::util::getYamlArg<int>(node,"nthreads",1,message_prefix,verbose,yamlargout);
    std::string cuts = saga::util::getYamlArg<std::string>(node,"cuts","",message_prefix,verbose,yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN ASYMMETRY INJECTION ARGUMENTS
    bool inject_asym = saga::util::getYamlArg<bool>(node, "inject_asym", false, message_prefix, verbose, yamlargout);
    int inject_seed = saga::util::getYamlArg<int>(node, "inject_seed", 2, message_prefix, verbose, yamlargout);
    std::string mc_cuts = saga::util::getYamlArg<std::string>(node,"mc_cuts","Q2>1",message_prefix,verbose,yamlargout); //NOTE: This may not be empty!
    std::vector<double> sgasyms = saga::util::getYamlArg<std::vector<double>>(node, "sgasyms", {}, message_prefix, verbose, yamlargout);
    std::vector<double> bgasyms = saga::util::getYamlArg<std::vector<double>>(node, "bgasyms", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> particle_suffixes = saga::util::getYamlArg<std::vector<std::string>>(node, "particle_suffixes", {}, message_prefix, verbose, yamlargout); // -> For MC matching with mc_sg_match cut
    std::string mc_sg_match_name = saga::util::getYamlArg<std::string>(node,"mc_sg_match_name","mc_sg_match",message_prefix,verbose,yamlargout); //NOTE: This may not be empty!
    std::string mc_sg_match_formula = saga::util::getYamlArg<std::string>(node,"mc_sg_match_formula","(bool)true",message_prefix,verbose,yamlargout); //NOTE: This may not be empty!
    std::string phi_s_original_name = saga::util::getYamlArg<std::string>(node, "phi_s_original_name", "phi_s_up", message_prefix, verbose, yamlargout);
    std::string phi_s_original_name_dn = saga::util::getYamlArg<std::string>(node, "phi_s_original_name_dn", "phi_s_dn", message_prefix, verbose, yamlargout);
    std::string phi_s_injected_name = saga::util::getYamlArg<std::string>(node, "phi_s_injected_name", "phi_s_injected", message_prefix, verbose, yamlargout);
    std::string fsgasyms_xs_uu_name = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_uu_name", "fsgasyms_xs_uu", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fsgasyms_xs_uu_formula = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_uu_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fsgasyms_xs_pu_name = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_pu_name", "fsgasyms_xs_pu", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fsgasyms_xs_pu_formula = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_pu_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fsgasyms_xs_up_name = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_up_name", "fsgasyms_xs_up", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fsgasyms_xs_up_formula = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_up_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fsgasyms_xs_pp_name = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_pp_name", "fsgasyms_xs_pp", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fsgasyms_xs_pp_formula = saga::util::getYamlArg<std::string>(node, "fsgasyms_xs_pp_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_uu_name = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_uu_name", "fbgasyms_xs_uu", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_uu_formula = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_uu_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_pu_name = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_pu_name", "fbgasyms_xs_pu", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_pu_formula = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_pu_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_up_name = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_up_name", "fbgasyms_xs_up", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_up_formula = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_up_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_pp_name = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_pp_name", "fbgasyms_xs_pp", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string fbgasyms_xs_pp_formula = saga::util::getYamlArg<std::string>(node, "fbgasyms_xs_pp_formula", "(float)0.0", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!
    std::string combined_spin_state = saga::util::getYamlArg<std::string>(node, "combined_spin_state", "ss", message_prefix, verbose, yamlargout); //NOTE: This may not be empty!

    //----------------------------------------------------------------------//
    // BEGIN RUN-DEPENDENT CSV VARIABLE ARGUMENTS
    std::vector<std::string> rdf_key_cols = saga::util::getYamlArg<std::vector<std::string>>(node, "rdf_key_cols", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> csv_paths = saga::util::getYamlArg<std::vector<std::string>>(node, "csv_paths", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> csv_key_cols = saga::util::getYamlArg<std::vector<std::string>>(node, "csv_key_cols", {}, message_prefix, verbose, yamlargout);
    std::vector<std::vector<std::string>> col_names = saga::util::getYamlArg<std::vector<std::vector<std::string>>>(node, "col_names", {}, message_prefix, verbose, yamlargout);
    std::map<std::string,std::string> col_aliases = saga::util::getYamlArg<std::map<std::string,std::string>>(node, "col_aliases", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // VAR_FORMULAS
    std::vector<std::vector<std::string>> var_formulas = saga::util::getYamlArg<std::vector<std::vector<std::string>>>(node, "var_formulas", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN HELICITY AND SPIN VARIABLES
    bool use_categories_as_float = saga::util::getYamlArg<bool>(node, "use_categories_as_float", false, message_prefix, verbose, yamlargout);
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
    std::vector<std::string> asymfitvars = saga::util::getYamlArg<std::vector<std::string>>(node, "asymfitvars", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> asymfitvar_titles = saga::util::getYamlArg<std::vector<std::string>>(node, "asymfitvar_titles", {}, message_prefix, verbose, yamlargout); //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    std::vector<std::vector<double>> asymfitvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "asymfitvar_lims", {}, message_prefix, verbose, yamlargout);
    std::vector<int> asymfitvar_bins = saga::util::getYamlArg<std::vector<int>>(node, "asymfitvar_bins", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN MASS FIT VARIABLES
    std::vector<std::string> massfitvars = saga::util::getYamlArg<std::vector<std::string>>(node, "massfitvars", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfitvar_titles = saga::util::getYamlArg<std::vector<std::string>>(node, "massfitvar_titles", {}, message_prefix, verbose, yamlargout); //NOTE: DEFAULT TO ACTUAL VARIABLE NAMES
    std::vector<std::vector<double>> massfitvar_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfitvar_lims", {}, message_prefix, verbose, yamlargout);
    std::vector<int> massfitvar_bins = saga::util::getYamlArg<std::vector<int>>(node, "massfitvar_bins", {}, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN ASYMMETRY FIT ARGUMENTS
    double bpol = saga::util::getYamlArg<double>(node, "bpol", 1.0, message_prefix, verbose, yamlargout); // Average beam polarization
    double tpol = saga::util::getYamlArg<double>(node, "tpol", 1.0, message_prefix, verbose, yamlargout); // Average target polarization
    std::string asymfit_formula_uu = saga::util::getYamlArg<std::string>(node, "asymfit_formula_uu", "", message_prefix, verbose, yamlargout);
    std::string asymfit_formula_pu = saga::util::getYamlArg<std::string>(node, "asymfit_formula_pu", "", message_prefix, verbose, yamlargout);
    std::string asymfit_formula_up = saga::util::getYamlArg<std::string>(node, "asymfit_formula_up", "", message_prefix, verbose, yamlargout);
    std::string asymfit_formula_pp = saga::util::getYamlArg<std::string>(node, "asymfit_formula_pp", "", message_prefix, verbose, yamlargout);
    std::vector<double> asymfitpar_inits = saga::util::getYamlArg<std::vector<double>>(node, "asymfitpar_inits", {}, message_prefix, verbose, yamlargout);
    std::vector<std::vector<double>> asymfitpar_initlims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "asymfitpar_initlims", {}, message_prefix, verbose, yamlargout);
    bool use_sumw2error = saga::util::getYamlArg<bool>(node, "use_sumw2error", true, message_prefix, verbose, yamlargout);
    bool use_average_depol = saga::util::getYamlArg<bool>(node, "use_average_depol", true, message_prefix, verbose, yamlargout);   
    bool use_extended_nll = saga::util::getYamlArg<bool>(node, "use_extended_nll", false, message_prefix, verbose, yamlargout);
    bool use_binned_fit = saga::util::getYamlArg<bool>(node, "use_binned_fit", false, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN MASS FIT ARGUMENTS
    std::map<std::string,std::string> massfit_yamlfile_map = saga::util::getYamlArg<std::map<std::string,std::string>>(node, "massfit_yamlfile_map", {}, message_prefix, verbose, yamlargout);
    std::string massfit_pdf_name = saga::util::getYamlArg<std::string>(node, "massfit_pdf_name", "", message_prefix, verbose, yamlargout); //NOTE: A mass fit and background correction will only be run if this is non-empty!
    std::string massfit_formula_sg = saga::util::getYamlArg<std::string>(node, "massfit_formula_sg", "gaus(x[0],x[1],x[2])", message_prefix, verbose, yamlargout); //NOTE: This is parsed by RooGenericPdf using TFormula
    std::string massfit_formula_bg = saga::util::getYamlArg<std::string>(node, "massfit_formula_bg", "cb2(x[0], x[1], x[2])", message_prefix, verbose, yamlargout); //NOTE: This is parsed by RooGenericPdf using TFormula
    std::string massfit_sgYield_name = saga::util::getYamlArg<std::string>(node, "massfit_sgYield_name", "sgYield", message_prefix, verbose, yamlargout);
    std::string massfit_bgYield_name = saga::util::getYamlArg<std::string>(node, "massfit_bgYield_name", "bgYield", message_prefix, verbose, yamlargout);
    double massfit_initsgfrac = saga::util::getYamlArg<double>(node, "massfit_initsgfrac", 0.1, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfit_parnames_sg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parnames_sg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfit_partitles_sg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_partitles_sg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfit_parunits_sg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parunits_sg", {}, message_prefix, verbose, yamlargout);
    std::vector<double> massfit_parinits_sg = saga::util::getYamlArg<std::vector<double>>(node, "massfit_parinits_sg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::vector<double>> massfit_parlims_sg = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfit_parlims_sg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfit_parnames_bg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parnames_bg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfit_partitles_bg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_partitles_bg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::string> massfit_parunits_bg = saga::util::getYamlArg<std::vector<std::string>>(node, "massfit_parunits_bg", {}, message_prefix, verbose, yamlargout);
    std::vector<double> massfit_parinits_bg = saga::util::getYamlArg<std::vector<double>>(node, "massfit_parinits_bg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::vector<double>> massfit_parlims_bg = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfit_parlims_bg", {}, message_prefix, verbose, yamlargout);
    std::vector<std::vector<double>> massfit_sgregion_lims = saga::util::getYamlArg<std::vector<std::vector<double>>>(node, "massfit_sgregion_lims", {}, message_prefix, verbose, yamlargout);
    double massfit_lg_text_size = saga::util::getYamlArg<double>(node, "massfit_lg_text_size", 0.04, message_prefix, verbose, yamlargout);
    double massfit_lg_margin = saga::util::getYamlArg<double>(node, "massfit_lg_margin", 0.1, message_prefix, verbose, yamlargout);
    double massfit_lg_ncols = saga::util::getYamlArg<double>(node, "massfit_lg_ncols", 1, message_prefix, verbose, yamlargout);
    bool massfit_plot_bg_pars = saga::util::getYamlArg<bool>(node, "massfit_plot_bg_pars", false, message_prefix, verbose, yamlargout);
    bool massfit_use_sumw2error = saga::util::getYamlArg<bool>(node, "massfit_use_sumw2error", false, message_prefix, verbose, yamlargout);
    bool massfit_use_extended_nll = saga::util::getYamlArg<bool>(node, "massfit_use_extended_nll", true, message_prefix, verbose, yamlargout);
    bool massfit_use_binned_fit = saga::util::getYamlArg<bool>(node, "massfit_use_binned_fit", false, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN SPLOT ARGUMENTS
    bool use_splot = saga::util::getYamlArg<bool>(node, "use_splot", true, message_prefix, verbose, yamlargout);

    //----------------------------------------------------------------------//
    // BEGIN SIDEBAND SUBTRACTION ARGUMENTS
    std::string massfit_sgcut = saga::util::getYamlArg<std::string>(node, "massfit_sgcut", "", message_prefix, verbose, yamlargout);
    std::string massfit_bgcut = saga::util::getYamlArg<std::string>(node, "massfit_bgcut", "", message_prefix, verbose, yamlargout);
    bool use_sb_subtraction = saga::util::getYamlArg<bool>(node, "use_sb_subtraction", false, message_prefix, verbose, yamlargout);
    bool use_binned_sb_bgfracs = saga::util::getYamlArg<bool>(node, "use_binned_sb_bgfracs", false, message_prefix, verbose, yamlargout);
    std::string bgfracvar = saga::util::getYamlArg<std::string>(node, "bgfracvar", "", message_prefix, verbose, yamlargout);
    std::vector<double> bgfracvar_lims = saga::util::getYamlArg<std::vector<double>>(node, "bgfracvar_lims", {}, message_prefix, verbose, yamlargout);
    int bgfrac_idx = saga::util::getYamlArg<int>(node, "bgfrac_idx", 0, message_prefix, verbose, yamlargout);

    // ASYMFITVAR_BINSCHEMES
    std::map<std::string,std::map<int,std::string>> asymfitvar_bincuts_map;
    if (node["asymfitvar_binschemes"] && node["asymfitvar_binschemes"].IsMap()) {

        // Get bin scheme node and get bin cuts maps
        asymfitvar_bincuts_map = saga::bins::getBinCutsMap(node["asymfitvar_binschemes"]);

        // Show the asymmetry fit variables bin cuts map
        if (verbose) {
            for (auto it = asymfitvar_bincuts_map.begin(); it != asymfitvar_bincuts_map.end(); ++it) {
                yamlargout << message_prefix.c_str() << "asymfitvar_bincuts_map["<<it->first<<"]: { \n";
                for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                    yamlargout <<"\t\t"<< it2->first<<": "<<it2->second<<", \n";
                }
                yamlargout << "}" << std::endl;
            }
        }
    }

    //----------------------------------------------------------------------------------------------------//
    // Additional arguments
    std::string logpath = saga::util::getYamlArg<std::string>(node, "logpath", "out.txt", message_prefix, verbose, yamlargout);
    bool dump_dataset = saga::util::getYamlArg<bool>(node, "dump_dataset", false, message_prefix, verbose, yamlargout);
    std::vector<std::string> dump_vars = saga::util::getYamlArg<std::vector<std::string>>(node, "dump_vars", {}, message_prefix, verbose, yamlargout); //NOTE: If empty ALL variables from dataset will be dumped if dump_dataset==true.

    //----------------------------------------------------------------------------------------------------//
    // Create list of categories to use as float  //NOTE: Ordering is important!
    std::vector<std::string> categories_as_float;
    if (use_categories_as_float) {
        if (asymfit_formula_pu!="" || asymfit_formula_pp!="") {
            categories_as_float.push_back(helicity_name);
        }
        if (asymfit_formula_up!="" || asymfit_formula_pp!="") {
            categories_as_float.push_back(tspin_name);
        }
    }

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
    yamlargout << message_prefix.c_str() << "cuts: "<<cuts.c_str() << std::endl;

    // Print actual inputs for debugging
    yamlargout << message_prefix.c_str() << "DEBUG: creating RDataFrame with tree='" << tree << "' inpath='" << inpath << "'" << std::endl;

    // Create RDataFrame
    ROOT::RDataFrame d(tree, inpath);

    // Create branch names for mc variables assuming they all append `_mc` to the corresponding data branch name
    std::vector<std::string> asymfitvars_mc;
    for (int idx=0; idx<asymfitvars.size(); idx++) {
        asymfitvars_mc.push_back(Form("%s_mc",asymfitvars[idx].c_str()));
        yamlargout << message_prefix.c_str() << "Defined MC variable : " << asymfitvars_mc[idx].c_str() << std::endl;
    }
    std::vector<std::string> depolvars_mc;
    for (int idx=0; idx<depolvars.size(); idx++) {
        depolvars_mc.push_back(Form("%s_mc",depolvars[idx].c_str()));
        yamlargout << message_prefix.c_str() << "Defined MC variable : " << depolvars_mc[idx].c_str() << std::endl;
    }

    // Define the tspin==+1 and tspin==-1 A_{UT} asymmetry names
    std::string fsgasyms_xs_uu_pos_name = Form("_%s_pos",fsgasyms_xs_uu_name.c_str());
    std::string fsgasyms_xs_uu_neg_name = Form("_%s_neg",fsgasyms_xs_uu_name.c_str());
    std::string fsgasyms_xs_pu_pos_name = Form("_%s_pos",fsgasyms_xs_pu_name.c_str());
    std::string fsgasyms_xs_pu_neg_name = Form("_%s_neg",fsgasyms_xs_pu_name.c_str());
    std::string fbgasyms_xs_uu_pos_name = Form("_%s_pos",fbgasyms_xs_uu_name.c_str());
    std::string fbgasyms_xs_uu_neg_name = Form("_%s_neg",fbgasyms_xs_uu_name.c_str());
    std::string fbgasyms_xs_pu_pos_name = Form("_%s_pos",fbgasyms_xs_pu_name.c_str());
    std::string fbgasyms_xs_pu_neg_name = Form("_%s_neg",fbgasyms_xs_pu_name.c_str());

    // Define the tspin==+1 and tspin==-1 A_{UT} asymmetry formulas
    std::string fsgasyms_xs_uu_pos_formula = fsgasyms_xs_uu_formula;
    std::string fsgasyms_xs_uu_neg_formula = fsgasyms_xs_uu_formula;
    std::string fsgasyms_xs_pu_pos_formula = fsgasyms_xs_pu_formula;
    std::string fsgasyms_xs_pu_neg_formula = fsgasyms_xs_pu_formula;
    std::string fbgasyms_xs_uu_pos_formula = fbgasyms_xs_uu_formula;
    std::string fbgasyms_xs_uu_neg_formula = fbgasyms_xs_uu_formula;
    std::string fbgasyms_xs_pu_pos_formula = fbgasyms_xs_pu_formula;
    std::string fbgasyms_xs_pu_neg_formula = fbgasyms_xs_pu_formula;

    // Replace variable names in XS formulas with their MC counter parts and asymmetries with their injected values
    if (inject_asym) {
        // Find and replace asymmetry names with injected values in fsgasyms_xs : example string fsgasyms_xs="0.747*depolvars_mc0*sgasym0*fitvar1_mc"
        for (int idx=0; idx<asymfitvars.size(); idx++) {
            saga::util::replaceAll(fsgasyms_xs_uu_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_pu_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_up_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_pp_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
        }
        for (int idx=0; idx<sgasyms.size(); idx++) {
            saga::util::replaceAll(fsgasyms_xs_uu_formula, Form("sgasyms[%d]",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fsgasyms_xs_pu_formula, Form("sgasyms[%d]",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fsgasyms_xs_up_formula, Form("sgasyms[%d]",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fsgasyms_xs_pp_formula, Form("sgasyms[%d]",idx), Form("%.8f",sgasyms[idx])); // Replace sgasyms[idx] with actual injected asymmetry value
        }
        for (int idx=0; idx<depolvars.size(); idx++) {
            saga::util::replaceAll(fsgasyms_xs_uu_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_pu_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_up_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fsgasyms_xs_pp_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
        }
        yamlargout << message_prefix.c_str() << "Updated " << fsgasyms_xs_uu_name.c_str() << " = " << fsgasyms_xs_uu_formula.c_str() << std::endl;
        yamlargout << message_prefix.c_str() << "Updated " << fsgasyms_xs_pu_name.c_str() << " = " << fsgasyms_xs_pu_formula.c_str() << std::endl;
        yamlargout << message_prefix.c_str() << "Updated " << fsgasyms_xs_up_name.c_str() << " = " << fsgasyms_xs_up_formula.c_str() << std::endl;
        yamlargout << message_prefix.c_str() << "Updated " << fsgasyms_xs_pp_name.c_str() << " = " << fsgasyms_xs_pp_formula.c_str() << std::endl;

        // Find and replace placeholder variable names with actual values in fbgasyms_xs : example string fbgasyms_xs="0.747*depolvars_mc0*bgasym0*fitvar1_mc"
        for (int idx=0; idx<asymfitvars.size(); idx++) {
            saga::util::replaceAll(fbgasyms_xs_uu_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_pu_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_up_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_pp_formula, asymfitvars[idx].c_str(), asymfitvars_mc[idx].c_str()); // Replace asymfitvars_mc[idx] with actual branch name
        }
        for (int idx=0; idx<bgasyms.size(); idx++) {
            saga::util::replaceAll(fbgasyms_xs_uu_formula, Form("bgasyms[%d]",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fbgasyms_xs_pu_formula, Form("bgasyms[%d]",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fbgasyms_xs_up_formula, Form("bgasyms[%d]",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
            saga::util::replaceAll(fbgasyms_xs_pp_formula, Form("bgasyms[%d]",idx), Form("%.8f",bgasyms[idx])); // Replace bgasyms[idx] with actual injected asymmetry value
        }
        for (int idx=0; idx<depolvars.size(); idx++) {
            saga::util::replaceAll(fbgasyms_xs_uu_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_pu_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_up_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
            saga::util::replaceAll(fbgasyms_xs_pp_formula, depolvars[idx].c_str(), depolvars_mc[idx].c_str()); // Replace depolvars_mc[idx] with actual branch name
        }
        yamlargout << message_prefix.c_str() << "Updated " << fbgasyms_xs_uu_name.c_str() << " = " << fbgasyms_xs_uu_formula.c_str() << std::endl;
        yamlargout << message_prefix.c_str() << "Updated " << fbgasyms_xs_pu_name.c_str() << " = " << fbgasyms_xs_pu_formula.c_str() << std::endl;
        yamlargout << message_prefix.c_str() << "Updated " << fbgasyms_xs_up_name.c_str() << " = " << fbgasyms_xs_up_formula.c_str() << std::endl;
        yamlargout << message_prefix.c_str() << "Updated " << fbgasyms_xs_pp_name.c_str() << " = " << fbgasyms_xs_pp_formula.c_str() << std::endl;

        // Reassign the tspin==+1 and tspin==-1 A_{UT} asymmetry formulas
        fsgasyms_xs_uu_pos_formula = fsgasyms_xs_uu_formula;
        fsgasyms_xs_uu_neg_formula = fsgasyms_xs_uu_formula;
        fsgasyms_xs_pu_pos_formula = fsgasyms_xs_pu_formula;
        fsgasyms_xs_pu_neg_formula = fsgasyms_xs_pu_formula;
        fbgasyms_xs_uu_pos_formula = fbgasyms_xs_uu_formula;
        fbgasyms_xs_uu_neg_formula = fbgasyms_xs_uu_formula;
        fbgasyms_xs_pu_pos_formula = fbgasyms_xs_pu_formula;
        fbgasyms_xs_pu_neg_formula = fbgasyms_xs_pu_formula;

        // Check if you have phi_s variable names to replace
        if (phi_s_original_name!="" && phi_s_original_name_dn!="") {

            // And put the appropriate sign on the phi_s variable for A_{UT} and A_{LT} asymmetries
            std::string phi_s_original_name_mc = Form("%s_mc", phi_s_original_name.c_str());
            std::string phi_s_original_name_mc_neg = Form("%s_mc",phi_s_original_name_dn.c_str());//NOTE: Need to multiply spin vector by -1 and then take phi in g*N CM frame which equates to flipping the angle vector, i.e., phi -> pi + phi = phi', if (phi'>2pi) phi' = phi' - 2pi
            saga::util::replaceAll(fsgasyms_xs_uu_neg_formula, phi_s_original_name_mc, phi_s_original_name_mc_neg);
            saga::util::replaceAll(fbgasyms_xs_uu_neg_formula, phi_s_original_name_mc, phi_s_original_name_mc_neg);
            saga::util::replaceAll(fsgasyms_xs_pu_neg_formula, phi_s_original_name_mc, phi_s_original_name_mc_neg);
            saga::util::replaceAll(fbgasyms_xs_pu_neg_formula, phi_s_original_name_mc, phi_s_original_name_mc_neg);
            yamlargout << message_prefix.c_str() << "Updated " << fsgasyms_xs_uu_neg_name.c_str() << " = " << fsgasyms_xs_uu_neg_formula.c_str() << std::endl;
            yamlargout << message_prefix.c_str() << "Updated " << fbgasyms_xs_uu_neg_name.c_str() << " = " << fbgasyms_xs_uu_neg_formula.c_str() << std::endl;
            yamlargout << message_prefix.c_str() << "Updated " << fsgasyms_xs_pu_neg_name.c_str() << " = " << fsgasyms_xs_pu_neg_formula.c_str() << std::endl;
            yamlargout << message_prefix.c_str() << "Updated " << fbgasyms_xs_pu_neg_name.c_str() << " = " << fbgasyms_xs_pu_neg_formula.c_str() << std::endl;
        }
    }

    // Define variables from formulas
    auto d2 = d.Define("__dummyvar__","(float)0.0"); //NOTE: Define a dummy variable to declare the data frame in this scope.
    for (int idx=0; idx<var_formulas.size(); idx++) {
        d2 = d2.Define(var_formulas[idx][0].c_str(),var_formulas[idx][1].c_str());
        yamlargout << message_prefix.c_str() << "Defined branch "<<var_formulas[idx][0].c_str()<<std::endl;
    }

    // Apply overall cuts AFTER defining depolarization and fit variables
    auto d2_filtered = (!inject_asym) ? d2.Filter(cuts.c_str()) :
                    d2.Filter(Form("(%s) && (%s)",cuts.c_str(),mc_cuts.c_str()));

    // Define angular difference variables
    if (inject_asym) d2_filtered = saga::data::defineAngularDiffVars(d2_filtered, particle_suffixes, "theta", "phi", "_mc");
    //TODO: Add output message about defined branches

    // Define signal matching condition, and XS values branches
    if (inject_asym) {
        d2_filtered = d2_filtered.Define(mc_sg_match_name.c_str(),mc_sg_match_formula.c_str()); //TODO: Throw error if formulas are empty
        d2_filtered = d2_filtered.Define(fsgasyms_xs_uu_pos_name.c_str(),fsgasyms_xs_uu_pos_formula.c_str());
        d2_filtered = d2_filtered.Define(fsgasyms_xs_uu_neg_name.c_str(),fsgasyms_xs_uu_neg_formula.c_str());
        d2_filtered = d2_filtered.Define(fsgasyms_xs_pu_pos_name.c_str(),fsgasyms_xs_pu_pos_formula.c_str());
        d2_filtered = d2_filtered.Define(fsgasyms_xs_pu_neg_name.c_str(),fsgasyms_xs_pu_neg_formula.c_str());
        d2_filtered = d2_filtered.Define(fsgasyms_xs_up_name.c_str(),fsgasyms_xs_up_formula.c_str());
        d2_filtered = d2_filtered.Define(fsgasyms_xs_pp_name.c_str(),fsgasyms_xs_pp_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_uu_pos_name.c_str(),fbgasyms_xs_uu_pos_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_uu_neg_name.c_str(),fbgasyms_xs_uu_neg_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_pu_pos_name.c_str(),fbgasyms_xs_pu_pos_formula.c_str());
        d2_filtered = d2_filtered.Define(fbgasyms_xs_pu_neg_name.c_str(),fbgasyms_xs_pu_neg_formula.c_str());
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
    std::string combined_spin_state_formula  = Form("(int)(10*(%s+1)+%s+1)",helicity_name.c_str(),tspin_name.c_str());
    auto frame = (!inject_asym) ?
                    d2_filtered.Define(helicity_name.c_str(), helicity_formula.c_str())
                                .Define(tspin_name.c_str(), tspin_formula.c_str())
                                .Define(combined_spin_state.c_str(), combined_spin_state_formula.c_str()) :
                    saga::data::injectAsym(
                        d2_filtered,
                        inject_seed,
                        bpol,
                        tpol,
                        mc_sg_match_name,
                        fsgasyms_xs_uu_pos_name,
                        fsgasyms_xs_uu_neg_name,
                        fsgasyms_xs_pu_pos_name,
                        fsgasyms_xs_pu_neg_name,
                        fsgasyms_xs_up_name,
                        fsgasyms_xs_pp_name,
                        fbgasyms_xs_uu_pos_name,
                        fbgasyms_xs_uu_neg_name,
                        fbgasyms_xs_pu_pos_name,
                        fbgasyms_xs_pu_neg_name,
                        fbgasyms_xs_up_name,
                        fbgasyms_xs_pp_name,
                        combined_spin_state,
                        helicity_name,
                        tspin_name,
                        phi_s_original_name,
                        phi_s_original_name_dn,
                        phi_s_injected_name
                    );
    //TODO: Add output message about defined branches

    // Reassign the phi_s fit variable name if present and injecting an asymmetry
    if (inject_asym && phi_s_original_name!="") {

        // Make sure the new name is initialized
        if (phi_s_injected_name=="") phi_s_injected_name = Form("%s_injected",phi_s_original_name.c_str());

        // Loop asymmetry fit variables and find the first match
        for (int idx=0; idx<asymfitvars.size(); idx++) {
            if (asymfitvars[idx]==phi_s_original_name) {
                asymfitvars[idx] = phi_s_injected_name;
                break;
            }
        }
    }

    // Make sure injection values are all computed before running analysis
    if (inject_asym) {
        double my_testvar  = (double)*frame.Mean(combined_spin_state.c_str());
        double my_testvar1 = (double)*frame.Mean(helicity_name.c_str());
        double my_testvar2 = (double)*frame.Mean(tspin_name.c_str());
    }

    // Dump dataset to ROOT file and exit
    if (dump_dataset) {
        std::string out_ds_path = Form("%sdataset.root", baseoutpath.c_str());
        yamlargout << message_prefix.c_str() << "Dumping dataset to: "<<out_ds_path.c_str()<<std::endl;
        if (dump_vars.size()==0) frame.Snapshot(tree.c_str(), out_ds_path.c_str());
        else frame.Snapshot(tree.c_str(), out_ds_path.c_str(), dump_vars);

        return;
    }

    // Create output dataset ROOT file with all variables
    yamlargout << message_prefix.c_str() << "Creating output dataset with all variables at: "<<out_ds_path.c_str()<<std::endl;
    frame = frame.Snapshot(tree.c_str(), out_ds_path.c_str(), frame.GetColumnNames());

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
        if (use_binned_sb_bgfracs) {
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
            scheme_name, // std::string                      scheme_name,
            frame, // RNode                            frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
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
            asymfitvars, // std::vector<std::string>         asymfitvars,
            asymfitvar_titles, // std::vector<std::string>         asymfitvar_titles,
            asymfitvar_lims, // std::vector<std::vector<double>> asymfitvar_lims,
            asymfitvar_bins, // std::vector<int>                 asymfitvar_bins,
            massfitvars, // std::vector<std::string>         massfitvars,
            massfitvar_titles, // std::vector<std::string>         massfitvar_titles,
            massfitvar_lims, // std::vector<std::vector<double>> massfitvar_lims,
            massfitvar_bins, // std::vector<int>                 massfitvar_bins,

            // // parameterss passed to analysis::fitAsym()
            bpol, // double                           bpol,
            tpol, // double                           tpol,
            asymfit_formula_uu, // std::string                      asymfit_formula_uu,
            asymfit_formula_pu, // std::string                      asymfit_formula_pu,
            asymfit_formula_up, // std::string                      asymfit_formula_up,
            asymfit_formula_pp, // std::string                      asymfit_formula_pp,
            asymfitpar_inits, // std::vector<double>              asymfitpar_inits,
            asymfitpar_initlims, // std::vector<std::vector<double>> asymfitpar_initlims,
            use_sumw2error, // bool                             use_sumw2error,
            use_average_depol, // bool                             use_average_depol,
            use_extended_nll, // bool                             use_extended_nll,
            use_binned_fit, // bool                             use_binned_fit,

            // // parameters passed to saga::signal::fitMass() //TODO: Add init fit parameter value and limits arguments here...assuming you always want a chebychev polynomial background...
            massfit_yamlfile_map, 
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

            // // Parameters passed to analysis::applySPlots()
            use_splot, // bool                             use_splot,

            // // Parameters used for sb subtraction
            massfit_sgcut, // std::string                      massfit_sgcut,
            massfit_bgcut, // std::string                      massfit_bgcut,
            use_sb_subtraction, // bool                             use_sb_subtraction,
            use_binned_sb_bgfracs, // bool                             use_binned_sb_bgfracs,
            asymfitvar_bincuts, // std::map<int,std::string>        asymfitvar_bincuts,
            bgfracvar, // std::string                      bgfracvar,
            bgfracvar_lims, // std::vector<double>              bgfracvar_lims,
            bgfrac_idx, // int                              bgfrac_idx               = 0,

            // // Parameters passed to signal::fitMass()
            massfit_lg_text_size, // double                           massfit_lg_text_size     = 0.04,
            massfit_lg_margin, // double                           massfit_lg_margin        = 0.1,
            massfit_lg_ncols, // int                              massfit_lg_ncols         = 1,
            massfit_plot_bg_pars, // bool                             massfit_plot_bg_pars     = false,
            massfit_use_sumw2error, // bool                             massfit_use_sumw2error   = false,
            massfit_use_extended_nll, // bool                             massfit_use_extended_nll = true,
            massfit_use_binned_fit, // bool                             massfit_use_binned_fit   = false,

            // // Ouput stream
            out // std::ostream                    &out                      = std::cout
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
