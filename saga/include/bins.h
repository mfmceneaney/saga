#include <iostream>
#include <memory>
#include <string>
#include <map>

// YAML Includes
#include <yaml-cpp/yaml.h>

// ROOT Includes
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include <TLatex.h>

// Local includes
#include <util.h>

#pragma once

/**
* @file
* @author Matthew McEneaney
* @date 16/Dec./24
* @version 0.0.0
* @brief Find ideal bin limits and compute bin migration fractions.
*/

namespace saga {

namespace bins {

/**
* @brief Find bin limits for equal bin statistics
*
* Given the desired number of bins in a distribution, find the bin limits
* that will ensure all bins have roughly equal statistics.
*
* @param frame ROOT RDataFrame with which to find bin limits
* @param varname Bin variable name
* @param nbins Number of bins
*
* @return Bin limits
*/
std::vector<double> findBinLims(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        std::string varname,
        const int nbins
    ) {

    // Initialize vectors
    std::vector<double> binmeans;
    std::vector<double> binmeans_err;
    std::vector<double> binlims;

    // Get overall limits and count and the necessary bin count with nbins
    double varmin = (double)*frame.Min(varname);
    double varmax = (double)*frame.Max(varname);
    int    count  = (int)*frame.Count();
    double targetbincount = (double)count/nbins;

    // Set initial step size and add initial bin limit
    double step = (double)(varmax-varmin)/nbins;
    binlims.push_back(varmin);

    // Loop each bin and find its upper limit
    int nsteps;
    double threshold = 0.01*targetbincount; //NOTE: This is a threshold in counts, which will be an int.
    for (int i=0; i<nbins; i++) {

        // Set the bin max cut starting at the bin min 
        double binmin = (i==0 ? (double)varmin : binlims.at(binlims.size()-1)); //NOTE: NEED TO SET OUTSIDE OR ADD.....
        double binmax = (i==nbins-1 ? (double)varmax: binmin);
        std::string bin_cut = Form("%s>=%.16f && %s<%.16f",varname.c_str(),binmin,varname.c_str(),binmax);
        int bincount = (int)*frame.Filter(bin_cut).Count();
        bool pass_flag = false;

        // Set the initial adjustment step
        step = (double)(varmax-varmin)/nbins;//IMPORTANT! NEED TO RESET IN CASE IT'S NEGATIVE BUT ALSO SO IT DOESN'T START REALLY SMALL.
        double delta = TMath::Abs(bincount-targetbincount);

        // Adjust the bin max cut until the statistics match within the threshold value
        if (i<nbins-1) {
            while (delta>threshold) {

                // Set the flags for whether or not you surpassed the bin count
                if (bincount>targetbincount)  { pass_flag = true;  }
                if (bincount<=targetbincount) { pass_flag = false; }

                // Reset upper bin limit
                binmax += step;

                // Reset bin cut with new bin lims
                bin_cut = Form("%s>=%.16f && %s<%.16f",varname.c_str(),binmin,varname.c_str(),binmax);

                // Reset bin count with new bin lims
                bincount = (int)*frame.Filter(bin_cut.c_str()).Count();

                // Reset step size if passed targetbincount and switch step direction
                if (pass_flag && bincount<=targetbincount) step *= -0.3; //NOTE: IF YOU DID PASS THE BIN COUNT ALREADY AND ARE NOW LESS 
                if (!pass_flag && bincount>targetbincount) step *= -0.3; //NOTE: IF YOU DIDN'T 

                // Reset delta of bin count to target
                delta = TMath::Abs(bincount-targetbincount);

            } // while(delta<threshold)
        } // if (i<nbins-1)

        // Compute the bin statistics
        double binmean = (double)*frame.Filter(bin_cut).Mean(varname.c_str());
        double binstd  = (double)*frame.Filter(bin_cut).StdDev(varname.c_str());

        // Add to bin lims vector
        binlims.push_back(binmax);

        // Show results message
        std::cout<<"------------------------------------------------------------"<<std::endl;
        std::cout<<" i        = "<<i<<std::endl;
        std::cout<<" varname  = "<<varname.c_str()<<std::endl;
        std::cout<<" bin_cut   = "<<bin_cut.c_str()<<std::endl;
        std::cout<<" varmax   = "<<varmax<<std::endl;
        std::cout<<" varmin   = "<<varmin<<std::endl;
        std::cout<<" mean     = "<<binmean<<std::endl;
        std::cout<<" stddev   = "<<binstd<<std::endl;
        std::cout<<" bincount = "<<bincount<<std::endl;
        std::cout<<" step                = "<<step<<std::endl;
        std::cout<<" bincount            = "<<bincount<<std::endl;
        std::cout<<" |bincount - target| = "<<delta<<std::endl;
        std::cout<<" \% diff             = "<<100*(delta/targetbincount)<<"\%"<<std::endl;
        std::cout<<" bin_cut             = "<<bin_cut<<std::endl;

    } // for (int i=0; i<nbins; i++) {

    // Print out bin limits
    std::cout<<varname<<" = [";
    for (int i=0; i<nbins; i++) {
        double binmin = binlims.at(i);
        std::string limstring = Form(" %.4f,",binmin);
        std::cout<<limstring.c_str();
    }
    std::cout<<" ]"<<std::endl;

    return binlims;

} // std::vector<double> findBinLims()

/**
* @brief Recursively set a map of bin scheme coordinates to bin variable limits for a nested bin scheme.
*
* Given a dataframe and a yaml node defining a nested binning scheme with the desired number
* of bins specified at each level, recursively set a map of bin scheme coordinates to bin variable limits
* and a list of bin variables encountered.
*
* @param frame ROOT RDataFrame with which to find bin limits
* @param node YAML node containing nested bin scheme definition
* @param node_name Name of YAML node
* @param nbins_key YAML key for number of bins at current depth
* @param lims_key YAML key for bin limits at current depth
* @param nested_key YAML key for nested binning
* @param bin_cuts List of bin cuts to apply to the dataframe
*/
void findNestedBinLims(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        YAML::Node node,
        std::string node_name             = "",
        std::string nbins_key             = "nbins",
        std::string lims_key              = "lims",
        std::string nested_key            = "nested",
        std::vector<std::string> bin_cuts = {}
    ) {
    
    // Check the YAML node
    if (node && node.IsMap()) {

        // Check if there is more depth and return if not
        if (node[nested_key] && node[nested_key].IsSequence()) {

            // Get nested YAML node
            auto node_nested = node[nested_key];

            // Loop bins
            for (int bin=0; bin<node_nested.size(); bin++) {

                // Get nested bin YAML node
                auto node_nested_bin = node_nested[bin];

                // Loop each nested bin looking for bin variables and select ONLY the first one found
                for (auto it_nested = node_nested_bin.begin(); it_nested != node_nested_bin.end(); ++it_nested) {

                    // Get bin variable name
                    std::string it_key = it_nested->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES OR THE NBINS_KEY

                    // Filter dataframe if a bin cut is available
                    auto df_filtered = (bin_cuts.size()==0) ? frame : frame.Filter(bin_cuts[bin].c_str());

                    // Check for bin limits and number of bins and find limits if not provided
                    std::vector<double> bin_lims;
                    if (it_nested->second[lims_key]) {
                        bin_lims = it_nested->second[lims_key].as<std::vector<double>>();
                    } else if (it_nested->second[nbins_key]) {
                        const int nbins = it_nested->second[nbins_key].as<int>();
                        bin_lims = findBinLims(df_filtered, it_key, nbins);
                        it_nested->second[lims_key] = bin_lims;
                    }

                    // Set bin cuts to carry to next depth of bin scheme
                    std::vector<std::string> new_bin_cuts;
                    for (int idx=0; idx<bin_lims.size()-1; idx++) {
                        std::string bincut = saga::util::addLimitCuts("",{it_key},{{bin_lims[idx], bin_lims[idx+1]}});
                        new_bin_cuts.push_back(bincut);
                    }

                    // Recursion call
                    findNestedBinLims(
                        df_filtered,
                        it_nested->second,
                        it_key,
                        nbins_key,
                        lims_key,
                        nested_key,
                        new_bin_cuts
                    );
                        
                    // Break on first nested bin variable found
                    break;

                } // for (auto it_nested = node_nested_bin.begin(); it_nested != node_nested_bin.end(); ++it_nested) {

            } // for (int bin=0; bin<node_nested.size(); bin++) {

        } else { return; } // if (node[nested_key] && node[nested_key].IsSequence()) {

    } // if (node && node.IsMap()) {
    else { return; }

} // void findNestedBinLims() {

/**
* @brief Compute bin limits on a regular interval between a minimum and maximum.
*
* @param nbins Number of bins
* @param xmin Minimum bound of bin variable
* @param xmax Maximum bound of bin variable
*
* @return Bin limits
*/
std::vector<double> getBinLims(
        const int nbins,
        double xmin,
        double xmax
    ) {

    std::vector<double> binlims;
    double step = (xmax-xmin)/nbins;
    for (int i=0; i<nbins+1; i++) {
        double lim = xmin + i*step; 
        binlims.push_back(lim);
    }

    return binlims;

} // std::vector<double> getBinLims(

/**
* @brief Set binning scheme cuts for a nested binning scheme.
*
* Recursively set a list of bin cuts given a YAML node defining a nested bin scheme.
* Note that this will set cuts for all bins within the nested bin scheme.
*
* @param cuts List of nested binning cuts
* @param node YAML node defining a nested bin scheme
* @param node_name Name of nested YAML node
* @param old_cuts Old list of cuts from previous recursion level
* @param lims_key YAML key for bin limits
* @param nested_key YAML key for nested binning
*
* @return List of bin cuts for nested binning scheme
*/
void setNestedBinCuts(
        std::vector<std::string> &cuts, //NOTE: List to set.
        YAML::Node               node,
        std::vector<std::string> &old_cuts, //NOTE: Modify this separately in each branch and then add cuts to the overall vector before returning
        std::string              node_name  = "",
        std::string              lims_key   = "lims",
        std::string              nested_key = "nested"
    ) {

    // Check the YAML node
    if (node && node.IsMap()) {

        // Check for bin limits
        std::vector<double> lims;
        if (node[lims_key] && node[lims_key].IsSequence()) {
            lims = node[lims_key].as<std::vector<double>>();
        }

        // Set nbins lower limit to 0 since you allow passing limits with length 0
        int nbins = lims.size()-1;
        if (nbins<0) nbins=0;

        // Loop bins and get bin cuts
        std::vector<std::string> varcuts;
        for (int bin=0; bin<nbins; bin++) {
            std::string cut = saga::util::addLimitCuts("",{node_name.c_str()},{{lims[bin],lims[bin+1]}});
            varcuts.push_back(cut);
        }

        // Loop previous cuts
        std::vector<std::string> newcuts;
        for (int idx=0; idx<old_cuts.size(); idx++) {
            
            // Loop this variable's cuts
            for (int bin=0; bin<varcuts.size(); bin++) {
                std::string newcut = Form("%s && %s",old_cuts[idx].c_str(),varcuts[bin].c_str());
                newcuts.push_back(newcut);
            }
        }

        // Reassign old_cuts
        if (old_cuts.size()==0) { old_cuts = varcuts;}
        else { old_cuts = newcuts; }

        // Check for nested binning
        if (node[nested_key] && node[nested_key].IsSequence()) {

            // Get nested YAML node
            auto node_nested = node[nested_key];

            // Loop nested bins
            for (int bin=0; bin<node_nested.size(); bin++) { //NOTE: These are not bin limits just maps to bin limits for each bin, so loop normally.

                // Get nested YAML node
                auto node_nested_bin = node_nested[bin];

                // Loop nested bin variables (only expect one!)
                for (auto it_nested = node_nested_bin.begin(); it_nested != node_nested_bin.end(); ++it_nested) {

                    // Get bin variable
                    std::string it_key = it_nested->first.as<std::string>();

                    // Create a new vector for uniqueness along different recursion branches
                    std::vector<std::string> new_old_cuts;
                    if (bin<old_cuts.size()) new_old_cuts = {old_cuts[bin]}; 

                    // Recursion call
                    setNestedBinCuts(
                        cuts,
                        it_nested->second,
                        new_old_cuts,
                        it_key,
                        lims_key,
                        nested_key
                    );

                    // Break on first nested variable found
                    break;
                }
            }
        } else {
            for (int bin=0; bin<old_cuts.size(); bin++) {
                cuts.push_back(old_cuts[bin]);
            }
            return;
        } // if (node[nested_key] && node[nested_key].IsSequence()) {
    } else {
        return;
    } // if (node && node.IsMap()) {
}

/**
* @brief Produce binning scheme cuts for a grid binning scheme.
*
* Produce a map of unique integer bin identifiers to bin cuts given a map of bin variables
* to their respective bin limits.  Note that this will produce cuts for all bins
* within the grid scheme and bin identifiers by default start at zero but can be made to
* start at any integer.
*
* @param binscheme Map of bin variable names to their respective bin limits
* @param start_bin_id Starting unique integer bin identifier
*
* @return Map of unique integer bin ids to bin cuts
*/
std::map<int,std::string> getBinCuts(
        std::map<std::string,std::vector<double>> binscheme,
        int                                       start_bin_id
    ) {

    std::vector<std::string> cuts;

    // Loop bin variables
    for (auto it = binscheme.begin(); it != binscheme.end(); ++it) {

        // Get bin variable name and limits
        std::string binvar = it->first;
        std::vector<double> lims = it->second;

        // Loop bin limits and get bin cuts
        std::vector<std::string> varcuts;
        for (int bin=0; bin<lims.size()-1; bin++) {
            double bin_min = lims[bin];
            double bin_max = lims[bin+1];
            std::string cut = Form("(%s>=%.8f && %s<%.8f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
            varcuts.push_back(cut);
        }

        // Loop previous cuts
        std::vector<std::string> newcuts;
        for (int idx=0; idx<cuts.size(); idx++) {
            
            // Loop this variable's cuts
            for (int bin=0; bin<varcuts.size(); bin++) {
                std::string newcut = Form("%s && %s",cuts[idx].c_str(),varcuts[bin].c_str());
                newcuts.push_back(newcut);
            }
        }

        // Reassign cuts
        if (cuts.size()==0) { cuts = varcuts; }
        else { cuts = newcuts; }
    }

    // Convert vector to map starting at given start index
    std::map<int,std::string> bincuts;
    for (int idx=0; idx<cuts.size(); idx++) {
        bincuts[start_bin_id+idx] = cuts[idx];
    }

    return bincuts;
}

/**
* @brief Read a YAML node and create a map of bin scheme names to maps of bin id to cuts.
*
* Produce a map of unique bin scheme names to maps of unique integer bin identifiers to bin cuts
* given a YAML node containing a map of bin variables to their respective bin limits.
* Note that this will produce cuts for all bins within the grid scheme and bin identifiers
* by default start at zero but can be made to start at any integer.
*
* @param node_binschemes YAML node containing bin scheme definitions
* @param start_bin_id Starting unique integer bin identifier
*
* @return Map of bin scheme names to maps of unique integer bin ids to bin cuts
*/
std::map<std::string,std::map<int,std::string>> getBinCutsMap(YAML::Node node_binschemes, int start_bin_id = 0) {

    // Set minimum allowed bin id
    int min_bin_id = start_bin_id;

    // Initialize bin cuts map
    std::map<std::string,std::map<int,std::string>> bincuts_map;

    // Loop bin schemes
    for (auto it_binschemes = node_binschemes.begin(); it_binschemes != node_binschemes.end(); ++it_binschemes) {

        // Get bin scheme name
        std::string binscheme_name = it_binschemes->first.as<std::string>();

        // Get bin scheme node
        auto node_binscheme = node_binschemes[binscheme_name];

        // Read bin scheme
        if (node_binscheme && node_binscheme.IsMap()) {

            // Recursively read nested bin scheme OR loop bin scheme yaml and create grid
            std::map<std::string,std::vector<double>> binscheme;
            std::map<int,std::string> bincuts;

            // Check if you have a nested bin scheme
            if (node_binscheme["nested"]) {

                try {
                    // Set nested bin cuts
                    std::vector<std::string> cuts;
                    std::vector<std::string> old_cuts;
                    setNestedBinCuts(cuts,node_binscheme,old_cuts,"");

                    // Convert vector to map starting at given start index
                    for (int idx=0; idx<cuts.size(); idx++) {
                        bincuts[start_bin_id+idx] = cuts[idx];
                    }
                } catch (std::exception& e) {
                    std::cerr<<"ERROR: Could not read nested bin limits for binscheme: "<<binscheme_name.c_str()<<std::endl;
                }
            }
            
            // Otherwise loop the yaml for a grid scheme
            else {
                for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

                    // Get bin variable name
                    std::string binvar = it->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES

                    // Compute bin limits or read directly from yaml 
                    int nbins = 0;
                    std::vector<double> binlims;
                    auto node_binvar = node_binscheme[binvar];
                    if (node_binvar.IsMap() && node_binvar["nbins"] && node_binvar["lims"]) {
                        try {
                            int nbins = node_binvar["nbins"].as<int>();
                            std::vector<double> lims = node_binvar["lims"].as<std::vector<double>>();
                            if (nbins>0 && lims.size()==2) {
                                binlims = getBinLims(nbins,lims[0],lims[1]);
                            }
                        } catch (std::exception& e) {
                            std::cerr<<"ERROR: Could not compute bin limits for binvar: "<<binvar.c_str()<<std::endl;
                        }
                    } else {
                        try {
                            binlims = node_binvar.as<std::vector<double>>();
                        } catch (std::exception& e) {
                            std::cerr<<"ERROR: Could not read bin limits for binvar: "<<binvar.c_str()<<std::endl;
                        }
                    }

                    // Add bin limits list to bin scheme
                    binscheme[binvar] = binlims;

                } // for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {
            }

            // Get bin cuts map and reset bin id minimum
            if (binscheme.size()>0) bincuts = getBinCuts(binscheme,min_bin_id);
            bincuts_map[binscheme_name] = bincuts;

        } // if (node_binscheme && node_binscheme.IsMap()) {

    } // for (auto it = node["binschemes"].begin(); it != node["binschemes"].end(); ++it) {

    return bincuts_map;
} // std::map<std::string,std::map<int,std::string>> getBinCutsMap(YAML::Node node_binschemes, int start_bin_id = 0) {

/**
* @brief Produce a list of lists of the bin variables used in each bin scheme defined in given a YAML node.
*
* @param node_binschemes YAML node containing bin scheme definitions
*
* @return Map of bin scheme names to lists of bin variable names used in each scheme
*/
std::map<std::string,std::vector<std::string>> getBinSchemesVars(YAML::Node node_binschemes) {

    // Initialize bin cuts map
    std::map<std::string,std::vector<std::string>> binschemes_vars;

    // Loop bin schemes
    for (auto it_binschemes = node_binschemes.begin(); it_binschemes != node_binschemes.end(); ++it_binschemes) {

        // Get bin scheme name
        std::string binscheme_name = it_binschemes->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN SCHEMES

        // Get bin scheme node
        auto node_binscheme = node_binschemes[binscheme_name];

        // Read bin scheme
        if (node_binscheme["nested"] && node_binscheme["nested"].IsSequence()) {

            // Follow nested yaml structure and create nested bin scheme
            YAML::Node node_nested = node_binscheme["nested"];
            std::vector<std::string> binvars;
            while (node_nested && node_nested.IsSequence()) {

                // Loop keys to find the bin variable (only expect one entry!)
                for (auto it = node_nested[0].begin(); it != node_nested[0].end(); ++it) {

                    // Get bin variable name and add to list
                    std::string binvar = it->first.as<std::string>();
                    binvars.push_back(binvar);
                    break;
                }

                // Reset the node
                node_nested = node_nested[0]["nested"];

            } // while (node_nested && node_nested.IsSequence()) {

            binschemes_vars[binscheme_name] = binvars;

        } else if (node_binscheme && node_binscheme.IsMap()) {

            // Loop bin scheme yaml and create grid bin scheme
            std::vector<std::string> binvars;
            for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

                // Get bin variable name and add to list
                std::string binvar = it->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES
                binvars.push_back(binvar);

            } // for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

            binschemes_vars[binscheme_name] = binvars;

        } // if (node_binscheme && node_binscheme.IsMap()) {

    } // for (auto it = node["binschemes"].begin(); it != node["binschemes"].end(); ++it) {

    return binschemes_vars;
} // std::vector<std::vector<std::string>> getBinSchemesVars(YAML::Node node_binschemes) {

/**
* @brief Compute bin migration fractions and save to a CSV file.
*
* Compute bin migration fraction and save to a CSV file.  Note that
* the truth bin cuts will be inferred from the provided cuts
* assuming they follow the form `(binvar>=binmin && binvar<=binmax)`.
*
* @param frame ROOT RDataframe from which to compute bin migration fraction
* @param scheme_name Bin scheme name
* @param bincuts Map of unique integer bin identifiers to bin cuts
* @param binvars List of bin variable names
* @param mc_suffix Suffix for forming the truth variable names
*/
void getBinMigration(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    std::string                                                   scheme_name,
    std::map<int,std::string>                                     bincuts,
    std::vector<std::string>                                      binvars,
    std::string                                                   mc_suffix = "_mc"
    ) {

    // Form MC Truth bin cut map
    std::map<int,std::string> bincuts_gen;
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {
        int binid = it->first;
        std::string bincut = it->second;
        for (int idx=0; idx<binvars.size(); idx++) {
            std::string binvar = binvars[idx];
            saga::util::replaceAll(bincut,Form("(%s>",binvar.c_str()),Form("(%s%s>",binvar.c_str(),mc_suffix.c_str())); //NOTE: ASSUME FORM OF BIN CUTS FROM getBinCuts() METHOD.
            saga::util::replaceAll(bincut,Form("& %s<",binvar.c_str()),Form("& %s%s<",binvar.c_str(),mc_suffix.c_str())); //NOTE: ASSUME FORM OF BIN CUTS FROM getBinCuts() METHOD.
        }
        bincuts_gen[binid] = bincut;
    }

    // Open output CSV
    std::string csvpath = Form("%s_bin_migration.csv",scheme_name.c_str());
    std::ofstream csvoutf; csvoutf.open(csvpath.c_str());
    std::ostream &csvout = csvoutf;
    std::string csv_separator = ",";

    // Set CSV column headers
    // COLS: binid_gen,binid_rec,migration_fraction=N_(REC && GEN)/N_GEN
    csvout << "binid_gen" << csv_separator.c_str();
    csvout << "binid_rec" << csv_separator.c_str();
    csvout << "mig" << std::endl;

    // Loop bin cut maps, compute bin migration, and write to CSV
    for (auto it_gen = bincuts_gen.begin(); it_gen != bincuts_gen.end(); ++it_gen) {

        // Get generated bin id and cut
        int binid_gen = it_gen->first;
        std::string bincut_gen = it_gen->second;

        // Get generated count
        auto frame_filtered = frame.Filter(bincut_gen.c_str());
        double count_gen = (double)*frame_filtered.Count();

        // Loop reconstructed bins
        for (auto it_rec = bincuts.begin(); it_rec != bincuts.end(); ++it_rec) {

            // Get generated bin id and cut
            int binid_rec = it_rec->first;
            std::string bincut_rec = it_rec->second;

            // Get reconstructed count
            double count_rec = (double)*frame_filtered.Filter(bincut_rec.c_str()).Count();

            // Compute bin migration fraction
            double mig = count_rec / count_gen; //NOTE: f[i->j] = [# generated in i AND reconstructed in j] / [# generated in bin i]

            // Set CSV column data
            // COLS: bin_id_gen,bin_id_rec,migration_fraction=N_(REC && GEN)/N_GEN
            csvout << binid_gen << csv_separator.c_str();
            csvout << binid_rec << csv_separator.c_str();
            csvout << mig << std::endl;

        } // for (auto it_rec = bincuts.begin(); it_rec != bincuts.end(); ++it_rec) {

    } // for (auto it_gen = bincuts.begin(); it_gen != bincuts.end(); ++it_gen) {

    // Close CSV file
    csvoutf.close();

} // void getBinMigration()

/**
* @brief Compute bin statistics and kinematics and save to a CSV file.
*
* @param frame ROOT RDataframe from which to compute bin migration fraction
* @param scheme_name Bin scheme name, csv file will be named `<scheme_name>_kinematics.csv`
* @param bincuts Map of unique integer bin identifiers to bin cuts
* @param kinvars List of kinematic variable names
*/
void getBinKinematics(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    std::string                                                   scheme_name,
    std::map<int,std::string>                                     bincuts,
    std::vector<std::string>                                      kinvars
    ) {

    // Open output CSV
    std::string csvpath = Form("%s_kinematics.csv",scheme_name.c_str());
    std::ofstream csvoutf; csvoutf.open(csvpath.c_str());
    std::ostream &csvout = csvoutf;
    std::string csv_separator = ",";

    // Set CSV column headers
    // COLS: bin,{kinvar,kinvar_err}
    csvout << "bin" << csv_separator.c_str();
    csvout << "count"; if (kinvars.size()>0) { csvout << csv_separator.c_str(); }
    for (int idx=0; idx<kinvars.size(); idx++) {
        csvout << kinvars[idx].c_str() << csv_separator.c_str();
        csvout << kinvars[idx].c_str() << "_err";
        if (idx<kinvars.size()-1) csvout << csv_separator.c_str();
    }
    csvout << std::endl;

    // Loop bin cut maps, compute count and kinematics, and write to CSV
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get generated bin id and cut
        int bin = it->first;
        std::string bincut = it->second;

        // Apply bin cut
        auto frame_filtered = frame.Filter(bincut.c_str());

        // Get bin count
        int count = (int)*frame_filtered.Count();

        // Set CSV column data
        // COLS: bin,{kinvar,kinvar_err}
        csvout << bin << csv_separator.c_str();
        csvout << count; if (kinvars.size()>0) { csvout << csv_separator.c_str(); }
        for (int idx=0; idx<kinvars.size(); idx++) {
            double binvar_mean = (double)*frame_filtered.Mean(kinvars[idx].c_str());
            double binvar_err  = (double)*frame_filtered.StdDev(kinvars[idx].c_str());
            csvout << binvar_mean << csv_separator.c_str();
            csvout << binvar_err;
            if (idx<kinvars.size()-1) csvout << csv_separator.c_str();
        }
        csvout << std::endl;

    } // for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

    // Close CSV file
    csvoutf.close();

} // void getBinKinematics()

/**
* @brief Create 1D kinematics histograms for each bin and save to a ROOT file.
*
* @param frame ROOT RDataframe from which to compute bin migration fraction
* @param scheme_name Bin scheme name, ROOT file will be named `<scheme_name>_kinematics.root`
* @param bincuts Map of unique integer bin identifiers to bin cuts
* @param kinvars List of kinematic variable names
* @param kinvar_lims List of outer bin limits for each kinematic variable
* @param kinvar_bins List of number of bins in each kinematic variable
*/
void getBinKinematicsTH1Ds(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        std::string                                                   scheme_name,
        std::map<int,std::string>                                     bincuts,
        std::vector<std::string>                                      kinvars,
        std::vector<std::vector<double>>                              kinvar_lims,
        std::vector<int>                                              kinvar_bins
    ) {

    // Open output ROOT file
    std::string path = Form("%s_kinematics.root",scheme_name.c_str());
    TFile *f = new TFile(path.c_str(),"RECREATE");

    // Loop bin cut maps, get kinematics histograms, and write to ROOT
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get generated bin id and cut
        int bin = it->first;
        std::string bincut = it->second;

        // Apply bin cut
        auto frame_filtered = frame.Filter(bincut.c_str());

        // Create 1D histograms for each kinematic variable
        for (int idx=0; idx<kinvars.size(); idx++) {
            std::string hist_name = Form("h1_bin%d_%s", bin, kinvars[idx].c_str());
            TH1D h1 = (TH1D) *frame_filtered.Histo1D({hist_name.c_str(),kinvars[idx].c_str(),kinvar_bins[idx],kinvar_lims[idx][0],kinvar_lims[idx][1]},kinvars[idx].c_str());
            f->WriteObject(&h1, hist_name.c_str());
        }

    } // for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

    // Close ROOT file
    f->Close();

} // void getBinKinematicsTH1Ds()

} // namespace bins {

} // namespace saga {
