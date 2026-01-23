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
// #include <TLegend.h>
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
// #include <TLatex.h>

// Local Includes
#include <log.h>
#include <util.h>
#include <data.h>

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

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::is_same;
using std::map;
using std::ofstream;
using std::ostream;
using std::runtime_error;
using std::string;
using std::unique_ptr;
using std::vector;

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
vector<double> findBinLims(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        string varname,
        const int nbins
    ) {

    string method_name = "findBinLims";

    // Initialize vectors
    vector<double> binmeans;
    vector<double> binmeans_err;
    vector<double> binlims;

    // Get overall limits and count and the necessary bin count with nbins
    LOG_DEBUG(Form("[%s]: Computing min, max, and count for variable %s", method_name.c_str(), varname.c_str()));
    double varmin = (double)*frame.Min(varname);
    double varmax = (double)*frame.Max(varname);
    int    count  = (int)*frame.Count();
    double targetbincount = (double)count/nbins;

    // Set initial step size and add initial bin limit
    double step = (double)(varmax-varmin)/nbins;
    LOG_DEBUG(Form("[%s]: Setting initial step size to %.8f", method_name.c_str(), step));
    binlims.push_back(varmin);
    LOG_DEBUG(Form("[%s]: Setting lowest bin limit to %.8f", method_name.c_str(), varmin));

    // Loop each bin and find its upper limit
    int nsteps;
    double threshold = 0.01*targetbincount; //NOTE: This is a threshold in counts, which will be an int.
    for (int i=0; i<nbins; i++) {

        LOG_DEBUG(Form("[%s]: Finding upper limit for bin %d", method_name.c_str(), i));

        // Set the bin max cut starting at the bin min 
        double binmin = (i==0 ? (double)varmin : binlims.at(binlims.size()-1)); //NOTE: NEED TO SET OUTSIDE OR ADD.....
        double binmax = (i==nbins-1 ? (double)varmax: binmin);
        string bin_cut = Form("%s>=%.16f && %s<%.16f",varname.c_str(),binmin,varname.c_str(),binmax);
        LOG_DEBUG(Form("[%s]: Filtering with bin cut %s", method_name.c_str(), bin_cut.c_str()));
        int bincount = (int)*frame.Filter(bin_cut).Count();
        bool pass_flag = false;

        // Set the initial adjustment step
        step = (double)(varmax-varmin)/nbins;//IMPORTANT! NEED TO RESET IN CASE IT'S NEGATIVE BUT ALSO SO IT DOESN'T START REALLY SMALL.
        double delta = TMath::Abs(bincount-targetbincount);
        LOG_DEBUG(Form("[%s]: Initializing with step=%.8f delta=%.8f", method_name.c_str(), step, delta));

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
                LOG_DEBUG(Form("[%s]: Filtering with updated bin cut %s", method_name.c_str(), bin_cut.c_str()));
                bincount = (int)*frame.Filter(bin_cut.c_str()).Count();

                // Reset step size if passed targetbincount and switch step direction
                if (pass_flag && bincount<=targetbincount) step *= -0.3; //NOTE: IF YOU DID PASS THE BIN COUNT ALREADY AND ARE NOW LESS 
                if (!pass_flag && bincount>targetbincount) step *= -0.3; //NOTE: IF YOU DIDN'T 

                // Reset delta of bin count to target
                delta = TMath::Abs(bincount-targetbincount);

                LOG_DEBUG(Form("[%s]: Updated parameters: step=%.8f delta=%.8f", method_name.c_str(), step, delta));

            } // while(delta<threshold)
        } // if (i<nbins-1)

        // Compute the bin statistics
        LOG_DEBUG(Form("[%s]: Initializing with step=%.8f delta=%.8f", method_name.c_str(), step, delta));
        double binmean = (double)*frame.Filter(bin_cut).Mean(varname.c_str());
        double binstd  = (double)*frame.Filter(bin_cut).StdDev(varname.c_str());

        // Add to bin lims vector
        binlims.push_back(binmax);

        // Show results message
        LOG_INFO("------------------------------------------------------------");
        LOG_INFO(Form(" i        = %d", i));
        LOG_INFO(Form(" varname  = %s", varname.c_str()));
        LOG_INFO(Form(" bin_cut  = %s", bin_cut.c_str()));
        LOG_INFO(Form(" varmax   = %.8f", varmax));
        LOG_INFO(Form(" varmin   = %.8f", varmin));
        LOG_INFO(Form(" mean     = %.8f", binmean));
        LOG_INFO(Form(" stddev   = %.8f", binstd));
        LOG_INFO(Form(" bincount = %d", bincount));
        LOG_INFO(Form(" step                = %.8f", step));
        LOG_INFO(Form(" bincount            = %d", bincount));
        LOG_INFO(Form(" |bincount - target| = %.8f", delta));
        LOG_INFO(Form(" %% diff             = %.3f%%", (double)100*(delta/targetbincount)));
        LOG_INFO(Form(" bin_cut             = %s", bin_cut.c_str()));

    } // for (int i=0; i<nbins; i++) {

    // Print out bin limits
    string msg = Form("%s = [", varname.c_str());
    for (int i=0; i<nbins; i++) {
        double binmin = binlims.at(i);
        msg = Form("%s %.4f,",msg.c_str(), binmin);
    }
    msg = Form("%s %.4f ]",msg.c_str(), varmax);
    LOG_INFO(msg);

    return binlims;

} // vector<double> findBinLims()

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
        string node_name             = "",
        string nbins_key             = "nbins",
        string lims_key              = "lims",
        string nested_key            = "nested",
        vector<string> bin_cuts = {}
    ) {

    string method_name = "findNestedBinLims";

    // Check the YAML node
    if (node && node.IsMap()) {

        // Check if there is more depth and return if not
        LOG_DEBUG(Form("[%s]: Checking yaml node for nested key %s", method_name.c_str(), nested_key.c_str()));
        if (node[nested_key] && node[nested_key].IsSequence()) {

            // Get nested YAML node
            auto node_nested = node[nested_key];

            // Loop bins
            for (int bin=0; bin<node_nested.size(); bin++) {

                // Get nested bin YAML node
                LOG_DEBUG(Form("[%s]: Grabbing element %d from yaml node[%s]...", method_name.c_str(), bin, nested_key.c_str()));
                auto node_nested_bin = node_nested[bin];

                // Loop each nested bin looking for bin variables and select ONLY the first one found
                for (auto it_nested = node_nested_bin.begin(); it_nested != node_nested_bin.end(); ++it_nested) {

                    // Get bin variable name
                    string it_key = it_nested->first.as<string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES OR THE NBINS_KEY
                    LOG_DEBUG(Form("[%s]: Found yaml key %s in node[%s][%d]", method_name.c_str(), it_key.c_str(), nested_key.c_str(), bin ));

                    // Filter dataframe if a bin cut is available
                    if (bin_cuts.size()!=0) LOG_DEBUG(Form("[%s]: Filtering frame with bin_cuts[%d] = %s", method_name.c_str(), bin, bin_cuts[bin].c_str()));
                    auto df_filtered = (bin_cuts.size()==0) ? frame : frame.Filter(bin_cuts[bin].c_str());

                    // Check for bin limits and number of bins and find limits if not provided
                    vector<double> bin_lims;
                    if (it_nested->second[lims_key]) {
                        LOG_DEBUG(Form("[%s]: Found bin limits key %s", method_name.c_str(), lims_key.c_str()));
                        bin_lims = it_nested->second[lims_key].as<vector<double>>();
                    } else if (it_nested->second[nbins_key]) {
                        LOG_DEBUG(Form("[%s]: Found nbins key %s", method_name.c_str(), nbins_key.c_str()));
                        const int nbins = it_nested->second[nbins_key].as<int>();
                        bin_lims = findBinLims(df_filtered, it_key, nbins);
                        it_nested->second[lims_key] = bin_lims;
                    }

                    // Set bin cuts to carry to next depth of bin scheme
                    LOG_DEBUG(Form("[%s]: Adding new bin limits cuts...", method_name.c_str()));
                    vector<string> new_bin_cuts;
                    for (int idx=0; idx<bin_lims.size()-1; idx++) {
                        string bincut = saga::util::addLimitCuts("",{it_key},{{bin_lims[idx], bin_lims[idx+1]}});
                        new_bin_cuts.push_back(bincut);
                    }

                    // Recursion call
                    LOG_DEBUG(Form("[%s]: Recursive call to findNestedBinLims...", method_name.c_str()));
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
vector<double> getBinLims(
        const int nbins,
        double xmin,
        double xmax
    ) {

    vector<double> binlims;
    double step = (xmax-xmin)/nbins;
    for (int i=0; i<nbins+1; i++) {
        double lim = xmin + i*step; 
        binlims.push_back(lim);
    }

    return binlims;

} // vector<double> getBinLims(

/**
* @brief Set binning scheme cuts for a nested binning scheme.
*
* Recursively set a list of bin cuts given a YAML node defining a nested bin scheme.
* Note that this will set cuts for all bins within the nested bin scheme.
*
* @param cuts List of nested binning cuts to set
* @param node YAML node defining a nested bin scheme
* @param node_name Name of nested YAML node
* @param old_cuts Old list of cuts from previous recursion level
* @param lims_key YAML key for bin limits
* @param nested_key YAML key for nested binning
*/
void setNestedBinCuts(
        vector<string> &cuts, //NOTE: List to set.
        YAML::Node               node,
        vector<string> &old_cuts, //NOTE: Modify this separately in each branch and then add cuts to the overall vector before returning
        string              node_name  = "",
        string              lims_key   = "lims",
        string              nested_key = "nested"
    ) {

    string method_name = "setNestedBinCuts";

    // Check the YAML node
    if (node && node.IsMap()) {

        // Check for bin limits
        vector<double> lims;
        if (node[lims_key] && node[lims_key].IsSequence()) {
            LOG_DEBUG(Form("[%s]: Found bin limits key %s", method_name.c_str(), lims_key.c_str()));
            lims = node[lims_key].as<vector<double>>();
        }

        // Set nbins lower limit to 0 since you allow passing limits with length 0
        int nbins = lims.size()-1;
        if (nbins<0) nbins=0;

        // Loop bins and get bin cuts
        LOG_DEBUG(Form("[%s]: Adding bin limits cuts...", method_name.c_str()));
        vector<string> varcuts;
        for (int bin=0; bin<nbins; bin++) {
            string cut = saga::util::addLimitCuts("",{node_name.c_str()},{{lims[bin],lims[bin+1]}});
            varcuts.push_back(cut);
        }

        // Loop previous cuts
        LOG_DEBUG(Form("[%s]: Looping previous bin cuts...", method_name.c_str()));
        vector<string> newcuts;
        for (int idx=0; idx<old_cuts.size(); idx++) {
            
            // Loop this variable's cuts
            for (int bin=0; bin<varcuts.size(); bin++) {
                string newcut = Form("%s && %s",old_cuts[idx].c_str(),varcuts[bin].c_str());
                newcuts.push_back(newcut);
            }
        }

        // Reassign old_cuts
        if (old_cuts.size()==0) { old_cuts = varcuts;}
        else { old_cuts = newcuts; }

        // Check for nested binning
        if (node[nested_key] && node[nested_key].IsSequence()) {

            // Get nested YAML node
            LOG_DEBUG(Form("[%s]: Found nested yaml node with key %s", method_name.c_str(), nested_key.c_str()));
            auto node_nested = node[nested_key];

            // Loop nested bins
            for (int bin=0; bin<node_nested.size(); bin++) { //NOTE: These are not bin limits just maps to bin limits for each bin, so loop normally.

                // Get nested YAML node
                LOG_DEBUG(Form("[%s]: Found bin node with index %d", method_name.c_str(), bin));
                auto node_nested_bin = node_nested[bin];

                // Loop nested bin variables (only expect one!)
                for (auto it_nested = node_nested_bin.begin(); it_nested != node_nested_bin.end(); ++it_nested) {

                    // Get bin variable
                    string it_key = it_nested->first.as<string>();
                    LOG_DEBUG(Form("[%s]: Found bin variable %s", method_name.c_str(), it_key.c_str()));

                    // Create a new vector for uniqueness along different recursion branches
                    vector<string> new_old_cuts;
                    if (bin<old_cuts.size()) new_old_cuts = {old_cuts[bin]}; 

                    // Recursion call
                    LOG_DEBUG(Form("[%s]: Recursive call to setNestedBinCuts...", method_name.c_str()));
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
            LOG_DEBUG(Form("[%s]: Looping old bin cuts...", method_name.c_str()));
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
map<int,string> getBinCuts(
        map<string,vector<double>> binscheme,
        int                                       start_bin_id
    ) {

    string method_name = "getBinCuts";

    vector<string> cuts;

    // Loop bin variables
    for (auto it = binscheme.begin(); it != binscheme.end(); ++it) {

        // Get bin variable name and limits
        string binvar = it->first;
        vector<double> lims = it->second;

        LOG_DEBUG(Form("[%s]: Getting bin cuts for variable %s", method_name.c_str(), binvar.c_str()));

        // Loop bin limits and get bin cuts
        LOG_DEBUG(Form("[%s]: Looping bin limits for cuts...", method_name.c_str()));
        vector<string> varcuts;
        for (int bin=0; bin<lims.size()-1; bin++) {
            double bin_min = lims[bin];
            double bin_max = lims[bin+1];
            string cut = Form("(%s>=%.8f && %s<%.8f)",binvar.c_str(),bin_min,binvar.c_str(),bin_max);
            varcuts.push_back(cut);
        }

        // Loop previous cuts
        LOG_DEBUG(Form("[%s]: Looping previous cuts...", method_name.c_str()));
        vector<string> newcuts;
        for (int idx=0; idx<cuts.size(); idx++) {
            
            // Loop this variable's cuts
            for (int bin=0; bin<varcuts.size(); bin++) {
                string newcut = Form("%s && %s",cuts[idx].c_str(),varcuts[bin].c_str());
                newcuts.push_back(newcut);
            }
        }

        // Reassign cuts
        if (cuts.size()==0) { cuts = varcuts; }
        else { cuts = newcuts; }
    }

    // Convert vector to map starting at given start index
    LOG_DEBUG(Form("[%s]: Converting bin cuts to map...", method_name.c_str()));
    map<int,string> bincuts;
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
*
* @throws Runtime error
*/
map<string,map<int,string>> getBinCutsMap(YAML::Node node_binschemes, int start_bin_id = 0) {

    string method_name = "getBinCutsMap";

    // Set minimum allowed bin id
    int min_bin_id = start_bin_id;

    // Initialize bin cuts map
    map<string,map<int,string>> bincuts_map;

    // Loop bin schemes
    for (auto it_binschemes = node_binschemes.begin(); it_binschemes != node_binschemes.end(); ++it_binschemes) {

        // Get bin scheme name
        string binscheme_name = it_binschemes->first.as<string>();
        LOG_DEBUG(Form("[%s]: Getting bin cuts for bin scheme %s", method_name.c_str(), binscheme_name.c_str()));

        // Get bin scheme node
        auto node_binscheme = node_binschemes[binscheme_name];

        // Read bin scheme
        LOG_DEBUG(Form("[%s]: Checking if yaml node[%s] is a map...", method_name.c_str(), binscheme_name.c_str()));
        if (node_binscheme && node_binscheme.IsMap()) {

            // Recursively read nested bin scheme OR loop bin scheme yaml and create grid
            map<string,vector<double>> binscheme;
            map<int,string> bincuts;

            // Check if you have a nested bin scheme
            if (node_binscheme["nested"]) {
                LOG_DEBUG(Form("[%s]: Found nested bin scheme...", method_name.c_str()));

                try {
                    // Set nested bin cuts
                    LOG_DEBUG(Form("[%s]: Setting nested bin cuts", method_name.c_str()));
                    vector<string> cuts;
                    vector<string> old_cuts;
                    setNestedBinCuts(cuts,node_binscheme,old_cuts,"");

                    // Convert vector to map starting at given start index
                    for (int idx=0; idx<cuts.size(); idx++) {
                        bincuts[start_bin_id+idx] = cuts[idx];
                    }
                } catch (exception& e) {
                    LOG_ERROR(Form("[%s]: Could not read nested bin limits for binscheme: %s", method_name.c_str(), binscheme_name.c_str()));
                    throw runtime_error(e.what());
                }
            }
            
            // Otherwise loop the yaml for a grid scheme
            else {
                LOG_DEBUG(Form("[%s]: Found grid bin scheme...", method_name.c_str()));
                for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

                    // Get bin variable name
                    string binvar = it->first.as<string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES
                    LOG_DEBUG(Form("[%s]: Found grid bin variable %s", method_name.c_str(), binvar.c_str()));

                    // Compute bin limits or read directly from yaml 
                    int nbins = 0;
                    vector<double> binlims;
                    auto node_binvar = node_binscheme[binvar];
                    if (node_binvar.IsMap() && node_binvar["nbins"] && node_binvar["lims"]) {
                        LOG_DEBUG(Form("[%s]: Found nbins and lims keys...", method_name.c_str()));
                        try {
                            int nbins = node_binvar["nbins"].as<int>();
                            vector<double> lims = node_binvar["lims"].as<vector<double>>();
                            if (nbins>0 && lims.size()==2) {
                                binlims = getBinLims(nbins,lims[0],lims[1]);
                            }
                        } catch (exception& e) {
                            LOG_ERROR(Form("[%s]: Could not compute bin limits for binvar: %s", method_name.c_str(), binvar.c_str()));
                            throw runtime_error(e.what());
                        }
                    } else {
                        LOG_DEBUG(Form("[%s]: Attempting to load bin limits...", method_name.c_str()));
                        try {
                            binlims = node_binvar.as<vector<double>>();
                        } catch (exception& e) {
                            LOG_ERROR(Form("[%s]: Could not read bin limits for binvar: %s", method_name.c_str(), binvar.c_str()));
                            throw runtime_error(e.what());
                        }
                    }

                    // Add bin limits list to bin scheme
                    binscheme[binvar] = binlims;

                } // for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {
            }

            // Get bin cuts map and reset bin id minimum
            if (binscheme.size()>0) {
                LOG_DEBUG(Form("[%s]: Getting bin cuts...", method_name.c_str()));
                bincuts = getBinCuts(binscheme,min_bin_id);
            }
            bincuts_map[binscheme_name] = bincuts;

        } // if (node_binscheme && node_binscheme.IsMap()) {

    } // for (auto it = node["binschemes"].begin(); it != node["binschemes"].end(); ++it) {

    return bincuts_map;
} // map<string,map<int,string>> getBinCutsMap(YAML::Node node_binschemes, int start_bin_id = 0) {

/**
* @brief Produce a list of lists of the bin variables used in each bin scheme defined in given a YAML node.
*
* @param node_binschemes YAML node containing bin scheme definitions
*
* @return Map of bin scheme names to lists of bin variable names used in each scheme
*/
map<string,vector<string>> getBinSchemesVars(YAML::Node node_binschemes) {

    string method_name = "getBinSchemesVars";

    // Initialize bin cuts map
    map<string,vector<string>> binschemes_vars;

    // Loop bin schemes
    for (auto it_binschemes = node_binschemes.begin(); it_binschemes != node_binschemes.end(); ++it_binschemes) {

        // Get bin scheme name
        string binscheme_name = it_binschemes->first.as<string>();//NOTE: THESE SHOULD BE NAMES OF BIN SCHEMES
        LOG_DEBUG(Form("[%s]: Found binscheme %s", method_name.c_str(), binscheme_name.c_str()));

        // Get bin scheme node
        auto node_binscheme = node_binschemes[binscheme_name];

        // Read bin scheme
        if (node_binscheme["nested"] && node_binscheme["nested"].IsSequence()) {

            LOG_DEBUG(Form("[%s]: Found nested bin scheme...", method_name.c_str()));

            // Follow nested yaml structure and create nested bin scheme
            YAML::Node node_nested = node_binscheme["nested"];
            vector<string> binvars;
            
            while (node_nested && node_nested.IsSequence()) {
                LOG_DEBUG(Form("[%s]: Found nested sequence...", method_name.c_str()));

                // Loop keys to find the bin variable (only expect one entry!)
                for (auto it = node_nested[0].begin(); it != node_nested[0].end(); ++it) {

                    // Get bin variable name and add to list
                    string binvar = it->first.as<string>();
                    binvars.push_back(binvar);
                    LOG_DEBUG(Form("[%s]: Found bin variable %s", method_name.c_str(), binvar.c_str()));

                    // Reset the node
                    node_nested = node_nested[0][binvar.c_str()]["nested"];
                    break;
                }

            } // while (node_nested && node_nested.IsSequence()) {

            binschemes_vars[binscheme_name] = binvars;

        } else if (node_binscheme && node_binscheme.IsMap()) {

            LOG_DEBUG(Form("[%s]: Found grid bin scheme...", method_name.c_str()));

            // Loop bin scheme yaml and create grid bin scheme
            vector<string> binvars;
            for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

                // Get bin variable name and add to list
                string binvar = it->first.as<string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES
                binvars.push_back(binvar);
                LOG_DEBUG(Form("[%s]: Found bin variable %s", method_name.c_str(), binvar.c_str()));

            } // for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

            binschemes_vars[binscheme_name] = binvars;

        } // if (node_binscheme && node_binscheme.IsMap()) {

    } // for (auto it = node["binschemes"].begin(); it != node["binschemes"].end(); ++it) {

    return binschemes_vars;
} // vector<vector<string>> getBinSchemesVars(YAML::Node node_binschemes) {


/**
* @brief Reduce a bin cuts map to a smaller batched version
*
* Reduce a bin cuts map to a smaller batched version given the
* total number of batches and the index of the batch.  This is useful for
* parallelizing results computed on a large bin cuts map.
*
* @param bincuts_map ROOT RDataframe from which to compute bin migration fraction
* @param nbatches Total number of batches
* @param ibatch Index of the batch \f$i\in[0,N_{batches}-1]\f$
*/
map<string,map<int,string>> getBinCutsMapBatch(
    map<string,map<int,string>> bincuts_map,
    int nbatches,
    int ibatch
    ) {

    string method_name = "getBinCutsMapBatch";

    // Initialize output map
    map<string,map<int,string>> bincuts_map_batch;

    // Loop bin schemes
    for (auto it = bincuts_map.begin(); it != bincuts_map.end(); ++it) {

        // Get bin scheme name and cuts
        string binscheme_name = it->first;
        LOG_DEBUG(Form("[%s]: Found binscheme %s", method_name.c_str(), binscheme_name.c_str()));
        map<int,string> bincuts = it->second;

        // Initialize output map
        map<int,string> bincuts_batch;

        // Loop bin cuts
        int n_bin_ids = bincuts.size();
        int idx = 0;
        int batch_size = n_bin_ids/nbatches;
        LOG_DEBUG(Form("[%s]: Looping bin cuts with batch size %d", method_name.c_str(), batch_size));
        for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

            // Get bin id and cut
            int binid = it->first;
            string bincut = it->second;

            LOG_DEBUG(Form("[%s]: Found index %d with bin cut %s", method_name.c_str(), binid, bincut.c_str()));

            // Check if bin id is in batch
            if ((                   idx>=batch_size*ibatch && idx<batch_size*(ibatch+1)) ||
                (idx==nbatches-1 && idx>=batch_size*ibatch && idx<batch_size*(ibatch+2))) { //NOTE: PUT REMAINDER IN LAST BATCH
                bincuts_batch[binid] = bincut;
            }
            idx++;
        }

        // Add to output map
        bincuts_map_batch[binscheme_name] = bincuts_batch;
    }

    return bincuts_map_batch;
}

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
* @param weight_name Name of the weight variable, ignored if empty
*/
void getBinMigration(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    string                                                   scheme_name,
    map<int,string>                                     bincuts,
    vector<string>                                      binvars,
    string                                                   mc_suffix = "_mc",
    string                                                   weight_name = ""
    ) {

    string method_name = "getBinMigration";

    // Form MC Truth bin cut map
    map<int,string> bincuts_gen;
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {
        int binid = it->first;
        string bincut = it->second;
        for (int idx=0; idx<binvars.size(); idx++) {
            string binvar = binvars[idx];
            saga::util::replaceAll(bincut,Form("(%s>",binvar.c_str()),Form("(%s%s>",binvar.c_str(),mc_suffix.c_str())); //NOTE: ASSUME FORM OF BIN CUTS FROM getBinCuts() METHOD.
            saga::util::replaceAll(bincut,Form("& %s<",binvar.c_str()),Form("& %s%s<",binvar.c_str(),mc_suffix.c_str())); //NOTE: ASSUME FORM OF BIN CUTS FROM getBinCuts() METHOD.
        }
        bincuts_gen[binid] = bincut;
    }

    // Open output CSV
    string csvpath = Form("%s_bin_migration.csv",scheme_name.c_str());
    LOG_DEBUG(Form("[%s]: Opening csv file %s", method_name.c_str(), csvpath.c_str()));
    ofstream csvoutf; csvoutf.open(csvpath.c_str());
    ostream &csvout = csvoutf;
    string csv_separator = ",";

    // Set CSV column headers
    // COLS: binid_gen,binid_rec,migration_fraction=N_(REC && GEN)/N_GEN
    LOG_DEBUG(Form("[%s]: Writing csv headers...", method_name.c_str()));
    csvout << "binid_gen" << csv_separator.c_str();
    csvout << "binid_rec" << csv_separator.c_str();
    csvout << "mig" << endl;

    // Loop bin cut maps, compute bin migration, and write to CSV
    for (auto it_gen = bincuts_gen.begin(); it_gen != bincuts_gen.end(); ++it_gen) {

        // Get generated bin id and cut
        int binid_gen = it_gen->first;
        string bincut_gen = it_gen->second;

        // Get generated count
        LOG_DEBUG(Form("[%s]: Filtering frame with generated bin cut %s", method_name.c_str(), bincut_gen.c_str()));
        auto frame_filtered = frame.Filter(bincut_gen.c_str());
        double count_gen = saga::data::get_weighted_count<double>(frame_filtered,weight_name);

        // Loop reconstructed bins
        for (auto it_rec = bincuts.begin(); it_rec != bincuts.end(); ++it_rec) {

            // Get generated bin id and cut
            int binid_rec = it_rec->first;
            string bincut_rec = it_rec->second;

            // Get reconstructed count
            LOG_DEBUG(Form("[%s]: Filtering frame with reconstructed bin cut %s", method_name.c_str(), bincut_rec.c_str()));
            auto frame_filtered_bincut_rec = frame_filtered.Filter(bincut_rec.c_str());
            double count_rec = saga::data::get_weighted_count<double>(frame_filtered_bincut_rec,weight_name);

            // Compute bin migration fraction
            double mig = count_rec / count_gen; //NOTE: f[i->j] = [# generated in i AND reconstructed in j] / [# generated in bin i]

            // Set CSV column data
            // COLS: bin_id_gen,bin_id_rec,migration_fraction=N_(REC && GEN)/N_GEN
            LOG_DEBUG(Form("[%s]: Writing csv data...", method_name.c_str()));
            csvout << binid_gen << csv_separator.c_str();
            csvout << binid_rec << csv_separator.c_str();
            csvout << mig << endl;

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
* @param weight_name Name of the weight variable, ignored if empty
*/
void getBinKinematics(
    ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
    string                                                   scheme_name,
    map<int,string>                                     bincuts,
    vector<string>                                      kinvars,
    string                                              weight_name
    ) {

    string method_name = "getBinKinematics";

    // Open output CSV
    string csvpath = Form("%s_kinematics.csv",scheme_name.c_str());
    LOG_DEBUG(Form("[%s]: Opening csv file %s", method_name.c_str(), csvpath.c_str()));
    ofstream csvoutf; csvoutf.open(csvpath.c_str());
    ostream &csvout = csvoutf;
    string csv_separator = ",";

    // Set CSV column headers
    // COLS: bin,{kinvar,kinvar_err}
    LOG_DEBUG(Form("[%s]: Writing csv headers...", method_name.c_str()));
    csvout << "bin" << csv_separator.c_str();
    csvout << "count"; if (kinvars.size()>0) { csvout << csv_separator.c_str(); }
    for (int idx=0; idx<kinvars.size(); idx++) {
        csvout << kinvars[idx].c_str() << csv_separator.c_str();
        csvout << kinvars[idx].c_str() << "_err";
        if (idx<kinvars.size()-1) csvout << csv_separator.c_str();
    }
    csvout << endl;

    // Loop bin cut maps, compute count and kinematics, and write to CSV
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get generated bin id and cut
        int bin = it->first;
        string bincut = it->second;

        // Apply bin cut
        LOG_DEBUG(Form("[%s]: Filtering with bin cut %s", method_name.c_str(), bincut.c_str()));
        auto frame_filtered = frame.Filter(bincut.c_str());

        // Get bin count
        int count = saga::data::get_weighted_count<int>(frame_filtered,weight_name);

        // Set CSV column data
        // COLS: bin,{kinvar,kinvar_err}
        LOG_DEBUG(Form("[%s]: Writing csv data...", method_name.c_str()));
        csvout << bin << csv_separator.c_str();
        csvout << count; if (kinvars.size()>0) { csvout << csv_separator.c_str(); }
        for (int idx=0; idx<kinvars.size(); idx++) {
            double binvar_mean = saga::data::get_weighted_mean<double>(frame_filtered,kinvars[idx],weight_name);
            double binvar_err  = saga::data::get_weighted_mean<double>(frame_filtered,kinvars[idx],weight_name,binvar_mean);
            csvout << binvar_mean << csv_separator.c_str();
            csvout << binvar_err;
            if (idx<kinvars.size()-1) csvout << csv_separator.c_str();
        }
        csvout << endl;

    } // for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

    // Close CSV file
    csvoutf.close();

} // void getBinKinematics()

/**
* @brief Save a TH1 or TH2 ROOT histogram to a CSV file.
*
* Save a TH1 or TH2 ROOT histogram to a CSV file.  The CSV file will be named
* have columns `bin`, `llimx`, `count` if it is 1D, or `binx`, `biny`, `llimx`,
* `llimy`, `count` if it is 2D.  Note that the lower bin limits are written for
* each bin so there are \f$N_{bins}+1\f$ rows in the CSV file for 1D histograms
* and \f$(N_{bins,x}+1)\times(N_{bins,y}+1)\f$ rows in the CSV file for 2D histograms.
*
* @param h1 ROOT histogram to save to CSV
* @param csv_name Path to CSV file
* 
* @throws Runtime error
*/
void saveTH1ToCSV(
    const TH1& h1,
    string csv_name
    ) {

    string method_name = "saveTH1ToCSV";

    // Check histogram dimensions
    LOG_DEBUG(Form("[%s]: Getting histogram bins...", method_name.c_str()));
    int nbinsx = h1.GetNbinsX();
    int nbinsy = h1.GetNbinsY();

    // Open output CSV
    LOG_DEBUG(Form("[%s]: Opening csv file %s", method_name.c_str(), csv_name.c_str()));
    ofstream csvoutf; csvoutf.open(csv_name.c_str());
    ostream &csvout = csvoutf;
    string csv_separator = ",";

    // Write CSV data
    if (nbinsx>0 && nbinsy<=1) { //TH1 case
        // Write column headers
        LOG_DEBUG(Form("[%s]: Writing csv headers for TH1...", method_name.c_str()));
        csvout << "bin" << csv_separator.c_str() << "llimx" << csv_separator.c_str() << "count" << endl;
        
        // Loop x bins
        LOG_DEBUG(Form("[%s]: Writing csv data for TH1...", method_name.c_str()));
        for (int idx=1; idx<=nbinsx+1; idx++) { //NOTE: ROOT HISTOGRAM INDICES BEGIN AT 1 AND YOU NEED TO WRITE ALL THE (N+1) BIN LIMITS
            
            // Get bin lower limit and count
            double llimx = h1.GetXaxis()->GetBinLowEdge(idx);
            int    count = (idx<nbinsx+1) ? h1.GetBinContent(idx) : 0;

            // Write data
            csvout << idx-1 << csv_separator.c_str(); //NOTE: WRITE PYTHON INDEX
            csvout << llimx << csv_separator.c_str();
            csvout << count << endl;
        }
    } else if (nbinsx>0 && nbinsy>1) { //TH2 case

        // Write column headers
        LOG_DEBUG(Form("[%s]: Writing csv headers for TH2...", method_name.c_str()));
        csvout << "binx" << csv_separator.c_str() << "biny" << csv_separator.c_str();
        csvout << "llimx" << csv_separator.c_str() << "llimy" << csv_separator.c_str() << "count" << endl;
        
        // Loop x bin
        LOG_DEBUG(Form("[%s]: Writing csv data for TH2...", method_name.c_str()));
        for (int idx=1; idx<=nbinsx+1; idx++) { //NOTE: ROOT HISTOGRAM INDICES BEGIN AT 1 AND YOU NEED TO WRITE ALL THE (N+1) BIN LIMITS
            
            // Get bin lower limit
            double llimx = h1.GetXaxis()->GetBinLowEdge(idx);

            // Loop y bins
            for (int idy=1; idy<=nbinsy+1; idy++) { //NOTE: ROOT HISTOGRAM INDICES BEGIN AT 1 AND YOU NEED TO WRITE ALL THE (N+1) BIN LIMITS
                
                // Get bin lower limit and count
                double llimy = h1.GetYaxis()->GetBinLowEdge(idy);
                int    count = (idx<nbinsx+1 && idy<nbinsy+1) ? h1.GetBinContent(idx,idy) : 0;

                // Write data
                csvout << idx-1 << csv_separator.c_str(); //NOTE: WRITE PYTHON INDEX
                csvout << idy-1 << csv_separator.c_str();
                csvout << llimx << csv_separator.c_str();
                csvout << llimy << csv_separator.c_str();
                csvout << count << endl;
            }
        }
    } else {
        string msg = Form("[%s]: Uknown histogram type! Must be TH1 or TH2.", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }
    
    // Close output CSV
    csvoutf.close();

} // void saveTH1ToCSV()

/**
* @brief Create 1D kinematics histograms for each bin and save to a ROOT file.
*
* @param frame ROOT RDataframe from which to compute bin migration fraction
* @param scheme_name Bin scheme name, ROOT file will be named `<scheme_name>_kinematics.root`
* @param bincuts Map of unique integer bin identifiers to bin cuts
* @param kinvars List of kinematic variable names
* @param kinvar_lims List of outer bin limits for each kinematic variable
* @param kinvar_bins List of number of bins in each kinematic variable
* @param save_pdfs Option to save 1D histograms as PDFs, files will be names `c1_<scheme_name>_bin<bin_id>_<kinvar>.pdf`
* @param save_csvs Option to save 1D histograms as CSVs, files will be names `<scheme_name>_bin<bin_id>_<kinvar>.csv`
*/
void getBinKinematicsTH1Ds(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        string                                                   scheme_name,
        map<int,string>                                     bincuts,
        vector<string>                                      kinvars,
        vector<vector<double>>                              kinvar_lims,
        vector<int>                                              kinvar_bins,
        bool                                                          save_pdfs = false,
        bool                                                          save_csvs = false
    ) {

    string method_name = "getBinKinematicsTH1Ds";

    // Open output ROOT file
    string path = Form("%s_kinematics.root",scheme_name.c_str());
    LOG_DEBUG(Form("[%s]: Opening TH1Ds ROOT file %s", method_name.c_str(), path.c_str()));
    TFile *f = new TFile(path.c_str(),"RECREATE");

    // Loop bin cut maps, get kinematics histograms, and write to ROOT
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get generated bin id and cut
        int bin = it->first;
        string bincut = it->second;

        // Apply bin cut
        LOG_DEBUG(Form("[%s]: Filtering frame with bin cut %s", method_name.c_str(), bincut.c_str()));
        auto frame_filtered = frame.Filter(bincut.c_str());

        // Create 1D histograms for each kinematic variable
        for (int idx=0; idx<kinvars.size(); idx++) {
            string hist_name = Form("h1_bin%d_%s", bin, kinvars[idx].c_str());
            LOG_DEBUG(Form("[%s]: Creating TH1D %s", method_name.c_str(), hist_name.c_str()));
            TH1D h1 = (TH1D) *frame_filtered.Histo1D({hist_name.c_str(),kinvars[idx].c_str(),kinvar_bins[idx],kinvar_lims[idx][0],kinvar_lims[idx][1]},kinvars[idx].c_str());
            f->WriteObject(&h1, hist_name.c_str());
            if (save_pdfs) {
                string canvas_name = Form("c1_%s_bin%d_%s", scheme_name.c_str(), bin, kinvars[idx].c_str());
                LOG_DEBUG(Form("[%s]: Creating TH1D pdf %s.pdf", method_name.c_str(), canvas_name.c_str()));
                TCanvas *c1 = new TCanvas(canvas_name.c_str());
                c1->cd();
                h1.Draw("COLZ");
                c1->Print(Form("%s.pdf", canvas_name.c_str()));
            }
            if (save_csvs) {
                string csv_name = Form("%s_bin%d_%s.csv", scheme_name.c_str(), bin, kinvars[idx].c_str());
                saveTH1ToCSV(h1, csv_name);
            }
        }

    } // for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

    // Close ROOT file
    f->Close();

} // void getBinKinematicsTH1Ds()

/**
* @brief Create 2D kinematics histograms for each bin and save to a ROOT file.
*
* @param frame ROOT RDataframe from which to compute bin migration fraction
* @param scheme_name Bin scheme name, ROOT file will be named `<scheme_name>_kinematics.root`
* @param bincuts Map of unique integer bin identifiers to bin cuts
* @param kinvars List of kinematic variable pairs (x-axis,y-axis) names
* @param kinvar_lims List of outer bin limits for each kinematic variable pair (x-axis,y-axis)
* @param kinvar_bins List of number of bins in each kinematic variable pair (x-axis,y-axis)
* @param save_pdfs Option to save 2D histograms as PDFs, files will be names `c2_<scheme_name>_bin<bin_id>_<kinvar_x>_<kinvar_y>.pdf`
* @param save_csvs Option to save 2D histograms as CSVs, files will be names `<scheme_name>_bin<bin_id>_<kinvar_x>_<kinvar_y>.csv`
*/
void getBinKinematicsTH2Ds(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        string                                                   scheme_name,
        map<int,string>                                     bincuts,
        vector<vector<string>>                         kinvars,
        vector<vector<vector<double>>>                 kinvar_lims,
        vector<vector<int>>                                 kinvar_bins,
        bool                                                          save_pdfs = false,
        bool                                                          save_csvs = false
    ) {

    string method_name = "getBinKinematicsTH2Ds";

    // Open output ROOT file
    string path = Form("%s_kinematics.root",scheme_name.c_str());
    LOG_DEBUG(Form("[%s]: Opening TH2Ds ROOT file %s", method_name.c_str(), path.c_str()));
    TFile *f = new TFile(path.c_str(),"RECREATE");

    // Loop bin cut maps, get kinematics histograms, and write to ROOT
    for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

        // Get generated bin id and cut
        int bin = it->first;
        string bincut = it->second;

        // Apply bin cut
        LOG_DEBUG(Form("[%s]: Filtering frame with bin cut %s", method_name.c_str(), bincut.c_str()));
        auto frame_filtered = frame.Filter(bincut.c_str());

        // Create 2D histograms for each kinematic variable
        for (int idx=0; idx<kinvars.size(); idx++) {
            string hist_name = Form("h2_bin%d_%s_%s", bin, kinvars[idx][0].c_str(), kinvars[idx][1].c_str());
            LOG_DEBUG(Form("[%s]: Creating TH2D %s", method_name.c_str(), hist_name.c_str()));
            string hist_title = Form("Bin %d : %s vs. %s", bin, kinvars[idx][0].c_str(), kinvars[idx][1].c_str());
            TH2D h2 = (TH2D) *frame_filtered.Histo2D({hist_name.c_str(),hist_title.c_str(),kinvar_bins[idx][0],kinvar_lims[idx][0][0],kinvar_lims[idx][0][1],kinvar_bins[idx][1],kinvar_lims[idx][1][0],kinvar_lims[idx][1][1]},kinvars[idx][0].c_str(),kinvars[idx][1].c_str());
            f->WriteObject(&h2, hist_name.c_str());
            if (save_pdfs) {
                string canvas_name = Form("c2_%s_bin%d_%s_%s", scheme_name.c_str(), bin, kinvars[idx][0].c_str(), kinvars[idx][1].c_str());
                LOG_DEBUG(Form("[%s]: Creating TH1D pdf %s.pdf", method_name.c_str(), canvas_name.c_str()));
                TCanvas *c1 = new TCanvas(canvas_name.c_str());
                c1->cd();
                h2.Draw("COLZ");
                c1->Print(Form("%s.pdf", canvas_name.c_str()));
            }
            if (save_csvs) {
                string csv_name = Form("%s_bin%d_%s_%s.csv", scheme_name.c_str(), bin, kinvars[idx][0].c_str(), kinvars[idx][1].c_str());
                saveTH1ToCSV(h2, csv_name);
            }
        }

    } // for (auto it = bincuts.begin(); it != bincuts.end(); ++it) {

    // Close ROOT file
    f->Close();

} // void getBinKinematicsTH2Ds()

} // namespace bins {

} // namespace saga {
