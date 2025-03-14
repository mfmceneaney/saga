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
* @param nested_bin_lims Map of bin scheme coordinates to bin variable limits
* @param binvars List of bin variables encountered (Note, the first bin variable will be set automatically from `node_nested_name`.)
* @param node_nested YAML node containing nested bin scheme definition
* @param node_nested_name Name of nested YAML node
* @param coordinates Current bin scheme coordinates
* @param nbins_key Key for number of bins at current depth
*/
void findNestedBinLims(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame,
        std::map<std::vector<int>,std::vector<double>> nested_bin_lims,
        std::vector<std::string> binvars,
        YAML::Node node_nested,
        std::string node_nested_name = "",
        std::vector<int> coordinates = {},
        std::string nbins_key = "nbins"
    ) {

    // Check if you were given the bin scheme node to start
    if (node_nested_name.size()==0) {  // NOTE: Assume you already check that this was a map.
        for (auto it_nested = node_nested.begin(); it_nested != node_nested.end(); ++it_nested) {

            // Get the entry key
            std::string it_key = it_nested->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES OR THE NBINS_KEY

            // Recursion call
            findNestedBinLims(
                frame,
                nested_bin_lims,
                binvars,
                it_nested->second,
                it_key,
                coordinates,
                nbins_key
            );
                
            // Break on first nested variable found
            break;
        }

    // Check the YAML node
    } else if (node_nested && node_nested.IsMap()) {

        // Add to bin variables
        if (binvars.size()==0) {
            binvars.push_back(node_nested_name);
        } else {
            if (binvars[binvars.size()-1]!=node_nested_name) binvars.push_back(node_nested_name);
        }

        // Check for number of bins
        std::vector<double> bin_lims;
        if (node_nested[nbins_key]) { //TODO: Also check if bins statically specified?
            const int nbins = node_nested[nbins_key].as<int>();
            bin_lims = findBinLims(frame, node_nested_name, nbins);
            nested_bin_lims[coordinates] = bin_lims;
        }

        // Check if there is more depth and return if not
        if (node_nested.size()<=1) return;

        // Loop the entries in your nested entry looking for bin variables and select ONLY the first one found
        for (auto it_nested = node_nested.begin(); it_nested != node_nested.end(); ++it_nested) {

            // Get the entry key
            std::string it_key = it_nested->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES OR THE NBINS_KEY

            // Check for requested bin number and produce bin limits if so
            if (it_key != nbins_key) {

                // Loop bins
                for (int bin=0; bin<bin_lims.size()-1; bin++) {

                    // Update coordinates
                    std::vector<int> new_coordinates = coordinates;
                    new_coordinates.push_back(bin);

                    // Create cut
                    std::string cut;
                    saga::util::addLimitCuts(cut,{it_key},{{bin_lims[bin],bin_lims[bin+1]}});

                    // Filter dataframe
                    auto df_filtered = frame.Filter(cut);

                    // Recursion call
                    findNestedBinLims(
                        df_filtered,
                        nested_bin_lims,
                        binvars,
                        it_nested->second,
                        it_key,
                        new_coordinates,
                        nbins_key
                    );

                }
                
                // Break on first nested variable found
                break;
            }
        }

        // Update yaml with bin limits
        node_nested["bins"] = bin_lims;

    } // if (node_nested && node_nested.IsMap()) {
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
        std::string binscheme_name = it_binschemes->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN SCHEMES

        // Get bin scheme node
        auto node_binscheme = node_binschemes[binscheme_name];

        // Read bin scheme
        if (node_binscheme && node_binscheme.IsMap()) {

            // Loop bin scheme yaml and create grid bin scheme
            std::map<std::string,std::vector<double>> binscheme;
            for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

                // Get bin variable name
                std::string binvar = it->first.as<std::string>();//NOTE: THESE SHOULD BE NAMES OF BIN VARIABLES

                // Compute bin limits or read directly from yaml 
                int nbins = 0;
                std::vector<double> binlims;
                auto node_binvar = node_binscheme[binvar];
                if (node_binvar.IsMap() && node_binvar["nbins"] && node_binvar["lims"]) {
                    int nbins = node_binvar["nbins"].as<int>();
                    std::vector<double> lims = node_binvar["lims"].as<std::vector<double>>();
                    if (nbins>0 && lims.size()==2) {
                        binlims = getBinLims(nbins,lims[0],lims[1]);
                    } else { std::cerr<<"ERROR: Could not read bins for binvar: "<<binvar.c_str()<<std::endl; }
                } else {
                    binlims = node_binvar.as<std::vector<double>>();
                }

                // Add bin limits list to bin scheme
                binscheme[binvar] = binlims;

            } // for (auto it = node_binscheme.begin(); it != node_binscheme.end(); ++it) {

            // Get bin cuts map and reset bin id minimum
            std::map<int,std::string> bincuts = getBinCuts(binscheme,min_bin_id);
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
        if (node_binscheme && node_binscheme.IsMap()) {

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
* the Monte Carlo truth bin cuts will be inferred from the provided cuts
* assuming they follow the form `(binvar>=binmin && binvar<=binmax)`.
*
* @param frame ROOT RDataframe from which to compute bin migration fraction
* @param scheme_name Bin scheme name
* @param bincuts Map of unique integer bin identifiers to bin cuts
* @param binvars List of bin variable names
* @param mc_suffix Suffix for forming Monte Carlo truth variable names
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
* @param scheme_name Bin scheme name
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
            csvout << binvar_err << "_err";
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
* @param scheme_name Bin scheme name
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
