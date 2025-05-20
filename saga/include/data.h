#include <iostream>
#include <memory>
#include <string>
#include <map>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RCsvDS.hxx>
#include <TRandom3.h>
#include <thread>
#include <functional>

// RooFit Includes
#include <RooCategory.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooAbsDataHelper.h>
#include <RooDataHist.h>
#include <RooArgList.h>
#include <RooWorkspace.h>

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 4/Feb./2025
* @version 0.0.0
* @brief Create asymmetry fit datasets.
*/

namespace saga {

namespace data {

using RNode = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;

/**
* @brief Create a dataset for an asymmetry fit.
*
* Create a RooFit dataset for an asymmetry fit from a ROOT RDataFrame,
* adding helicity, target spin, binning, depolarization, asymmetry fit, and invariant mass fit variables.
* Store all variables and RooDataSet in a RooWorkspace.
*
* @param frame ROOT RDataframe from which to create a RooDataSet
* @param w RooWorkspace in which to work
* @param name Dataset name
* @param title Dataset title
* @param helicity Name of helicity variable
* @param helicity_states Map of state names to helicity values
* @param tspin Name of target spin variable
* @param tspin_states Map of state names to target spin values
* @param htspin Name of helicity times target spin variable
* @param htspin_states Map of state names to helicity times target spin values
* @param combined_spin_state Name of combined spin state variable
* @param binvars List of kinematic binning variables names
* @param binvar_titles List of kinematic binning variables titles
* @param binvar_lims List kinematic binning variable minimum and maximum bounds 
* @param binvar_bins List of kinematic binning variables bins
* @param depolvars List of depolarization variables names
* @param depolvar_titles List of depolarization variables titles
* @param depolvar_lims List depolarization variable minimum and maximum bounds 
* @param depolvar_bins List of depolarization variables bins
* @param asymfitvars List of asymmetry fit variables names
* @param asymfitvar_titles List of asymmetry fit variables titles
* @param asymfitvar_lims List asymmetry fit variable minimum and maximum bounds 
* @param asymfitvar_bins List of asymmetry fit variables bins
* @param massfitvars List of invariant mass fit variables names
* @param massfitvar_titles List of invariant mass fit variables titles
* @param massfitvar_lims List invariant mass fit variable minimum and maximum bounds 
* @param massfitvar_bins List of invariant mass fit variables bins
*/
void createDataset(
        RNode frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        RooWorkspace *w,
        std::string name,
        std::string title,
        std::string helicity,
        std::map<std::string,int> helicity_states,
        std::string tspin,
        std::map<std::string,int> tspin_states,
        std::string htspin,
        std::map<std::string,int> htspin_states,
        std::string combined_spin_state,
        std::vector<std::string> binvars,
        std::vector<std::string> binvar_titles,
        std::vector<std::vector<double>> binvar_lims,
        std::vector<int> binvar_bins,
        std::vector<std::string> depolvars,
        std::vector<std::string> depolvar_titles,
        std::vector<std::vector<double>> depolvar_lims,
        std::vector<int> depolvar_bins,
        std::vector<std::string> asymfitvars,
        std::vector<std::string> asymfitvar_titles,
        std::vector<std::vector<double>> asymfitvar_lims,
        std::vector<int> asymfitvar_bins,
        std::vector<std::string> massfitvars,
        std::vector<std::string> massfitvar_titles,
        std::vector<std::vector<double>> massfitvar_lims,
        std::vector<int> massfitvar_bins
    ) {

    // Define the helicity variable
    RooCategory h(helicity.c_str(), helicity.c_str());
    for (auto it = helicity_states.begin(); it != helicity_states.end(); it++) {
        h.defineType(it->first.c_str(), it->second);
    }

    // Define the target spin variable
    RooCategory t(tspin.c_str(), tspin.c_str());
    for (auto it = tspin_states.begin(); it != tspin_states.end(); it++) {
        t.defineType(it->first.c_str(), it->second);
    }

    // Define the target spin variable
    RooCategory ht(htspin.c_str(), htspin.c_str());
    for (auto it = htspin_states.begin(); it != htspin_states.end(); it++) {
        ht.defineType(it->first.c_str(), it->second);
    }

    // Define the combined spin state variable
    RooCategory ss(combined_spin_state.c_str(), combined_spin_state.c_str());
    for (auto it_h = helicity_states.begin(); it_h != helicity_states.end(); it_h++) {
        for (auto it_t = htspin_states.begin(); it_t != htspin_states.end(); it_t++) {
            std::string state_name = Form("%s_%s",it_h->first.c_str(),it_t->first.c_str());
            int state_value = (it_h->second + 1) * 10 + (it_t->second + 1);
            ss.defineType(state_name.c_str(), state_value);
        }
    }

    // Define the full variables and limits lists
    std::vector<std::string> vars;
    for (int idx=0; idx<binvars.size();     idx++) vars.push_back(binvars[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) vars.push_back(depolvars[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) vars.push_back(asymfitvars[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) vars.push_back(massfitvars[idx]);
    std::vector<std::string> var_titles;
    for (int idx=0; idx<binvars.size();     idx++) var_titles.push_back(binvar_titles[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_titles.push_back(depolvar_titles[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_titles.push_back(asymfitvar_titles[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_titles.push_back(massfitvar_titles[idx]);
    std::vector<std::vector<double>> var_lims;
    for (int idx=0; idx<binvars.size();     idx++) var_lims.push_back(binvar_lims[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_lims.push_back(depolvar_lims[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_lims.push_back(asymfitvar_lims[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_lims.push_back(massfitvar_lims[idx]);
    std::vector<int> var_bins;
    for (int idx=0; idx<binvars.size();     idx++) var_bins.push_back(binvar_bins[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_bins.push_back(depolvar_bins[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_bins.push_back(asymfitvar_bins[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_bins.push_back(massfitvar_bins[idx]);
    int nvars = vars.size();

    // Define RooRealVar variables
    RooRealVar *rrvars[nvars];
    for (int rr=0; rr<nvars; rr++) {
        rrvars[rr] = new RooRealVar(vars[rr].c_str(), var_titles[rr].c_str(), var_lims[rr][0], var_lims[rr][1]);
        rrvars[rr]->setBins(var_bins[rr]);
    }

    // Define variable list for RooDataSetHelper
    RooArgSet *argset = new RooArgSet();
    for (int rr=0; rr<nvars; rr++) {
        argset->add(*rrvars[rr]);
    }

    // Create RDataFrame to RooDataSet pointer
    ROOT::RDF::RResultPtr<RooDataSet> rooDataSetResult;
    switch (nvars) {
        case 2: //NOTE: Need at least one fit variable and one bin variable.
            rooDataSetResult = frame.Book<float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 3:
            rooDataSetResult = frame.Book<float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 4:
            rooDataSetResult = frame.Book<float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 5:
            rooDataSetResult = frame.Book<float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 6:
            rooDataSetResult = frame.Book<float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 7:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 8:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 9:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 10:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 11:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 12:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 13:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 14:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 15:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 16:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 17:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        case 18:
            rooDataSetResult = frame.Book<float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float>(
                RooDataSetHelper(name.c_str(),title.c_str(),*argset),
                vars
            );
            break;
        default:
            std::cerr<<"ERROR: nvars="<<nvars<<" is outside the allowed range [2,18]"<<std::endl;
            return;
    }

    // Manually create dataset containing helicity as a RooCategory variable
    RooDataSet *ds_h = new RooDataSet("ds_h","ds_h", RooArgSet(h,t,ht,ss));

    // Set cuts for variable limits so that new dataset will have same length as old dataset
    std::string varlims_cuts = "";
    for (int vv=0; vv<nvars; vv++) {
        if (varlims_cuts.size()==0) {
            varlims_cuts = Form("(%s>=%.8f && %s<=%.8f)", vars[vv].c_str(), rrvars[vv]->getMin(), vars[vv].c_str(), rrvars[vv]->getMax());
        } else {
            varlims_cuts = Form("%s && (%s>=%.8f && %s<=%.8f)", varlims_cuts.c_str(), vars[vv].c_str(), rrvars[vv]->getMin(), vars[vv].c_str(), rrvars[vv]->getMax());
        }
    }

    // Loop RDataFrame and fill helicity dataset
    frame.Filter(varlims_cuts.c_str()).Foreach(
                  [&h,&t,&ht,&ss,&ds_h](float h_val, float tspin_val, int ss_val){

                    // Assign the state for the beam helicity (assuming 1 or -1 as the states)
                    if (h_val > 0) {
                      h.setIndex(1);
                    } else if (h_val < 0) {
                      h.setIndex(-1);
                    } else {
                      h.setIndex(0);
                    }

                    // Assign the state for the target spin (assuming 1 or -1 as the states)
                    if (tspin_val > 0) {
                      t.setIndex(1);
                    } else if (tspin_val < 0) {
                      t.setIndex(-1);
                    } else {
                      t.setIndex(0);
                    }

                    // Assign the state for the helicity times target spin (assuming 1 or -1 as the states)
                    float htspin_val = h_val * tspin_val;
                    if (htspin_val > 0) {
                      ht.setIndex(1);
                    } else if (htspin_val < 0) {
                      ht.setIndex(-1);
                    } else {
                      ht.setIndex(0);
                    }

                    // Assign the state for the combined spin state
                    ss.setIndex(ss_val);

                    // Add the event to the RooDataSet
                    ds_h->add(RooArgSet(h,t,ht,ss));
                  },
                  {h.GetName(),t.GetName(),ss.GetName()}
                  );

    // Merge datasets
    static_cast<RooDataSet&>(*rooDataSetResult).merge(&static_cast<RooDataSet&>(*ds_h));

    // Import variables into workspace
    w->import(h);
    w->import(t);
    w->import(ht);
    for (int rr=0; rr<nvars; rr++) { w->import(*rrvars[rr]); }

    // Import data into the workspace
    w->import(*rooDataSetResult);

    return;
}

/**
* @brief Map values from a CSV file into an existing RDataFrame.
*
* Load a CSV file containing, e.g., run-dependent values, with `ROOT::RDataFrame::FromCSV`.
* Then, add the data from the requested column names to an existing RDataFrame
* by matching entries for `csv_key_col` in the CSV to entries for `rdf_key_col`
* in the RDataFrame.  Note that column values will automatically be cast to float in the RDataFrame.
*
* @param filtered_df RDataFrame in which to load data from CSV
* @param rdf_key_col Name of the key column in the RDataFrame
* @param csv_path Path to the CSV file
* @param csv_key_col Name of the key column in the CSV
* @param col_names List of column names for values to map from the CSV file
* @param col_aliases Map of column names to aliases for defining branches in the RDataFrame
* @param readHeaders Option to read the headers from the CSV file
* @param delimiter Delimiter used in the CSV file
*
* @return RDataFrame with run dependent values loaded from the CSV file
*/
template<typename CsvKeyType, typename CsvValueType>
RNode mapDataFromCSV(RNode filtered_df,
                            std::string rdf_key_col,
                            std::string csv_path,
                            std::string csv_key_col,
                            std::vector<std::string> col_names,
                            std::map<std::string,std::string> col_aliases,
                            bool readHeaders=true,
                            char delimiter=','
    ) {

    // Read CSV once
    ROOT::RDataFrame csv_df = ROOT::RDF::FromCSV(csv_path, readHeaders, delimiter);

    // Get keys from csv
    auto keys = csv_df.Take<CsvKeyType>(csv_key_col);

    // Loop the column names and define variables from a CSV map
    auto df_with_new_column = filtered_df;
    for (int cc=0; cc<col_names.size(); cc++) {
        
        // Set column name using alias if available
        const std::string& csv_value_col = col_names[cc];
        std::string new_column_name = col_names[cc];
        for (auto it = col_aliases.begin(); it != col_aliases.end(); it++) {
            if (it->first == new_column_name) {
                new_column_name = it->second;
                break;
            }
        }

        // Get values from csv
        auto values = csv_df.Take<CsvValueType>(csv_value_col);

        // Capture data as simple vectors to share across threads
        std::vector<CsvKeyType> keys_vec = *keys;
        std::vector<CsvValueType> values_vec = *values;

        // Define column
        df_with_new_column = filtered_df.Define(new_column_name,
            [keys_vec, values_vec](float key_in) -> float {
                // Build map thread-local, initialized on first use
                static thread_local std::map<float, float> map;
                static thread_local bool initialized = false;

                if (!initialized) {
                    for (size_t i = 0; i < keys_vec.size(); ++i) {
                        float k = static_cast<float>(keys_vec[i]);
                        float v = static_cast<float>(values_vec[i]);
                        map[k] = v;
                    }
                    initialized = true;
                }

                auto it = map.find(key_in);
                return it != map.end() ? it->second : 0.0f;
            },
            {rdf_key_col});

    } // for (int cc=0; cc<col_names.size(); cc++) {

    return df_with_new_column;
}

/**
* @brief Inject an asymmetry into an existing RDataFrame.
*
* Inject an asymmetry into an existing `ROOT::RDataFrame` given a random seed,
* beam and target polarizations, and the relevant signal and background 
* asymmetry formulas separated into terms dependent on beam helicity, target spin, or both.
* The algorithm proceeds as follows.
* For each event, a random number \f$r\in[0,1)\f$, beam helicity \f$h_b\in(-1,0,1)\f$, and target spin \f$h_t\in(-1,0,1)\f$ are all randomly generated.
* A non-zero \f$h_b\f$ and \f$h_t\f$ are generated with probabilities taken from the beam and target polarizations respectively:
* \f$P(h_b\neq0) = \overline{P^2_b}\f$ and
* \f$P(h_t\neq0) = \overline{P^2_t}\f$.
* The event weight \f$w\f$ is computed as:
* \f[
*   w = \frac{1}{2} (1 + h_b \cdot A_{PU} + h_b \cdot A_{UP} + h_b \cdot h_t \cdot A_{PP}),
* \f]
* where \f$A_{PU}\f$, \f$A_{UP}\f$, and \f$A_{PP}\f$ are the asymmetry terms
* dependent on beam helicity, target spin, or both.
* Note that the asymmetry terms will taken from either the signal or background asymmetries
* depending on whether the event has been marked as signal or background.
* If \f$w>r\f$ the beam helicity and target spin values for that event are accepted,
* otherwise all random values are regenerated and the process repeats until \f$w>r\f$.
*
* @param df `ROOT::RDataFrame` in which to inject asymmetry
* @param seed Seed for random number generator
* @param bpol Average beam polarization
* @param tpol Average target polarization
* @param mc_sg_match_name Name of boolean column indicating signal events
* @param asyms_sg_pu_name Name of column containing the true signal asymmetries dependent on beam helicity
* @param asyms_sg_up_name Name of column containing the true signal asymmetries dependent on target spin
* @param asyms_sg_pp_name Name of column containing the true signal asymmetries dependent on beam helicity and target spin
* @param asyms_bg_pu_name Name of column containing the true background asymmetries dependent on beam helicity
* @param asyms_bg_up_name Name of column containing the true background asymmetries dependent on target spin
* @param asyms_bg_pp_name Name of column containing the true background asymmetries dependent on beam helicity and target spin
* @param randvar_name Name of column containing randomly generated values used to determine beam helicity and target spin
* @param helicity_name Name of column containing the beam helicity
* @param tspin_name Name of column containing the target spin
*
* @return `ROOT::RDataFrame` with helicity and target spin values injected
*/
RNode injectAsym(
    RNode df,
    int seed,
    double bpol,
    double tpol,
    std::string mc_sg_match_name,
    std::string asyms_sg_pu_name,
    std::string asyms_sg_up_name,
    std::string asyms_sg_pp_name,
    std::string asyms_bg_pu_name,
    std::string asyms_bg_up_name,
    std::string asyms_bg_pp_name,
    std::string randvar_name,
    std::string helicity_name,
    std::string tspin_name
    ) {

    // Define a lambda to inject an asymmetry for each rdf entry
    auto getEntrySlot = [seed,bpol,tpol](
                        ULong64_t iEntry,
                        bool mc_sg_match,
                        float asyms_sg_pu, float asyms_sg_up, float asyms_sg_pp,
                        float asyms_bg_pu, float asyms_bg_up, float asyms_bg_pp) -> int {

        // Combine global seed and row index for determinism
        TRandom3 rng(seed + static_cast<UInt_t>(iEntry));

        // Initialize beam helicity, target spin, random variable and cross-section asymmetry value
        int bhelicity    = 0;
        int tspin        = 0;
        double random_var = 0.0;
        double xs_val     = 0.0;
        
        // Generate random variables
        double b_rand_var = 0.0;
        double t_rand_var = 0.0;

        // Reassign helicity and target spin based on XS value
        while (random_var<=xs_val) {

            // Generate random variables
            b_rand_var = rng.Uniform();
            t_rand_var = rng.Uniform();
            random_var = rng.Uniform();

            // Assign helicity with probabilities given by polarization
            if ((b_rand_var>bpol && bpol>0.0) || (t_rand_var>tpol && tpol>0.0)) break;

            // Assign beam helicity and target spin
            if (bpol>0.0) bhelicity = (b_rand_var<=bpol/2.0) ? 1 : -1;
            if (tpol>0.0) tspin     = (t_rand_var<=tpol/2.0) ? 1 : -1;

            // Compute the XS value
            if (mc_sg_match) {
                xs_val = 0.5*(1.0 + bhelicity*asyms_sg_pu + tspin*asyms_sg_up + bhelicity*tspin*asyms_sg_pp);
            } else {
                xs_val = 0.5*(1.0 + bhelicity*asyms_bg_pu + tspin*asyms_bg_up + bhelicity*tspin*asyms_bg_pp);
            }

        } //  while (bhelicity==0.0 || random_var<=xs_val) {

        // Once helicity and target spin satisfy XS(h,tspin)>random_var
        return (int)((bhelicity+1)*10 + (tspin+1)); //NOTE: ENCODE AS A 2 DIGIT NUMBER ASSUMING 3 STATES FOR BEAM HELICITY AND TARGET SPIN EACH.
    };

    // Define random variable
    auto frame = df.Define(
                        randvar_name.c_str(),
                        getEntrySlot,
                        {
                            "rdfentry_", mc_sg_match_name.c_str(),
                            asyms_sg_pu_name.c_str(), asyms_sg_up_name.c_str(), asyms_sg_pp_name.c_str(),
                            asyms_bg_pu_name.c_str(), asyms_bg_up_name.c_str(), asyms_bg_pp_name.c_str()
                        }
                    )
                    .Define(helicity_name.c_str(), [](int my_rand_var) -> float {
                        if (my_rand_var / 10 == 2) return (float) 1.0;
                        if (my_rand_var / 10 == 1) return (float) 0.0;
                        if (my_rand_var / 10 == 0) return (float)-1.0;
                        return (float)0.0;
                    },
                    {randvar_name.c_str()})
                    .Define(tspin_name.c_str(), [](int my_rand_var) -> float {
                        if (my_rand_var % 10 == 2) return (float) 1.0;
                        if (my_rand_var % 10 == 1) return (float) 0.0;
                        if (my_rand_var % 10 == 0) return (float)-1.0;
                        return (float)0.0;
                    },
                    {randvar_name.c_str()});

    return frame;

}// RNode injectAsym()

} // namespace data

} // namespace saga {
