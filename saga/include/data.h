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

// Local Includes
#include <log.h>

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
* @param categories_as_float List of category variables to include as floats named `<category>_as_float` in dataset
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
        string name,
        string title,
        vector<string> categories_as_float,
        string helicity,
        map<string,int> helicity_states,
        string tspin,
        map<string,int> tspin_states,
        string htspin,
        map<string,int> htspin_states,
        string combined_spin_state,
        vector<string> binvars,
        vector<string> binvar_titles,
        vector<vector<double>> binvar_lims,
        vector<int> binvar_bins,
        vector<string> depolvars,
        vector<string> depolvar_titles,
        vector<vector<double>> depolvar_lims,
        vector<int> depolvar_bins,
        vector<string> &asymfitvars,
        vector<string> &asymfitvar_titles,
        vector<vector<double>> &asymfitvar_lims,
        vector<int> &asymfitvar_bins,
        vector<string> massfitvars,
        vector<string> massfitvar_titles,
        vector<vector<double>> massfitvar_lims,
        vector<int> massfitvar_bins
    ) {

    // Define the helicity variable
    LOG_DEBUG(Form("Defining RooCategory for helicity variable: %s", helicity.c_str()));
    RooCategory h(helicity.c_str(), helicity.c_str());
    for (auto it = helicity_states.begin(); it != helicity_states.end(); it++) {
        h.defineType(it->first.c_str(), it->second);
    }

    // Define the target spin variable
    LOG_DEBUG(Form("Defining RooCategory for target spin variable: %s", tspin.c_str()));
    RooCategory t(tspin.c_str(), tspin.c_str());
    for (auto it = tspin_states.begin(); it != tspin_states.end(); it++) {
        t.defineType(it->first.c_str(), it->second);
    }

    // Define the helicity timestarget spin variable
    LOG_DEBUG(Form("Defining RooCategory for helicity times target spin variable: %s", htspin.c_str()));
    RooCategory ht(htspin.c_str(), htspin.c_str());
    for (auto it = htspin_states.begin(); it != htspin_states.end(); it++) {
        ht.defineType(it->first.c_str(), it->second);
    }

    // Define the combined spin state variable
    LOG_DEBUG(Form("Defining RooCategory for combined spin state variable: %s", combined_spin_state.c_str()));
    RooCategory ss(combined_spin_state.c_str(), combined_spin_state.c_str());
    for (auto it_h = helicity_states.begin(); it_h != helicity_states.end(); it_h++) {
        for (auto it_t = htspin_states.begin(); it_t != htspin_states.end(); it_t++) {
            string state_name = Form("%s_%s",it_h->first.c_str(),it_t->first.c_str());
            int state_value = (it_h->second + 1) * 10 + (it_t->second + 1);
            ss.defineType(state_name.c_str(), state_value);
        }
    }

    // Loop categories to convert include as floats
    for (int idx=0; idx<categories_as_float.size(); idx++) {

        // Check for beam helicity variable
        if (categories_as_float[idx]==helicity) {
            LOG_DEBUG(Form("Converting category to float: %s", helicity.c_str()));

            // Set helicity variable name and formula
            string h_as_float = Form("%s_as_float",helicity.c_str());
            string h_as_float_formula = Form("(float)(%s)",helicity.c_str());

            // Define the new branch in the frame
            frame = frame.Define(h_as_float.c_str(),h_as_float_formula.c_str());

            // Check if it's already been added before
            bool contains_var = false;
            for (int ff = 0; ff<asymfitvars.size(); ff++) {
                if (asymfitvars[ff]==h_as_float) contains_var = true;
            }
            if (!contains_var) {
                string h_as_float_title = helicity;
                int h_as_float_bins=0;
                double h_as_float_min = 0.0;
                double h_as_float_max = 0.0;
                for (auto it = helicity_states.begin(); it != helicity_states.end(); it++) {
                    if (it->second<h_as_float_min) h_as_float_min = it->second;
                    if (it->second>h_as_float_max) h_as_float_max = it->second;
                    h_as_float_bins++;
                }
                vector<double> h_as_float_lims = {h_as_float_min,h_as_float_max};

                // Add helicity variable as a fit variable
                asymfitvars.insert(asymfitvars.begin(),h_as_float);
                asymfitvar_titles.insert(asymfitvar_titles.begin(),h_as_float_title);
                asymfitvar_lims.insert(asymfitvar_lims.begin(),h_as_float_lims);
                asymfitvar_bins.insert(asymfitvar_bins.begin(),h_as_float_bins);
            }
        }
    }

    // Loop categories to convert include as floats
    for (int idx=0; idx<categories_as_float.size(); idx++) {

        // Check for target spin variable
        if (categories_as_float[idx]==tspin) {
            LOG_DEBUG(Form("Converting category to float: %s", tspin.c_str()));

            // Set helicity variable name and formula
            string h_as_float = Form("%s_as_float",tspin.c_str());
            string h_as_float_formula = Form("(float)(%s)",tspin.c_str());

            // Define the new branch in the frame
            frame = frame.Define(h_as_float.c_str(),h_as_float_formula.c_str());

            // Check if it's already been added before
            bool contains_var = false;
            for (int ff = 0; ff<asymfitvars.size(); ff++) {
                if (asymfitvars[ff]==h_as_float) contains_var = true;
            }
            if (!contains_var) {
                string h_as_float_title = tspin;
                int h_as_float_bins=0;
                double h_as_float_min = 0.0;
                double h_as_float_max = 0.0;
                for (auto it = tspin_states.begin(); it != tspin_states.end(); it++) {
                    if (it->second<h_as_float_min) h_as_float_min = it->second;
                    if (it->second>h_as_float_max) h_as_float_max = it->second;
                    h_as_float_bins++;
                }
                vector<double> h_as_float_lims = {h_as_float_min,h_as_float_max};

                // Add helicity variable as a fit variable
                asymfitvars.insert(asymfitvars.begin(),h_as_float);
                asymfitvar_titles.insert(asymfitvar_titles.begin(),h_as_float_title);
                asymfitvar_lims.insert(asymfitvar_lims.begin(),h_as_float_lims);
                asymfitvar_bins.insert(asymfitvar_bins.begin(),h_as_float_bins);
            }
        }
    }

    // Define the full variables and limits lists
    LOG_DEBUG("Creating list of variables...");
    vector<string> vars;
    for (int idx=0; idx<binvars.size();     idx++) vars.push_back(binvars[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) vars.push_back(depolvars[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) vars.push_back(asymfitvars[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) vars.push_back(massfitvars[idx]);
    LOG_DEBUG("Creating list of variable titles...");
    vector<string> var_titles;
    for (int idx=0; idx<binvars.size();     idx++) var_titles.push_back(binvar_titles[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_titles.push_back(depolvar_titles[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_titles.push_back(asymfitvar_titles[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_titles.push_back(massfitvar_titles[idx]);
    LOG_DEBUG("Creating list of variable limits...");
    vector<vector<double>> var_lims;
    for (int idx=0; idx<binvars.size();     idx++) var_lims.push_back(binvar_lims[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_lims.push_back(depolvar_lims[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_lims.push_back(asymfitvar_lims[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_lims.push_back(massfitvar_lims[idx]);
    LOG_DEBUG("Creating list of variable bins...");
    vector<int> var_bins;
    for (int idx=0; idx<binvars.size();     idx++) var_bins.push_back(binvar_bins[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_bins.push_back(depolvar_bins[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_bins.push_back(asymfitvar_bins[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_bins.push_back(massfitvar_bins[idx]);
    int nvars = vars.size();

    // Define RooRealVar variables
    LOG_DEBUG("Creating RooRealVar variables...");
    RooRealVar *rrvars[nvars];
    for (int rr=0; rr<nvars; rr++) {
        rrvars[rr] = new RooRealVar(vars[rr].c_str(), var_titles[rr].c_str(), var_lims[rr][0], var_lims[rr][1]);
        rrvars[rr]->setBins(var_bins[rr]);
    }

    // Define variable list for RooDataSetHelper
    LOG_DEBUG("Creating RooArgSet of variables...");
    RooArgSet *argset = new RooArgSet();
    for (int rr=0; rr<nvars; rr++) {
        argset->add(*rrvars[rr]);
    }

    // Create RDataFrame to RooDataSet pointer
    LOG_DEBUG(Form("Creating dataset: %s , %s, with %d variables...", name.c_str(), title.c_str(), nvars));
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
            string msg = Form("Number of variables %d is outside the allowed range [2,18]", nvars);
            LOG_ERROR(msg);
            throw runtime_error(msg);
    }

    // Manually create dataset containing helicity as a RooCategory variable
    LOG_DEBUG("Creating RooDataSet with RooCategories...");
    RooDataSet *ds_h = new RooDataSet("ds_h","ds_h", RooArgSet(h,t,ht,ss));

    // Set cuts for variable limits so that new dataset will have same length as old dataset
    LOG_DEBUG("Creating variable limits cuts...");
    string varlims_cuts = "";
    for (int vv=0; vv<nvars; vv++) {
        if (varlims_cuts.size()==0) {
            varlims_cuts = Form("(%s>=%.8f && %s<=%.8f)", vars[vv].c_str(), rrvars[vv]->getMin(), vars[vv].c_str(), rrvars[vv]->getMax());
        } else {
            varlims_cuts = Form("%s && (%s>=%.8f && %s<=%.8f)", varlims_cuts.c_str(), vars[vv].c_str(), rrvars[vv]->getMin(), vars[vv].c_str(), rrvars[vv]->getMax());
        }
    }

    // Loop RDataFrame and fill helicity dataset
    LOG_DEBUG("Filling RooDataSet with RooCategories...");
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
    LOG_DEBUG("Merging RooDataSets...");
    static_cast<RooDataSet&>(*rooDataSetResult).merge(&static_cast<RooDataSet&>(*ds_h));

    // Import variables into workspace
    LOG_DEBUG("Importing variables into RooWorkspace...");
    w->import(h);
    w->import(t);
    w->import(ht);
    w->import(ss);
    for (int rr=0; rr<nvars; rr++) { w->import(*rrvars[rr]); }

    // Import data into the workspace
    LOG_DEBUG("Importing RooDataSet into RooWorkspace...");
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
                            string rdf_key_col,
                            string csv_path,
                            string csv_key_col,
                            vector<string> col_names,
                            map<string,string> col_aliases,
                            bool readHeaders=true,
                            char delimiter=','
    ) {

    // Read CSV once
    LOG_DEBUG(Form("Loading CSV from: %s", csv_path.c_str()));
    ROOT::RDataFrame csv_df = ROOT::RDF::FromCSV(csv_path, readHeaders, delimiter);

    // Get keys from csv
    LOG_DEBUG(Form("Grabing keys from column: %s", csv_key_col.c_str()));
    auto keys = csv_df.Take<CsvKeyType>(csv_key_col);

    // Loop the column names and define variables from a CSV map
    LOG_DEBUG("Looping CSV columns...");
    auto df_with_new_column = filtered_df;
    for (int cc=0; cc<col_names.size(); cc++) {
        
        // Set column name using alias if available
        const string& csv_value_col = col_names[cc];
        string new_column_name = col_names[cc];
        for (auto it = col_aliases.begin(); it != col_aliases.end(); it++) {
            if (it->first == new_column_name) {
                new_column_name = it->second;
                break;
            }
        }

        // Get values from csv
        auto values = csv_df.Take<CsvValueType>(csv_value_col);

        // Capture data as simple vectors to share across threads
        vector<CsvKeyType> keys_vec = *keys;
        vector<CsvValueType> values_vec = *values;

        // Define column
        df_with_new_column = filtered_df.Define(new_column_name,
            [keys_vec, values_vec](float key_in) -> float {
                // Build map thread-local, initialized on first use
                static thread_local map<float, float> map;
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
* asymmetry formulas separated into unpolarized and (only) transverse target spin, i.e., \f$\phi_{S}\f$, dependent asymmetry terms,
* as well as asymmetry terms dependent on beam helicity, target spin, or both.
* In the case of an asymmetry term depending on transverse target spin, the \f$\phi_{S}\f$ variable can be injected into the dataset if a variable name is supplied.
* The injection algorithm proceeds as follows.
* For each event, a random number \f$r\in[0,1)\f$, beam helicity \f$\lambda_{\ell}\in(-1,0,1)\f$, and target spin \f$S\in(-1,0,1)\f$ are all randomly generated.
* A non-zero \f$\lambda_{\ell}\f$ and \f$S\f$ are generated with probabilities taken from the beam and target polarizations respectively:
* \f$P(\lambda_{\ell}\neq0) = \overline{\lambda_{\ell}^2}\f$ and
* \f$P(S\neq0) = \overline{S^2}\f$.
* Otherwise, positive and negative helicity and spin values are generated with equal probability.
* The probability \f$w\f$ of accepting the proposed \f$(\lambda_{\ell},S)\f$ pair is:
* \f[
*   w = \frac{1}{N} (1 + A_{UU} + A_{UT}(\phi_{S}) + \lambda_{\ell} \cdot A_{PU}(\phi_{S}) + \lambda_{\ell} \cdot A_{UP} + \lambda_{\ell} \cdot S_{||} \cdot A_{PP}),
* \f]
* where \f$N\f$ is the number of possible combinations of \f$(\lambda_{\ell},S)\f$, given whether either has already been set to \f$0\f$.
* For example, if \f$(\lambda_{\ell},S)=(0,\pm1)\f$ or \f$(\lambda_{\ell},S)=(\pm1,0)\f$ then \f$N=2\f$,
* but if \f$(\lambda_{\ell},S)=(\pm1,\pm1)\f$ then \f$N=4\f$.
* \f$A_{UU}\f$, \f$A_{UT}\f$ are the unpolarized and (only) transverse target spin dependent asymmetry terms
* and \f$A_{PU}\f$, \f$A_{UP}\f$, and \f$A_{PP}\f$ are the asymmetry terms
* dependent on beam helicity, target spin, or both.
* The asymmetry terms will be taken from either the signal or background asymmetries
* according to the boolean variable `mc_sg_match_name` indicating signal events.
* If \f$r<w\f$ the beam helicity and target spin values for that event are accepted,
* otherwise all random values are regenerated and the process repeats until \f$r<w\f$.
*
* @param df `ROOT::RDataFrame` in which to inject asymmetry
* @param seed Seed for random number generator
* @param bpol Average beam polarization
* @param tpol Average target polarization
* @param mc_sg_match_name Name of boolean column indicating signal events
* @param asyms_sg_uu_pos_name Name of column containing the true signal unpolarized asymmetries and asymmetries only dependent on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=+1\f$
* @param asyms_sg_uu_neg_name Name of column containing the true signal unpolarized asymmetries and asymmetries only dependent on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=-1\f$
* @param asyms_sg_pu_pos_name Name of column containing the true signal asymmetries dependent on beam helicity, for \f$S_{\perp}=+1\f$
* @param asyms_sg_pu_neg_name Name of column containing the true signal asymmetries dependent on beam helicity, for \f$S_{\perp}=-1\f$
* @param asyms_sg_up_name Name of column containing the true signal asymmetries dependent on target spin
* @param asyms_sg_pp_name Name of column containing the true signal asymmetries dependent on beam helicity and target spin
* @param asyms_bg_uu_pos_name Name of column containing the true background unpolarized asymmetries and asymmetries only dependent on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=+1\f$
* @param asyms_bg_uu_neg_name Name of column containing the true background unpolarized asymmetries and asymmetries only dependent on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=-1\f$
* @param asyms_bg_pu_pos_name Name of column containing the true background asymmetries dependent on beam helicity, for \f$S_{\perp}=+1\f$
* @param asyms_bg_pu_neg_name Name of column containing the true background asymmetries dependent on beam helicity, for \f$S_{\perp}=-1\f$
* @param asyms_bg_up_name Name of column containing the true background asymmetries dependent on target spin
* @param asyms_bg_pp_name Name of column containing the true background asymmetries dependent on beam helicity and target spin
* @param combined_spin_state_name Name of column containing combined beam helicity and target spin state encoded as \f$ss = (\lambda_{\ell}+1)\cdot10 + (S+1)\f$
* @param helicity_name Name of column containing the beam helicity
* @param tspin_name Name of column containing the target spin
* @param phi_s_up_name Name of column containing the injected \f$\phi_{S}\f$ variable for \f$S_{\perp}=+1\f$ events
* @param phi_s_dn_name Name of column containing the injected \f$\phi_{S}\f$ variable for \f$S_{\perp}=-1\f$ events
* @param phi_s_name_injected Name of column to contain the injected \f$\phi_{S}\f$ variable
*
* @return `ROOT::RDataFrame` with helicity and target spin values injected
*/
RNode injectAsym(
    RNode df,
    int seed,
    double bpol,
    double tpol,
    string mc_sg_match_name,
    string asyms_sg_uu_pos_name,
    string asyms_sg_uu_neg_name,
    string asyms_sg_pu_pos_name,
    string asyms_sg_pu_neg_name,
    string asyms_sg_up_name,
    string asyms_sg_pp_name,
    string asyms_bg_uu_pos_name,
    string asyms_bg_uu_neg_name,
    string asyms_bg_pu_pos_name,
    string asyms_bg_pu_neg_name,
    string asyms_bg_up_name,
    string asyms_bg_pp_name,
    string combined_spin_state_name,
    string helicity_name,
    string tspin_name,
    string phi_s_up_name,
    string phi_s_dn_name,
    string phi_s_name_injected
    ) {

    // Define a lambda to inject an asymmetry for each rdf entry
    LOG_DEBUG("Defining lambda function for injected spin state variable...");
    auto getEntrySlot = [seed,bpol,tpol](
                        ULong64_t iEntry,
                        bool mc_sg_match,
                        float asyms_sg_uu_pos, float asyms_sg_uu_neg,
                        float asyms_sg_pu_pos, float asyms_sg_pu_neg,
                        float asyms_sg_up, float asyms_sg_pp,
                        float asyms_bg_uu_pos, float asyms_bg_uu_neg,
                        float asyms_bg_pu_pos, float asyms_bg_pu_neg,
                        float asyms_bg_up, float asyms_bg_pp) -> int {

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
        while (random_var>=xs_val) {

            // Generate random variables
            b_rand_var = rng.Uniform();
            t_rand_var = rng.Uniform();
            random_var = rng.Uniform();

            // Assign helicity with probabilities given by polarization
            if (b_rand_var>bpol && bpol>0.0 && tpol==0.0) break; //NOTE: Case you inject only beam helicity and beam helicity==0.
            if (t_rand_var>tpol && tpol>0.0 && bpol==0.0) break; //NOTE: Case you inject only target spin and target spin==0.
            if (b_rand_var>bpol && bpol>0.0 && t_rand_var>tpol && tpol>0.0) break; //NOTE: Case you inject both beam helicity and target spin and both are 0.

            // Assign beam helicity and target spin
            if (bpol>0.0 && b_rand_var<=bpol) bhelicity = (b_rand_var<=bpol/2.0) ? 1 : -1;
            if (tpol>0.0 && t_rand_var<=tpol) tspin     = (t_rand_var<=tpol/2.0) ? 1 : -1;

            // Set the phi_s dependent asymmetries
            float asyms_sg_uu = (tspin>0) ? asyms_sg_uu_pos : asyms_sg_uu_neg; //NOTE: Double check this...what to do if tspin==0? -> Need to add separate arguments for UT LT XSs so that you can drop those in case tspin=0.
            float asyms_bg_uu = (tspin>0) ? asyms_bg_uu_pos : asyms_bg_uu_neg;
            float asyms_sg_pu = (tspin>0) ? asyms_sg_pu_pos : asyms_sg_pu_neg;
            float asyms_bg_pu = (tspin>0) ? asyms_bg_pu_pos : asyms_bg_pu_neg;

            // Set the number of pdfs
            double npdfs     = (bhelicity!=0 && bpol>0.0 && tspin!=0 && tpol>0.0) ? 4.0 : 2.0;
            double prefactor = 1.0/npdfs;

            // Compute the XS value
            if (mc_sg_match) {
                xs_val = prefactor*(1.0 + asyms_sg_uu + bhelicity*asyms_sg_pu + tspin*asyms_sg_up + bhelicity*tspin*asyms_sg_pp);
            } else {
                xs_val = prefactor*(1.0 + asyms_bg_uu + bhelicity*asyms_bg_pu + tspin*asyms_bg_up + bhelicity*tspin*asyms_bg_pp);
            }

        } // while (random_var>=xs_val) {

        // Once helicity and target spin satisfy XS(h,tspin)>random_var
        return (int)((bhelicity+1)*10 + (tspin+1)); //NOTE: ENCODE AS A 2 DIGIT NUMBER ASSUMING 3 STATES FOR BEAM HELICITY AND TARGET SPIN EACH.
    };

    // Define random variable
    LOG_DEBUG(Form("Defining injected spin state variable %s", combined_spin_state_name.c_str()));
    auto frame = df.Define(
                        combined_spin_state_name.c_str(),
                        getEntrySlot,
                        {
                            "rdfentry_", mc_sg_match_name.c_str(),
                            asyms_sg_uu_pos_name.c_str(), asyms_sg_uu_neg_name.c_str(),
                            asyms_sg_pu_pos_name.c_str(), asyms_sg_pu_neg_name.c_str(),
                            asyms_sg_up_name.c_str(), asyms_sg_pp_name.c_str(),
                            asyms_bg_uu_pos_name.c_str(), asyms_bg_uu_neg_name.c_str(),
                            asyms_bg_pu_pos_name.c_str(), asyms_bg_pu_neg_name.c_str(),
                            asyms_bg_up_name.c_str(), asyms_bg_pp_name.c_str()
                        }
                    )
                    .Define(helicity_name.c_str(), [](int my_rand_var) -> float {
                        if (my_rand_var / 10 == 2) return (float) 1.0;
                        if (my_rand_var / 10 == 1) return (float) 0.0;
                        if (my_rand_var / 10 == 0) return (float)-1.0;
                        return (float)0.0;
                    },
                    {combined_spin_state_name.c_str()})
                    .Define(tspin_name.c_str(), [](int my_rand_var) -> float {
                        if (my_rand_var % 10 == 2) return (float) 1.0;
                        if (my_rand_var % 10 == 1) return (float) 0.0;
                        if (my_rand_var % 10 == 0) return (float)-1.0;
                        return (float)0.0;
                    },
                    {combined_spin_state_name.c_str()});

    // Define tspin dependent phi_s AFTER injecting the asymmetry
    if (phi_s_up_name!="" && phi_s_dn_name!="") {
        string _phi_s_name_injected = phi_s_name_injected;
        if (_phi_s_name_injected=="") _phi_s_name_injected = Form("%s_injected",phi_s_up_name.c_str());
        LOG_DEBUG(Form("Defining injected phi_s variable %s", _phi_s_name_injected.c_str()));
        frame = frame.Define(_phi_s_name_injected.c_str(), [](float tspin, float phi_s_up, float phi_s_dn) -> float { return (float)(tspin>0 ? phi_s_up : phi_s_dn);},{tspin_name.c_str(),phi_s_up_name.c_str(),phi_s_dn_name.c_str()});
    }

    return frame;

}// RNode injectAsym()

/**
* @brief Define Monte Carlo (MC) simulation angular difference variables
*
* Define the angular difference variables by taking the difference of reconstructed and MC values:
* - \f$\Delta\theta = |\theta_{Rec} - \theta_{MC}|\f$
* - \f$\Delta\phi = |\phi_{Rec} - \phi_{MC}|\f$.
*
* Note that \f$\phi\f$ is a cyclic variable on \f$2\pi\f$, so if \f$|\phi_{Rec} - \phi_{MC}|>\pi\f$ then:
* - \f$\Delta\phi = 2\pi - |\phi_{Rec} - \phi_{MC}|\f$.
*
* @param frame `ROOT::RDataFrame` in which to define angular difference variables
* @param particle_suffixes Suffixes of particle variables to define angular difference variables for
* @param theta_name Name of theta variable
* @param phi_name Name of phi variable
* @param mc_suffix Suffix of MC variables
*
* @return `ROOT::RDataFrame` with angular difference variables defined
*/
RNode defineAngularDiffVars(
        RNode frame,
        vector<string> particle_suffixes,
        string theta_name = "theta",
        string phi_name   = "phi",
        string mc_suffix  = "_mc"
    ) {
    
    // Define angular difference variable names
    vector<string> theta_vars;
    vector<string> phi_vars;
    vector<string> theta_mc_vars;
    vector<string> phi_mc_vars;
    vector<string> dtheta_vars;
    vector<string> dphi_vars;
    for (int idx=0; idx<particle_suffixes.size(); idx++) {
        theta_vars.push_back(Form("%s%s",theta_name.c_str(),particle_suffixes[idx].c_str()));
        phi_vars.push_back(Form("%s%s",phi_name.c_str(),particle_suffixes[idx].c_str()));
        theta_mc_vars.push_back(Form("%s%s%s",theta_name.c_str(),particle_suffixes[idx].c_str(),mc_suffix.c_str()));
        phi_mc_vars.push_back(Form("%s%s%s",phi_name.c_str(),particle_suffixes[idx].c_str(),mc_suffix.c_str()));
        dtheta_vars.push_back(Form("d%s%s",theta_name.c_str(),particle_suffixes[idx].c_str()));
        dphi_vars.push_back(Form("d%s%s",phi_name.c_str(),particle_suffixes[idx].c_str()));
    }

    // Define angular difference variable branches
    if (particle_suffixes.size()==0) {
        return frame;
    }
    LOG_DEBUG(Form("Defining angular difference variables %s and %s", dtheta_vars[0].c_str(), dphi_vars[0].c_str()));
    auto newframe = frame.Define(dtheta_vars[0].c_str(),[](float theta, float theta_mc){ return TMath::Abs(theta-theta_mc); },{theta_vars[0].c_str(),theta_mc_vars[0].c_str()})
        .Define(dphi_vars[0].c_str(),[](float phi, float phi_mc){
            return (float) (TMath::Abs(phi-phi_mc)<TMath::Pi()
            ? TMath::Abs(phi-phi_mc) : 2*TMath::Pi() - TMath::Abs(phi-phi_mc));
            },{phi_vars[0].c_str(),phi_mc_vars[0].c_str()});
    for (int idx=1; idx<particle_suffixes.size(); idx++) {
        LOG_DEBUG(Form("Defining angular difference variables %s and %s", dtheta_vars[idx].c_str(), dphi_vars[idx].c_str()));
        newframe = newframe.Define(dtheta_vars[idx].c_str(),[](float theta, float theta_mc){ return TMath::Abs(theta-theta_mc); },{theta_vars[idx].c_str(),theta_mc_vars[idx].c_str()})
            .Define(dphi_vars[idx].c_str(),[](float phi, float phi_mc){
                return (float) (TMath::Abs(phi-phi_mc)<TMath::Pi()
                ? TMath::Abs(phi-phi_mc) : 2*TMath::Pi() - TMath::Abs(phi-phi_mc));
                },{phi_vars[idx].c_str(),phi_mc_vars[idx].c_str()});
    }
    return newframe;
}

} // namespace data

} // namespace saga {
