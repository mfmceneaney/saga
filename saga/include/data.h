#include <iostream>
#include <memory>
#include <string>
#include <map>
#include <thread>
#include <functional>
#include <unordered_set>
#include <vector>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RCsvDS.hxx>
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandomGen.h>

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
using std::unordered_multiset;
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
* @param weight_name Name of weight variable, ignored if empty
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
        string weight_name,
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

    string method_name = "createDataset";

    // Define weight variable
    LOG_DEBUG(Form("[%s]: Defining RooRealVar for event weight: %s", method_name.c_str(), weight_name.c_str()));
    RooRealVar *weightvar = new RooRealVar(
        weight_name.c_str(),
        "Event weight",
        1.0,        // initial value (irrelevant for datasets)
        -1e6, 1e6   // range (important!)
    );

    // Define the helicity variable
    LOG_DEBUG(Form("[%s]: Defining RooCategory for helicity variable: %s", method_name.c_str(), helicity.c_str()));
    RooCategory h(helicity.c_str(), helicity.c_str());
    for (auto it = helicity_states.begin(); it != helicity_states.end(); it++) {
        h.defineType(it->first.c_str(), it->second);
    }

    // Define the target spin variable
    LOG_DEBUG(Form("[%s]: Defining RooCategory for target spin variable: %s", method_name.c_str(), tspin.c_str()));
    RooCategory t(tspin.c_str(), tspin.c_str());
    for (auto it = tspin_states.begin(); it != tspin_states.end(); it++) {
        t.defineType(it->first.c_str(), it->second);
    }

    // Define the helicity timestarget spin variable
    LOG_DEBUG(Form("[%s]: Defining RooCategory for helicity times target spin variable: %s", method_name.c_str(), htspin.c_str()));
    RooCategory ht(htspin.c_str(), htspin.c_str());
    for (auto it = htspin_states.begin(); it != htspin_states.end(); it++) {
        ht.defineType(it->first.c_str(), it->second);
    }

    // Define the combined spin state variable
    LOG_DEBUG(Form("[%s]: Defining RooCategory for combined spin state variable: %s", method_name.c_str(), combined_spin_state.c_str()));
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
            LOG_DEBUG(Form("[%s]: Converting category to float: %s", method_name.c_str(), helicity.c_str()));

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
            LOG_DEBUG(Form("[%s]: Converting category to float: %s", method_name.c_str(), tspin.c_str()));

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
    LOG_DEBUG(Form("[%s]: Creating list of variables...", method_name.c_str()));
    vector<string> vars;
    for (int idx=0; idx<binvars.size();     idx++) vars.push_back(binvars[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) vars.push_back(depolvars[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) vars.push_back(asymfitvars[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) vars.push_back(massfitvars[idx]);
    LOG_DEBUG(Form("[%s]: Creating list of variable titles...", method_name.c_str()));
    vector<string> var_titles;
    for (int idx=0; idx<binvars.size();     idx++) var_titles.push_back(binvar_titles[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_titles.push_back(depolvar_titles[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_titles.push_back(asymfitvar_titles[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_titles.push_back(massfitvar_titles[idx]);
    LOG_DEBUG(Form("[%s]: Creating list of variable limits...", method_name.c_str()));
    vector<vector<double>> var_lims;
    for (int idx=0; idx<binvars.size();     idx++) var_lims.push_back(binvar_lims[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_lims.push_back(depolvar_lims[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_lims.push_back(asymfitvar_lims[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_lims.push_back(massfitvar_lims[idx]);
    LOG_DEBUG(Form("[%s]: Creating list of variable bins...", method_name.c_str()));
    vector<int> var_bins;
    for (int idx=0; idx<binvars.size();     idx++) var_bins.push_back(binvar_bins[idx]);
    for (int idx=0; idx<depolvars.size();   idx++) var_bins.push_back(depolvar_bins[idx]);
    for (int idx=0; idx<asymfitvars.size(); idx++) var_bins.push_back(asymfitvar_bins[idx]);
    for (int idx=0; idx<massfitvars.size(); idx++) var_bins.push_back(massfitvar_bins[idx]);
    int nvars = vars.size();

    // Define RooRealVar variables
    LOG_DEBUG(Form("[%s]: Creating RooRealVar variables...", method_name.c_str()));
    RooRealVar *rrvars[nvars];
    for (int rr=0; rr<nvars; rr++) {
        rrvars[rr] = new RooRealVar(vars[rr].c_str(), var_titles[rr].c_str(), var_lims[rr][0], var_lims[rr][1]);
        rrvars[rr]->setBins(var_bins[rr]);
    }

    // Define variable list for RooDataSetHelper
    LOG_DEBUG(Form("[%s]: Creating RooArgSet of variables...", method_name.c_str()));
    RooArgSet *argset = new RooArgSet();
    for (int rr=0; rr<nvars; rr++) {
        argset->add(*rrvars[rr]);
    }
    argset->add(*weightvar);

    // Create RDataFrame to RooDataSet pointer
    LOG_DEBUG(Form("[%s]: Creating dataset: %s , %s, with %d variables...", method_name.c_str(), name.c_str(), title.c_str(), nvars));
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
            string msg = Form("[%s]: Number of variables %d is outside the allowed range [2,18]", method_name.c_str(), nvars);
            LOG_ERROR(msg);
            throw runtime_error(msg);
    }

    // Manually create dataset containing helicity as a RooCategory variable
    LOG_DEBUG(Form("[%s]: Creating RooDataSet with RooCategories...", method_name.c_str()));
    RooDataSet *ds_h = new RooDataSet("ds_h","ds_h", RooArgSet(h,t,ht,ss));

    // Set cuts for variable limits so that new dataset will have same length as old dataset
    LOG_DEBUG(Form("[%s]: Creating variable limits cuts...", method_name.c_str()));
    string varlims_cuts = "";
    for (int vv=0; vv<nvars; vv++) {
        if (varlims_cuts.size()==0) {
            varlims_cuts = Form("(%s>=%.8f && %s<=%.8f)", vars[vv].c_str(), rrvars[vv]->getMin(), vars[vv].c_str(), rrvars[vv]->getMax());
        } else {
            varlims_cuts = Form("%s && (%s>=%.8f && %s<=%.8f)", varlims_cuts.c_str(), vars[vv].c_str(), rrvars[vv]->getMin(), vars[vv].c_str(), rrvars[vv]->getMax());
        }
    }

    // Loop RDataFrame and fill helicity dataset
    LOG_DEBUG(Form("[%s]: Filling RooDataSet with RooCategories...", method_name.c_str()));
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
    LOG_DEBUG(Form("[%s]: Merging RooDataSets...", method_name.c_str()));
    static_cast<RooDataSet&>(*rooDataSetResult).merge(&static_cast<RooDataSet&>(*ds_h));

    // Import variables into workspace
    LOG_DEBUG(Form("[%s]: Importing variables into RooWorkspace...", method_name.c_str()));
    w->import(h);
    w->import(t);
    w->import(ht);
    w->import(ss);
    for (int rr=0; rr<nvars; rr++) { w->import(*rrvars[rr]); }
    if (weight_name!="") w->import(*weightvar);

    // Import data into the workspace
    if (weight_name=="") {
        LOG_DEBUG(Form("[%s]: Importing RooDataSet into RooWorkspace...", method_name.c_str()));
        w->import(*rooDataSetResult);
    } else {
        LOG_DEBUG(Form("[%s]: Setting weight variable...", method_name.c_str()));
        auto& data = static_cast<RooDataSet&>(*rooDataSetResult);
        RooDataSet * rooDataSetResult_weighted = new RooDataSet(name.c_str(), title.c_str(), &data, *data.get(), nullptr, weight_name.c_str());
        LOG_DEBUG(Form("[%s]: Importing RooDataSet into RooWorkspace...", method_name.c_str()));
        w->import(*rooDataSetResult_weighted);
    }

    return;

}// createDataset()

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

    string method_name = "mapDataFromCSV";

    // Read CSV once
    LOG_DEBUG(Form("[%s]: Loading CSV from: %s", method_name.c_str(), csv_path.c_str()));
    ROOT::RDataFrame csv_df = ROOT::RDF::FromCSV(csv_path, readHeaders, delimiter);

    // Get keys from csv
    LOG_DEBUG(Form("[%s]: Grabbing keys from column: %s", method_name.c_str(), csv_key_col.c_str()));
    auto keys = csv_df.Take<CsvKeyType>(csv_key_col);

    // Loop the column names and define variables from a CSV map
    LOG_DEBUG(Form("[%s]: Looping CSV columns...", method_name.c_str()));
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
* @brief Initialize a `TRandom` generator
*
* Initialize a `TRandom` generator from the available algorithms
* provided by ROOT.
*
* @param seed Seed for random number generator
* @param trandom_type Type name of ROOT TRandom number generator
*
* @return `TRandom` generator of given type initialized with the given seed
*
* @throws Runtime Error
*/
TRandom * initializeTRandom(
    UInt_t seed,
    string trandom_type
    ) {

    string method_name = "initializeTRandom";

    // Initialize the generator
    LOG_DEBUG(Form("[%s]: Initializing random number generator with type %s", method_name.c_str(), trandom_type.c_str()));
    TRandom * rng;

    // Select the generator
    if (trandom_type=="TRandom3") { rng = new TRandom3(seed); }
    else if (trandom_type=="TRandomRanluxpp") { rng = new TRandomRanluxpp(seed); }
    else if (trandom_type=="TRandomMixMax") { rng = new TRandomMixMax(seed); }
    else if (trandom_type=="TRandomMixMax17") { rng = new TRandomMixMax17(seed); }
    else if (trandom_type=="TRandomMixMax256") { rng = new TRandomMixMax256(seed); }
    else if (trandom_type=="TRandomMT64") { rng = new TRandomMT64(seed); }
    else if (trandom_type=="TRandom1") { rng = new TRandom1(seed); }
    else if (trandom_type=="TRandomRanlux48") { rng = new TRandomRanlux48(seed); }
    else if (trandom_type=="TRandom2") { rng = new TRandom2(seed); }
    else {
        string msg = Form("[%s]: trandom_type %s is not recognized.  Please consult the ROOT documentation for allowed types", method_name.c_str(), trandom_type.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }

    return rng;

}// initializeTRandom()

/**
* @brief Inject an asymmetry into an existing RDataFrame.
*
* Inject an asymmetry into an existing `ROOT::RDataFrame` given a random seed,
* beam and target polarizations, and the relevant signal and background 
* asymmetry formulas separated into unpolarized modulations and modulations even under transverse target spin flips,
* i.e., modulations even under a flip of \f$\phi_{S}\f$,
* as well as asymmetry terms dependent on beam helicity, target spin, or both.
*
* In almost **all** scenarios, the unpolarized and even \f$\phi_{S}\f$ dependent modulations will **not** be needed.
* However, in the case of a term with an even dependence on \f$\phi_{S}\f$,
* the \f$\phi_{S}\f$ dependence can be injected into the dataset if a variable name is supplied
* for \f$\phi_{S}\f$ in both spin states via the arguments `phi_s_up_name` and `phi_s_dn_name`.
*
* The injection algorithm proceeds as follows.
* For each event, a random number \f$r\in[0,1)\f$, beam helicity \f$\lambda_{\ell}\in(-1,0,1)\f$, and target spin \f$S\in(-1,0,1)\f$ are all randomly generated.
* A non-zero \f$\lambda_{\ell}\f$ and \f$S\f$ are generated with probabilities taken from the beam and target polarizations respectively:
* \f$P(\lambda_{\ell}\neq0) = \overline{\lambda_{\ell}^2}\f$ and
* \f$P(S\neq0) = \overline{S^2}\f$.
* Otherwise, positive and negative helicity and spin values are generated with equal probability.
* The probability \f$w\f$ of accepting the proposed \f$(\lambda_{\ell},S)\f$ pair is:
* \f[
*   w &= \frac{1}{N} \bigg{\{} 1 + A_{UU} + \, S_{||} \, A_{UL} + A_{UT}(\phi^{True}_{S}) \\
    & \quad + \,\lambda_{\ell} \, \big{[}A_{LU} + S_{||} \, A_{LL} + A_{LT}(\phi^{True}_{S})\big{]} \bigg{\}}\,,
* \f]
* where \f$N\f$ is the number of possible combinations of \f$(\lambda_{\ell},S)\f$, given whether either has already been set to \f$0\f$.
* For example, if \f$(\lambda_{\ell},S)=(0,\pm1)\f$ or \f$(\lambda_{\ell},S)=(\pm1,0)\f$ then \f$N=2\f$,
* but if \f$(\lambda_{\ell},S)=(\pm1,\pm1)\f$ then \f$N=4\f$.
* Note that since we rely on the fact that the \f$A_{UT}\f$ terms are odd under a transverse target spin flip,
* this formulation is equivalent to the following
* \f[
*   w &= \frac{1}{N} \bigg{\{} 1 + A_{UU} + \, S \, A_{UP} \\
    & \quad + \,\lambda_{\ell} \, \big{[}A_{PU} + S \, A_{PP} \big{]} \bigg{\}}\,,
* \f]
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
* @param asyms_sg_uu_pos_name Name of column containing the true signal unpolarized modulations and modulations with an even dependence on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=+1\f$
* @param asyms_sg_uu_neg_name Name of column containing the true signal unpolarized modulations and modulations with an even dependence on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=-1\f$
* @param asyms_sg_pu_pos_name Name of column containing the true signal asymmetries dependent on beam helicity, for \f$S_{\perp}=+1\f$ in the case of a modulation even under transverse target spin flips
* @param asyms_sg_pu_neg_name Name of column containing the true signal asymmetries dependent on beam helicity, for \f$S_{\perp}=-1\f$ in the case of a modulation even under transverse target spin flips
* @param asyms_sg_up_name Name of column containing the true signal asymmetries dependent on target spin
* @param asyms_sg_pp_name Name of column containing the true signal asymmetries dependent on beam helicity and target spin
* @param asyms_bg_uu_pos_name Name of column containing the true background unpolarized modulations and modulations with an even dependence on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=+1\f$
* @param asyms_bg_uu_neg_name Name of column containing the true background unpolarized modulations and modulations with an even dependence on transverse target spin, i.e., \f$\phi_{S}\f$, for \f$S_{\perp}=-1\f$
* @param asyms_bg_pu_pos_name Name of column containing the true background asymmetries dependent on beam helicity, for \f$S_{\perp}=+1\f$ in the case of a modulation even under transverse target spin flips
* @param asyms_bg_pu_neg_name Name of column containing the true background asymmetries dependent on beam helicity, for \f$S_{\perp}=-1\f$ in the case of a modulation even under transverse target spin flips
* @param asyms_bg_up_name Name of column containing the true background asymmetries dependent on target spin
* @param asyms_bg_pp_name Name of column containing the true background asymmetries dependent on beam helicity and target spin
* @param combined_spin_state_name Name of column containing combined beam helicity and target spin state encoded as \f$ss = (\lambda_{\ell}+1)\cdot10 + (S+1)\f$
* @param helicity_name Name of column containing the beam helicity
* @param tspin_name Name of column containing the target spin
* @param phi_s_up_name Name of column containing the injected \f$\phi_{S}\f$ variable for \f$S_{\perp}=+1\f$ events
* @param phi_s_dn_name Name of column containing the injected \f$\phi_{S}\f$ variable for \f$S_{\perp}=-1\f$ events
* @param phi_s_name_injected Name of column to contain the injected \f$\phi_{S}\f$ variable
* @param trandom_type Type name of ROOT TRandom number generator
*
* @return `ROOT::RDataFrame` with helicity and target spin values injected
*
* @throws Runtime Error
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
    string phi_s_name_injected,
    string trandom_type
    ) {

    string method_name = "injectAsym";

    // Define a lambda to inject an asymmetry for each rdf entry
    LOG_DEBUG(Form("[%s]: Defining lambda function for injected spin state variable...", method_name.c_str()));
    auto getEntrySlot = [seed,bpol,tpol,trandom_type](
                        ULong64_t iEntry,
                        bool mc_sg_match,
                        float asyms_sg_uu_pos, float asyms_sg_uu_neg,
                        float asyms_sg_pu_pos, float asyms_sg_pu_neg,
                        float asyms_sg_up, float asyms_sg_pp,
                        float asyms_bg_uu_pos, float asyms_bg_uu_neg,
                        float asyms_bg_pu_pos, float asyms_bg_pu_neg,
                        float asyms_bg_up, float asyms_bg_pp) -> int {

        // Combine global seed and row index for determinism
        UInt_t seed_iEntry = seed + static_cast<UInt_t>(iEntry);
        TRandom * rng = initializeTRandom(seed_iEntry,trandom_type);

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
            b_rand_var = rng->Uniform();
            t_rand_var = rng->Uniform();
            random_var = rng->Uniform();

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
    LOG_DEBUG(Form("[%s]: Defining injected spin state variable %s", method_name.c_str(), combined_spin_state_name.c_str()));
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
        LOG_DEBUG(Form("[%s]: Defining injected phi_s variable %s", method_name.c_str(), _phi_s_name_injected.c_str()));
        frame = frame.Define(_phi_s_name_injected.c_str(), [](float tspin, float phi_s_up, float phi_s_dn) -> float { return (float)(tspin>0 ? phi_s_up : phi_s_dn);},{tspin_name.c_str(),phi_s_up_name.c_str(),phi_s_dn_name.c_str()});
    }

    return frame;

}// RNode injectAsym()

/**
* @brief Weight a dataframe by resampling with Poissonian statistics
*
* Weight an existing `ROOT::RDataFrame` following the Poissonian bootstrapping method
* of resampling each event randomly from a Poissonian distribution with mean \f$\lambda=1\f$.
*
* @param df `ROOT::RDataFrame` to weight
* @param seed Seed for random number generator
* @param weight_name Name of column containing the event weights
* @param trandom_type Type name of ROOT TRandom number generator
*
* @return `ROOT::RDataFrame` filtered for non-zero resampling weights
*/
RNode bootstrapPoisson(
    RNode df,
    int seed,
    string weight_name,
    string trandom_type
    ) {

    string method_name = "bootstrapPoisson";

    // Define a lambda to inject an asymmetry for each rdf entry
    LOG_DEBUG(Form("[%s]: Defining lambda function for Poisson resampling...", method_name.c_str()));
    auto getEntrySlot = [seed,trandom_type](ULong64_t iEntry) -> int {

        // Combine global seed and row index for determinism
        UInt_t seed_iEntry = seed + static_cast<UInt_t>(iEntry);
        TRandom * rng = initializeTRandom(seed_iEntry,trandom_type);

        // Draw from Poisson distribution with mean lambda=1
        return rng->Poisson(1.0);
    };

    // Define the event weights and filter for non-zero weights
    return df.Define(
                weight_name.c_str(),
                getEntrySlot,
                {"rdfentry_"}
            ).Filter(Form("%s > 0",weight_name.c_str()));

}// RNode bootstrapPoisson()

/**
* @brief Weight a dataframe by resampling with replacement
*
* Weight an existing `ROOT::RDataFrame` following the classical bootstrapping method
* of resampling with replacement.
*
* @param df `ROOT::RDataFrame` to weight
* @param n Sample size
* @param seed Seed for random number generator
* @param weight_name Name of column containing the event weights
* @param trandom_type Type name of ROOT TRandom number generator
*
* @return `ROOT::RDataFrame` filtered for non-zero resampling weights
*
* @throws Runtime Error
*/
RNode bootstrapClassical(
        RNode df,
        int n,
        int seed,
        string weight_name,
        string trandom_type
    ) {

    string method_name = "bootstrapClassical";

    // Grab indices
    LOG_DEBUG(Form("[%s]: Grabbing event indices...", method_name.c_str()));
    auto entries = *df.Take<ULong64_t>("rdfentry_");

    // Check dataframe and resampling size
    const ULong64_t N = entries.size();
    if (N==0) {
        string msg = Form("[%s]: Dataframe size is 0.  Cannot resample.", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }
    if (n==0) {
        string msg = Form("[%s]: Resampling size 0 is not allowed.", method_name.c_str());
        LOG_ERROR(msg);
        throw runtime_error(msg);
    }

    // Initialize the generator
    LOG_DEBUG(Form("[%s]: Initialize random number generator with type, seed = %s, %d ...", method_name.c_str(), trandom_type.c_str(), seed));
    TRandom * rng = initializeTRandom(seed,trandom_type);

    // Resample the dataset
    LOG_DEBUG(Form("[%s]: Resampling dataset...", method_name.c_str()));
    unordered_multiset<ULong64_t> selected;
    for (ULong64_t i = 0; i < n; ++i) {
        selected.insert(rng->Integer(N));
    }

    // Set event weights and filter for non-zero weights
    LOG_DEBUG(Form("[%s]: Defining event weight %s...", method_name.c_str(), weight_name.c_str()));
    return df.Define(
                weight_name.c_str(),
                [selected](ULong64_t entry) {
                    return selected.count(entry);
                },
                {"rdfentry_"}
            ).Filter(Form("%s > 0",weight_name.c_str()));

}// RNode bootstrapClassical()

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
    
    string method_name = "defineAngularDiffVars";
    
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
    LOG_DEBUG(Form("[%s]: Defining angular difference variables %s and %s", method_name.c_str(), dtheta_vars[0].c_str(), dphi_vars[0].c_str()));
    auto newframe = frame.Define(dtheta_vars[0].c_str(),[](float theta, float theta_mc){ return TMath::Abs(theta-theta_mc); },{theta_vars[0].c_str(),theta_mc_vars[0].c_str()})
        .Define(dphi_vars[0].c_str(),[](float phi, float phi_mc){
            return (float) (TMath::Abs(phi-phi_mc)<TMath::Pi()
            ? TMath::Abs(phi-phi_mc) : 2*TMath::Pi() - TMath::Abs(phi-phi_mc));
            },{phi_vars[0].c_str(),phi_mc_vars[0].c_str()});
    for (int idx=1; idx<particle_suffixes.size(); idx++) {
        LOG_DEBUG(Form("[%s]: Defining angular difference variables %s and %s", method_name.c_str(), dtheta_vars[idx].c_str(), dphi_vars[idx].c_str()));
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
