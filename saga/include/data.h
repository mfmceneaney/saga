#include <iostream>
#include <memory>
#include <string>

// ROOT Includes
#include <ROOT/RDataFrame.hxx>

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

namespace data {

/**
* @brief Create a dataset for an asymmetry fit.
*
* Create a RooFit dataset for an asymmetry fit from a ROOT RDataFrame,
* adding helicity, fit, invariant mass fit, binning, and depolarization variables.
* Store RooFit dataset and corresponding variables in a RooWorkspace.
*
* @param frame ROOT RDataframe from which to create a RooDataSet
* @param w RooWorkspace in which to work
* @param name Dataset name
* @param title Dataset title
* @param helicity Name of helicity variable
* @param helicity_states Map of state names to helicity values
* @param fitvars List of fit variables
* @param fitvarlims List of minimum and maximum bounds for each fit variable
* @param binvars List of kinematic binning variables (up to 4)
* @param binvarlims List of minimum and maximum bounds for each kinematic binning variable
* @param depolvars List of depolarization variables (up to 5)
* @param depolvarlims List of minimum and maximum bounds for each depolarization variable
*/
void createDataset(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> frame, //NOTE: FRAME SHOULD ALREADY BE FILTERED
        RooWorkspace *w,
        std::string name,
        std::string title,
        std::string helicity,
        std::map<std::string,int> helicity_states,
        std::vector<std::string> fitvars,
        std::vector<std::vector<double>> fitvarlims,
        std::vector<std::string> binvars,
        std::vector<std::vector<double>> binvarlims, //NOTE: THAT THESE SHOULD JUST BE THE OUTERMOST LIMITS!!!
        std::vector<std::string> depolvars,
        std::vector<std::vector<double>> depolvarlims //NOTE: THAT THESE SHOULD JUST BE THE OUTERMOST LIMITS!!!
    ) {

    // Define the helicity variable
    RooCategory h(helicity.c_str(), helicity.c_str());
    for (auto it = helicity_states.begin(); it != helicity_states.end(); it++) {
        h.defineType(it->first.c_str(), it->second);
    }

    // Define the full variables and limits lists
    std::vector<std::string> vars;
    for (int idx=0; idx<fitvars.size();   idx++) vars.push_back(fitvars[idx]);
    for (int idx=0; idx<binvars.size();   idx++) vars.push_back(binvars[idx]);
    for (int idx=0; idx<depolvars.size(); idx++) vars.push_back(depolvars[idx]);
    std::vector<std::vector<double>> varlims;
    for (int idx=0; idx<fitvarlims.size();   idx++) varlims.push_back(fitvarlims[idx]);
    for (int idx=0; idx<binvarlims.size();   idx++) varlims.push_back(binvarlims[idx]);
    for (int idx=0; idx<depolvarlims.size(); idx++) varlims.push_back(depolvarlims[idx]);
    int nvars = vars.size();

    // Define RooRealVar variables
    RooRealVar *rrvars[nvars];
    for (int rr=0; rr<nvars; rr++) {
        rrvars[rr] = new RooRealVar(vars[rr].c_str(), vars[rr].c_str(), varlims[rr][0], varlims[rr][1]);
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
    RooDataSet *ds_h = new RooDataSet("ds_h","ds_h", RooArgSet(h));

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
                  [&h,&ds_h](float val){

                    // Assign the state for the categorical variable (assuming 1 or -1 as the states)
                    if (val > 0) {
                      h.setIndex(1);
                    } else if (val < 0) {
                      h.setIndex(-1);
                    } else {
                      h.setIndex(0);
                    }
                    // Add the event to the RooDataSet
                    ds_h->add(RooArgSet(h));
                  },
                  {h.GetName()}
                  );

    // Merge datasets
    static_cast<RooDataSet&>(*rooDataSetResult).merge(&static_cast<RooDataSet&>(*ds_h));

    // Import variables into workspace
    w->import(h);
    for (int rr=0; rr<nvars; rr++) { w->import(*rrvars[rr]); }

    // Import data into the workspace
    w->import(*rooDataSetResult);

    return;
}

} // namespace data
