#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <map>
#include <TString.h>

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 7/Feb./2024
* @version 0.0.0
* @brief Utility functions for asymmetry fits.
*/

namespace saga {

namespace util {

/**
* @brief Find and replace function for std::string.
*
* Find all occurences of a substring and replace them with another.
* Note that this is an in place operation. REGEX is NOT supported.
*
* @param s String to search
* @param to_replace Substring to find and replace
* @param replace_with Substring to insert
*/
void replaceAll(
        std::string& s,
        std::string const& to_replace,
        std::string const& replace_with
    ) {
    std::string buf;
    std::size_t pos = 0;
    std::size_t prev_pos;

    // Reserves rough estimate of final size of string.
    buf.reserve(s.size());

    while (true) {
        prev_pos = pos;
        pos = s.find(to_replace, pos);
        if (pos == std::string::npos)
            break;
        buf.append(s, prev_pos, pos - prev_pos);
        buf += replace_with;
        pos += to_replace.size();
    }

    buf.append(s, prev_pos, s.size() - prev_pos);
    s.swap(buf);
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
* @brief Add additional variable limit cuts to an existing cut string.
*
* Add additional variable limit cuts to an existing cut string.
* Note that this is an in place operation.
*
* @param cuts Substring to insert
* @param vars String to search
* @param varlims Substring to find and replace
*/
void addLimitCuts(
        std::string                      cuts,
        std::vector<std::string>         vars,
        std::vector<std::vector<double>> varlims
    ) {
    
    for (int idx=0; idx<vars.size(); idx++) {

        // Get variable name and limits
        std::string var = vars[idx];
        double varmin   = varlims[idx][0];
        double varmax   = varlims[idx][1];

        // Add variable limit cuts to overall cuts
        if (cuts.size()>0) {
            cuts = Form("%s && %s>=%.8f && %s<%.8f",cuts.c_str(),var.c_str(),varmin,var.c_str(),varmax);
        } else {
            cuts = Form("%s>=%.8f && %s<%.8f",var.c_str(),varmin,var.c_str(),varmax);
        }

    } // for (int idx=0; idx<vars.size(); idx++) {
}

} // using namespace util;

} // namespace saga {
