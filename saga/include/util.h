#include <iostream>
#include <fstream>
#include <memory>
#include <string>
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
* @brief Add additional variable limit cuts to an existing cut string.
*
* Add additional variable limit cuts to an existing cut string.
* Note that this is an in place operation.
*
* @param cuts Substring to insert
* @param vars String to search
* @param varlims Substring to find and replace
*/
std::string addLimitCuts(
        std::string                      cuts,
        std::vector<std::string>         vars,
        std::vector<std::vector<double>> varlims
    ) {
    
    std::string newcuts = cuts;
    for (int idx=0; idx<vars.size(); idx++) {

        // Get variable name and limits
        std::string var = vars[idx];
        double varmin   = varlims[idx][0];
        double varmax   = varlims[idx][1];

        // Add variable limit cuts to overall cuts
        if (newcuts.size()>0) {
            newcuts = Form("%s && %s>=%.8f && %s<%.8f",newcuts.c_str(),var.c_str(),varmin,var.c_str(),varmax);
        } else {
            newcuts = Form("%s>=%.8f && %s<%.8f",var.c_str(),varmin,var.c_str(),varmax);
        }

    } // for (int idx=0; idx<vars.size(); idx++) {

    return newcuts;
}

} // using namespace util;

} // namespace saga {
