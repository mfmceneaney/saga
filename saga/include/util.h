#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <TString.h>

#include <yaml-cpp/yaml.h>

#pragma once

/**
* @file
* @author Matthew F. McEneaney
* @date 7/Feb./2025
* @version 0.0.0
* @brief Utility functions for asymmetry fits.
*/

namespace saga {

namespace util {

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
* @brief Load an argument from YAML file.
*
* @tparam T Type of argument to get

* @param node Yaml node from which to load argument
* @param argname Name of argument to load
* @param defaultval Default value to use if argument is not found
* @param message_prefix Prefix to output message
* @param verbose Option to print out argument name and value
* @param out Output stream to print to
*
* @return Argument of type T
*/
template<typename T>
T getYamlArg(
        const YAML::Node &node,
        string       argname,
        T                 defaultval,
        string       message_prefix,
        bool              verbose,
        ostream     &out = cout
    ){
    T arg = defaultval;
    if (node[argname]) {

        // Load argument
        arg = node[argname].as<T>();

    } // if (node[argname]) {

    // Print argument
    if (verbose) {

        // First, set prefix
        out << message_prefix.c_str();
        if (!node[argname]) { out << "USE DEFAULT: "; }

        // Single value numeric
        if constexpr (is_same<T, bool>::value || is_same<T, int>::value || is_same<T, float>::value || is_same<T, double>::value) {
            out << argname.c_str() << ": " << arg << endl;
        }
        // Single value string
        else if constexpr (is_same<T, string>::value) {
            out << argname.c_str() << ": " << arg.c_str() << endl;
        }
        // Vector numeric
        else if constexpr (is_same<T, vector<bool>>::value || is_same<T, vector<int>>::value || is_same<T, vector<float>>::value || is_same<T, vector<double>>::value) {
            out << argname.c_str() << ": [ ";
            for (int idx=0; idx<arg.size(); idx++) {
                if (idx!=arg.size()-1) { out << arg[idx] <<", "; }
                else { out << arg[idx]; }
            }
            out << " ]" << endl;
        }
        // Vector string
        else if constexpr (is_same<T, vector<string>>::value) {
            out << argname.c_str() << ": [ ";
            for (int idx=0; idx<arg.size(); idx++) {
                if (idx!=arg.size()-1) { out << arg[idx].c_str() <<", "; }
                else { out << arg[idx].c_str(); }
            }
            out << " ]" << endl;
        }
        // Vector vector numeric
        else if constexpr (is_same<T, vector<vector<bool>>>::value || is_same<T, vector<vector<int>>>::value || is_same<T, vector<vector<float>>>::value || is_same<T, vector<vector<double>>>::value) {
            out << argname.c_str() << ": [ ";
            for (int idx=0; idx<arg.size(); idx++) {
                out << "\n\t [ ";
                for (int idx2=0; idx2<arg[idx].size(); idx2++) {
                    if (idx2!=arg[idx].size()-1) { out << arg[idx][idx2] <<", "; }
                    else { out << arg[idx][idx2]; }
                }
                if (idx!=arg.size()-1) { out <<"],"; }
                else { out <<"]\n"; }
            }
            out << " ]" << endl;
        }
        // Vector vector string
        else if constexpr (is_same<T, vector<vector<string>>>::value) {
            out << argname.c_str() << ": [ ";
            for (int idx=0; idx<arg.size(); idx++) {
                out << "\n\t [ ";
                for (int idx2=0; idx2<arg[idx].size(); idx2++) {
                    if (idx2!=arg[idx].size()-1) { out << arg[idx][idx2].c_str() <<", "; }
                    else { out << arg[idx][idx2].c_str(); }
                }
                if (idx!=arg.size()-1) { out <<"],"; }
                else { out <<"]\n"; }
            }
            out << " ]" << endl;
        }
        // Map string numeric
        else if constexpr (is_same<T, map<string,bool>>::value || is_same<T, map<string,int>>::value || is_same<T, map<string,float>>::value || is_same<T, map<string,double>>::value) {
            out << argname.c_str() << ": { ";
            for (auto it = arg.begin(); it != arg.end(); ++it) {
                out << it->first<<" : "<<it->second<<", ";
            }
            out << " }" << endl;
        }
        // Map string string
        else if constexpr (is_same<T, map<string,string>>::value) {
            out << argname.c_str() << ": { ";
            for (auto it = arg.begin(); it != arg.end(); ++it) {
                out << it->first.c_str()<<" : "<<it->second.c_str()<<", ";
            }
            out << " }" << endl;
        }
        else { out << argname.c_str() << ": STRING CONVERSION NOT IMPLEMENTED FOR ARGUMENT TYPE" << endl; }
    } // if (verbose) {

    return arg;

}

/**
* @brief Find and replace function for string.
*
* Find all occurences of a substring and replace them with another.
* Note that this is an in place operation. REGEX is NOT supported.
*
* @param s String to search
* @param to_replace Substring to find and replace
* @param replace_with Substring to insert
*/
void replaceAll(
        string& s,
        string const& to_replace,
        string const& replace_with
    ) {
    string buf;
    size_t pos = 0;
    size_t prev_pos;

    // Reserves rough estimate of final size of string.
    buf.reserve(s.size());

    while (true) {
        prev_pos = pos;
        pos = s.find(to_replace, pos);
        if (pos == string::npos)
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
*
* @param cuts Base cut string
* @param vars Variables for which to add limits cuts
* @param varlims Limits of provided variables
*
* @return Updated cut string
*/
string addLimitCuts(
        string                      cuts,
        vector<string>         vars,
        vector<vector<double>> varlims
    ) {
    
    string newcuts = cuts;
    for (int idx=0; idx<vars.size(); idx++) {

        // Get variable name and limits
        string var = vars[idx];
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
