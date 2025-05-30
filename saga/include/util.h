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
        std::string       argname,
        T                 defaultval,
        std::string       message_prefix,
        bool              verbose,
        std::ostream     &out = std::cout
    ){
    T arg = defaultval;
    if (node[argname]) {

        // Load argument
        arg = node[argname].as<T>();

        // Print argument
        if (verbose) {

            // Single value numeric
            if constexpr (std::is_same<T, bool>::value || std::is_same<T, int>::value || std::is_same<T, float>::value || std::is_same<T, double>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": " << arg << std::endl;
            }
            // Single value string
            else if constexpr (std::is_same<T, std::string>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": " << arg.c_str() << std::endl;
            }
            // Vector numeric
            else if constexpr (std::is_same<T, std::vector<bool>>::value || std::is_same<T, std::vector<int>>::value || std::is_same<T, std::vector<float>>::value || std::is_same<T, std::vector<double>>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": [ ";
                for (int idx=0; idx<arg.size(); idx++) {
                    if (idx!=arg.size()-1) { out << arg[idx] <<", "; }
                    else { out << arg[idx]; }
                }
                out << " ]" << std::endl;
            }
            // Vector string
            else if constexpr (std::is_same<T, std::vector<std::string>>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": [ ";
                for (int idx=0; idx<arg.size(); idx++) {
                    if (idx!=arg.size()-1) { out << arg[idx].c_str() <<", "; }
                    else { out << arg[idx].c_str(); }
                }
                out << " ]" << std::endl;
            }
            // Vector vector numeric
            else if constexpr (std::is_same<T, std::vector<std::vector<bool>>>::value || std::is_same<T, std::vector<std::vector<int>>>::value || std::is_same<T, std::vector<std::vector<float>>>::value || std::is_same<T, std::vector<std::vector<double>>>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": [ ";
                for (int idx=0; idx<arg.size(); idx++) {
                    out << "\n\t [ ";
                    for (int idx2=0; idx2<arg[idx].size(); idx2++) {
                        if (idx2!=arg[idx].size()-1) { out << arg[idx][idx2] <<", "; }
                        else { out << arg[idx][idx2]; }
                    }
                    if (idx!=arg.size()-1) { out <<"],"; }
                    else { out <<"]\n"; }
                }
                out << " ]" << std::endl;
            }
            // Vector vector string
            else if constexpr (std::is_same<T, std::vector<std::vector<std::string>>>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": [ ";
                for (int idx=0; idx<arg.size(); idx++) {
                    out << "\n\t [ ";
                    for (int idx2=0; idx2<arg[idx].size(); idx2++) {
                        if (idx2!=arg[idx].size()-1) { out << arg[idx][idx2].c_str() <<", "; }
                        else { out << arg[idx][idx2].c_str(); }
                    }
                    if (idx!=arg.size()-1) { out <<"],"; }
                    else { out <<"]\n"; }
                }
                out << " ]" << std::endl;
            }
            // Map string numeric
            else if constexpr (std::is_same<T, std::map<std::string,bool>>::value || std::is_same<T, std::map<std::string,int>>::value || std::is_same<T, std::map<std::string,float>>::value || std::is_same<T, std::map<std::string,double>>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": { ";
                for (auto it = arg.begin(); it != arg.end(); ++it) {
                    out << it->first<<" : "<<it->second<<", ";
                }
                out << " }" << std::endl;
            }
            // Map string string
            else if constexpr (std::is_same<T, std::map<std::string,std::string>>::value) {
                out << message_prefix.c_str() << argname.c_str() << ": { ";
                for (auto it = arg.begin(); it != arg.end(); ++it) {
                    out << it->first.c_str()<<" : "<<it->second.c_str()<<", ";
                }
                out << " }" << std::endl;
            }
            else { out << message_prefix.c_str() << argname.c_str() << ": STRING CONVERSION NOT IMPLEMENTED FOR ARGUMENT TYPE" << std::endl; }
        } // if (verbose) {

    } // if (node[argname]) {
    return arg;

}

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
*
* @param cuts Base cut string
* @param vars Variables for which to add limits cuts
* @param varlims Limits of provided variables
*
* @return Updated cut string
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
