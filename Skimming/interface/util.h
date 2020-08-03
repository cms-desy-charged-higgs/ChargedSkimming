#ifndef UTIL_H
#define UTIL_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Util{
    template <typename T>
    std::vector<T> GetVector(boost::property_tree::ptree& tree, const std::string& key){
        std::vector<T> vec;
    
        for(std::pair<std::string, boost::property_tree::ptree> values : tree.get_child(key)){
            vec.push_back(values.second.get_value<T>());
        }
    
        return vec;
    };
};

#endif
