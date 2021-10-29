#ifndef UTIL_H
#define UTIL_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Util{
    template <typename T>
    std::vector<T> GetVector(const boost::property_tree::ptree& tree, const std::string& key){
        std::vector<T> vec;
    
        for(std::pair<std::string, boost::property_tree::ptree> values : tree.get_child(key)){
            vec.push_back(values.second.get_value<T>());
        }
    
        return vec;
    };

    std::vector<std::string> GetKeys(const pt::ptree& tree, const std::string path){
        std::vector<std::string> keys;
        pt::ptree node = tree.get_child(path);

        for(const std::pair<const std::string, pt::ptree>& p : node){
            keys.push_back(p.first);
        }

        return keys;
    }

    float DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2){
        return std::sqrt(std::pow(eta1 - eta2, 2) + std::pow(phi1 - phi2, 2));
    }
};

#endif
