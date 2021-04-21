#ifndef BTagCSVReader_H
#define BTagCSVReader_H

#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>

#include <TF1.h>

template <typename T>
std::vector<T> SplitString(const std::string& splitString, const std::string& delimeter){
    T value;
    std::vector<T> values;

    std::size_t current, previous = 0;
    current = splitString.find(delimeter);

    while (current != std::string::npos){
        std::istringstream iss(splitString.substr(previous, current - previous));
        iss >> value;
        values.push_back(value);

        previous = current + 1;
        current = splitString.find(delimeter, previous);
    }

    std::istringstream iss(splitString.substr(previous, current - previous));
    iss >> value;
    values.push_back(value);

    return values;
}

class BTagCSVReader{
    private:        
        std::map<std::tuple<int, int, int>, std::vector<std::shared_ptr<TF1>>> SF;
        std::map<std::tuple<int, int, int>, std::vector<std::pair<float, float>>> ptRange;

    public: 
        BTagCSVReader(){};
        BTagCSVReader(const std::string& fileName);

        double Get(const float& pt, const int& wp, const int& flavour);
        double GetUp(const float& pt, const int& wp, const int& flavour);
        double GetDown(const float& pt, const int& wp, const int& flavour);
};

#endif
