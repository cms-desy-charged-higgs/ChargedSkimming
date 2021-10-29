#ifndef BASEANALYZER_H
#define BASEANALYZER_H

#include <memory>

#include <TFile.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <ChargedSkimming/Core/interface/output.h>

namespace pt = boost::property_tree;

/// Virtual base class for analyzer, templated for nano/mini AOD input

template <typename T>
class BaseAnalyzer {
    protected:
        std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/";

    public:
        virtual ~BaseAnalyzer() = default;

        virtual void BeginJob(const pt::ptree& skim, const pt::ptree& sf) = 0;
        virtual void Analyze(T& input, Output& out) = 0;
        virtual void EndJob(const std::shared_ptr<TFile>& outFile) = 0;
};

#endif
