#ifndef METFILTERANALYZER_H
#define METFILTERANALYZER_H

#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/util.h>

#include <numeric>

class MetFilterAnalyzer : public BaseAnalyzer {
    private:
        //Era
        std::string era;
        std::shared_ptr<Token> tokens;

        //Vector with filter names
        std::vector<std::string> filterNames;

        //Vector with TTreeReaderValues
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> filterValues;

    public:
        MetFilterAnalyzer(const std::string& era, TTreeReader& reader);
        MetFilterAnalyzer(const std::string& era, const std::shared_ptr<Token>& tokens);

        void BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
