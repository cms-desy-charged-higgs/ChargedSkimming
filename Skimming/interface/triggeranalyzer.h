#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <numeric>

#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/tokens.h>

class TriggerAnalyzer : public BaseAnalyzer {
    private:
        std::string era;
        std::vector<std::string> channels;
        std::shared_ptr<Token> tokens;

        //Trigger strings and vector with values
        std::map<std::string, std::vector<std::string>> triggerPaths;
        std::map<std::string, std::vector<std::shared_ptr<TTreeReaderValue<bool>>>> trigger;

        //EDM std::shared_ptr<Token> for MINIAOD analysis
        std::map<std::string, std::vector<int>> triggerIndex;
        
        //Vector with triger results
        std::map<std::string, std::vector<int>> triggerResults;

    public:
        TriggerAnalyzer(const std::string& era, const std::vector<std::string>& channels, TTreeReader& reader);
        TriggerAnalyzer(const std::string& era, const std::vector<std::string>& channels, const std::shared_ptr<Token>& tokens);

        void BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
