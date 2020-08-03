#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <numeric>

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

class TriggerAnalyzer : public BaseAnalyzer {
    private:
        //Trigger strings and vector with values
        std::map<std::string, std::vector<std::string>> triggerPaths;
        std::map<std::string, std::vector<std::shared_ptr<TTreeReaderValue<bool>>>> trigger;

        //EDM Token for MINIAOD analysis
        trigToken triggerToken;
        std::map<std::string, std::vector<int>> triggerIndex;
        
        //Vector with triger results
        std::map<std::string, std::vector<int>> triggerResults;

    public:
        TriggerAnalyzer(const std::map<std::string, std::vector<std::string>> triggerPaths, TTreeReader &reader);
        TriggerAnalyzer(const std::map<std::string, std::vector<std::string>> triggerPaths, trigToken& triggerToken);
        void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
