#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <numeric>

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

class TriggerAnalyzer : public BaseAnalyzer {
    private:
        //Trigger strings and vector with values
        std::vector<std::string> muPaths, elePaths;
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> triggerEle, triggerMu;

        //EDM Token for MINIAOD analysis
        trigToken triggerToken;
        
        //Vector with triger results
        std::vector<int> muResults, eleResults;

    public:
        TriggerAnalyzer(const std::vector<std::string> &muPaths, const std::vector<std::string> &elePaths, TTreeReader &reader);
        TriggerAnalyzer(const std::vector<std::string> &muPaths, const std::vector<std::string> &elePaths, trigToken& triggerToken);
        void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
