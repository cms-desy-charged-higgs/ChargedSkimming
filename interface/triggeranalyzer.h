#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <ChargedAnalysis/Skimming/interface/baseanalyzer.h>

#include <numeric>

class TriggerAnalyzer : public BaseAnalyzer {
    private:
        //Trigger strings and vector with values
        std::vector<std::string> muPaths;
        std::vector<std::string> elePaths;
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> triggerEle;
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> triggerMu;

        //EDM Token for MINIAOD analysis
        trigToken triggerToken;
        
        //Vector with triger results
        std::vector<int> muResults;
        std::vector<int> eleResults;

    public:
        TriggerAnalyzer(const std::vector<std::string> &muPaths, const std::vector<std::string> &elePaths, TTreeReader &reader);
        TriggerAnalyzer(const std::vector<std::string> &muPaths, const std::vector<std::string> &elePaths, trigToken& triggerToken);
        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
