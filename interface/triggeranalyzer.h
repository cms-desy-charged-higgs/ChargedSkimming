#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <ChargedHiggs/Skimming/interface/baseanalyzer.h>

#include <numeric>

class TriggerAnalyzer : public BaseAnalyzer {
    private:
        //Trigger strings and vector with values
        std::vector<std::string> triggerPaths;
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> triggerValues;

        //EDM Token for MINIAOD analysis
        trigToken triggerToken;
        
        //Vector with triger results
        std::vector<int> triggerResults;

    public:
        TriggerAnalyzer(const std::vector<std::string> &triggerPaths, TTreeReader &reader);
        TriggerAnalyzer(const std::vector<std::string> &triggerPaths, trigToken& triggerToken);
        void BeginJob(TTree *tree, bool &isData);
        bool Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
