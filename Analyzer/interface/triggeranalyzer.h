#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class TriggerAnalyzer : public BaseAnalyzer<T> {
    public:
        TriggerAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){}

        void Analyze(T& input, Output& out){
            input.ReadTrigger();
            input.GetTrigger();
        
            for(std::size_t i = 0; i < input.triggers.size(); ++i){
                out.triggers[i] = input.triggers[i];
            }
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
