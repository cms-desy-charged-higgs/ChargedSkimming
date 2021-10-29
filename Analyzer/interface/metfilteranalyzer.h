#ifndef METFILTERANALYZER_H
#define METFILTERANALYZER_H

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class METFilterAnalyzer : public BaseAnalyzer<T> {
    public:
        METFilterAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){}

        void Analyze(T& input, Output& out){
            input.ReadMETFilter();
            input.GetMETFilter();
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
