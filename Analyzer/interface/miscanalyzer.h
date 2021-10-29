#ifndef MISCANALYZER_H
#define MISCANALYZER_H

#include <TParameter.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class MiscAnalyzer : public BaseAnalyzer<T> {
    private: 
        bool isData;

    public:
        MiscAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            isData = skim.get<std::string>("run") != "MC";
        }

        void Analyze(T& input, Output& out){
            input.ReadMiscEntry(isData);
            input.GetMisc();

            out.evNr = input.evNr;
            out.nParton = input.nParton;
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){
        };
};

#endif
