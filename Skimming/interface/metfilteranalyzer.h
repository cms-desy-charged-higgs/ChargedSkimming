#ifndef METFILTERANALYZER_H
#define METFILTERANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/util.h>

#include <numeric>

class MetFilterAnalyzer : public BaseAnalyzer {
    private:
        //Era
        int era;

        //Vector with filter names
        std::vector<std::string> filterNames;

        //EDM Token for MINIAOD analysis
        trigToken triggerToken;

        //Vector with TTreeReaderValues
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> filterValues;

    public:
        MetFilterAnalyzer(const int &era, TTreeReader &reader);
        MetFilterAnalyzer(const int &era, trigToken& triggerToken);
        void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
