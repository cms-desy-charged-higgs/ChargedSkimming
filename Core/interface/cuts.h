#ifndef CUTS_H
#define CUTS_H

#include <vector>
#include <string>
#include <functional>
#include <memory>		

#include <TFile.h>
#include <TH1F.h>

#include <ChargedSkimming/Core/interface/output.h>

class Cuts{
    private:
        std::shared_ptr<TH1F> cutFlow;
        std::vector<std::function<bool()>> cuts;
        std::vector<std::string> cutNames;

        std::function<bool()> ConstructCut(short& value, const std::string& op, const short& threshold);

    public:
        Cuts(){}
        Cuts(const std::shared_ptr<TFile>& outFile, const std::string& channel);
        void AddCut(const std::string& part, Output& out, const std::string& op, const short& threshold);

        template <typename T>
        void AddTrigger(const std::vector<int>& triggerIdx, T& input){
            if(triggerIdx.size() != 0){
                cuts.insert(cuts.begin(), [&input, triggerIdx](){for(const int& idx : triggerIdx){if(input.triggers[idx]) return true;} return false;});
                cutNames.insert(cutNames.begin(), "Trigger");
            }

            else{
                cuts.insert(cuts.begin(), [&input](){for(const bool& passed : input.METFilter){if(!passed) return false;} return true;});
                cutNames.insert(cutNames.begin(), "MET Filter");
            }
        }

        bool Passed();
        void Count(){cutFlow->Fill("No cuts", 1);};
        void FillCutflow();
        void WriteOutput(){cutFlow->Write();};
};

#endif
