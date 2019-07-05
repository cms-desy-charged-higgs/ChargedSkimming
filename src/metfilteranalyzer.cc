#include <ChargedHiggs/Skimming/interface/metfilteranalyzer.h>

MetFilterAnalyzer::MetFilterAnalyzer(const int &era, TTreeReader &reader):
    BaseAnalyzer(&reader),
    era(era){}

void MetFilterAnalyzer::BeginJob(TTree *tree, bool &isData){
    //Set Filter names for each era
    filterNames = {
                {2017, {"Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        "Flag_eeBadScFilter",
                        }
                 },
    };

    //Set TTreeReaderValues
    for(std::string filterName: filterNames[era]){
        filterValues.push_back(std::make_unique<TTreeReaderValue<bool>>(*reader, filterName.c_str()));
    }

}

bool MetFilterAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow){
    bool passedFilter = true;

    for(std::unique_ptr<TTreeReaderValue<bool>> &filter: filterValues){
        passedFilter *= *filter->Get();
    }

    if(passedFilter){
        cutflow.first->Fill("Met filter", cutflow.second);
       return true;
    }

    return false;
}

void MetFilterAnalyzer::EndJob(TFile* file){}
