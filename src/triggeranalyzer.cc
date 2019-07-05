#include <ChargedHiggs/Skimming/interface/triggeranalyzer.h>

TriggerAnalyzer::TriggerAnalyzer(const std::vector<std::string> &triggerPaths, TTreeReader& reader):
    BaseAnalyzer(&reader),
    triggerPaths(triggerPaths){}

void TriggerAnalyzer::BeginJob(TTree *tree, bool &isData){
    //TTreeReader Values
    for(std::string triggerPath: triggerPaths){
        triggerValues.push_back(std::make_unique<TTreeReaderValue<bool>>(*reader, triggerPath.c_str()));
    }

    triggerResults = std::vector<int>(triggerPaths.size(), 0);

    //Set Branches of output tree
    for(unsigned int i=0; i < triggerPaths.size(); i++){
        tree->Branch(triggerPaths[i].c_str(), &triggerResults[i]);
    }
}

bool TriggerAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow){
    //Clear result vector
    triggerResults.clear();

    for(std::unique_ptr<TTreeReaderValue<bool>> &triggerValue: triggerValues){
        triggerResults.push_back(*triggerValue->Get());
    }

    if(std::find(triggerResults.begin(), triggerResults.end(), 1) != triggerResults.end()){
        std::string cutName = triggerPaths[0];
 
        cutflow.first->Fill(cutName.c_str(), cutflow.second);        
        return true;
    }

    return false;
}

void TriggerAnalyzer::EndJob(TFile* file){}
