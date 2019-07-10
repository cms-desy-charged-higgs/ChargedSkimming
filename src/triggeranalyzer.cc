#include <ChargedHiggs/Skimming/interface/triggeranalyzer.h>

TriggerAnalyzer::TriggerAnalyzer(const std::vector<std::string> &triggerPaths, TTreeReader& reader):
    BaseAnalyzer(&reader),
    triggerPaths(triggerPaths){}

TriggerAnalyzer::TriggerAnalyzer(const std::vector<std::string> &triggerPaths, trigToken& triggerToken):
    BaseAnalyzer(),
    triggerPaths(triggerPaths),
    triggerToken(triggerToken){}

void TriggerAnalyzer::BeginJob(TTree *tree, bool &isData){
    if(isNANO){
        //TTreeReader Values
        for(std::string triggerPath: triggerPaths){
            triggerValues.push_back(std::make_unique<TTreeReaderValue<bool>>(*reader, triggerPath.c_str()));
        }
    }

    triggerResults = std::vector<int>(triggerPaths.size(), 0);

    //Set Branches of output tree
    for(unsigned int i=0; i < triggerPaths.size(); i++){
        tree->Branch(triggerPaths[i].c_str(), &triggerResults[i]);
    }
}

bool TriggerAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event){
    //Clear result vector
    triggerResults.clear();

    //Get Event info is using MINIAOD
    edm::Handle<edm::TriggerResults> triggers;

    if(!isNANO){
        event->getByToken(triggerToken, triggers);
        const edm::TriggerNames &names = event->triggerNames(*triggers);

        //Find trigger result with given trigger paths with MINIAOD
        for(std::string pathName: triggerPaths){
            std::vector<int> triggerVersion; 

            for(unsigned int i = 0; i < names.size(); i++){
                if(names.triggerName(i).find(pathName) != std::string::npos){
                    triggerVersion.push_back(triggers->accept(i));
                }
            }    
    
            if(triggerVersion.size() != 0){
                triggerResults.push_back(triggerVersion.back());
            }    
        }
    }

    else{
        //Trigger result in NANOAOD
        for(std::unique_ptr<TTreeReaderValue<bool>> &triggerValue: triggerValues){
            triggerResults.push_back(*triggerValue->Get());
        }
    }

    if(std::find(triggerResults.begin(), triggerResults.end(), 1) != triggerResults.end()){
        std::string cutName = triggerPaths[0];
 
        cutflow.first->Fill(cutName.c_str(), cutflow.second);        
        return true;
    }

    return false;
}

void TriggerAnalyzer::EndJob(TFile* file){}
