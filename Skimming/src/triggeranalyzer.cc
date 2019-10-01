#include <ChargedSkimming/Skimming/interface/triggeranalyzer.h>

TriggerAnalyzer::TriggerAnalyzer(const std::vector<std::string> &muPaths, const std::vector<std::string> &elePaths, TTreeReader& reader):
    BaseAnalyzer(&reader),
    muPaths(muPaths),
    elePaths(elePaths){}

TriggerAnalyzer::TriggerAnalyzer(const std::vector<std::string> &muPaths, const std::vector<std::string> &elePaths, trigToken& triggerToken):
    BaseAnalyzer(),
    muPaths(muPaths),
    elePaths(elePaths),
    triggerToken(triggerToken){}

void TriggerAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData){
    if(isNANO){
        //TTreeReader Values
        for(std::string triggerPath: muPaths){
            triggerMu.push_back(std::make_unique<TTreeReaderValue<bool>>(*reader, triggerPath.c_str()));
        }

        for(std::string triggerPath: elePaths){
            triggerEle.push_back(std::make_unique<TTreeReaderValue<bool>>(*reader, triggerPath.c_str()));
        }
    }

    eleResults = std::vector<int>(elePaths.size(), 0);
    muResults = std::vector<int>(muPaths.size(), 0);

    for(TTree* tree: trees){
        std::string treeName(tree->GetName());

        if(treeName.find("mu") != std::string::npos){
            //Set Branches of output tree
            for(unsigned int i=0; i < muPaths.size(); i++){
                tree->Branch(muPaths[i].c_str(), &muResults[i]);
            }
        }

        if(treeName.find("e") != std::string::npos){
            //Set Branches of output tree
            for(unsigned int i=0; i < elePaths.size(); i++){
                tree->Branch(elePaths[i].c_str(), &eleResults[i]);
            }
        }
    }
}

void TriggerAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear result vector
    eleResults.clear();
    muResults.clear();

    //Get Event info is using MINIAOD
    edm::Handle<edm::TriggerResults> triggers;

    if(!isNANO){
        event->getByToken(triggerToken, triggers);
        const edm::TriggerNames &names = event->triggerNames(*triggers);

        for(unsigned int i = 0; i < names.size(); i++){
            //Find trigger result with given trigger paths with MINIAOD
            for(std::string& pathName: muPaths){
                if(names.triggerName(i).find(pathName + "_v") != std::string::npos){
                    muResults.push_back(triggers->accept(i));
                }
            }

            for(std::string& pathName: elePaths){
                if(names.triggerName(i).find(pathName + "_v") != std::string::npos){
                    eleResults.push_back(triggers->accept(i));
                }
            }
        }   
    }

    else{
        //Trigger result in NANOAOD
        for(std::unique_ptr<TTreeReaderValue<bool>> &triggerValue: triggerEle){
            eleResults.push_back(*triggerValue->Get());
        }

        for(std::unique_ptr<TTreeReaderValue<bool>> &triggerValue: triggerMu){
            muResults.push_back(*triggerValue->Get());
        }
    }

    for(CutFlow& cutflow: cutflows){
        if(cutflow.nMinMu>=1){
            if(std::find(muResults.begin(), muResults.end(), 1) != muResults.end()){
                if(cutflow.passed){
                    std::string cutName = muPaths[0];
                    cutflow.hist->Fill(cutName.c_str(), cutflow.weight);                
                }    
            }

            else cutflow.passed = false;
        }

        if(cutflow.nMinEle>=1){
            if(std::find(eleResults.begin(), eleResults.end(), 1) != eleResults.end()){
                if(cutflow.passed){
                    std::string cutName = elePaths[0];
                    cutflow.hist->Fill(cutName.c_str(), cutflow.weight);        
                }
            }

            else cutflow.passed = false;
        }
    }
}

void TriggerAnalyzer::EndJob(TFile* file){}
