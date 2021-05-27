#include <ChargedSkimming/Skimming/interface/triggeranalyzer.h>

TriggerAnalyzer::TriggerAnalyzer(const std::string& era, const std::vector<std::string>& channels, TTreeReader& reader):
    BaseAnalyzer(&reader),
    era(era),
    channels(channels){}

TriggerAnalyzer::TriggerAnalyzer(const std::string& era, const std::vector<std::string>& channels, const std::shared_ptr<Token>& tokens):
    BaseAnalyzer(),
    era(era),
    channels(channels),
    tokens(tokens){}

void TriggerAnalyzer::BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf){
    for(const std::string& channel: channels){
        triggerPaths[channel] = Util::GetVector<std::string>(skim, "Channel." + channel + ".Trigger." + era);
    }

    if(isNANO){
        //TTreeReader Values
        for(const std::pair<std::string, std::vector<std::string>>& t : triggerPaths){
            for(const std::string& triggerName : t.second){
                trigger[t.first].push_back(std::make_shared<TTreeReaderValue<bool>>(*reader, triggerName.c_str()));
            }
        }
    }

    for(const std::pair<std::string, std::vector<std::string>>& t : triggerPaths){
        triggerResults[t.first] = std::vector<int>(t.second.size(), 0);
    }

    for(TTree* tree: trees){
        std::string treeName(tree->GetName());

        for(const std::pair<std::string, std::vector<std::string>>& t : triggerPaths){
            if(treeName == t.first){
                for(unsigned int i = 0; i < t.second.size(); i++){
                    tree->Branch(t.second[i].c_str(), &triggerResults[t.first][i]);
                }
            }
        }
    }
}

void TriggerAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear result vector
    for(const std::pair<std::string, std::vector<std::string>>& t : triggerPaths){
        for(unsigned int i = 0; i < triggerResults[t.first].size(); i++){
            triggerResults[t.first][i] = 0;
        }
    }

    //Get Event info is using MINIAOD
    edm::TriggerResults triggers;

    if(!isNANO){
        triggers = Token::GetTokenValue<edm::TriggerResults>(event, tokens->triggerToken);
        const edm::TriggerNames &names = event->triggerNames(triggers);

        //Find position of trigger name in trigger collection to avoid looping everytime over all trigger names
        if(triggerIndex.empty()){
            for(unsigned int i = 0; i < names.size(); i++){
                for(const std::pair<std::string, std::vector<std::string>>& t : triggerPaths){
                    for(unsigned int j = 0; j < t.second.size(); j++){
                        if(names.triggerName(i).find(t.second[j] + "_v") != std::string::npos){
                            triggerIndex[t.first].push_back(i);
                        }
                    }
                }
            }
        }

        //Check trigger result
        for(const std::pair<std::string, std::vector<std::string>>& t : triggerPaths){
            for(unsigned int i = 0; i < t.second.size(); i++){
                triggerResults[t.first][i] = triggers.accept(triggerIndex[t.first][i]);
            }
        }
    }

    else{
        //Trigger result in NANOAOD
        for(const std::pair<std::string, std::vector<std::shared_ptr<TTreeReaderValue<bool>>>>& t : trigger){
            for(unsigned int i = 0; i < t.second.size(); i++){
                triggerResults[t.first][i] = *t.second[i]->Get();
            }
        }
    }

    for(CutFlow& cutflow: cutflows){
        if(triggerResults.count(cutflow.channel)){
            bool passed = false; 
            for(int& value : triggerResults[cutflow.channel]) passed = passed or value;
    
            if(passed){
                std::string cutName = triggerPaths[cutflow.channel][0];
                cutflow.hist->Fill(cutName.c_str(), cutflow.weight);                   
            }

            else cutflow.passed = false;
        }
    }
}

void TriggerAnalyzer::EndJob(TFile* file){}
