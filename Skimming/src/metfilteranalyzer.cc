#include <ChargedSkimming/Skimming/interface/metfilteranalyzer.h>

MetFilterAnalyzer::MetFilterAnalyzer(const std::string& era, TTreeReader& reader):
    BaseAnalyzer(&reader),
    era(era){}

MetFilterAnalyzer::MetFilterAnalyzer(const std::string& era, const std::shared_ptr<Token>& tokens):
    BaseAnalyzer(),
    era(era),
    tokens(tokens)
    {}

void MetFilterAnalyzer::BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf){
    //Set Filter names for each era
    filterNames = Util::GetVector<std::string>(skim, "Analyzer.METFilter." + era);

    if(isNANO){
        //Set TTreeReaderValues
        for(std::string filterName: filterNames){
            filterValues.push_back(std::make_unique<TTreeReaderValue<bool>>(*reader, filterName.c_str()));
        }
    }
}

void MetFilterAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    bool passedFilter = true;

    //Get Event info is using MINIAOD
    edm::TriggerResults triggers;

    if(!isNANO){
        triggers = Token::GetTokenValue(event, tokens->triggerToken);
        const edm::TriggerNames &names = event->triggerNames(triggers);

        //Find result with given filter name with MINIAOD
        for(std::string filterName: filterNames){
            std::vector<int> filterVersion; 

            for(unsigned int i = 0; i < names.size(); i++){
                if(names.triggerName(i).find(filterName) != std::string::npos){
                    filterVersion.push_back(triggers.accept(i));
                }
            }    
    
            if(filterVersion.size() != 0){
                passedFilter *= filterVersion.back();
            }    
        }
    }

    else{
        //Filter result in NANOAOD
        for(std::unique_ptr<TTreeReaderValue<bool>> &filter: filterValues){
            passedFilter *= *filter->Get();
        }
    }

    if(passedFilter){
        for(CutFlow& cutflow: cutflows){
            if(cutflow.passed){
                cutflow.hist->Fill("Met filter", cutflow.weight);        
            }
        }
    }

    else{
        for(CutFlow& cutflow: cutflows){     
            cutflow.passed = false;
        }
    }
}

void MetFilterAnalyzer::EndJob(TFile* file){}
