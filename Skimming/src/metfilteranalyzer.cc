#include <ChargedSkimming/Skimming/interface/metfilteranalyzer.h>

MetFilterAnalyzer::MetFilterAnalyzer(const int &era, TTreeReader &reader):
    BaseAnalyzer(&reader),
    era(era){}

MetFilterAnalyzer::MetFilterAnalyzer(const int &era, trigToken& triggerToken):
    BaseAnalyzer(),
    era(era),
    triggerToken(triggerToken)
    {}

void MetFilterAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    //Read in json config with sf files
    boost::property_tree::ptree sf; 
    boost::property_tree::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/config/skim.json", sf);

    //Set Filter names for each era
    filterNames = Util::GetVector<std::string>(sf, "Analyzer.METFilter." + std::to_string(era));

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
    edm::Handle<edm::TriggerResults> triggers;

    if(!isNANO){
        event->getByToken(triggerToken, triggers);
        const edm::TriggerNames &names = event->triggerNames(*triggers);

        //Find result with given filter name with MINIAOD
        for(std::string filterName: filterNames){
            std::vector<int> filterVersion; 

            for(unsigned int i = 0; i < names.size(); i++){
                if(names.triggerName(i).find(filterName) != std::string::npos){
                    filterVersion.push_back(triggers->accept(i));
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
