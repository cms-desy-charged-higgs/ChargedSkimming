#include <ChargedSkimming/Skimming/interface/genpartanalyzer.h>

GenPartAnalyzer::GenPartAnalyzer(genPartToken& genParticleToken, const std::vector<std::string>& particleNames):
    BaseAnalyzer(),
    genParticleToken(genParticleToken),
    particleNames(particleNames)
    {}

GenPartAnalyzer::GenPartAnalyzer(TTreeReader& reader, const std::vector<std::string>& particleNames):
    BaseAnalyzer(&reader),
    particleNames(particleNames){}

void GenPartAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    //Set data bool
    this->isData = isData;
    
    //Set TTreeReader for genpart and trigger obj from baseanalyzer
    if(isNANO) SetCollection(this->isData);
    
    //Convert particle names to ID
    std::map<std::string, int> convert = {
        {"HPlus", 37}, 
        {"higgs", 25}, 
        {"top", 6}, 
        {"b", 5}, 
        {"W", 24},
        {"Z", 23}, 
    };

    for(std::string& name : particleNames) particleIDs.push_back(convert.at(name));

    floatVar = {
            {"Pt", Pt},
            {"Eta", Eta},
            {"Phi", Phi},
            {"Mass", Mass},
    };

    intVar = {
            {"ParticleID", partID},
            {"MotherID", mothID},
    };
    //Set Branches of output tree
    for(TTree* tree: trees){
        for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
            tree->Branch(("GenParticle_" + var.first).c_str(), &var.second);
        }

        for(std::pair<const std::string, std::vector<char>&>& var: intVar){
            tree->Branch(("GenParticle_" + var.first).c_str(), &var.second);
        }
    }
}

void GenPartAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
        var.second.clear();
    }

    for(std::pair<const std::string, std::vector<char>&>& var: intVar){
        var.second.clear();
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    if(!isData){
        if(!isNANO){
            event->getByToken(genParticleToken, genParts);
        }
    
        std::vector<int> alreadySeenNANO;
        std::vector<const reco::Candidate*> alreadySeenMINI;
        int size = isNANO ? genID->GetSize() : genParts->size();

        //Fill 4 four vectors
        for(int i = 0; i < size; i++){
            int ID = isNANO ? abs(genID->At(i)) : abs(genParts->at(i).pdgId());  

            if(std::find(particleIDs.begin(), particleIDs.end(), ID) != particleIDs.end()){
                const reco::Candidate* part=nullptr;
                int index=0;

                if(isNANO) index = FirstCopy(i, ID);
                else part = FirstCopy(genParts->at(i), ID);

                if(isNANO){
                    if(std::find(alreadySeenNANO.begin(), alreadySeenNANO.end(), index) != alreadySeenNANO.end()) continue; 
                    else alreadySeenNANO.push_back(index);
                }

                else{
                    if(std::find(alreadySeenMINI.begin(), alreadySeenMINI.end(), part) != alreadySeenMINI.end()) continue; 
                    else alreadySeenMINI.push_back(part);
                }

                int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(part->mother()->pdgId()); 

                float pt, eta, phi, m;
                pt = isNANO ? genPt->At(index) : part->pt();
                phi = isNANO ? genPhi->At(index) : part->phi();
                eta = isNANO ? genEta->At(index) : part->eta();
                m = isNANO ? genMass->At(index) : part->mass();

                Pt.push_back(pt);
                Eta.push_back(eta);
                Phi.push_back(phi);
                Mass.push_back(m);
            
                partID.push_back(ID);
                mothID.push_back(motherID);
            }
        }
    }
}

void GenPartAnalyzer::EndJob(TFile* file){
}
