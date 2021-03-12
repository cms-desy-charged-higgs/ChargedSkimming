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
        {"e", 11},
        {"nu_e", 12},
        {"muon", 13},
        {"nu_muon", 14},
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
        tree->Branch("GenParticle_Size", &nGenPart, "GenParticle_Size/S");

        for(std::pair<const std::string, float*>& var: floatVar){
            tree->Branch(("GenParticle_" + var.first).c_str(), var.second, ("GenParticle_" + var.first + "[GenParticle_Size]/F").c_str());
        }

        for(std::pair<const std::string, short*>& var: intVar){
            tree->Branch(("GenParticle_" + var.first).c_str(), var.second, ("GenParticle_" + var.first + "[GenParticle_Size]/S").c_str());
        }
    }
}

void GenPartAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Get Event info is using MINIAOD
    edm::Handle<std::vector<reco::GenParticle>> genParts;
    nGenPart = 0;

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
                if(motherID > 100) continue;

                float pt, eta, phi, m;
                pt = isNANO ? genPt->At(index) : part->pt();
                phi = isNANO ? genPhi->At(index) : part->phi();
                eta = isNANO ? genEta->At(index) : part->eta();
                m = isNANO ? genMass->At(index) : part->mass();

                Pt[nGenPart] = pt;
                Eta[nGenPart] = eta;
                Phi[nGenPart] = phi;
                Mass[nGenPart] = m;
            
                partID[nGenPart] = ID;
                mothID[nGenPart] = motherID;

                ++nGenPart;
            }
        }
    }
}

void GenPartAnalyzer::EndJob(TFile* file){
}
