#include <ChargedHiggs/Skimming/interface/genpartanalyzer.h>

GenPartAnalyzer::GenPartAnalyzer(genPartToken& genParticleToken):
    BaseAnalyzer(),
    genParticleToken(genParticleToken)
    {}

GenPartAnalyzer::GenPartAnalyzer(TTreeReader& reader):
    BaseAnalyzer(&reader){}

void GenPartAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData){
    //Set data bool
    this->isData = isData;
    
    //Set TTreeReader for genpart and trigger obj from baseanalyzer
    if(isNANO) SetCollection(this->isData);

    floatNames = {"E", "Px", "Py", "Pz"};
    leptonVariables = std::vector<std::vector<float>>(floatNames.size(), std::vector<float>(2, -999.));
    h1Variables = std::vector<std::vector<float>>(floatNames.size(), std::vector<float>());
    h2Variables = std::vector<std::vector<float>>(floatNames.size(), std::vector<float>());

    //Set Branches of output tree
    for(TTree* tree: trees){
        for(unsigned int i=0; i<floatNames.size(); i++){
            tree->Branch(("GenLepton_" + floatNames[i]).c_str(), &leptonVariables[i]);
        }

        for(unsigned int i=0; i<floatNames.size(); i++){
            tree->Branch(("GenPartFromh1_" + floatNames[i]).c_str(), &h1Variables[i]);
        }

       for(unsigned int i=0; i<floatNames.size(); i++){
            tree->Branch(("GenPartFromh2_" + floatNames[i]).c_str(), &h2Variables[i]);
        }
    }
}


void GenPartAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear vectors
    for(std::vector<float>& variable: leptonVariables){
        variable = std::vector<float>(2, -999.);
    }

    for(std::vector<float>& variable: h1Variables){
        variable.clear();
    }

    for(std::vector<float>& variable: h2Variables){
        variable.clear();
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    if(!isData){
        if(!isNANO){
            event->getByToken(genParticleToken, genParts);
        }

        int size = isNANO ? genID->GetSize() : genParts->size();

        //Fill 4 four vectors
        for(int i = 0; i < size; i++){
            int ID = isNANO ? abs(genID->At(i)) : abs(genParts->at(i).pdgId());  

            if(ID == 11 or ID == 12 or ID == 13 or ID == 14){
                const reco::Candidate* part=NULL;
                int index=0;

                if(isNANO) index = LastCopy(i, ID);
                else part = LastCopy(genParts->at(i), ID);

                int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(part->mother()->pdgId()); 

                float pt, eta, phi, m;
                pt = isNANO ? genPt->At(index) : part->pt();
                phi = isNANO ? genPhi->At(index) : part->phi();
                eta = isNANO ? genEta->At(index) : part->eta();
                m = isNANO ? genMass->At(index) : part->mass();

                TLorentzVector lVec;
                lVec.SetPtEtaPhiM(pt, eta, phi, m);

                if(motherID == 24){
                    const reco::Candidate* part=NULL;
                    int index=0;
                
                    if(isNANO) index = LastCopy(i, 24);
                    else part = LastCopy(genParts->at(i), 24);

                    int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(part->mother()->pdgId()); 

                    if(motherID == 37){
                        //Lepton four momentum components
                        leptonVariables[0][ID%2==0 ? 0 : 1] = lVec.E();   //Energy
                        leptonVariables[1][ID%2==0 ? 0 : 1] = lVec.Px();  //Px
                        leptonVariables[2][ID%2==0 ? 0 : 1] = lVec.Py();  //Py
                        leptonVariables[3][ID%2==0 ? 0 : 1] = lVec.Pz();  //Pz
                    }
                }

                if(ID == 25){
                    const reco::Candidate* part=NULL;
                    int index=0;
                
                    if(isNANO) index = LastCopy(i, 25);
                    else part = LastCopy(genParts->at(i), 25);

                    int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(part->mother()->pdgId()); 

                    float pt, eta, phi, m;
                    pt = isNANO ? genPt->At(index) : part->pt();
                    phi = isNANO ? genPhi->At(index) : part->phi();
                    eta = isNANO ? genEta->At(index) : part->eta();
                    m = isNANO ? genMass->At(index) : part->mass();

                    if(motherID == 37){
                        genColl.h1.SetPtEtaPhiM(pt, eta, phi, m);
                    }

                    else{
                        genColl.h2.SetPtEtaPhiM(pt, eta, phi, m);
                    }
                }
            } 

            if(ID == 5){
                const reco::Candidate* part=NULL;
                int index=0;

                if(isNANO) index = LastCopy(i, ID);
                else part = LastCopy(genParts->at(i), ID);

                int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(part->mother()->pdgId()); 

                float pt, eta, phi, m;
                pt = isNANO ? genPt->At(index) : part->pt();
                phi = isNANO ? genPhi->At(index) : part->phi();
                eta = isNANO ? genEta->At(index) : part->eta();
                m = isNANO ? genMass->At(index) : part->mass();

                TLorentzVector lVec;
                lVec.SetPtEtaPhiM(pt, eta, phi, m);

                if(motherID == 25){
                    const reco::Candidate* part=NULL;
                    int index=0;
                
                    if(isNANO) index = LastCopy(i, 25);
                    else part = LastCopy(genParts->at(i), 25);

                    int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(part->mother()->pdgId()); 

                    if(motherID == 37 and h1Variables[0].size() < 2){
                        //Quark four momentum components
                        h1Variables[0].push_back(lVec.E());   //Energy
                        h1Variables[1].push_back(lVec.Px());  //Px
                        h1Variables[2].push_back(lVec.Py());  //Py
                        h1Variables[3].push_back(lVec.Pz());  //Pz
                    }
            
                    else if(h2Variables[0].size() < 2){
                        //Quark four momentum components
                        h2Variables[0].push_back(lVec.E());   //Energy
                        h2Variables[1].push_back(lVec.Px());  //Px
                        h2Variables[2].push_back(lVec.Py());  //Py
                        h2Variables[3].push_back(lVec.Pz());  //Pz
                    }
                }
            }
        }

        if(!h1Variables[0].empty() and !h2Variables[0].empty()){
            //Sort b quarks for pt
            std::function<bool(std::vector<float>, std::vector<float>)> sortFunc = [&](std::vector<float> b1, std::vector<float> b2){return std::pow(b1[1], 2) + std::pow(b1[2], 2) > std::pow(b2[1], 2) + std::pow(b2[2], 2);};

            std::sort(h1Variables.begin(), h1Variables.end(), sortFunc);
            std::sort(h2Variables.begin(), h2Variables.end(), sortFunc);
        }
    }
}


void GenPartAnalyzer::EndJob(TFile* file){
}
