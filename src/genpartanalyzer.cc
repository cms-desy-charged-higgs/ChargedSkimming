#include <ChargedHiggs/Skimming/interface/genpartanalyzer.h>

GenPartAnalyzer::GenPartAnalyzer(genPartToken& genParticleToken):
    BaseAnalyzer(),
    genParticleToken(genParticleToken)
    {}

GenPartAnalyzer::GenPartAnalyzer(TTreeReader& reader):
    BaseAnalyzer(&reader){}

void GenPartAnalyzer::BeginJob(TTree* tree, bool &isData){
    //Set data bool
    this->isData = isData;
    
    //Set TTreeReader for genpart and trigger obj from baseanalyzer
    if(isNANO) SetCollection(this->isData);

    //Set Branches of output tree
    tree->Branch("genPart", &genColl);
}


bool GenPartAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event){
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

            if(ID == 24){
                const reco::Candidate* part=NULL;
                int index=0;
            
                if(isNANO) index = LastCopy(i, 24);
                else part = LastCopy(genParts->at(i), 24);

                int motherID = isNANO ? abs(genID->At(genMotherIdx->At(i))) : abs(part->mother()->pdgId()); 

                float pt, eta, phi, m;
                pt = isNANO ? genPt->At(index) : part->pt();
                phi = isNANO ? genPhi->At(index) : part->phi();
                eta = isNANO ? genEta->At(index) : part->eta();
                m = isNANO ? genMass->At(index) : part->mass();

                if(motherID == 37){
                    genColl.W.SetPtEtaPhiM(pt, eta, phi, m);
                }
            }

            if(ID == 37){
                const reco::Candidate* part=NULL;
                int index=0;
            
                if(isNANO) index = LastCopy(i, 37);
                else part = LastCopy(genParts->at(i), 37);

                float pt, eta, phi, m;
                pt = isNANO ? genPt->At(index) : part->pt();
                phi = isNANO ? genPhi->At(index) : part->phi();
                eta = isNANO ? genEta->At(index) : part->eta();
                m = isNANO ? genMass->At(index) : part->mass();

                genColl.Hc.SetPtEtaPhiM(pt, eta, phi, m);
            }
        }
    }
    
    return true;
}


void GenPartAnalyzer::EndJob(TFile* file){
}
