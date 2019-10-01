#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

BaseAnalyzer::BaseAnalyzer(): isNANO(false){}
BaseAnalyzer::BaseAnalyzer(TTreeReader* reader): reader(reader), isNANO(true){}

void BaseAnalyzer::SetCollection(bool &isData){
    if(!isData){
        genPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_pt");
        genEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_eta");
        genPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_phi");
        genMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenPart_mass");
        genID = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_pdgId");
        genMotherIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_genPartIdxMother");
        genStatus = std::make_unique<TTreeReaderArray<int>>(*reader, "GenPart_statusFlags");
        eleGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_genPartIdx");
        muonGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_genPartIdx");
    }

    run = std::make_unique<TTreeReaderValue<unsigned int>>(*reader, "run");

    trigObjPt = std::make_unique<TTreeReaderArray<float>>(*reader, "TrigObj_pt");
    trigObjEta = std::make_unique<TTreeReaderArray<float>>(*reader, "TrigObj_eta");
    trigObjPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "TrigObj_phi");
    trigObjID = std::make_unique<TTreeReaderArray<int>>(*reader, "TrigObj_id");
    trigObjFilterBit = std::make_unique<TTreeReaderArray<int>>(*reader, "TrigObj_filterBits");
}


const reco::Candidate* BaseAnalyzer::LastCopy(const reco::GenParticle& part, const int& pdgID){
    const reco::Candidate* cand = &part;

    while(abs(cand->mother()->pdgId()) == pdgID){
        cand = cand->mother();
    }

    return cand;
}

const reco::Candidate* BaseAnalyzer::LastCopy(const reco::Candidate* part, const int& pdgID){
    const reco::Candidate* cand = part;

    while(abs(cand->mother()->pdgId()) == pdgID){
        cand = cand->mother();
    }

    return cand;
}

int BaseAnalyzer::LastCopy(const int& index, const int& pdgID){
    int idx = genMotherIdx->At(index);

    while(abs(genID->At(idx)) == pdgID){
        idx = genMotherIdx->At(idx);
    }

    return idx; 
}

bool BaseAnalyzer::triggerMatching(const TLorentzVector &particle, const std::vector<pat::TriggerObjectStandAlone> trigObj){
    int size = isNANO ? trigObjEta->GetSize() : trigObj.size();

    for(int i=0; i < size; i++){
        //Check if trig obj is electron
        float pt, eta, phi;

        pt = isNANO ? trigObjPt->At(i) : trigObj[i].pt();
        eta = isNANO ? trigObjEta->At(i) : trigObj[i].eta();
        phi = isNANO ? trigObjPhi->At(i) : trigObj[i].phi();

        if(5e-2 > std::sqrt(std::pow(eta - particle.Eta(), 2) + std::pow(phi - particle.Phi(), 2)) and abs(particle.Pt()-pt)/particle.Pt() < 0.02){
            return true;
        }
    }

    return false;
}

bool BaseAnalyzer::SetGenParticles(TLorentzVector &validLepton, const int &i, const int &pdgID, const std::vector<reco::GenParticle>& genParticle){
    const reco::GenParticle* matchedLep=NULL;

    if(!isNANO){
        for(const reco::GenParticle &part: genParticle){
            if(part.isPromptFinalState()){
                if(0.3 > std::sqrt(std::pow(part.eta() - validLepton.Eta(), 2) + std::pow(part.phi() - validLepton.Phi(), 2)) and abs(validLepton.Pt()-part.pt())/validLepton.Pt() < 0.05){
                    matchedLep = &part;
                }
            }
        }
    }   

    bool isgenMatched = isNANO ? eleGenIdx->At(i) != -1 : matchedLep!=NULL;

    //Check if gen matched particle exist
    if(isgenMatched){
        const reco::Candidate* lepton=NULL;
        int index=0;
            
        if(isNANO) index = LastCopy(pdgID == 11 ? eleGenIdx->At(i) : muonGenIdx->At(i), pdgID);
        else lepton = LastCopy(*matchedLep, pdgID);

        int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(lepton->mother()->pdgId()); 

        if(motherID == 24){
            const reco::Candidate* WBoson=NULL;
            int index=0;
                
            if(isNANO) index = LastCopy(eleGenIdx->At(i), 24);
            else WBoson = LastCopy(lepton->mother(), 24);

            int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(WBoson->mother()->pdgId());     

            if(motherID == 37){
                return true;
            }

            return false;
        }

        return false;
    }

    return false;
}
