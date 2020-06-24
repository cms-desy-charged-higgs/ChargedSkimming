#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

BaseAnalyzer::BaseAnalyzer(): isNANO(false){}
BaseAnalyzer::BaseAnalyzer(TTreeReader* reader): reader(reader), isNANO(true){}

float BaseAnalyzer::DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2){
    return std::sqrt(std::pow(phi1 - phi2, 2) + std::pow(eta1 - eta2, 2)); 
}

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

const reco::Candidate* BaseAnalyzer::FirstCopy(const reco::GenParticle& part, const int& pdgID){
    const reco::Candidate* daughter = &part;
    const reco::Candidate* mother = daughter->mother();

    if(mother == nullptr) return daughter;

    while(abs(mother->pdgId()) == pdgID){
        daughter = mother;
        mother = daughter->mother();
    }

    return daughter;
}

const reco::Candidate* BaseAnalyzer::FirstCopy(const reco::Candidate* part, const int& pdgID){
    const reco::Candidate* daughter = part;
    const reco::Candidate* mother = daughter->mother();
        
    if(mother == nullptr) return daughter;

    while(abs(mother->pdgId()) == pdgID){
        daughter = mother;
        mother = daughter->mother();
    }

    return daughter;
}

int BaseAnalyzer::FirstCopy(const int& index, const int& pdgID){
    int daughterIdx = index;
    int motherIdx = genMotherIdx->At(index);

    while(abs(genID->At(motherIdx)) == pdgID){
        daughterIdx = motherIdx;
        motherIdx = genMotherIdx->At(daughterIdx);
    }

    return daughterIdx; 
}

std::tuple<int, int, int> BaseAnalyzer::SetGenParticles(const float& pt, const float& eta, const float& phi, const int &i, const int &pdgID, const std::vector<reco::GenParticle>& genParticle){
    const reco::GenParticle* matchedLep=nullptr;
    std::tuple<int, int, int> IDs = std::make_tuple(-99., -99., -99.);

    if(!isNANO){
        for(const reco::GenParticle& part: genParticle){
            if(part.isPromptFinalState()){
                if(0.5 > BaseAnalyzer::DeltaR(part.eta(), part.phi(), eta, phi) and abs(pt-part.pt())/pt < 0.5){
                    matchedLep = &part;
                }
            }
        }
    }  

    std::map<int, int> genIndex;
    if(isNANO) genIndex = {{11, eleGenIdx->At(i)}, {13, muonGenIdx->At(i)}};

    bool isgenMatched = isNANO ? genIndex[pdgID] != -1 : matchedLep!=nullptr;

    //Check if gen matched particle exist
    if(isgenMatched){
        const reco::Candidate* lepton=nullptr;
        const reco::Candidate* leptonMother=nullptr;
        int index=0, motherIdx=0;
            
        if(isNANO) index = FirstCopy(genIndex[pdgID], pdgID);
        else lepton = FirstCopy(matchedLep, pdgID);

        const int& motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(lepton->mother()->pdgId());

        if(isNANO) motherIdx = FirstCopy(genMotherIdx->At(index), motherID);
        else leptonMother = FirstCopy(lepton->mother(), motherID);

        const int& grandMotherID = isNANO ? abs(genID->At(genMotherIdx->At(motherIdx))) : leptonMother->mother() != nullptr ? abs(leptonMother->mother()->pdgId()) : -99.;  
    
        IDs = std::make_tuple(pdgID, motherID, grandMotherID);
    }

    return IDs;
}
