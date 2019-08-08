#ifndef BASEANALYZER_H
#define BASEANALYZER_H

#include <memory>
#include <map>
#include <vector>
#include <cmath>
#include <bitset>

#include <TFile.h>
#include <TH2F.h>
#include <Rtypes.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <FWCore/Framework/interface/Event.h>

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include "DataFormats/Candidate/interface/Candidate.h"
#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>
#include <FWCore/Common/interface/TriggerNames.h>

//Struct for saving gen particle lorentz vectors
struct GenParticles{
    TLorentzVector W;
    TLorentzVector h1;
    TLorentzVector h2;
    TLorentzVector Hc;
};

//Struct for cutflow
struct CutFlow {
    TH1F* hist;
    Float_t weight = 1.;

    unsigned int nMinEle=0;
    unsigned int nMinMu=0;
    unsigned int nMinJet=0;
    unsigned int nMinFatjet=0;
    
    bool passed = true;

};

typedef edm::EDGetTokenT<std::vector<pat::Jet>> jToken;
typedef edm::EDGetTokenT<std::vector<reco::GenJet>> genjToken;
typedef edm::EDGetTokenT<std::vector<pat::Electron>> eToken;
typedef edm::EDGetTokenT<std::vector<pat::Muon>> muToken;
typedef edm::EDGetTokenT<edm::TriggerResults> trigToken;
typedef edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken;
typedef edm::EDGetTokenT<GenEventInfoProduct> genToken;
typedef edm::EDGetTokenT<std::vector<pat::MET>> mToken;
typedef edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> trigObjToken;
typedef edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken;
typedef edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken;
typedef edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> secvtxToken;

class BaseAnalyzer {
    protected:
        //File path for SF etc.
        std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedHiggs/Skimming/data/";

        //Collection which are used in several analyzers if NANO AOD is analyzed
        TTreeReader* reader = NULL;
        bool isNANO;

        std::map<int, std::map<std::string, std::pair<int, int>>> runEras = {
              {2017, {
                        {"B", {297046, 299329}}, 
                        {"C", {299368, 302029}}, 
                        {"DE", {302030, 304797}},
                        {"F", {305040, 306462}}, 
                     }
              },
        };

        std::unique_ptr<TTreeReaderValue<unsigned int>> run;

        std::unique_ptr<TTreeReaderArray<float>> trigObjPt;
        std::unique_ptr<TTreeReaderArray<float>> trigObjPhi;
        std::unique_ptr<TTreeReaderArray<float>> trigObjEta;
        std::unique_ptr<TTreeReaderArray<int>> trigObjID;
        std::unique_ptr<TTreeReaderArray<int>> trigObjFilterBit;
    
        std::unique_ptr<TTreeReaderArray<float>> genPhi;
        std::unique_ptr<TTreeReaderArray<float>> genEta;
        std::unique_ptr<TTreeReaderArray<float>> genPt;
        std::unique_ptr<TTreeReaderArray<float>> genMass;
        std::unique_ptr<TTreeReaderArray<int>> genID;
        std::unique_ptr<TTreeReaderArray<int>> genMotherIdx;
        std::unique_ptr<TTreeReaderArray<int>> genStatus;
        std::unique_ptr<TTreeReaderArray<int>> eleGenIdx;
        std::unique_ptr<TTreeReaderArray<int>> muonGenIdx;

        //Set trihObj and Gen particle collection
        void SetCollection(bool &isData);

        //Check for gen particle if it is last copy
        const reco::Candidate* LastCopy(const reco::GenParticle& part, const int& pdgID);
        const reco::Candidate* LastCopy(const reco::Candidate* part, const int& pdgID); //MINIAOD
        int LastCopy(const int& index, const int& pdgID); //NANOAOD

        //Match Reco to gen particles
        template<typename Lepton>
        void SetGenParticles(Lepton &lepton, const int &i, const int &pdgID, const std::vector<reco::GenParticle>& genParticle={});

        //Trigger matching
        bool triggerMatching(const TLorentzVector &particle, const std::vector<pat::TriggerObjectStandAlone> trigObj = {});


    public:
        virtual ~BaseAnalyzer(){};
        BaseAnalyzer();
        BaseAnalyzer(TTreeReader* reader);
        virtual void BeginJob(std::vector<TTree*>& trees, bool &isData) = 0;
        virtual void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event = NULL) = 0;
        virtual void EndJob(TFile* file) = 0;
};


template<typename Lepton>
void BaseAnalyzer::SetGenParticles(Lepton &validLepton, const int &i, const int &pdgID, const std::vector<reco::GenParticle>& genParticle){
    const reco::GenParticle* matchedLep=NULL;

    if(!isNANO){
        for(const reco::GenParticle &part: genParticle){
            if(part.isPromptFinalState()){
                if(0.3 > std::sqrt(std::pow(part.eta() - validLepton.fourVec.Eta(), 2) + std::pow(part.phi() - validLepton.fourVec.Phi(), 2)) and abs(validLepton.fourVec.Pt()-part.pt())/validLepton.fourVec.Pt() < 0.05){
                    matchedLep = &part;
                }
            }
        }
    }   

    validLepton.isgenMatched = isNANO ? eleGenIdx->At(i) != -1 : matchedLep!=NULL;

    //Check if gen matched particle exist
    if(validLepton.isgenMatched){
        const reco::Candidate* lepton=NULL;
        int index=0;
            
        if(isNANO) index = LastCopy(pdgID == 11 ? eleGenIdx->At(i) : muonGenIdx->At(i), pdgID);
        else lepton = LastCopy(*matchedLep, pdgID);

        float pt, eta, phi, m;
        pt = isNANO ? genPt->At(index) : lepton->pt();
        phi = isNANO ? genPhi->At(index) : lepton->phi();
        eta = isNANO ? genEta->At(index) : lepton->eta();
        m = isNANO ? genMass->At(index) : lepton->mass();

        validLepton.genVec.SetPtEtaPhiM(pt, eta, phi, m);

        int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(lepton->mother()->pdgId()); 

        if(motherID == 24){
            const reco::Candidate* WBoson=NULL;
            int index=0;
                
            if(isNANO) index = LastCopy(eleGenIdx->At(i), 24);
            else WBoson = LastCopy(lepton->mother(), 24);

            int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(WBoson->mother()->pdgId());     

            if(motherID == 37){
                validLepton.isFromHc = true;
            }
        }
    }
}

#endif
