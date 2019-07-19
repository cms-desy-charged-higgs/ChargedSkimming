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

class BaseAnalyzer {
    protected:
        //File path for SF etc.
        std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedHiggs/Skimming/data/";

        //Collection which are used in several analyzers if NANO AOD is analyzed
        TTreeReader* reader = NULL;
        bool isNANO;

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

        //Set trihObj and Gen particle collection
        void SetCollection(bool &isData);

        //Check for gen particle if it is last copy
        const reco::Candidate* LastCopy(const reco::GenParticle& part, const int& pdgID); //MINIAOD
        int LastCopy(const int& index, const int& pdgID); //NANOAOD

        //Trigger matching
        bool triggerMatching(const TLorentzVector &particle, const std::vector<pat::TriggerObjectStandAlone> trigObj = {});

    public:
        virtual ~BaseAnalyzer(){};
        BaseAnalyzer();
        BaseAnalyzer(TTreeReader* reader);
        virtual void BeginJob(TTree *tree, bool &isData) = 0;
        virtual bool Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event = NULL) = 0;
        virtual void EndJob(TFile* file) = 0;
};

#endif
