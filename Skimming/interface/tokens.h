#ifndef TOKENS_H
#define TOKENS_H

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <DataFormats/PatCandidates/interface/IsolatedTrack.h>
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

struct Token {
    public:
        //Jet related edm token
        edm::EDGetTokenT<std::vector<pat::Jet>> jetToken, fatJetToken;
        edm::EDGetTokenT<std::vector<reco::GenJet>> genJetToken, genFatJetToken;
        edm::EDGetTokenT<std::vector<pat::MET>> METToken;
        edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>> secVertexToken; 

        //Electron related token
        edm::EDGetTokenT<std::vector<pat::Electron>> eleToken;

        //Muon related token
        edm::EDGetTokenT<std::vector<pat::Muon>> muonToken;

        //Gen related token
        edm::EDGetTokenT<GenEventInfoProduct> genToken;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken;
        edm::EDGetTokenT<LHEEventProduct> lheToken;

        //PDF related token
        edm::EDGetTokenT<std::vector<float>> pdfToken, scaleToken;

        //Weight token
        edm::EDGetTokenT<double> prefireToken, prefireTokenUp, prefireTokenDown;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileUpToken;

        //Miscellaneous
        edm::EDGetTokenT<std::vector<pat::IsolatedTrack>> isoTrackToken;
        edm::EDGetTokenT<edm::TriggerResults> triggerToken;
        edm::EDGetTokenT<double> rhoToken;

        template <typename T>
        static const T GetTokenValue(const edm::Event* event, const edm::EDGetTokenT<T>& token){
            edm::Handle<T> handle;
            event->getByToken(token, handle);

            return *handle;
        }
};

#endif
