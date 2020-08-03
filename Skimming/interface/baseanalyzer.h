#ifndef BASEANALYZER_H
#define BASEANALYZER_H

#include <memory>
#include <map>
#include <vector>
#include <cmath>
#include <bitset>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <Rtypes.h>
#include <Math/Vector4Dfwd.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <FWCore/Framework/interface/Event.h>

#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

typedef edm::EDGetTokenT<edm::TriggerResults> trigToken;
typedef edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken;

//Struct for cutflow
struct CutFlow {
    TH1F* hist;
    Float_t weight = 1.;
    std::string channel;

    unsigned char nMinEle=0;
    unsigned char nMinMu=0;
    unsigned char nMinJet=0;
    unsigned char nMinFatjet=0;
    
    bool passed = true;
};

class BaseAnalyzer {
    protected:
        //File path for SF etc.
        std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/";

        std::map<int, std::map<std::string, std::pair<int, int>>> runEras = {
              {2016, {
                        {"BCD", {272007, 276811}}, 
                        {"EF", {276831, 278808}}, 
                        {"GH", {278820, 284044}},
                     }
              },

              {2017, {
                        {"B", {297046, 299329}}, 
                        {"C", {299368, 302029}}, 
                        {"DE", {302030, 304797}},
                        {"F", {305040, 306462}}, 
                     }
              },

              {2018, {
                        {"A", {315252, 316995}}, 
                        {"B", {316998, 319312}}, 
                        {"C", {319313, 320393}},
                        {"D", {320394, 325273}}, 
                     }
              },
        };

        //Collection which are used in several analyzers if NANO AOD is analyzed
        TTreeReader* reader = nullptr;
        bool isNANO;

        std::unique_ptr<TTreeReaderValue<unsigned int>> run;

        std::unique_ptr<TTreeReaderArray<float>> trigObjPt, trigObjPhi, trigObjEta;
        std::unique_ptr<TTreeReaderArray<int>> trigObjID, trigObjFilterBit;
    
        std::unique_ptr<TTreeReaderArray<float>> genPhi, genEta, genPt, genMass;
        std::unique_ptr<TTreeReaderArray<int>> genID, genMotherIdx, genStatus, eleGenIdx, muonGenIdx;

        //Set trihObj and Gen particle collection
        void SetCollection(bool &isData);
        bool isSyst=false;

        //Check for gen particle if it is last copy
        const reco::Candidate* FirstCopy(const reco::GenParticle& part, const int& pdgID);
        const reco::Candidate* FirstCopy(const reco::Candidate* part, const int& pdgID); //MINIAOD
        int FirstCopy(const int& index, const int& pdgID); //NANOAOD

        //Match Reco to gen particles
        std::tuple<int, int, int> SetGenParticles(const float& Pt, const float& Eta, const float& Phi, const int &i, const int& pdgID, const std::vector<reco::GenParticle>& genParticle={});

    public:
        virtual ~BaseAnalyzer(){};
        BaseAnalyzer();
        BaseAnalyzer(TTreeReader* reader);
        virtual void BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst=false) = 0;
        virtual void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event = NULL) = 0;
        virtual void EndJob(TFile* file) = 0;

        static float DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2);
};
#endif
