#include <memory>
#include <vector>
#include <string>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <ChargedHiggs/Skimming/interface/baseanalyzer.h>
#include <ChargedHiggs/Skimming/interface/jetanalyzer.h>
#include <ChargedHiggs/Skimming/interface/electronanalyzer.h>
#include <ChargedHiggs/Skimming/interface/muonanalyzer.h>
#include <ChargedHiggs/Skimming/interface/triggeranalyzer.h>
#include <ChargedHiggs/Skimming/interface/metfilteranalyzer.h>
#include <ChargedHiggs/Skimming/interface/weightanalyzer.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

class MiniSkimmer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit MiniSkimmer(const edm::ParameterSet&);
        ~MiniSkimmer();

    private:
        //Measure execution time
        std::chrono::steady_clock::time_point start;
        std::chrono::steady_clock::time_point end;

        //Output file
        TFile* outputFile;
        std::vector<TTree*> outputTrees;
        std::vector<std::pair<TH1F*, float>> cutflows; 

        //Map with vector with analyzers
        std::map<std::string, std::vector<std::shared_ptr<BaseAnalyzer>>> analyzerMap;

        //Vector with wished analyzers
        std::vector<std::vector<std::shared_ptr<BaseAnalyzer>>> analyzers;

        //EDM tokens
        jToken jetToken; 
        jToken fatjetToken; 
        genjToken genjetToken; 
        genjToken genfatjetToken; 
        mToken metToken; 
        eToken eleToken; 
        muToken muonToken;
        trigToken triggerToken; 
        puToken pileupToken;
        genToken geninfoToken;
        edm::EDGetTokenT<double> rhoToken;

        //Channel
        std::vector<std::string> channels;
        float xSec;
        std::string outFile;
        bool isData;            

        //Number of analyzed events
        int nEvents=0;

        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
};
