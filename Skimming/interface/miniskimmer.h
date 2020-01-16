#ifndef MINISKIMMER_H
#define MINISKIMMER_H

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

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/jetanalyzer.h>
#include <ChargedSkimming/Skimming/interface/electronanalyzer.h>
#include <ChargedSkimming/Skimming/interface/muonanalyzer.h>
#include <ChargedSkimming/Skimming/interface/triggeranalyzer.h>
#include <ChargedSkimming/Skimming/interface/metfilteranalyzer.h>
#include <ChargedSkimming/Skimming/interface/weightanalyzer.h>
#include <ChargedSkimming/Skimming/interface/genpartanalyzer.h>

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

        //Trees for each systematic analysis
        std::vector<TFile*> outputFiles;
        std::vector<std::vector<TTree*>> outputTrees;
        std::vector<std::vector<CutFlow>> cutflows; 

        //Vector with wished analyzers for each systematic
        std::vector<std::vector<std::shared_ptr<BaseAnalyzer>>> analyzers;

        //Map with systematic uncertanties
        std::vector<std::pair<std::string, std::string>> systNames;

        //EDM tokens
        jToken jetToken; 
        jToken fatjetToken; 
        genjToken genjetToken; 
        genjToken genfatjetToken; 
        mToken metToken; 
        eToken eleToken; 
        muToken muonToken;
        trigToken triggerToken; 
        trigObjToken triggerObjToken;
        puToken pileupToken;
        genToken geninfoToken;
        genPartToken genParticleToken;
        edm::EDGetTokenT<double> rhoToken;
        vtxToken vertexToken;
        secvtxToken secVertexToken;
        edm::EDGetTokenT<double> prefireToken;
        edm::EDGetTokenT<double> prefireTokenUp;
        edm::EDGetTokenT<double> prefireTokenDown;
        wgtToken pdfToken;
        wgtToken scaleToken;

        //Channel
        std::vector<std::string> channels;
        float xSec;
        std::string outFile;
        bool isData;           

        std::map<std::string, std::vector<unsigned int>> nMin;

        //Number of analyzed events
        int nEvents=0;

        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
};

#endif
