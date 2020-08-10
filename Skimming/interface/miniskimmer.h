#ifndef MINISKIMMER_H
#define MINISKIMMER_H

#include <memory>
#include <vector>
#include <string>
#include <map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/one/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>

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
        std::chrono::steady_clock::time_point start, end;

        //Output file
        std::vector<TFile*> outputFiles;

        //Trees for each systematic analysis
        std::vector<std::vector<TTree*>> outputTrees;
        std::vector<std::vector<CutFlow>> cutflows; 

        //Vector with wished analyzers for each systematic
        std::vector<std::vector<std::shared_ptr<BaseAnalyzer>>> analyzers;

        //EDM tokens
        jToken jetToken, fatjetToken; 
        genjToken genjetToken, genfatjetToken; 
        mToken metToken; 
        eToken eleToken; 
        muToken muonToken;
        trigToken triggerToken; 
        puToken pileupToken;
        genToken geninfoToken;
        genPartToken genParticleToken;
        edm::EDGetTokenT<double> rhoToken;
        secvtxToken secVertexToken;
        edm::EDGetTokenT<double> prefireToken, prefireTokenUp, prefireTokenDown;
        wgtToken pdfToken, scaleToken;

        //Channel
        std::vector<std::string> channels;
        float xSec;
        std::string outFile;
        int era;
        bool isData;     

        //Number of analyzed events
        int nEvents=0;

        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
};

#endif
