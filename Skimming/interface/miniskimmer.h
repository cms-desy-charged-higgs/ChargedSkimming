#ifndef MINISKIMMER_H
#define MINISKIMMER_H

#include <memory>
#include <vector>
#include <string>
#include <map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/one/EDAnalyzer.h>

#include <ChargedSkimming/Skimming/interface/tokens.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

namespace pt = boost::property_tree;

class MiniSkimmer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit MiniSkimmer(const edm::ParameterSet&);
        ~MiniSkimmer();

    private:
        friend class Token;

        //Measure execution time
        std::chrono::steady_clock::time_point start, end;

        //Output file
        std::vector<TFile*> outputFiles;

        //Trees for each systematic analysis
        std::vector<std::vector<TTree*>> outputTrees;

        //Vector with wished analyzers for each systematic

        //Channel
        std::shared_ptr<Token> tokens;
        std::vector<std::string> channels;
        float xSec;
        std::string outFile;
        std::string era;
        bool isData;     

        //Number of analyzed events
        int nEvents=0;

        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
};

#endif
