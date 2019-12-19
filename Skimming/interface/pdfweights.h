#ifndef PDFWEIGHTS_H
#define PDFWEIGHTS_H

#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include <FWCore/Utilities/interface/InputTag.h>

#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h>

class PDFWeights : public edm::EDProducer {
    public:
        explicit PDFWeights(const edm::ParameterSet&);
        ~PDFWeights();

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

        //Member functions
        std::vector<std::string> SplitString(const std::string& splitString, const std::string& delimeter);

        //Input tags
        edm::EDGetTokenT<LHEEventProduct> LHEToken;
        edm::EDPutTokenT<std::vector<float>> PDFToken;
        edm::EDPutTokenT<std::vector<float>> ScaleToken;
        edm::InputTag LHETag;
        int LHAID;

        //
        std::vector<std::string> variationIDs; 
        std::vector<std::string> scaleIDs; 
};

#endif
