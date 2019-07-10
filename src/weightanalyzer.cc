#include <ChargedHiggs/Skimming/interface/weightanalyzer.h>

WeightAnalyzer::WeightAnalyzer(const float era, const float xSec, TTreeReader &reader):
    BaseAnalyzer(&reader),
    era(era),
    xSec(xSec)
    {}

WeightAnalyzer::WeightAnalyzer(const float era, const float xSec, puToken &pileupToken, genToken &geninfoToken):
    BaseAnalyzer(),
    era(era),
    xSec(xSec),
    pileupToken(pileupToken),
    geninfoToken(geninfoToken)
    {}


void WeightAnalyzer::BeginJob(TTree *tree, bool &isData){
    //Set lumi map
    lumis = {{2016, 35.92*1e3}, {2017, 41.53*1e3}};

    //Set data bool
    this->isData = isData;

    if(!this->isData){
        if(isNANO){
            //Initiliaze TTreeReaderValues
            genWeightValue = std::make_unique<TTreeReaderValue<float>>(*reader, "Generator_weight");
            nPU = std::make_unique<TTreeReaderValue<float>>(*reader, "Pileup_nTrueInt");
        }        

        puMC = new TH1F("puMC", "puMC", 100, 0, 100);
        nGenHist = new TH1F("nGen", "nGen", 100, 0, 2);
    }

    evtNumber = std::make_unique<TTreeReaderValue<ULong64_t>>(*reader, "event");

    //Branches for output tree
    tree->Branch("lumi", &lumi);
    tree->Branch("xsec", &xSec);
    tree->Branch("genWeight", &genWeight);
    tree->Branch("nTrueInt", &nTrueInt);
    tree->Branch("eventNumber", &eventNumber);
}

bool WeightAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event){
    edm::Handle<std::vector<PileupSummaryInfo>> pileUp; 
    edm::Handle<GenEventInfoProduct> genInfo;
    
    if(!isNANO){
        event->getByToken(pileupToken, pileUp);
        event->getByToken(geninfoToken, genInfo);
    }

    //Set values if not data
    if(!this->isData){
        lumi = lumis[era];
        nTrueInt = isNANO ? *nPU->Get() : pileUp->at(1).getTrueNumInteractions(); 
        genWeight = isNANO ? *genWeightValue->Get() : genInfo->weight();

        nGenHist->Fill(1);
        puMC->Fill(nTrueInt);
    }

    //eventNumber = *evtNumber->Get();

    cutflow.second = xSec*lumi;

    cutflow.first->Fill("No cuts", cutflow.second);
    return true;
}

void WeightAnalyzer::EndJob(TFile* file){
    if(!this->isData){
        if(!file->GetListOfKeys()->Contains("nGen")){
            nGenHist->Write();
            puMC->Write();
        }
    }
}
