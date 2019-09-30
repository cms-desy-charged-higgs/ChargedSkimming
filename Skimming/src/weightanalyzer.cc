#include <ChargedAnalysis/Skimming/interface/weightanalyzer.h>

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


void WeightAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData){
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
    for(TTree* tree: trees){
        tree->Branch("Weight_lumi", &lumi);
        tree->Branch("Weight_xsec", &xSec);
        tree->Branch("Weight_genWeight", &genWeight);
        tree->Branch("Misc_TrueInteraction", &nTrueInt);
        tree->Branch("Misc_eventNumber", &eventNumber);
    }
}

void WeightAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
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

    eventNumber = isNANO ? *evtNumber->Get() : event->eventAuxiliary().id().event();

    for(CutFlow& cutflow: cutflows){
        cutflow.weight = xSec*lumi;
        cutflow.hist->Fill("No cuts", cutflow.weight);
        cutflow.passed *= true;
    }
}

void WeightAnalyzer::EndJob(TFile* file){
    if(!this->isData){
        nGenHist->Write();
        puMC->Write();
    }

    delete puMC;
    delete nGenHist;
}
