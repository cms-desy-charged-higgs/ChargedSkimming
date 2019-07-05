#include <ChargedHiggs/Skimming/interface/weightanalyzer.h>

WeightAnalyzer::WeightAnalyzer(const float era, const float xSec, TTreeReader &reader):
    BaseAnalyzer(&reader),
    era(era),
    xSec(xSec)
    {}

void WeightAnalyzer::BeginJob(TTree *tree, bool &isData){
    //Set lumi map
    lumis = {{2016, 35.92*1e3}, {2017, 41.53*1e3}};

    pileUpFiles = {
            {2017, filePath + "pileUp/data_pileUp2017.root"}, 
    };

    //Set data bool
    this->isData = isData;

    if(!this->isData){
        //Calculate pile up weights
        TFile* puFile = TFile::Open(pileUpFiles[era].c_str());
        TH1F* puData = (TH1F*)puFile->Get("pileup");

        TH1F* puMC = new TH1F("puMC", "puMC", 100, 0, 100);
        reader->GetTree()->Draw("Pileup_nTrueInt>>puMC");

        puData->Scale(1./puData->Integral(), "width");
        puMC->Scale(1./puMC->Integral(), "width");

        pileUpWeights = new TH1F("pileUpWeights", "pileUpWeights", 100, 0, 100);
        pileUpWeights->Divide(puData, puMC);

        //Initiliaze TTreeReaderValues
        genWeightValue = std::make_unique<TTreeReaderValue<float>>(*reader, "Generator_weight");
        nPU = std::make_unique<TTreeReaderValue<float>>(*reader, "Pileup_nTrueInt");
        
        nGenHist = new TH1F("nGen", "nGen", 100, 0, 2);
    }

    evtNumber = std::make_unique<TTreeReaderValue<ULong64_t>>(*reader, "event");

    //Branches for output tree
    tree->Branch("lumi", &lumi);
    tree->Branch("xsec", &xSec);
    tree->Branch("genWeight", &genWeight);
    tree->Branch("puWeight", &puWeight);
    tree->Branch("eventNumber", &eventNumber);
}

bool WeightAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow){
    //Set values if not data
    if(!this->isData){
        lumi = lumis[era];
        puWeight = pileUpWeights->GetBinContent(pileUpWeights->FindBin(*nPU->Get())) != 0 ? pileUpWeights->GetBinContent(pileUpWeights->FindBin(*nPU->Get())) : 1.;
        genWeight = *genWeightValue->Get();

        nGenHist->Fill(1);
    }

    eventNumber = *evtNumber->Get();

    cutflow.second = xSec*lumi*puWeight;

    cutflow.first->Fill("No cuts", cutflow.second);
    return true;
}

void WeightAnalyzer::EndJob(TFile* file){
    if(!this->isData){
        if(!file->GetListOfKeys()->Contains("nGen")){
            nGenHist->Write();
        }
    }
}
