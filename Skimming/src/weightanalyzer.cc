#include <ChargedSkimming/Skimming/interface/weightanalyzer.h>

WeightAnalyzer::WeightAnalyzer(const int& era, const float& xSec, TTreeReader& reader):
    BaseAnalyzer(&reader),
    era(era),
    xSec(xSec)
    {}

WeightAnalyzer::WeightAnalyzer(const int& era, const float& xSec, puToken& pileupToken, genToken& geninfoToken, const std::vector<edm::EDGetTokenT<double>>& prefireTokens, wgtToken& pdfToken, wgtToken& scaleToken):
    BaseAnalyzer(),
    era(era),
    xSec(xSec),
    pileupToken(pileupToken),
    geninfoToken(geninfoToken),
    prefireTokens(prefireTokens),
    scaleToken(scaleToken),
    pdfToken(pdfToken)
    {}

void WeightAnalyzer::BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst){
    //Set data bool
    this->isData = isData;
    this->isSyst = isSyst;

    //Read infos from skim.json
    boost::property_tree::ptree parser; 
    boost::property_tree::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/skim.json", parser);
    
    pileUpFile = filePath + parser.get<std::string>("Analyzer.Weight.PileUp." + std::to_string(era));
    lumi = parser.get<float>("Analyzer.Weight.Lumi." + std::to_string(era))*1e3;

    if(!this->isData){
        if(isNANO){
            //Initiliaze TTreeReaderValues
            genWeightValue = std::make_unique<TTreeReaderValue<float>>(*reader, "Generator_weight");
            nPU = std::make_unique<TTreeReaderValue<float>>(*reader, "Pileup_nTrueInt");
        }        

        puMC = std::make_shared<TH1F>("puMC", "puMC", 100, 0, 100);
        nGenHist = std::make_shared<TH1F>("nGen", "nGen", 100, 0, 2);
        nGenWeightedHist = std::make_shared<TH1F>("nGenWeighted", "nGenWeighted", 100, -1e7, 1e7);
    }

    evtNumber = std::make_unique<TTreeReaderValue<ULong64_t>>(*reader, "event");

    prefireWeights = {1., 1., 1.};

    //Branches for output tree
    for(TTree* tree: trees){
        tree->Branch("Weight_genWeight", &genWeight);
        tree->Branch("Weight_prefireWeight", &prefireWeights[0]);
        tree->Branch("Misc_TrueInteraction", &nTrueInt);
        tree->Branch("Misc_eventNumber", &eventNumber);

        std::fill_n(scaleWeights, 8, 1);
        std::fill_n(pdfWeights, 102, 1);

        if(!isSyst){
            tree->Branch("Weight_pdfVariations", pdfWeights, "Weight_pdfVariations[102]/F");
            tree->Branch("Weight_scaleVariations", scaleWeights, "Weight_scaleVariations[8]/F");
            tree->Branch("Weight_prefireWeightUp", &prefireWeights[1]);
            tree->Branch("Weight_prefireWeightDown", &prefireWeights[2]);
        }
    }
}

void WeightAnalyzer::Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event){
    edm::Handle<std::vector<PileupSummaryInfo>> pileUp; 
    edm::Handle<GenEventInfoProduct> genInfo;
    std::vector<edm::Handle<double>> prefire;
    edm::Handle<std::vector<float>> pdfVariations, scaleVariations;
    prefireWeights.clear();
    
    if(!isNANO){
        event->getByToken(pileupToken, pileUp);
        event->getByToken(geninfoToken, genInfo);
        event->getByToken(pdfToken, pdfVariations);
        event->getByToken(scaleToken, scaleVariations);

        for(unsigned int i = 0; i < prefireTokens.size(); i ++){
            edm::Handle<double> fireHandle;
            if(era != 2018) event->getByToken(prefireTokens[i], fireHandle);
            prefire.push_back(fireHandle);
        }
    }

    //Set values if not data
    if(!this->isData){
        nTrueInt = isNANO ? *nPU->Get() : 1.;
        genWeight = isNANO ? *genWeightValue->Get() : genInfo->weight();

        if(!isNANO){
            if(!isSyst){
                //PDF uncertainties
                if(pdfVariations->size() == 102 and scaleVariations->size() == 8){
                    std::copy(scaleVariations->begin(), scaleVariations->begin()+8, scaleWeights);
                    std::copy(pdfVariations->begin(), pdfVariations->begin()+102, pdfWeights);
                }

                //Fill prefire weight
                for(unsigned int i = 0; i < 3; i ++){
                    prefireWeights.push_back(era != 2018 ? *(prefire[i]) : 1.);
                }
            }

            else{
                prefireWeights.push_back(era != 2018 ? *(prefire[0]) : 1.);
            }

            //Get true number of interaction https://twiki.cern.ch/twiki/bin/view/CMS/PileupSystematicErrors
            for(PileupSummaryInfo puInfo: *pileUp) {
                const int& BX = puInfo.getBunchCrossing();

                if(BX == 0) { 
                    nTrueInt = puInfo.getTrueNumInteractions();
                    continue;
                }
            }
        }

        nGenHist->Fill(1);
        nGenWeightedHist->Fill(genWeight);
        puMC->Fill(nTrueInt);
    }

    eventNumber = isNANO ? *evtNumber->Get() : event->id().event();

    for(CutFlow& cutflow: cutflows){
        cutflow.weight = isData ? 1. : xSec*lumi;
        cutflow.hist->Fill("No cuts", cutflow.weight);
    }
}

void WeightAnalyzer::EndJob(TFile* file){
    if(!this->isData){
        if(!file->GetListOfKeys()->Contains("nGen")){
            nGenHist->Write();
            nGenWeightedHist->Write();
            puMC->Write();

            //Also get measured PU distributions https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
            for(const std::string& syst: {"", "Up", "Down"}){
                std::string name = pileUpFile;
                name.replace(name.find("@"), 1, syst);

                std::shared_ptr<TFile> pileFile = std::make_shared<TFile>(name.c_str(), "READ");
                std::shared_ptr<TH1F> realPile(static_cast<TH1F*>(pileFile->Get("pileup")));
                realPile->SetName(("pileUp" + syst).c_str());
                file->cd();
                realPile->Write();
            }
        }

        std::shared_ptr<TH1F> forXSec = std::make_shared<TH1F>("xSec", "xSec", 1, 0, 1);
        forXSec->SetBinContent(1, xSec);
        forXSec->Write();

        std::shared_ptr<TH1F> forLumi = std::make_shared<TH1F>("Lumi", "Lumi", 1, 0, 1);
        forLumi->SetBinContent(1, lumi);
        forLumi->Write();
    }
}
