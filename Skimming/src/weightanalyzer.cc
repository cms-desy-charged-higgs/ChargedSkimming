#include <ChargedSkimming/Skimming/interface/weightanalyzer.h>

WeightAnalyzer::WeightAnalyzer(const float era, const float xSec, TTreeReader &reader):
    BaseAnalyzer(&reader),
    era(era),
    xSec(xSec)
    {}

WeightAnalyzer::WeightAnalyzer(const float era, const float xSec, puToken &pileupToken, genToken &geninfoToken, const std::vector<edm::EDGetTokenT<double>> prefireTokens, wgtToken& pdfToken, wgtToken& scaleToken):
    BaseAnalyzer(),
    era(era),
    xSec(xSec),
    pileupToken(pileupToken),
    geninfoToken(geninfoToken),
    prefireTokens(prefireTokens),
    scaleToken(scaleToken),
    pdfToken(pdfToken)
    {}

void WeightAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    //Set lumi map
    lumis = {{2016, 35.92*1e3}, {2017, 41.53*1e3}};

    //Set data bool
    this->isData = isData;
    this->isSyst = isSyst;

    if(!this->isData){
        if(isNANO){
            //Initiliaze TTreeReaderValues
            genWeightValue = std::make_unique<TTreeReaderValue<float>>(*reader, "Generator_weight");
            nPU = std::make_unique<TTreeReaderValue<float>>(*reader, "Pileup_nTrueInt");
        }        

        puMC = new TH1F("puMC", "puMC", 100, 0, 100);
        nGenHist = new TH1F("nGen", "nGen", 100, 0, 2);
        nGenWeightedHist = new TH1F("nGenWeighted", "nGenWeighted", 100, -1e7, 1e7);
    }

    evtNumber = std::make_unique<TTreeReaderValue<ULong64_t>>(*reader, "event");

    prefireWeights = {1., 1., 1.};

    //Branches for output tree
    for(TTree* tree: trees){
        tree->Branch("Weight_lumi", &lumi);
        tree->Branch("Weight_xsec", &xSec);
        tree->Branch("Weight_genWeight", &genWeight);
        tree->Branch("Weight_prefireWeight", &prefireWeights[0]);
        tree->Branch("Misc_TrueInteraction", &nTrueInt);
        tree->Branch("Misc_eventNumber", &eventNumber);

        if(!isSyst){
            tree->Branch("Weight_pdfVariations", &pdfWeights);
            tree->Branch("Weight_scaleVariations", &scaleWeights);
            tree->Branch("Weight_prefireWeightUp", &prefireWeights[1]);
            tree->Branch("Weight_prefireWeightDown", &prefireWeights[2]);
        }
    }
}

void WeightAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
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
            event->getByToken(prefireTokens[i], fireHandle);
            prefire.push_back(fireHandle);
        }
    }

    //Set values if not data
    if(!this->isData){
        lumi = lumis[era];
        nTrueInt = isNANO ? *nPU->Get() : 1.;
        genWeight = isNANO ? *genWeightValue->Get() : genInfo->weight();

        if(!isNANO){
            if(!isSyst){
                //PDF uncertainties
                scaleWeights = *scaleVariations;
                pdfWeights = *pdfVariations;

                //Fill prefire weight
                for(unsigned int i = 0; i < prefire.size(); i ++){
                    prefireWeights.push_back(*(prefire[i]));
                }
            }

            else{
                prefireWeights.push_back(*(prefire[0]));
            }

            //Get true number of interaction https://twiki.cern.ch/twiki/bin/view/CMS/PileupSystematicErrors
            for(PileupSummaryInfo puInfo: *pileUp) {
                int BX = puInfo.getBunchCrossing();

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
        cutflow.weight = xSec*lumi;
        cutflow.hist->Fill("No cuts", cutflow.weight);
        cutflow.passed *= true;
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
                std::string filename = filePath + "/pileUp/pileUp@!.root";
                filename.replace(filename.find("@"), 1, std::to_string(int(era)));
                filename.replace(filename.find("!"), 1, syst);

                TFile* pileFile = TFile::Open(filename.c_str(), "READ");
                TH1F* realPile = (TH1F*)pileFile->Get("pileup");
                realPile->SetName(("pileUp" + syst).c_str());
                file->cd();
                realPile->Write();

                delete realPile;
                delete pileFile;
            }
        }

        delete puMC;
        delete nGenHist;
        delete nGenWeightedHist;
    }
}
