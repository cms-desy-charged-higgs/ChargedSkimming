#include <ChargedSkimming/Skimming/interface/weightanalyzer.h>

WeightAnalyzer::WeightAnalyzer(const std::string& era, const float& xSec, TTreeReader &reader):
    BaseAnalyzer(&reader),
    era(era),
    xSec(xSec)
    {}

WeightAnalyzer::WeightAnalyzer(const std::string& era, const float& xSec, const std::shared_ptr<Token>& tokens):
    BaseAnalyzer(),
    era(era),
    xSec(xSec),
    tokens(tokens)
    {}

void WeightAnalyzer::BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf){
    //Set data bool
    this->isData = skim.get<bool>("isData");
    this->isSyst = skim.get<bool>("isSyst");
    
    pileUpFile = filePath + skim.get<std::string>("Analyzer.Weight.PileUp." + era);
    lumi = skim.get<float>("Analyzer.Weight.Lumi.nominal." + era)*1e3;
    lumiUp = skim.get<float>("Analyzer.Weight.Lumi.Up." + era)*1e3;
    lumiDown = skim.get<float>("Analyzer.Weight.Lumi.Down." + era)*1e3;

    if(!this->isData){
        if(isNANO){
            //Initiliaze TTreeReaderValues
            genWeightValue = std::make_unique<TTreeReaderValue<float>>(*reader, "Generator_weight");
            nPU = std::make_unique<TTreeReaderValue<float>>(*reader, "Pileup_nTrueInt");
        }        

        puMC = std::make_shared<TH1F>("puMC", "puMC", 100, 0, 100);
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
    std::vector<PileupSummaryInfo> pileUp; 
    GenEventInfoProduct genInfo;
    std::vector<double> prefire;
    std::vector<float> pdfVariations, scaleVariations;
    prefireWeights.clear();
    
    if(!isNANO and !this->isData){
        pileUp = Token::GetTokenValue(event, tokens->pileUpToken);
        genInfo = Token::GetTokenValue(event, tokens->genToken);
        pdfVariations = Token::GetTokenValue(event, tokens->pdfToken);
        scaleVariations = Token::GetTokenValue(event, tokens->scaleToken);

        for(const edm::EDGetTokenT<double>& token : {tokens->prefireToken, tokens->prefireTokenUp, tokens->prefireTokenDown}){
            if(era != "2018") prefire.push_back(Token::GetTokenValue(event, token));
        }
    }

    //Set values if not data
    if(!this->isData){
        nTrueInt = isNANO ? *nPU->Get() : 1.;
        genWeight = isNANO ? *genWeightValue->Get() : genInfo.weight();

        if(!isNANO){
            if(!isSyst){
                //PDF uncertainties
                if(pdfVariations.size() == 102 and scaleVariations.size() == 8){
                    std::copy(scaleVariations.begin(), scaleVariations.begin()+8, scaleWeights);
                    std::copy(pdfVariations.begin(), pdfVariations.begin()+102, pdfWeights);
                }

                //Fill prefire weight
                for(unsigned int i = 0; i < 3; i ++){
                    prefireWeights.push_back(era != "2018" ? prefire[i] : 1.);
                }
            }

            else{
                prefireWeights.push_back(era != "2018" ? prefire[0] : 1.);
            }

            //Get true number of interaction https://twiki.cern.ch/twiki/bin/view/CMS/PileupSystematicErrors
            for(PileupSummaryInfo puInfo: pileUp) {
                const int& BX = puInfo.getBunchCrossing();

                if(BX == 0) { 
                    nTrueInt = puInfo.getTrueNumInteractions();
                    continue;
                }
            }
        }

        ++nGen;
        nGenWeighted += genWeight;
        puMC->Fill(nTrueInt);
    }

    eventNumber = isNANO ? *evtNumber->Get() : event->id().event();

    for(CutFlow& cutflow: cutflows){
        cutflow.weight = 1.;
        cutflow.hist->Fill("No cuts", cutflow.weight);
    }
}

void WeightAnalyzer::EndJob(TFile* file){
    if(!this->isData){
        if(!file->GetListOfKeys()->Contains("nGen")){
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

        TParameter<float> xSecV("xSec", xSec);
        xSecV.Write();

        TParameter<float> lumiV("Lumi", lumi);
        lumiV.Write();

        TParameter<float> lumiUpV("LumiUp", lumiUp);
        lumiUpV.Write();

        TParameter<float> lumiDownV("LumiDown", lumiDown);
        lumiDownV.Write();

        TParameter<float> nGenV("nGen", nGen);
        nGenV.Write();

        TParameter<float> nGenWeightedV("nGenWeighted", nGenWeighted);
        nGenWeightedV.Write();
    }
}
