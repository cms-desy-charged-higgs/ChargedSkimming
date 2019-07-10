#include <ChargedHiggs/Skimming/interface/miniskimmer.h>

MiniSkimmer::MiniSkimmer(const edm::ParameterSet& iConfig):
      //Tokens
      jetToken(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      fatjetToken(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatjets"))),
      genjetToken(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"))),
      genfatjetToken(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genfatjets"))),
      metToken(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
      eleToken(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      muonToken(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigger"))),
      pileupToken(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileUp"))),
      geninfoToken(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
      rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),

      //Other stuff
      channels(iConfig.getParameter<std::vector<std::string>>("channels")),
      xSec(iConfig.getParameter<double>("xSec")),
      outFile(iConfig.getParameter<std::string>("outFile")),
      isData(iConfig.getParameter<bool>("isData")){

        start = std::chrono::steady_clock::now();
}

MiniSkimmer::~MiniSkimmer(){
    end = std::chrono::steady_clock::now();
    std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;

    //Print stats
    for(TTree* tree: outputTrees){
        std::cout << tree->GetName() << " analysis: Selected " << tree->GetEntries() << " events of " << nEvents << " (" << 100*(float)tree->GetEntries()/nEvents << "%)" << std::endl;
    }
}

void MiniSkimmer::beginJob(){
    //Set analyzer modules for each final state
    std::vector<jToken> jetTokens = {jetToken, fatjetToken};
    std::vector<genjToken> genjetTokens = {genjetToken, genfatjetToken};

    analyzerMap = {
        {"e4j", 
                {
                    std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, pileupToken, geninfoToken)),
                    std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"}, triggerToken)),
                    std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, triggerToken)),
                    std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, {{4,0}}, jetTokens, genjetTokens, metToken, rhoToken)),
                    std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 25., 2.4, 1, eleToken)),
                    std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 20., 2.4, 0, muonToken)),
                }
        },
    }; 

    for(const std::string &channel: channels){
        //Create output trees
        TTree* tree = new TTree();
        tree->SetName(channel.c_str());
        outputTrees.push_back(tree);

        //Create cutflow histograms
        TH1F* cutflow = new TH1F();
        cutflow->SetName(("cutflow_" + channel).c_str());
        cutflow->SetName(("cutflow_" + channel).c_str());
        cutflow->GetYaxis()->SetName("Events");
        cutflows.push_back(std::make_pair(cutflow, 1.)); 

        //Push back analyzers for each final state 
        analyzers.push_back(analyzerMap[channel]);
    }

    //Begin jobs for all analyzers
    for(unsigned int i = 0; i < analyzers.size(); i++){
        for(std::shared_ptr<BaseAnalyzer> analyzer: analyzers[i]){
            analyzer->BeginJob(outputTrees[i], isData);
        }
    }
}

void MiniSkimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    nEvents++;
    bool eventPassed = true;

    for(unsigned int i = 0; i < analyzers.size(); i++){
        //Call each analyzer and break if one analyzer reject the event
        for(std::shared_ptr<BaseAnalyzer> analyzer: analyzers[i]){
            if(!analyzer->Analyze(cutflows[i], &iEvent)){
                eventPassed = false;
                break;
            }
        }
            
        //Fill trees is event passed all analyzers
        if(eventPassed){
            outputTrees[i]->Fill();
        }
    }
}

void MiniSkimmer::endJob(){
    TFile* file = TFile::Open(outFile.c_str(), "RECREATE");
    
    for(TTree* tree: outputTrees){
        tree->Write();
    }

    //End jobs for all analyzers
    for(unsigned int i = 0; i < analyzers.size(); i++){
        for(std::shared_ptr<BaseAnalyzer> analyzer: analyzers[i]){
            analyzer->EndJob(file);
        }

        cutflows[i].first->Write();
    }

    file->Write();
    file->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniSkimmer);
