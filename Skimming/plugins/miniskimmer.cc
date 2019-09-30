#include <ChargedAnalysis/Skimming/interface/miniskimmer.h>

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
      triggerObjToken(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
      pileupToken(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileUp"))),
      geninfoToken(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
      genParticleToken(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPart"))),
      rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
      vertexToken(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vtx"))),
      secVertexToken(consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("svtx"))),

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

    nMin = {
            {"mu4j", {1, 0, 4, 0}},
            {"e4j", {0, 1, 4, 0}},
            {"mu2j1f", {1, 0, 2, 1}},
            {"e2j1f", {0, 1, 2, 1}},
            {"mu2f", {1, 0, 0, 2}},
            {"e2f", {0, 1, 0, 2}},
    };

    for(const std::string &channel: channels){
        //Create output trees
        TTree* tree = new TTree();
        tree->SetName(channel.c_str());
        outputTrees.push_back(tree);

        //Create cutflow histograms
        CutFlow cutflow;

        cutflow.hist = new TH1F();
        cutflow.hist->SetName(("cutflow_" + channel).c_str());
        cutflow.hist->SetName(("cutflow_" + channel).c_str());
        cutflow.hist->GetYaxis()->SetName("Events");

        cutflow.nMinMu=nMin[channel][0];
        cutflow.nMinEle=nMin[channel][1];
        cutflow.nMinJet=nMin[channel][2];
        cutflow.nMinFatjet=nMin[channel][3];
        
        cutflow.weight = 1;    

        cutflows.push_back(cutflow); 
    }

    analyzers = {
        std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, pileupToken, geninfoToken)),
        std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_IsoMu27"}, {"HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"}, triggerToken)),
        std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, triggerToken)),
        std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, jetTokens, genjetTokens, metToken, rhoToken, genParticleToken, secVertexToken)),
        std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 20., 2.4, muonToken, triggerObjToken, genParticleToken)),
        std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 20., 2.4, eleToken, triggerObjToken, genParticleToken)),
        std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(genParticleToken)),
    };

    //Begin jobs for all analyzers
    for(std::shared_ptr<BaseAnalyzer> analyzer: analyzers){
        analyzer->BeginJob(outputTrees, isData);
    }
}

void MiniSkimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    nEvents++;
    unsigned int nFailed = 0;

    //Call each analyzer
    for(unsigned int i = 0; i < analyzers.size(); i++){
        nFailed = 0;
        analyzers[i]->Analyze(cutflows, &iEvent);

        for(CutFlow &cutflow: cutflows){
            if(!cutflow.passed) nFailed++;
        }

        //If for all channels one analyzer fails, reject event
        if(nFailed == cutflows.size()){
            break;
        }        
    }

    //Check individual for each channel, if event should be filled
    for(unsigned int i = 0; i < outputTrees.size(); i++){
        if(cutflows[i].passed){
            outputTrees[i]->Fill();
        }

        cutflows[i].passed = true;
    }
}

void MiniSkimmer::endJob(){
    TFile* file = TFile::Open(outFile.c_str(), "RECREATE");
    
    for(TTree* tree: outputTrees){
        tree->Write();
    }

    //End jobs for all analyzers
    for(unsigned int i = 0; i < analyzers.size(); i++){
        analyzers[i]->EndJob(file);
    }

    for(CutFlow& cutflow: cutflows){
        cutflow.hist->Write();
        delete cutflow.hist;
    }

    file->Write();
    file->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniSkimmer);
