#include <ChargedSkimming/Skimming/interface/miniskimmer.h>

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
      genParticleToken(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPart"))),
      rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
      secVertexToken(consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("svtx"))), 
      prefireToken(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
      prefireTokenUp(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
      prefireTokenDown(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
      pdfToken(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("pdf"))),
      scaleToken(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("scale"))),

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
    for(unsigned int i = 0; i < analyzers.size(); i++){
        for(TTree* tree: outputTrees[i]){
            std::cout << tree->GetName() << " analysis: Selected " << tree->GetEntries() << " events of " << nEvents << " (" << 100*(float)tree->GetEntries()/nEvents << "%)" << std::endl;
        }
    }
}

void MiniSkimmer::beginJob(){
    //Set analyzer modules for each final state
    std::vector<jToken> jetTokens = {jetToken, fatjetToken};
    std::vector<genjToken> genjetTokens = {genjetToken, genfatjetToken};

    nMin = {
            {"MuonIncl", {1, 0, 0, 0}},
            {"EleIncl", {0, 1, 0, 0}},
    };

    //Name of systematic uncertainties
    systNames = {
                {"", ""}, 
            //        {"energyScale", "e"}, 
            //        {"energySigma", "e"},
            //{"JECTotal", "incl"},
             //       {"JER", "incl"},
    };

    for(const std::pair<std::string, std::string>& systInfo: systNames){
        //Dont do systematics for data
        if(isData and systInfo.first != "") continue;

        for(const std::string& shift: {"Up", "Down"}){
            //Get systematic name and relevant particle
            std::string systematic, particle;
            std::tie(systematic, particle) = systInfo;

            //Check if nominal analysis
            bool isNomi = systematic == "" ? true : false;
            std::string systName = isNomi ? "" : systematic + shift;

            //Skip Up Down loop for nominal for Down variation
            if(systematic == "" and shift == "Down") continue;

            //Open output file
            std::string outname = outFile;

            if(!isNomi){
                outname.insert(outname.find("."), "_" + systName);            
            } 

            TFile* outFile = TFile::Open(outname.c_str(), "RECREATE"); 

            //Vector of objects (One object for each channel)
            std::vector<TTree*> treesPerSyst;
            std::vector<CutFlow> flowPerSyst;
            std::vector<std::shared_ptr<BaseAnalyzer>> analyzerPerSyst;

            for(const std::string &channel: channels){
                //Check channel has relevant particle for the systematic
                if(channel.find(particle) == std::string::npos) continue;

                //Create output trees
                TTree* tree = new TTree();
                tree->SetName(channel.c_str());
                tree->SetAutoSave(0);
                tree->SetAutoFlush(0);
                treesPerSyst.push_back(tree);

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

                flowPerSyst.push_back(cutflow); 
            }

            analyzerPerSyst = {
                std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, pileupToken, geninfoToken, {prefireToken, prefireTokenUp, prefireTokenDown}, pdfToken, scaleToken)),
                std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_IsoMu27"}, {"HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"}, triggerToken)),
                std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, triggerToken)),
                std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, jetTokens, genjetTokens, metToken, rhoToken, genParticleToken, secVertexToken, particle == "j" ? systName : "")),
                std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 15., 2.4, muonToken, genParticleToken)),
                std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 15., 2.4, eleToken, genParticleToken, particle == "e" ? systName : "")),
                std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(genParticleToken)),
            };

            //Begin jobs for all analyzers
            for(std::shared_ptr<BaseAnalyzer> analyzer: analyzerPerSyst){
                analyzer->BeginJob(treesPerSyst, isData, !isNomi);
            }

            outputTrees.push_back(treesPerSyst);
            cutflows.push_back(flowPerSyst);
            analyzers.push_back(analyzerPerSyst);
            outputFiles.push_back(outFile);
        }
    }
}

void MiniSkimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    nEvents++;

    for(unsigned int i = 0; i < analyzers.size(); i++){
        //Call each analyzer
        for(unsigned int j = 0; j < analyzers[i].size(); j++){
            analyzers[i][j]->Analyze(cutflows[i], &iEvent);
        }

        //Check individual for each channel, if event should be filled
        for(unsigned int j = 0; j < outputTrees[i].size(); j++){
            if(cutflows[i][j].passed){
                outputTrees[i][j]->Fill();
            }

            cutflows[i][j].passed = true;
        }
    }
}

void MiniSkimmer::endJob(){    
    for(unsigned int i = 0; i < analyzers.size(); i++){
        outputFiles[i]->cd();

        for(TTree* tree: outputTrees[i]){
            tree->Write();
        }

        //End jobs for all analyzers
        for(unsigned int j = 0; j < analyzers[i].size(); j++){
            analyzers[i][j]->EndJob(outputFiles[i]);
        }

        for(CutFlow& cutflow: cutflows[i]){
            cutflow.hist->Write();
            delete cutflow.hist;
        }

        outputFiles[i]->Write();
        outputFiles[i]->Close();
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniSkimmer);
