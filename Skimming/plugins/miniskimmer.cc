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
      isoTrackToken(consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("isoTrack"))),
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
      era(iConfig.getParameter<int>("era")),
      isData(iConfig.getParameter<bool>("isData")){

        start = std::chrono::steady_clock::now();
}

MiniSkimmer::~MiniSkimmer(){
    end = std::chrono::steady_clock::now();
    std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;
}

void MiniSkimmer::beginJob(){
    //Read in json config
    boost::property_tree::ptree parser; 
    boost::property_tree::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/skim.json", parser);

    //Set analyzer modules for each final state
    std::vector<jToken> jetTokens = {jetToken, fatjetToken};
    std::vector<genjToken> genjetTokens = {genjetToken, genfatjetToken};
    std::vector<edm::EDGetTokenT<double>> prefireTokens = {prefireToken, prefireTokenUp, prefireTokenDown};

    for(std::pair<std::string, boost::property_tree::ptree> syst : parser.get_child("Systematics")){
        std::string systematic = syst.first;

        //Dont do systematics for data
        if(isData and systematic != "") continue;

        for(const std::string& shift: {"Up", "Down"}){
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
            std::cout << "Open output file: " << outname << std::endl;

            //Vector of objects (One object for each channel)
            std::vector<TTree*> treesPerSyst;
            std::vector<CutFlow> flowPerSyst;
            std::vector<std::shared_ptr<BaseAnalyzer>> analyzerPerSyst;

            for(const std::string &channel: channels){
                //Check channel has relevant particle for the systematic
                bool useSystematic = false;
                for(std::pair<std::string, boost::property_tree::ptree> channelName: syst.second){
                    if(channel == channelName.second.get_value<std::string>()) useSystematic = true;
                }
    
                if(!useSystematic) continue;

                //Create output trees
                TTree* tree = new TTree(channel.c_str(), channel.c_str());
                tree->SetAutoFlush(1000);
                treesPerSyst.push_back(tree);

                //Create cutflow histograms
                CutFlow cutflow;

                cutflow.channel = channel;
                cutflow.hist = std::make_shared<TH1F>();
                cutflow.hist->SetName(("cutflow_" + channel).c_str());
                cutflow.hist->SetName(("cutflow_" + channel).c_str());
                cutflow.hist->GetYaxis()->SetName("Events");
    
                cutflow.nMinMu = parser.get<int>("Channel." + channel + ".Selection.nMinMuon");
                cutflow.nMinEle = parser.get<int>("Channel." + channel + ".Selection.nMinElectron");
                cutflow.nMinJet = parser.get<int>("Channel." + channel + ".Selection.nMinJet");
                cutflow.nMinFatjet =parser.get<int>("Channel." + channel + ".Selection.nMinFatJet");
                
                cutflow.weight = 1;    

                flowPerSyst.push_back(cutflow);
            }

            //Get information for analyzers
            float jetPt = parser.get<float>("Analyzer.Jet.pt." + std::to_string(era));
            float jetEta = parser.get<float>("Analyzer.Jet.eta." + std::to_string(era));
            float elePt = parser.get<float>("Analyzer.Electron.pt." + std::to_string(era));
            float eleEta = parser.get<float>("Analyzer.Electron.eta." + std::to_string(era));
            float muonPt = parser.get<float>("Analyzer.Muon.pt." + std::to_string(era));
            float muonEta = parser.get<float>("Analyzer.Muon.eta." + std::to_string(era));

            std::map<std::string, std::vector<std::string>> triggers;

            for(const std::string &channel: channels){
                triggers[channel] = Util::GetVector<std::string>(parser, "Channel." + channel + ".Trigger." + std::to_string(era));
            }

            std::vector<std::string> genParts = Util::GetVector<std::string>(parser, "Analyzer.Gen");

            //Set Analyzer
            analyzerPerSyst = {
                std::make_shared<WeightAnalyzer>(era, xSec, pileupToken, geninfoToken, prefireTokens, pdfToken, scaleToken),
                std::make_shared<TriggerAnalyzer>(triggers, triggerToken),
                std::make_shared<MetFilterAnalyzer>(era, triggerToken),
                std::make_shared<JetAnalyzer>(era, jetPt, jetEta, jetTokens, genjetTokens, metToken, rhoToken, genParticleToken, secVertexToken, systName),
                std::make_shared<MuonAnalyzer>(era, muonPt, muonEta, muonToken, genParticleToken),
                std::make_shared<ElectronAnalyzer>(era, elePt, eleEta, eleToken, genParticleToken, systName),
                std::make_shared<MiscAnalyzer>(era, std::max({jetEta, muonEta, eleEta}), isoTrackToken),
                std::make_shared<GenPartAnalyzer>(genParticleToken, genParts),
            };

            //Begin jobs for all analyzers
            for(std::shared_ptr<BaseAnalyzer>& analyzer: analyzerPerSyst){
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

        std::cout << "Overview for closed output file: " << outputFiles[i]->GetName() << std::endl;

        for(TTree* tree: outputTrees[i]){
            tree->Write("", TObject::kWriteDelete);
            std::cout << "Output tree '" << tree->GetName() << "': Selected " << tree->GetEntries() << " events of " << nEvents << " (" << 100*(float)tree->GetEntries()/nEvents << "%)" << std::endl;
        }

        //End jobs for all analyzers
        for(unsigned int j = 0; j < analyzers[i].size(); j++){
            analyzers[i][j]->EndJob(outputFiles[i]);
        }

        for(CutFlow& cutflow: cutflows[i]){
            cutflow.hist->Write();
        }

        outputFiles[i]->Close();
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniSkimmer);
