#include <ChargedSkimming/Skimming/interface/miniskimmer.h>

MiniSkimmer::MiniSkimmer(const edm::ParameterSet& iConfig){
    //Tokens
    tokens = std::make_shared<Token>();

    //Jet related edm token
    tokens->jetToken = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
    tokens->fatJetToken = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatjets"));
    tokens->genJetToken = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"));
    tokens->genFatJetToken = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genfatjets"));
    tokens->METToken = consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"));
    tokens->secVertexToken = consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("svtx")); 

    //Electron related token
    tokens->eleToken = consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"));

    //Muon related token
    tokens->muonToken = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));

    //Gen related token
    tokens->genToken = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"));
    tokens->genPartToken = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPart"));
    tokens->lheToken = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lhe"));

    //PDF related token
    tokens->pdfToken = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("pdf"));
    tokens->scaleToken = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("scale"));

    //Weight token
    tokens->prefireToken = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"));
    tokens->prefireTokenUp = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
    tokens->prefireTokenDown = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
    tokens->pileUpToken = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileUp"));

    //Miscellaneous
    tokens->isoTrackToken = consumes<std::vector<pat::IsolatedTrack>>(iConfig.getParameter<edm::InputTag>("isoTrack"));
    tokens->triggerToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigger"));
    tokens->rhoToken = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));

    //Other stuff
    channels = iConfig.getParameter<std::vector<std::string>>("channels");
    xSec = iConfig.getParameter<double>("xSec");
    outFile = iConfig.getParameter<std::string>("outFile");
    era = iConfig.getParameter<std::string>("era");
    isData = iConfig.getParameter<bool>("isData");

    start = std::chrono::steady_clock::now();
}

MiniSkimmer::~MiniSkimmer(){
    end = std::chrono::steady_clock::now();
    std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;
}

void MiniSkimmer::beginJob(){
    //Read in json configs
    pt::ptree sf, skim; 
    pt::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/UL/skim.json", skim);
    pt::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/UL/sf.json", sf);

    //Set analyzer modules for each final state
    /*for(std::pair<std::string, boost::property_tree::ptree> syst : skim.get_child("Systematics")){
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
    
                cutflow.nMinMu = skim.get<int>("Channel." + channel + ".Selection.nMinMuon");
                cutflow.nMinEle = skim.get<int>("Channel." + channel + ".Selection.nMinElectron");
                cutflow.nMinJet = skim.get<int>("Channel." + channel + ".Selection.nMinJet");
                cutflow.nMinFatjet = skim.get<int>("Channel." + channel + ".Selection.nMinFatJet");
                
                cutflow.weight = 1;    

                flowPerSyst.push_back(cutflow);
            }

            skim.put<bool>("isData", isData);
            skim.put<bool>("isSyst", !isNomi);

            //Set Analyzer
            analyzerPerSyst = {
                std::make_shared<WeightAnalyzer>(era, xSec, tokens),
                std::make_shared<TriggerAnalyzer>(era, channels, tokens),
                std::make_shared<MetFilterAnalyzer>(era, tokens),
                std::make_shared<JetAnalyzer>(era, tokens, systName),
                std::make_shared<MuonAnalyzer>(era, tokens),
                std::make_shared<ElectronAnalyzer>(era, tokens, systName),
                std::make_shared<MiscAnalyzer>(era, tokens),
                std::make_shared<GenPartAnalyzer>(tokens),
            };

            //Begin jobs for all analyzers
            for(std::shared_ptr<BaseAnalyzer>& analyzer: analyzerPerSyst){
                analyzer->BeginJob(treesPerSyst, skim, sf);
            }

            outputTrees.push_back(treesPerSyst);
            cutflows.push_back(flowPerSyst);
            analyzers.push_back(analyzerPerSyst);
            outputFiles.push_back(outFile);
        }
    }*/
}

void MiniSkimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    ++nEvents;

   /* for(unsigned int i = 0; i < analyzers.size(); ++i){
        //Call each analyzer
        for(unsigned int j = 0; j < analyzers[i].size(); ++j){
            analyzers[i][j]->Analyze(cutflows[i], &iEvent);
        }

        //Check individual for each channel, if event should be filled
        for(unsigned int j = 0; j < outputTrees[i].size(); ++j){
            if(cutflows[i][j].passed){
                outputTrees[i][j]->Fill();
            }

            cutflows[i][j].passed = true;
        }
    }*/
}

void MiniSkimmer::endJob(){    
    /*for(unsigned int i = 0; i < analyzers.size(); i++){
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
    }*/
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniSkimmer);
