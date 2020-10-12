#include <ChargedSkimming/Skimming/interface/nanoskimmer.h>

#include <ChargedSkimming/Skimming/interface/triggeranalyzer.h>
#include <ChargedSkimming/Skimming/interface/metfilteranalyzer.h>
#include <ChargedSkimming/Skimming/interface/electronanalyzer.h>
#include <ChargedSkimming/Skimming/interface/muonanalyzer.h>
#include <ChargedSkimming/Skimming/interface/jetanalyzer.h>
#include <ChargedSkimming/Skimming/interface/genpartanalyzer.h>
#include <ChargedSkimming/Skimming/interface/weightanalyzer.h>

NanoSkimmer::NanoSkimmer(){}

NanoSkimmer::NanoSkimmer(const std::string &inFile, const std::string &outFile, const std::vector<std::string>& channels, const float& xSec, const int& era, const bool &isData):
    inFile(inFile),
    outFile(outFile),
    channels(channels),
    xSec(xSec),
    era(era), 
    isData(isData)
    {    
        start = std::chrono::steady_clock::now();
        std::cout << "Input file for analysis: " + inFile << std::endl;
    }

void NanoSkimmer::ProgressBar(const int &progress){
    std::string progressBar = "["; 

    for(int i = 0; i < progress; i++){
        if(i%2 == 0) progressBar += "#";
    }

    for(int i = 0; i < 100 - progress; i++){
        if(i%2 == 0) progressBar += " ";
    }

    progressBar = progressBar + "] " + std::to_string(progress) + "% of Events processed";
    std::cout << "\r" << progressBar << std::flush;

    if(progress == 100) std::cout << std::endl;

}

void NanoSkimmer::Configure(TTreeReader& reader){
    //Read in json config
    boost::property_tree::ptree parser; 
    boost::property_tree::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/skim.json", parser);

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
                tree->SetAutoFlush(400000);
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
            float jetEta = parser.get<float>("Analyzer.Jet.pt." + std::to_string(era));
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
                  std::make_shared<WeightAnalyzer>(era, xSec, reader),
                  std::make_shared<TriggerAnalyzer>(triggers, reader),
                  std::make_shared<MetFilterAnalyzer>(era, reader),
                  std::make_shared<JetAnalyzer>(era, jetPt, jetEta, reader),
                  std::make_shared<MuonAnalyzer>(era, muonPt, muonEta, reader),
                  std::make_shared<ElectronAnalyzer>(era, elePt, eleEta, reader),
                  std::make_shared<GenPartAnalyzer>(reader, genParts),
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

void NanoSkimmer::EventLoop(){
    //TTreeReader preperation
    TFile* inputFile = TFile::Open(inFile.c_str(), "READ");
    TTree* eventTree = (TTree*)inputFile->Get("Events");
    TTreeReader reader(eventTree);

    Configure(reader);

    //Progress bar at 0%
    ProgressBar(0.);
    nEvents = 0;

    while(reader.Next()){
        for(unsigned int i = 0; i < analyzers.size(); i++){
            //Call each analyzer
            for(unsigned int j = 0; j < analyzers[i].size(); j++){
                analyzers[i][j]->Analyze(cutflows[i]);
            }

            //Check individual for each channel, if event should be filled
            for(unsigned int j = 0; j < outputTrees[i].size(); j++){
                if(cutflows[i][j].passed){
                    outputTrees[i][j]->Fill();
                }

                cutflows[i][j].passed = true;
 
            }
        }

        //progress bar
        nEvents++;
        if(nEvents % 10000 == 0){
            int progress = 100*(float)nEvents/eventTree->GetEntries();
            ProgressBar(progress);        
        }
    }

    ProgressBar(100);
}

void NanoSkimmer::WriteOutput(){
    for(unsigned int i = 0; i < analyzers.size(); i++){
        outputFiles[i]->cd();

        std::cout << "Overview for closed output file: " << outputFiles[i]->GetName() << std::endl;

        for(TTree* tree: outputTrees[i]){
            std::cout << "Output tree '" << tree->GetName() << "': Selected " << tree->GetEntries() << " events of " << nEvents << " (" << 100*(float)tree->GetEntries()/nEvents << "%)" << std::endl;
        }

        //End jobs for all analyzers
        for(unsigned int j = 0; j < analyzers[i].size(); j++){
            analyzers[i][j]->EndJob(outputFiles[i]);
        }

        outputFiles[i]->Write(0, TObject::kOverwrite);
        outputFiles[i]->Close();
    }

    end = std::chrono::steady_clock::now();
    std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;
}
