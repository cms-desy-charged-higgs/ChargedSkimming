#include <ChargedSkimming/Skimming/interface/nanoskimmer.h>

#include <ChargedSkimming/Skimming/interface/triggeranalyzer.h>
#include <ChargedSkimming/Skimming/interface/metfilteranalyzer.h>
#include <ChargedSkimming/Skimming/interface/electronanalyzer.h>
#include <ChargedSkimming/Skimming/interface/muonanalyzer.h>
#include <ChargedSkimming/Skimming/interface/jetanalyzer.h>
#include <ChargedSkimming/Skimming/interface/genpartanalyzer.h>
#include <ChargedSkimming/Skimming/interface/weightanalyzer.h>

NanoSkimmer::NanoSkimmer(){}

NanoSkimmer::NanoSkimmer(const std::string &inFile, const std::string &outFile, const std::vector<std::string>& channels, const float& xSec, const std::string& era, const bool &isData):
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
    pt::ptree sf, skim; 
    pt::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/skim.json", skim);
    pt::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/sf.json", sf);

    for(std::pair<std::string, boost::property_tree::ptree> syst : skim.get_child("Systematics")){
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
    
                cutflow.nMinMu = skim.get<int>("Channel." + channel + ".Selection.nMinMuon");
                cutflow.nMinEle = skim.get<int>("Channel." + channel + ".Selection.nMinElectron");
                cutflow.nMinJet = skim.get<int>("Channel." + channel + ".Selection.nMinJet");
                cutflow.nMinFatjet = skim.get<int>("Channel." + channel + ".Selection.nMinFatJet");
                
                cutflow.weight = 1;    

                flowPerSyst.push_back(cutflow);
            }

            //Set Analyzer
            analyzerPerSyst = {
                std::make_shared<WeightAnalyzer>(era, xSec, reader),
                std::make_shared<TriggerAnalyzer>(era, channels, reader),
                std::make_shared<MetFilterAnalyzer>(era, reader),
                std::make_shared<JetAnalyzer>(era, reader),
                std::make_shared<MuonAnalyzer>(era, reader),
                std::make_shared<ElectronAnalyzer>(era, reader),
            //    std::make_shared<MiscAnalyzer>(era, reader),
                std::make_shared<GenPartAnalyzer>(reader),
            };

            //Begin jobs for all analyzers
            for(std::shared_ptr<BaseAnalyzer>& analyzer: analyzerPerSyst){
                analyzer->BeginJob(treesPerSyst, sf, skim);
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
