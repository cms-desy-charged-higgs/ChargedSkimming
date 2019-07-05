#include <ChargedHiggs/Skimming/interface/nanoskimmer.h>

#include <ChargedHiggs/Skimming/interface/triggeranalyzer.h>
#include <ChargedHiggs/Skimming/interface/metfilteranalyzer.h>
#include <ChargedHiggs/Skimming/interface/electronanalyzer.h>
#include <ChargedHiggs/Skimming/interface/muonanalyzer.h>
#include <ChargedHiggs/Skimming/interface/jetanalyzer.h>
#include <ChargedHiggs/Skimming/interface/genpartanalyzer.h>
#include <ChargedHiggs/Skimming/interface/weightanalyzer.h>

NanoSkimmer::NanoSkimmer(){}

NanoSkimmer::NanoSkimmer(const std::string &inFile, const bool &isData):
    inFile(inFile),
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

void NanoSkimmer::Configure(const float &xSec, TTreeReader& reader){
    analyzerMap = {
        {"e4j", 
                {
                    std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, reader)),
                    std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"}, reader)),
                    std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, reader)),
                    std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, {{4,0}}, reader)),
                    std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 25., 2.4, 1, reader)),
                    std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 20., 2.4, 0, reader)),
                    std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(reader)),
                }
        },

        {"mu4j", 
                {
                    std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, reader)),
                    std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_IsoMu27"}, reader)),
                    std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, reader)),
                    std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, {{4,0}}, reader)),
                    std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 25., 2.4, 1, reader)),
                    std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 20., 2.4, 0, reader)),
                    std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(reader)),
                }
        },

        {"e2j1f", 
                {
                    std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, reader)),
                    std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"}, reader)),
                    std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, reader)),
                    std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, {{2,1}}, reader)),
                    std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 25., 2.4, 1, reader)),
                    std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 20., 2.4, 0, reader)),
                    std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(reader)),
                }
        },

        {"mu2j1f", 
                {
                    std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, reader)),
                    std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_IsoMu27"}, reader)),
                    std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, reader)),
                    std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, {{2,1}}, reader)),
                    std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 25., 2.4, 1, reader)),
                    std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 20., 2.4, 0, reader)),
                    std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(reader)),
                }
        },

        {"e2f", 
                {
                    std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, reader)),
                    std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"}, reader)),
                    std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, reader)),
                    std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, {{0,2}}, reader)),
                    std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 25., 2.4, 1, reader)),
                    std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 20., 2.4, 0, reader)),
                    std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(reader)),
                }
        },

        {"mu2f", 
                {
                    std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, reader)),
                    std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_IsoMu27"}, reader)),
                    std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, reader)),
                    std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, {{0,2}}, reader)),
                    std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 25., 2.4, 1, reader)),
                    std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 20., 2.4, 0, reader)),
                    std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(reader)),
                }
        },
    }; 
}

void NanoSkimmer::EventLoop(const std::vector<std::string> &channels, const float &xSec){

    //TTreeReader preperation
    TFile* file = TFile::Open(inFile.c_str(), "READ");
    TTree* eventTree = (TTree*)file->Get("Events");
    TTreeReader reader(eventTree);

    Configure(xSec, reader);

    for(const std::string &channel: channels){
        //Create output trees
        TTree* tree = new TTree();
        tree->SetName(channel.c_str());
        trees.push_back(tree);

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
            analyzer->BeginJob(trees[i], isData);
        }
    }

    //Progress bar at 0%
    bool eventPassed = true;
    int processed = 0;
    ProgressBar(0.);

    while(reader.Next()){
        for(unsigned int i = 0; i < analyzers.size(); i++){
            //Call each analyzer and break if one analyzer reject the event
            for(std::shared_ptr<BaseAnalyzer> analyzer: analyzers[i]){
                if(!analyzer->Analyze(cutflows[i])){
                    eventPassed = false;
                    break;
                }
            }
            
            //Fill trees is event passed all analyzers
            if(eventPassed){
                 trees[i]->Fill();
            }

            eventPassed = true;
        }
        
        //progress bar
        processed++;
        if(processed % 10000 == 0){
            int progress = 100*(float)processed/eventTree->GetEntries();
            ProgressBar(progress);        
        }
    }

    ProgressBar(100);

    //Print stats
    for(TTree* tree: trees){
        std::cout << tree->GetName() << " analysis: Selected " << tree->GetEntries() << " events of " << eventTree->GetEntries() << " (" << 100*(float)tree->GetEntries()/eventTree->GetEntries() << "%)" << std::endl;
    }

    file->Close();
}

void NanoSkimmer::WriteOutput(const std::string &outFile){
    TFile* file = TFile::Open(outFile.c_str(), "RECREATE");
    
    for(TTree* tree: trees){
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

    end = std::chrono::steady_clock::now();
    std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;

    std::cout << "Output file created: " + outFile << std::endl;
}
