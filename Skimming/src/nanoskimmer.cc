#include <ChargedAnalysis/Skimming/interface/nanoskimmer.h>

#include <ChargedAnalysis/Skimming/interface/triggeranalyzer.h>
#include <ChargedAnalysis/Skimming/interface/metfilteranalyzer.h>
#include <ChargedAnalysis/Skimming/interface/electronanalyzer.h>
#include <ChargedAnalysis/Skimming/interface/muonanalyzer.h>
#include <ChargedAnalysis/Skimming/interface/jetanalyzer.h>
#include <ChargedAnalysis/Skimming/interface/genpartanalyzer.h>
#include <ChargedAnalysis/Skimming/interface/weightanalyzer.h>

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
    analyzers = {
        std::shared_ptr<WeightAnalyzer>(new WeightAnalyzer(2017, xSec, reader)),
//        std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_IsoMu27"}, triggerToken)),
  //      std::shared_ptr<TriggerAnalyzer>(new TriggerAnalyzer({"HLT_Ele35_WPTight_Gsf", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150", "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"}, triggerToken)),
        std::shared_ptr<MetFilterAnalyzer>(new MetFilterAnalyzer(2017, reader)),
        std::shared_ptr<JetAnalyzer>(new JetAnalyzer(2017, 30., 2.4, reader)),
        std::shared_ptr<MuonAnalyzer>(new MuonAnalyzer(2017, 25., 2.4, reader)),
        std::shared_ptr<ElectronAnalyzer>(new ElectronAnalyzer(2017, 20., 2.4, reader)),
        std::shared_ptr<GenPartAnalyzer>(new GenPartAnalyzer(reader))
    };
}

void NanoSkimmer::EventLoop(const std::vector<std::string> &channels, const float &xSec){

    //TTreeReader preperation
    TFile* inputFile = TFile::Open(inFile.c_str(), "READ");
    TTree* eventTree = (TTree*)inputFile->Get("Events");
    TTreeReader reader(eventTree);

    Configure(xSec, reader);

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

    //Begin jobs for all analyzers
    for(std::shared_ptr<BaseAnalyzer> analyzer: analyzers){
        analyzer->BeginJob(outputTrees, isData);
    }

    //Progress bar at 0%
    int processed = 0;
    ProgressBar(0.);

    while(reader.Next()){
        //Call each analyzer
        for(unsigned int i = 0; i < analyzers.size(); i++){
            unsigned int nFailed = 0;
            analyzers[i]->Analyze(cutflows);

            //If for all channels one analyzer failes, reject event
        for(CutFlow &cutflow: cutflows){
                if(!cutflow.passed) nFailed++;
            }

            //If for all channels one analyzer failes, reject event
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
        
        //progress bar
        processed++;
        if(processed % 10000 == 0){
            int progress = 100*(float)processed/eventTree->GetEntries();
            ProgressBar(progress);        
        }
    }

    ProgressBar(100);

    //Print stats
    for(TTree* tree: outputTrees){
        std::cout << tree->GetName() << " analysis: Selected " << tree->GetEntries() << " events of " << eventTree->GetEntries() << " (" << 100*(float)tree->GetEntries()/eventTree->GetEntries() << "%)" << std::endl;
    }
}

void NanoSkimmer::WriteOutput(const std::string &outFile){
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

    end = std::chrono::steady_clock::now();
    std::cout << "Finished event loop (in seconds): " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << std::endl;

    std::cout << "Output file created: " + outFile << std::endl;
}
