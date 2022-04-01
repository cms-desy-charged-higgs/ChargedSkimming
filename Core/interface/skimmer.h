#ifndef SKIMMER_H
#define SKIMMER_H

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <experimental/filesystem>

#include <TFile.h>
#include <TTree.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <ChargedSkimming/Core/interface/output.h>
#include <ChargedSkimming/Core/interface/cuts.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>
#include <ChargedSkimming/Analyzer/interface/triggeranalyzer.h>
#include <ChargedSkimming/Analyzer/interface/metfilteranalyzer.h>
#include <ChargedSkimming/Analyzer/interface/electronanalyzer.h>
#include <ChargedSkimming/Analyzer/interface/muonanalyzer.h>
#include <ChargedSkimming/Analyzer/interface/jetanalyzer.h>
#include <ChargedSkimming/Analyzer/interface/weightanalyzer.h>
#include <ChargedSkimming/Analyzer/interface/isotrkanalyzer.h>
#include <ChargedSkimming/Analyzer/interface/miscanalyzer.h>
#include <ChargedSkimming/Analyzer/interface/sfanalyzer.h>

namespace pt = boost::property_tree;

template <typename T>
class Skimmer{
    private:
        //Output files/trees
        std::vector<std::shared_ptr<TFile>> outFiles;
        std::vector<std::shared_ptr<TTree>> outTrees;

        //Core classes used for skimming
        std::vector<Cuts> cuts;
        std::vector<std::shared_ptr<BaseAnalyzer<T>>> analyzer;

        //Input information
        std::vector<std::string> channels;
        std::string xSec, xSecUnc, era, run, systematic, shift;
        int nEvents;

    public:
        Skimmer(const std::vector<std::string>& channels, const std::string& xSec, const std::string& xSecUnc, const std::string& era, const std::string& run) : channels(channels), xSec(xSec), xSecUnc(xSecUnc), era(era), run(run) {}

        void Configure(T& input, Output& output, const std::string& outDir, const std::string& outFile){
            //Read in json config
            pt::ptree sf, skim; 
            pt::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/UL/skim.json", skim);
            pt::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/config/UL/sf.json", sf);

            //Trigger/METFilter names
            std::vector<std::string> triggerNames;

            for(const std::string& channel : channels){
                //Open output file and create trees
                std::string outD =  outDir;
                outD.replace(outDir.find("[C]"), 3, channel);
                std::experimental::filesystem::create_directories(outD);

                outFiles.push_back(std::make_shared<TFile>((outD + "/" + outFile).c_str(), "RECREATE")); 
                std::cout << "Open output file: " <<  outFiles.back()->GetName() << std::endl;

                outTrees.push_back(std::make_shared<TTree>(channel.c_str(), channel.c_str()));
                outTrees.back()->SetDirectory(outFiles.back().get());
                outTrees.back()->SetAutoFlush(10000);

                //Read out trigger needed and register in cut class
                for(const std::string& name : Util::GetVector<std::string>(skim, "Channel." + channel + ".Trigger." + era)){
                    if(std::find(triggerNames.begin(), triggerNames.end(), name) == triggerNames.end()) triggerNames.push_back(name);
                }

                //Define cut class
                cuts.push_back(Cuts(outFiles.back(), channel));

                //Register cut requirements
                std::string path = "Channel." + channel + ".Selection";

                for(const std::string part : Util::GetKeys(skim, path)){
                    cuts.back().AddCut(part, output, skim.get<std::string>(path + "." + part + ".operator"), skim.get<short>(path + "." + part + ".threshold"));
                }
            }

            //Register trigger to input class (Only thing which is needed to register for input class)
            input.SetTrigger(triggerNames, false);
            input.SetTrigger(Util::GetVector<std::string>(skim, "Analyzer.METFilter." + era), true);

            //Add trigger/METFilter to cuts
            for(std::size_t i = 0; i < channels.size(); ++i){
                std::vector<int> triggerIdx;

                for(std::size_t j = 0; j < triggerNames.size(); ++j){
                    for(const std::string& name : Util::GetVector<std::string>(skim, "Channel." + channels[i] + ".Trigger." + era)){
                        if(name == triggerNames[j]) triggerIdx.push_back(j);
                    }
                }
       
                //Input class instead output class is used to register cut for trigger!
                cuts[i].AddTrigger<T>({}, input);
                cuts[i].AddTrigger<T>(triggerIdx, input);
            }

            //Add. information for analyzer
            bool isData = run != "MC";

            skim.put<std::string>("xSec", xSec);
            skim.put<std::string>("xSecUnc", xSecUnc);
            skim.put<std::string>("run", run);
            skim.put<std::string>("era", era);
            skim.put<bool>("isData", isData);

            //Register branches to output trees
            output.RegisterTrigger(triggerNames, outTrees);
            output.Register("Weight", outTrees, skim, isData);
            output.Register("Electron", outTrees, skim, isData);
            output.Register("Muon", outTrees, skim, isData);
            output.Register("Jet", outTrees, skim, isData);
            output.Register("Isotrack", outTrees, skim, isData);
            output.Register("Misc", outTrees, skim, isData);

            //List of analyzer
            analyzer = {
                std::make_shared<TriggerAnalyzer<T>>(),
                std::make_shared<METFilterAnalyzer<T>>(),
                std::make_shared<JetAnalyzer<T>>(),
                std::make_shared<ElectronAnalyzer<T>>(),
                std::make_shared<MuonAnalyzer<T>>(),
                std::make_shared<IsotrkAnalyzer<T>>(),
                std::make_shared<WeightAnalyzer<T>>(),
                std::make_shared<MiscAnalyzer<T>>(),
                std::make_shared<SFAnalyzer<T>>(),
            };

            //Initialize analyzers
            for(std::shared_ptr<BaseAnalyzer<T>>& a : analyzer){
                a->BeginJob(skim, sf);
            }

            nEvents = input.GetEntries();
        };

        void Loop(T& input, Output& output){
            for(std::shared_ptr<BaseAnalyzer<T>>& a : analyzer){
                a->Analyze(input, output);
            }

            for(std::size_t i = 0; i < outTrees.size(); ++i){
                cuts[i].Count();
                cuts[i].FillCutflow();
                if(cuts[i].Passed()){
                    outTrees[i]->Fill();
                }
            }
        }

        void WriteOutput(){
            for(std::size_t i = 0; i < outTrees.size(); ++i){
                for(std::shared_ptr<BaseAnalyzer<T>>& a : analyzer){
                    a->EndJob(outFiles[i]);
                }

                outFiles[i]->cd();
                outTrees[i]->Write();
                cuts[i].WriteOutput();

                std::cout << "Close output file: " << outFiles[i]->GetName() << std::endl;
                std::cout << "Written Tree: " << outTrees[i]->GetName() <<  " with " << outTrees[i]->GetEntries() << " of " << 
                              nEvents << " (" << outTrees[i]->GetEntries()/float(nEvents)*100 << " %) events selected" << std::endl;
            }
        }
};

#endif
