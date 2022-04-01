#ifndef WEIGHTANALYZER_H
#define WEIGHTANALYZER_H

#include <TParameter.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class WeightAnalyzer : public BaseAnalyzer<T> {
    private:
        std::string run, era, pileUpFile;
        TParameter<float> lumi, lumiUp, lumiDown, xSec, xSecUp, xSecDown;
        float nGen;

        std::shared_ptr<TH1F> puMC;

    public:
        WeightAnalyzer() {}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            //Read information needed
            run = skim.get<std::string>("run");
            era = skim.get<std::string>("era");

            pileUpFile = this->filePath + skim.get<std::string>("Analyzer.Weight.PileUp." + era);

            //Define TParameter
            lumi = TParameter<float>("Lumi", skim.get<float>("Analyzer.Weight.Lumi.nominal." + era)*1e3);
            lumiUp = TParameter<float>("LumiUp", skim.get<float>("Analyzer.Weight.Lumi.Up." + era)*1e3);
            lumiDown = TParameter<float>("LumiDown", skim.get<float>("Analyzer.Weight.Lumi.Down." + era)*1e3);
            xSec  = TParameter<float>("xSec", skim.get<float>("xSec"));
            xSecUp  = TParameter<float>("xSecUp", skim.get<float>("xSec") + skim.get<float>("xSecUnc"));
            xSecDown  = TParameter<float>("xSecDown", skim.get<float>("xSec") - skim.get<float>("xSecUnc"));

            nGen = 0;

            //Histogram for MC PileUp distribution
            puMC = std::make_shared<TH1F>("puMC", "puMC", 100, 0, 100);
        }

        void Analyze(T& input, Output& out){
            if(run == "MC") {
                input.SetWeight();
                input.GetWeightEntry();

                std::copy(std::begin(input.pdfWeight), std::end(input.pdfWeight), std::begin(out.pdfWeight));
                std::copy(std::begin(input.scaleWeight), std::end(input.scaleWeight), std::begin(out.scaleWeight));

                out.nTrueInt = input.nTrueInt;
                nGen += 1;

                puMC->Fill(out.nTrueInt);
            }
            
            out.preFire = input.preFire;
            out.preFireUp = input.preFireUp;
            out.preFireDown = input.preFireDown;
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){
            outFile->cd();

            //Measured PU distributions https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
            if(run == "MC") {
                puMC->Write();

                for(const std::string& syst: {"", "Up", "Down"}){
                    std::string name = pileUpFile;
                    name.replace(name.find("@"), 1, syst);

                    std::shared_ptr<TFile> pileFile = std::make_shared<TFile>(name.c_str(), "READ");
                    std::shared_ptr<TH1F> realPile(static_cast<TH1F*>(pileFile->Get("pileup")));
                    realPile->SetName(("pileUp" + syst).c_str());
                    outFile->cd();
                    realPile->Write();
                }

                TParameter<float>("nGen", nGen).Write();
                xSec.Write();
                xSecUp.Write();
                xSecDown.Write();
                lumi.Write();
                lumiUp.Write();
                lumiDown.Write();
            }
        };
};

#endif
