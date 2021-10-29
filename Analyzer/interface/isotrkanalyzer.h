#ifndef ISOTRKANALYZER_H
#define ISOTRKANALYZER_H

#include <TParameter.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class IsotrkAnalyzer : public BaseAnalyzer<T> {
    private:
        std::string era;
        float etaCut;

    public:
        IsotrkAnalyzer() {}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            //Read information needed
            era = skim.get<std::string>("era");

            etaCut = std::min(skim.get<float>("Analyzer.Electron.eta." + era), skim.get<float>("Analyzer.Electron.eta." + era));
        }

        void Analyze(T& input, Output& out){
            out.isotrkSize = 0;
            input.ReadIsotrkEntry();
        
            //Loop over all electrons
            for(int i = 0; i < input.isotrkSize; ++i){
                if(out.isotrkSize >= isotrkMax) break;

                input.GetIsotrk(i);

                if(std::abs(input.isotrkEta) < etaCut and input.isotrkPt > 20 and std::abs(input.isotrkDz) < 0.1){
                    out.isotrkPt[out.isotrkSize] = input.isotrkPt;
                    out.isotrkEta[out.isotrkSize] = input.isotrkEta;
                    out.isotrkPhi[out.isotrkSize] = input.isotrkPhi;

                    out.isotrkIso[out.isotrkSize] = input.isotrkIso;
                    out.isotrkMiniIso[out.isotrkSize] = input.isotrkMiniIso;
                    out.isotrkDxy[out.isotrkSize] = input.isotrkDxy;
                    out.isotrkDz[out.isotrkSize] = input.isotrkDz;

                    //Charge
                    out.isotrkCharge[out.isotrkSize] = input.isotrkPDG > 0 ? 1 : -1;
                    out.isotrkPDG[out.isotrkSize] = input.isotrkPDG;

                    ++out.isotrkSize;
                }
            }
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
