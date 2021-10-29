#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class MuonAnalyzer : public BaseAnalyzer<T> {
    private:
        //Era information
        std::string era;

        //Kinematic cut criteria
        float ptCut, etaCut;

    public:
        MuonAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            //Read information needed
            era = skim.get<std::string>("era");

            ptCut = skim.get<float>("Analyzer.Muon.pt." + era);
            etaCut = skim.get<float>("Analyzer.Muon.eta." + era);
        };

        void Analyze(T& input, Output& out){
            out.nMuons = 0;
            input.ReadMuEntry();
        
            //Loop over all electrons
            for(int i = 0; i < input.muSize; ++i){
                if(out.nMuons >= muMax) break;
                input.GetMuon(i);

                if(input.muPt > ptCut && std::abs(input.muEta) < etaCut){
                    //Electron four momentum components
                    out.muPt[out.nMuons] = input.muPt;
                    out.muPhi[out.nMuons] =  input.muPhi;
                    out.muEta[out.nMuons] =  input.muEta;
                    out.muIso[out.nMuons] =  input.muIso;
                    out.muMiniIso[out.nMuons] =  input.muMiniIso;
                    out.muDxy[out.nMuons] =  input.muDxy;
                    out.muDz[out.nMuons] =  input.muDz;
                    out.muRelJetIso[out.nMuons] =  input.muRelJetIso;

                    out.muCharge[out.nMuons] = input.muCharge;
                    out.muCutID[out.nMuons] = input.muCutID;
                    out.muMVAID[out.nMuons] =  input.muMVAID;

                    ++out.nMuons;
                };
            }
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
