#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>
#include <RoccoR/RoccoR.cc>

template <typename T>
class MuonAnalyzer : public BaseAnalyzer<T> {
    private:
        //Era information
        bool isData;
        std::string era;

        //Kinematic cut criteria
        float ptCut, etaCut;

        //Muon scale corrector
        RoccoR rc; 


    public:
        MuonAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            isData = skim.get<std::string>("run") != "MC";

            //Read information needed
            era = skim.get<std::string>("era");

            ptCut = skim.get<float>("Analyzer.Muon.pt." + era);
            etaCut = skim.get<float>("Analyzer.Muon.eta." + era);

            rc.init(std::string(std::getenv("CMSSW_BASE")) + "/src/" + sf.get<std::string>("Muon.Scale." + era));
        };

        void Analyze(T& input, Output& out){
            out.nMuons = 0;
            input.ReadMuEntry();
            if(!isData) input.ReadGenEntry();
        
            //Loop over all electrons
            for(int i = 0; i < input.muSize; ++i){
                if(out.nMuons >= muMax) break;
                input.GetMuon(i);

                if(input.muPt > ptCut && std::abs(input.muEta) < etaCut){
                    //Rochester pt correction
                    int genIdx = -1;
                    double dtSF = 1., mcSF = 1., unc = 0.;

                    if(!isData){
                        genIdx = input.GenMatch(input, input.muPt, input.muPhi, input.muEta, 13, 0.4, 0.4);

                        if(genIdx != -1){
                            input.alreadyMatchedIdx.push_back(genIdx);
                            input.GetGenPart(genIdx);

                            mcSF = rc.kSpreadMC(input.muCharge, input.muPt, input.muEta, input.muPhi, input.genPt, 0, 0);
                            unc = rc.kSpreadMCerror(input.muCharge, input.muPt, input.muEta, input.muPhi, input.genPt); 
                        }

                        else{
                            mcSF = rc.kSmearMC(input.muCharge, input.muPt, input.muEta, input.muPhi, input.muNTrackerLayers, input.muRandomNumber, 0, 0);
                            unc = rc.kSmearMCerror(input.muCharge, input.muPt, input.muEta, input.muPhi, input.muNTrackerLayers, input.muRandomNumber); 
                        }                       
                    }

                    else dtSF = rc.kScaleDT(input.muCharge, input.muPt, input.muEta, input.muPhi, 0, 0);

                    if(!(mcSF < 2) or !(mcSF > 0)) mcSF = 1;
                    if(!(dtSF < 2) or !(dtSF > 0)) dtSF = 1;
                    if(!(unc < 2) or !(unc > 0)) unc = 0;

                    float muPt = isData ? input.muPt*dtSF : input.muPt*mcSF;
                    float muPtUp = isData ? input.muPt*dtSF : input.muPt*(mcSF + unc);
                    float muPtDown = isData ? input.muPt*dtSF : input.muPt*(mcSF - unc);

                    //Check if one pt value passed pt criteria
                    bool ptCriteria = false;
            
                    for(const float pt : {muPt, muPtUp, muPtDown}){
                        if(pt > ptCut){
                            ptCriteria = true;
                            break;
                        }
                    }
        
                    if(ptCriteria){
                        //Electron four momentum components
                        out.muPt[out.nMuons] = muPt;
                        out.muPtUp[out.nMuons] = muPtUp;
                        out.muPtDown[out.nMuons] = muPtDown;
                        out.muPhi[out.nMuons] =  input.muPhi;
                        out.muEta[out.nMuons] =  input.muEta;
                        out.muIso03[out.nMuons] =  input.muIso03;
                        out.muIso04[out.nMuons] =  input.muIso04;
                        out.muMiniIso[out.nMuons] =  input.muMiniIso;
                        out.muDxy[out.nMuons] =  input.muDxy;
                        out.muDz[out.nMuons] =  input.muDz;
                        out.muRelJetIso[out.nMuons] =  input.muRelJetIso;

                        out.muCharge[out.nMuons] = input.muCharge;
                        out.muCutID[out.nMuons] = input.muCutID;
                        out.muMVAID[out.nMuons] =  input.muMVAID;

                        if(genIdx != -1){
                            input.GetGenPart(genIdx);
                            out.muGenPt[out.nMuons] = input.genPt;
                            out.muGenPhi[out.nMuons] = input.genPhi;
                            out.muGenEta[out.nMuons] = input.genEta;
                            out.muGenID[out.nMuons] = input.genPDG;

                            input.GetGenPart(input.genMotherIdx);
                            out.muGenMotherID[out.nMuons] = input.genPDG;

                            input.GetGenPart(input.genMotherIdx);
                            out.muGenGrandMotherID[out.nMuons] = input.genPDG;
                        }

                        else{
                            out.muGenPt[out.nMuons] = -999.;
                            out.muGenPhi[out.nMuons] = -999.;
                            out.muGenEta[out.nMuons] = -999.;
                            out.muGenID[out.nMuons] = -999;
                            out.muGenMotherID[out.nMuons] = -999;
                            out.muGenGrandMotherID[out.nMuons] = -999;
                        }

                        ++out.nMuons;
                    }
                }
            }
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
