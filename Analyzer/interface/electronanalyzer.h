#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <TH2F.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class ElectronAnalyzer : public BaseAnalyzer<T> {
    private:
        //Era information
        bool isData;
        std::string era;

        //Kinematic cut criteria
        float ptCut, etaCut;

        //Name of the energy correction (dependent on systematic study)
        bool isScaleSyst, isSigmaSyst;
        std::string shift;

    public:
        ElectronAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            //Read information needed
            isData = skim.get<std::string>("run") != "MC";
            era = skim.get<std::string>("era");

            ptCut = skim.get<float>("Analyzer.Electron.pt." + era);
            etaCut = skim.get<float>("Analyzer.Electron.eta." + era);
        };

        void Analyze(T& input, Output& out){
            out.nElectrons = 0;
            input.ReadEleEntry();
            if(!isData) input.ReadGenEntry();
        
            //Loop over all electrons
            for(int i = 0; i < input.eleSize; ++i){
                if(out.nElectrons >= eleMax) break;
                input.GetElectron(i);

                //Check if one pt value passed pt criteria
                bool ptCriteria = false;
        
                for(const float pt : {input.elePt, input.elePtScaleDown, input.elePtSigmaDown, input.elePtScaleUp, input.elePtSigmaUp}){
                    if(pt > ptCut){
                        ptCriteria = true;
                        break;
                    }
                }
              
                if(ptCriteria && std::abs(input.eleEta) < etaCut && input.eleConvVeto){
                    if(!isData){
                        int genIdx = input.GenMatch(input, input.elePt, input.elePhi, input.eleEta, 11, 0.4, 0.4);

                        if(genIdx != -1){
                            input.alreadyMatchedIdx.push_back(genIdx);

                            input.GetGenPart(genIdx);
                            out.eleGenPt[out.nElectrons] = input.genPt;
                            out.eleGenPhi[out.nElectrons] = input.genPhi;
                            out.eleGenEta[out.nElectrons] = input.genEta;
                            out.eleGenID[out.nElectrons] = input.genPDG;

                            input.GetGenPart(input.genMotherIdx);
                            out.eleGenMotherID[out.nElectrons] = input.genPDG;

                            input.GetGenPart(input.genMotherIdx);
                            out.eleGenGrandMotherID[out.nElectrons] = input.genPDG;
                        }

                        else{
                            out.eleGenPt[out.nElectrons] = -999.;
                            out.eleGenPhi[out.nElectrons] = -999.;
                            out.eleGenEta[out.nElectrons] = -999.;
                            out.eleGenID[out.nElectrons] = -999;
                            out.eleGenMotherID[out.nElectrons] = -999;
                            out.eleGenGrandMotherID[out.nElectrons] = -999;
                        }
                    }

                    out.elePt[out.nElectrons] = input.elePt;
                    out.elePtEnergyScaleUp[out.nElectrons] = input.elePtScaleUp;
                    out.elePtEnergyScaleDown[out.nElectrons] = input.elePtScaleDown;
                    out.elePtEnergySigmaUp[out.nElectrons] = input.elePtSigmaUp;
                    out.elePtEnergySigmaDown[out.nElectrons] = input.elePtSigmaDown;
                    out.elePt[out.nElectrons] = input.elePt;
                    out.elePt[out.nElectrons] = input.elePt;
                    out.elePhi[out.nElectrons] =  input.elePhi;
                    out.eleEta[out.nElectrons] =  input.eleEta;
                    out.eleIso03[out.nElectrons] =  input.eleIso03;
                    out.eleMiniIso[out.nElectrons] =  input.eleMiniIso;
                    out.eleDxy[out.nElectrons] =  input.eleDxy;
                    out.eleDz[out.nElectrons] =  input.eleDz;
                    out.eleRelJetIso[out.nElectrons] =  input.eleRelJetIso;

                    out.eleCharge[out.nElectrons] = input.eleCharge;
                    out.eleCutID[out.nElectrons] = input.eleCutID;
                    out.eleMVAID[out.nElectrons] =  input.eleMVAID;

                    ++out.nElectrons;
                }
            }
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
