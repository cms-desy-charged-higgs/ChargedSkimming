#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <TH2F.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class ElectronAnalyzer : public BaseAnalyzer<T> {
    private:
        //Era information
        std::string era;

        //Kinematic cut criteria
        float ptCut, etaCut;

        //Name of the energy correction (dependent on systematic study)
        bool isScaleSyst, isSigmaSyst;
        std::string shift;

    public:
        ElectronAnalyzer(const std::string& systematic, const std::string& shift) : shift(shift){
            isScaleSyst = systematic.find("eleEnergyScale") != std::string::npos;
            isSigmaSyst = systematic.find("eleEnergySigma") != std::string::npos;
        }

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            //Read information needed
            era = skim.get<std::string>("era");

            ptCut = skim.get<float>("Analyzer.Electron.pt." + era);
            etaCut = skim.get<float>("Analyzer.Electron.eta." + era);
        };

        void Analyze(T& input, Output& out){
            out.nElectrons = 0;
            input.ReadEleEntry();
        
            //Loop over all electrons
            for(int i = 0; i < input.eleSize; ++i){
                if(out.nElectrons >= eleMax) break;
                input.GetElectron(i);

                float elePt;

                if(isScaleSyst){
                    elePt = shift == "Up" ? input.elePtScaleUp : input.elePtScaleDown;
                }

                else if(isSigmaSyst){
                    elePt = shift == "Up" ? input.elePtSigmaUp : input.elePtSigmaDown;
                }

                else{
                    elePt = input.elePt;
                }
              
                if(elePt > ptCut && std::abs(input.eleEta) < etaCut && input.eleConvVeto){
                    out.elePt[out.nElectrons] = elePt;
                    out.elePhi[out.nElectrons] =  input.elePhi;
                    out.eleEta[out.nElectrons] =  input.eleEta;
                    out.eleIso[out.nElectrons] =  input.eleIso;
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
