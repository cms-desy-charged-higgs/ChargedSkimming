#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <vector>

#include <ChargedSkimming/Skimming/interface/util.h>

#include <TRandom.h>

struct Input{
    public:
        //Weight related
        float pdfWeight[102], scaleWeight[8];
        short nTrueInt;

        //Trigger related
        std::vector<short> triggers;
        std::vector<short> METFilter;

        //Electron related
        short eleSize, eleCharge, eleCutID, eleMVAID, eleConvVeto;
        float elePt, elePtSigmaUp, elePtSigmaDown, elePtScaleUp, elePtScaleDown, eleEta, elePhi, eleDxy, eleDz, eleRelJetIso, eleIso03, eleIso04, eleMiniIso;

        //Muon related
        short muSize, muCharge, muCutID, muMVAID, muNTrackerLayers;
        float muPt, muEta, muPhi, muDxy, muDz, muRelJetIso, muIso03, muIso04, muMiniIso, muRandomNumber;

        //Jet related
        short jetSize, fatJetSize, genJetSize, genFatJetSize,
              jetGenIdx, fatJetGenIdx, fatJetDAK8ID,
              jetID, jetPUID;
        float metPt, metPhi, metDeltaUnClustX, metDeltaUnClustY,
              jetPt, jetMass, jetPtRaw, jetMassRaw,
              jetEta, jetPhi, jetArea,
              jetDeepJet, jetDeepCSV, jetPartFlav,
              genJetEta, genJetPhi, genJetPt,
              fatJetPt, fatJetMass, fatJetPtRaw, fatJetMassRaw,
              fatJetPhi, fatJetEta, fatJetArea,
              fatJetTau1, fatJetTau2, fatJetTau3,
              genFatJetEta, genFatJetPhi, genFatJetPt,
              rho;

        //Iso. track related
        short isotrkSize, isotrkPDG;
        float isotrkPt, isotrkEta, isotrkPhi, isotrkDxy, isotrkDz, isotrkIso03, isotrkIso04, isotrkMiniIso;

        //Misc related
        short nParton;
        long evNr;
        float preFire, preFireUp, preFireDown;

        //Gen part related
        short genSize, genPDG, genMotherIdx;
        float genPt, genPhi, genEta, genMass;

        std::vector<int> alreadyMatchedIdx;

        //Function to get last copy of particle in gen collection
        template <typename T>
        int LastGenCopy(T& input, const int& idx){
            input.GetGenPart(idx);
        
            int partIdx = idx, motherIdx = input.genMotherIdx;
            int partPDG = input.genPDG;

            while(true){
                input.GetGenPart(motherIdx);

                if(partPDG == input.genPDG){
                    partIdx = motherIdx;
                    motherIdx = input.genMotherIdx;
                }

                else break;
            }

            return partIdx;
        }

        //Matching function
        template <typename T>
        int GenMatch(T& input, const float& pt, const float& phi, const float& eta, const int& PDG, const float& dRthr, const float& dPTthr){
            int genIdx = -1;
            float dR, dPT, 
            dRmin = std::numeric_limits<float>::max(), 
            dPTmin = std::numeric_limits<float>::max();

            for(int i = 0; i < input.genSize; ++i){
                input.GetGenPart(i);

                dR = Util::DeltaR(eta, phi, input.genEta, input.genPhi);
                dPT = std::abs(pt - input.genPt)/pt;

                if(dR > dRthr or dPT > dPTthr) continue;
            
                if(dR < dRmin and dPT < dPTmin and PDG == std::abs(input.genPDG)){
                    int idx = LastGenCopy(input, i); 
                    if(std::find(alreadyMatchedIdx.begin(), alreadyMatchedIdx.end(), idx) != alreadyMatchedIdx.end()) continue;

                    genIdx = idx;
                    dRmin = dR;
                    dPTmin = dPT;
                }
            }

            return genIdx;
        }
};

#endif
