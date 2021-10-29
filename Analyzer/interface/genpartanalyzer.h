#ifndef GENPARTANALYZER_H
#define GENPARTANALYZER_H

#include <limits>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>

template <typename T>
class GenPartAnalyzer : public BaseAnalyzer<T> {
    private:
        bool isData;
        std::vector<int> alreadyMatchedIdx;

        //Function to get last copy of particle in gen collection
        int LastCopy(T& input, const int& idx){
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
        int Match(T& input, const float& pt, const float& phi, const float& eta, const int& PDG, const float& dRthr, const float& dPTthr){
            int genIdx = -1;
            float dR, dPT, 
            dRmin = std::numeric_limits<float>::max(), 
            dPTmin = std::numeric_limits<float>::max();

            for(int i = 0; i < input.genSize; ++i){
                input.GetGenPart(i);

                dR = Util2::DeltaR(eta, phi, input.genEta, input.genPhi);
                dPT = std::abs(pt - input.genPt)/pt;

                if(dR > dRthr or dPT > dPTthr) continue;
            
                if(dR < dRmin and dPT < dPTmin and PDG == std::abs(input.genPDG)){
                    int idx = LastCopy(input, i); 
                    if(std::find(alreadyMatchedIdx.begin(), alreadyMatchedIdx.end(), idx) != alreadyMatchedIdx.end()) continue;

                    genIdx = idx;
                    dRmin = dR;
                    dPTmin = dPT;
                }
            }

            return genIdx;
        }
    
    public:
        GenPartAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            isData = skim.get<std::string>("run") != "MC";
        }

        void Analyze(T& input, Output& out){
            if(isData) return;
    
            input.ReadGenEntry();
            int genIdx = -1;
            alreadyMatchedIdx.clear();

            for(int i = 0; i < out.nElectrons; ++i){
                genIdx = Match(input, out.elePt[i], out.elePhi[i], out.eleEta[i], 11, 0.4, 0.4);

                if(genIdx != -1){
                    alreadyMatchedIdx.push_back(genIdx);

                    input.GetGenPart(genIdx);
                    out.eleGenPt[i] = input.genPt;
                    out.eleGenPhi[i] = input.genPhi;
                    out.eleGenEta[i] = input.genEta;
                    out.eleGenID[i] = input.genPDG;

                    input.GetGenPart(LastCopy(input, input.genMotherIdx));
                    out.eleGenMotherID[i] = input.genPDG;

                    input.GetGenPart(input.genMotherIdx);
                    out.eleGenGrandMotherID[i] = input.genPDG;
                }

                else{
                    out.eleGenPt[i] = -999.;
                    out.eleGenPhi[i] = -999.;
                    out.eleGenEta[i] = -999.;
                    out.eleGenID[i] = -999;
                    out.eleGenMotherID[i] = -999;
                    out.eleGenGrandMotherID[i] = -999;
                }
            }

            for(int i = 0; i < out.nMuons; ++i){
                genIdx = Match(input, out.muPt[i], out.muPhi[i], out.muEta[i], 13, 0.4, 0.4);

                if(genIdx != -1){
                    alreadyMatchedIdx.push_back(genIdx);

                    input.GetGenPart(genIdx);
                    out.muGenPt[i] = input.genPt;
                    out.muGenPhi[i] = input.genPhi;
                    out.muGenEta[i] = input.genEta;
                    out.muGenID[i] = input.genPDG;

                    input.GetGenPart(input.genMotherIdx);
                    out.muGenMotherID[i] = input.genPDG;

                    input.GetGenPart(input.genMotherIdx);
                    out.muGenGrandMotherID[i] = input.genPDG;
                }

                else{
                    out.muGenPt[i] = -999.;
                    out.muGenPhi[i] = -999.;
                    out.muGenEta[i] = -999.;
                    out.muGenID[i] = -999;
                    out.muGenMotherID[i] = -999;
                    out.muGenGrandMotherID[i] = -999;
                }
            }

            for(int i = 0; i < out.nJets; ++i){
                genIdx = Match(input, out.jetPt[i], out.jetPhi[i], out.jetEta[i], 5, 0.4, 3.);

                if(genIdx != -1){
                    alreadyMatchedIdx.push_back(genIdx);

                    input.GetGenPart(genIdx);
                    out.jetGenPt[i] = input.genPt;
                    out.jetGenPhi[i] = input.genPhi;
                    out.jetGenEta[i] = input.genEta;
                    out.jetGenID[i] = input.genPDG;

                    input.GetGenPart(input.genMotherIdx);
                    out.jetGenMotherID[i] = input.genPDG;

                    input.GetGenPart(input.genMotherIdx);
                    out.jetGenGrandMotherID[i] = input.genPDG;
                }

                else{
                    out.jetGenPt[i] = -999.;
                    out.jetGenPhi[i] = -999.;
                    out.jetGenEta[i] = -999.;
                    out.jetGenID[i] = -999;
                    out.jetGenMotherID[i] = -999;
                    out.jetGenGrandMotherID[i] = -999;
                }
            }

            for(int i = 0; i < out.nSubJets; ++i){
                genIdx = Match(input, out.subJetPt[i], out.subJetPhi[i], out.subJetEta[i], 5, 0.4, 3.);

                if(genIdx != -1){
                    alreadyMatchedIdx.push_back(genIdx);

                    input.GetGenPart(genIdx);
                    out.subJetGenPt[i] = input.genPt;
                    out.subJetGenPhi[i] = input.genPhi;
                    out.subJetGenEta[i] = input.genEta;
                    out.subJetGenID[i] = input.genPDG;

                    input.GetGenPart(input.genMotherIdx);
                    out.subJetGenMotherID[i] = input.genPDG;

                    input.GetGenPart(input.genMotherIdx);
                    out.subJetGenGrandMotherID[i] = input.genPDG;
                }

                else{
                    out.subJetGenPt[i] = -999.;
                    out.subJetGenPhi[i] = -999.;
                    out.subJetGenEta[i] = -999.;
                    out.subJetGenID[i]= -999;
                    out.subJetGenMotherID[i] = -999;
                    out.subJetGenGrandMotherID[i] = -999;
                }
            }
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){
        }
};

#endif
