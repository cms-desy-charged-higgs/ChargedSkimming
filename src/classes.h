#include <ChargedAnalysis/Skimming/interface/nanoskimmer.h>
#include <ChargedAnalysis/Skimming/interface/jetanalyzer.h>
#include <ChargedAnalysis/Skimming/interface/electronanalyzer.h>
#include <ChargedAnalysis/Skimming/interface/muonanalyzer.h>
#include <ChargedAnalysis/Skimming/interface/genpartanalyzer.h>

namespace{
    namespace{
        Jet jet;
        FatJet fatjet;
        Electron electron;
        Muon muon;
        GenPart genPart;
        Particle particle;

        std::vector<Jet> jets;
        std::vector<Particle> Particles;
        std::vector<FatJet> fatjets;
        std::vector<Electron> electrons;
        std::vector<Muon> muons;

        NanoSkimmer skimmer;
    }
}
