#include <ChargedHiggs/Skimming/interface/nanoskimmer.h>
#include <ChargedHiggs/Skimming/interface/jetanalyzer.h>
#include <ChargedHiggs/Skimming/interface/electronanalyzer.h>
#include <ChargedHiggs/Skimming/interface/muonanalyzer.h>
#include <ChargedHiggs/Skimming/interface/genpartanalyzer.h>

namespace{
    namespace{
        Jet jet;
        FatJet fatjet;
        Electron electron;
        Muon muon;
        GenPart genPart;

        std::vector<Jet> jets;
        std::vector<FatJet> fatjets;
        std::vector<Electron> electrons;
        std::vector<Muon> muons;

        NanoSkimmer skimmer;
    }
}
