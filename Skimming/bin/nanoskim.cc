#include <ChargedSkimming/Core/interface/skimmer.h>
#include <ChargedSkimming/Core/interface/nanoinput.h>
#include <ChargedSkimming/Core/interface/output.h>

#include <vector>
#include <string>

std::string ParseLine(int argc, char* argv[], const std::string& name){
    std::string result;

    for(int i = 0; i < argc; ++i){
        std::string arg(argv[i]);

        if(arg.find(name) != std::string::npos){
            result = std::string(argv[i+1]);
        }
    }

    return result;
}

std::vector<std::string> SplitString(const std::string& splitString, const std::string& delimeter){
    //Function which handles splitting of string input
    std::vector<std::string> splittedString;
    std::string string;
    std::istringstream splittedStream(splitString);
    while (std::getline(splittedStream, string, delimeter.c_str()[0])){
        splittedString.push_back(string);
    }

    return splittedString;
}

int main(int argc, char* argv[]){
    //Extract informations of command line
    std::string fileName = ParseLine(argc, argv, "file-name");
    std::string outDirTmp = ParseLine(argc, argv, "out-dir");
    std::string outFile = ParseLine(argc, argv, "out-file");
    std::string run = ParseLine(argc, argv, "run");
    std::string era = ParseLine(argc, argv, "era");
    std::string xSec = ParseLine(argc, argv, "xSec");
    std::vector<std::string> channels = SplitString(ParseLine(argc, argv, "channels"), " ");
    std::vector<std::string> systematics = SplitString(ParseLine(argc, argv, "systematics"), " ");   

    NanoInput input(fileName, "Events");
    Output output;
    std::vector<Skimmer<NanoInput>> skimmers;

    for(const std::string& syst : systematics){
        for(const std::string shift : {"Up", "Down"}){
            if(syst == "Nominal" and shift == "Down") continue;
            const std::string systName = syst == "Nominal" ? syst : syst + shift;

            std::string outDir = outDirTmp;
            outDir.replace(outDir.find("[SYS]"), 5, systName);

            skimmers.push_back(Skimmer<NanoInput>(channels, xSec, era, run, syst, shift));
            skimmers.back().Configure(input, output, outDir, outFile);
        }
    }

    for(std::size_t entry = 0; entry < input.GetEntries(); ++entry){
        if(entry % 10000 == 0 and entry != 0){
            std::cout << "Events analyzed: " << entry  << " of " << input.GetEntries() << " (" << float(entry)/input.GetEntries()*100 << " %)"<< std::endl;
        }

        input.SetEntry(entry);

        for(Skimmer<NanoInput>& skimmer : skimmers) skimmer.Loop(input, output);
    }
   
    for(Skimmer<NanoInput>& skimmer : skimmers) skimmer.WriteOutput();
}


