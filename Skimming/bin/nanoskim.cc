#include <ChargedSkimming/Core/interface/skimmer.h>
#include <ChargedSkimming/Core/interface/nanoinput.h>
#include <ChargedSkimming/Core/interface/output.h>

#include <vector>
#include <string>
#include <chrono>

std::string ParseLine(int argc, char* argv[], const std::string& name){
    std::string result;

    for(int i = 0; i < argc; ++i){
        std::string arg = std::string(argv[i]);

        if(arg.size() > 2 and name == arg.substr(2)){
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
    std::string outDir = ParseLine(argc, argv, "out-dir");
    std::string outFile = ParseLine(argc, argv, "out-file");
    std::string run = ParseLine(argc, argv, "run");
    std::string era = ParseLine(argc, argv, "era");
    std::string xSec = ParseLine(argc, argv, "xSec");
    std::string xSecUnc = ParseLine(argc, argv, "xSecUnc");
    std::vector<std::string> channels = SplitString(ParseLine(argc, argv, "channels"), " ");

    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

    NanoInput input(fileName, "Events");
    Output output;

    Skimmer<NanoInput> skimmer(channels, xSec, xSecUnc, era, run);
    skimmer.Configure(input, output, outDir, outFile);
  
    for(std::size_t entry = 0; entry < input.GetEntries(); ++entry){
        if(entry % 10000 == 0 and entry != 0){
            std::cout << "Events analyzed: " << entry  << " of " << input.GetEntries() << " (" << float(entry)/input.GetEntries()*100 << " %) [" << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count() << " seconds]" << std::endl;
        }

        input.SetEntry(entry);
        skimmer.Loop(input, output);
    }
   
    skimmer.WriteOutput();
}


