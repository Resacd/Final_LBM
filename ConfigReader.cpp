#include "ConfigReader.h"

ConfigReader::ConfigReader(const std::string& configFile) : configFile(configFile) {}

bool ConfigReader::readConfigData(std::map<std::string, double>& configData) {
    std::ifstream file(configFile);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open " << configFile << " for reading." << std::endl;
        return false;
    }

    configData.clear();

    std::string line;
    while (std::getline(file, line)) {
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            double value = std::stod(line.substr(pos + 1));

            configData[key] = value;
        }
    }

    file.close();
    return true;
}
