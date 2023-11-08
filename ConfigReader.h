#pragma once

#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <map>
#include <string>
#include <iostream>
#include <fstream>

class ConfigReader {
public:
    ConfigReader(const std::string& configFile);
    bool readConfigData(std::map<std::string, double>& configData);

private:
    std::string configFile;
};

#endif

