#include "yaml_parser.hh"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

SimpleYamlParser::SimpleYamlParser() {}

std::unordered_map<std::string, YamlValue> SimpleYamlParser::parseFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open YAML file: " + filename);
    }
    
    std::unordered_map<std::string, YamlValue> data;
    std::string line;
    int lineNumber = 0;
    
    while (std::getline(file, line)) {
        ++lineNumber;
        line = trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Find the colon separator
        size_t colonPos = line.find(':');
        if (colonPos == std::string::npos) {
            continue; // Skip invalid lines
        }
        
        std::string key = trim(line.substr(0, colonPos));
        std::string valuePart = line.substr(colonPos + 1);
        
        // Remove inline comments
        size_t commentPos = valuePart.find('#');
        if (commentPos != std::string::npos) {
            valuePart = valuePart.substr(0, commentPos);
        }
        
        std::string value = trim(valuePart);
        
        if (key.empty()) {
            throw std::runtime_error("Empty key at line " + std::to_string(lineNumber));
        }
        
        try {
            data[key] = parseValue(value);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing line " + std::to_string(lineNumber) + ": " + e.what());
        }
    }
    
    file.close();
    return data;
}

std::string SimpleYamlParser::getString(const std::unordered_map<std::string, YamlValue>& data, const std::string& key) {
    auto it = data.find(key);
    if (it == data.end()) {
        throw std::runtime_error("Key not found: " + key);
    }
    
    if (std::holds_alternative<std::string>(it->second)) {
        return std::get<std::string>(it->second);
    }
    throw std::runtime_error("Key '" + key + "' is not a string");
}

double SimpleYamlParser::getDouble(const std::unordered_map<std::string, YamlValue>& data, const std::string& key) {
    auto it = data.find(key);
    if (it == data.end()) {
        throw std::runtime_error("Key not found: " + key);
    }
    
    if (std::holds_alternative<double>(it->second)) {
        return std::get<double>(it->second);
    } else if (std::holds_alternative<int>(it->second)) {
        return static_cast<double>(std::get<int>(it->second));
    }
    throw std::runtime_error("Key '" + key + "' is not a number");
}

int SimpleYamlParser::getInt(const std::unordered_map<std::string, YamlValue>& data, const std::string& key) {
    auto it = data.find(key);
    if (it == data.end()) {
        throw std::runtime_error("Key not found: " + key);
    }
    
    if (std::holds_alternative<int>(it->second)) {
        return std::get<int>(it->second);
    }
    throw std::runtime_error("Key '" + key + "' is not an integer");
}

std::vector<double> SimpleYamlParser::getDoubleArray(const std::unordered_map<std::string, YamlValue>& data, const std::string& key) {
    auto it = data.find(key);
    if (it == data.end()) {
        throw std::runtime_error("Key not found: " + key);
    }
    
    if (std::holds_alternative<std::vector<double>>(it->second)) {
        return std::get<std::vector<double>>(it->second);
    } else if (std::holds_alternative<std::vector<int>>(it->second)) {
        auto intVec = std::get<std::vector<int>>(it->second);
        std::vector<double> doubleVec;
        for (int val : intVec) {
            doubleVec.push_back(static_cast<double>(val));
        }
        return doubleVec;
    }
    throw std::runtime_error("Key '" + key + "' is not an array");
}

std::vector<int> SimpleYamlParser::getIntArray(const std::unordered_map<std::string, YamlValue>& data, const std::string& key) {
    auto it = data.find(key);
    if (it == data.end()) {
        throw std::runtime_error("Key not found: " + key);
    }
    
    if (std::holds_alternative<std::vector<int>>(it->second)) {
        return std::get<std::vector<int>>(it->second);
    }
    throw std::runtime_error("Key '" + key + "' is not an integer array");
}

bool SimpleYamlParser::hasKey(const std::unordered_map<std::string, YamlValue>& data, const std::string& key) {
    return data.find(key) != data.end();
}

std::string SimpleYamlParser::trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    
    size_t end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

std::vector<std::string> SimpleYamlParser::split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(trim(token));
    }
    return tokens;
}

YamlValue SimpleYamlParser::parseValue(const std::string& value) {
    if (value.empty()) {
        return std::string("");
    }
    
    // Handle arrays [1, 2, 3] or [1.0, 2.0, 3.0]
    if (value.front() == '[' && value.back() == ']') {
        std::string arrayContent = value.substr(1, value.length() - 2);
        auto tokens = split(arrayContent, ',');
        
        if (tokens.empty()) {
            return std::vector<double>();
        }
        
        // Check if all elements are integers
        bool allInts = true;
        for (const auto& token : tokens) {
            if (!isInteger(token)) {
                allInts = false;
                break;
            }
        }
        
        if (allInts) {
            std::vector<int> intArray;
            for (const auto& token : tokens) {
                intArray.push_back(std::stoi(token));
            }
            return intArray;
        } else {
            std::vector<double> doubleArray;
            for (const auto& token : tokens) {
                if (isNumber(token)) {
                    doubleArray.push_back(std::stod(token));
                } else {
                    throw std::runtime_error("Invalid number in array: " + token);
                }
            }
            return doubleArray;
        }
    }
    
    // Handle strings with quotes
    if ((value.front() == '"' && value.back() == '"') || 
        (value.front() == '\'' && value.back() == '\'')) {
        return value.substr(1, value.length() - 2);
    }
    
    // Handle numbers
    if (isInteger(value)) {
        return std::stoi(value);
    } else if (isNumber(value)) {
        return std::stod(value);
    }
    
    // Default to string
    return value;
}

bool SimpleYamlParser::isNumber(const std::string& str) {
    if (str.empty()) return false;
    
    size_t start = 0;
    if (str[0] == '-' || str[0] == '+') start = 1;
    if (start >= str.length()) return false;
    
    bool hasDecimal = false;
    for (size_t i = start; i < str.length(); ++i) {
        if (str[i] == '.') {
            if (hasDecimal) return false;
            hasDecimal = true;
        } else if (!std::isdigit(str[i])) {
            return false;
        }
    }
    return true;
}

bool SimpleYamlParser::isInteger(const std::string& str) {
    if (str.empty()) return false;
    
    size_t start = 0;
    if (str[0] == '-' || str[0] == '+') start = 1;
    if (start >= str.length()) return false;
    
    for (size_t i = start; i < str.length(); ++i) {
        if (!std::isdigit(str[i])) {
            return false;
        }
    }
    return true;
}