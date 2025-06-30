#ifndef YAML_PARSER_HH
#define YAML_PARSER_HH

#include <string>
#include <unordered_map>
#include <vector>
#include <variant>

// Simple YAML value type
using YamlValue = std::variant<std::string, double, int, std::vector<double>, std::vector<int>>;

class SimpleYamlParser {
public:
    SimpleYamlParser();
    
    // Parse YAML file and return key-value map
    std::unordered_map<std::string, YamlValue> parseFile(const std::string& filename);
    
    // Helper functions to extract typed values
    std::string getString(const std::unordered_map<std::string, YamlValue>& data, const std::string& key);
    double getDouble(const std::unordered_map<std::string, YamlValue>& data, const std::string& key);
    int getInt(const std::unordered_map<std::string, YamlValue>& data, const std::string& key);
    std::vector<double> getDoubleArray(const std::unordered_map<std::string, YamlValue>& data, const std::string& key);
    std::vector<int> getIntArray(const std::unordered_map<std::string, YamlValue>& data, const std::string& key);
    
    // Check if key exists
    bool hasKey(const std::unordered_map<std::string, YamlValue>& data, const std::string& key);

private:
    std::string trim(const std::string& str);
    std::vector<std::string> split(const std::string& str, char delimiter);
    YamlValue parseValue(const std::string& value);
    bool isNumber(const std::string& str);
    bool isInteger(const std::string& str);
};

#endif