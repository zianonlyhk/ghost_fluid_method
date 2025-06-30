#ifndef CONFIG_MANAGER_HH
#define CONFIG_MANAGER_HH

#include <string>
#include <array>
#include "yaml_parser.hh"

struct SimulationConfig {
    // Grid parameters
    int nCells_x;
    int nCells_y;
    
    // Domain boundaries
    double x0, x1, y0, y1;
    
    // Time parameters
    double c;        // CFL number
    double tStop;    // Simulation end time
    
    // Output parameters
    int loggingFactor;
    int reinitFactor;
    
    // Rigid body parameters
    int rigidBodyType;
    int acceleratingRigidBody;
    std::array<double, 2> initRigidBodyVel;
    std::array<double, 2> initRigidBodyLoc;
    double rigidBodyRadiusOrLength;
    double rigidBodyAdditionalFactor;
    
    // Initial conditions
    std::array<double, 4> leftState;
    std::array<double, 4> rightState;
    double leftRightStateBoundary;
    
    // Run parameters
    std::string runName;
    std::string repoDir;
};

class ConfigManager {
public:
    ConfigManager();
    
    // Load configuration from YAML file
    SimulationConfig loadConfig(const std::string& configPath = "./config.yaml");
    
    // Validate configuration
    void validateConfig(const SimulationConfig& config);
    
private:
    void setDefaults(SimulationConfig& config);
    void parseYamlConfig(const std::string& configPath, SimulationConfig& config);
};

#endif