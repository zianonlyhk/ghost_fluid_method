#include "config_manager.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <unordered_set>
#include <unordered_map>
#include <functional>

ConfigManager::ConfigManager() {}

SimulationConfig ConfigManager::loadConfig(const std::string& configPath) {
    SimulationConfig config;
    setDefaults(config);
    parseYamlConfig(configPath, config);
    validateConfig(config);
    return config;
}

void ConfigManager::setDefaults(SimulationConfig& config) {
    config.runName = "test";
    config.repoDir = std::filesystem::current_path().string();
    config.reinitFactor = 1;
    config.rigidBodyAdditionalFactor = 0.0;
    config.acceleratingRigidBody = 0;
}


void ConfigManager::validateConfig(const SimulationConfig& config) {
    // Validate domain boundaries
    if (config.x0 >= config.x1) {
        throw std::runtime_error("Invalid domain: x0 (" + std::to_string(config.x0) + 
                                ") must be less than x1 (" + std::to_string(config.x1) + ")");
    }
    if (config.y0 >= config.y1) {
        throw std::runtime_error("Invalid domain: y0 (" + std::to_string(config.y0) + 
                                ") must be less than y1 (" + std::to_string(config.y1) + ")");
    }
    
    // Validate grid resolution
    if (config.nCells_x < 10 || config.nCells_x > 2000) {
        throw std::runtime_error("Grid resolution nCells_x (" + std::to_string(config.nCells_x) + 
                                ") should be between 10 and 2000 for reasonable performance");
    }
    if (config.nCells_y < 10 || config.nCells_y > 2000) {
        throw std::runtime_error("Grid resolution nCells_y (" + std::to_string(config.nCells_y) + 
                                ") should be between 10 and 2000 for reasonable performance");
    }
    
    // Validate leftRightStateBoundary
    if (config.leftRightStateBoundary <= config.x0 || config.leftRightStateBoundary >= config.x1) {
        throw std::runtime_error("leftRightStateBoundary (" + std::to_string(config.leftRightStateBoundary) + 
                                ") must be between x0 (" + std::to_string(config.x0) + 
                                ") and x1 (" + std::to_string(config.x1) + ")");
    }
    
    // Validate CFL number
    if (config.c <= 0.0 || config.c > 1.0) {
        throw std::runtime_error("CFL number (" + std::to_string(config.c) + 
                                ") must be between 0 and 1 for numerical stability");
    }
    
    // Validate rigid body position
    const double margin = 0.01;  // Small margin from boundaries
    if (config.initRigidBodyLoc[0] < config.x0 + margin || 
        config.initRigidBodyLoc[0] > config.x1 - margin) {
        throw std::runtime_error("Rigid body x-position (" + std::to_string(config.initRigidBodyLoc[0]) + 
                                ") should be within domain [" + std::to_string(config.x0 + margin) + 
                                ", " + std::to_string(config.x1 - margin) + "]");
    }
    if (config.initRigidBodyLoc[1] < config.y0 + margin || 
        config.initRigidBodyLoc[1] > config.y1 - margin) {
        throw std::runtime_error("Rigid body y-position (" + std::to_string(config.initRigidBodyLoc[1]) + 
                                ") should be within domain [" + std::to_string(config.y0 + margin) + 
                                ", " + std::to_string(config.y1 - margin) + "]");
    }
    
    // Validate rigid body size
    const double domain_size_x = config.x1 - config.x0;
    const double domain_size_y = config.y1 - config.y0;
    const double max_size = std::min(domain_size_x, domain_size_y) * 0.4;  // Max 40% of smaller dimension
    
    if (config.rigidBodyRadiusOrLength <= 0.0) {
        throw std::runtime_error("Rigid body size (rigidBodyRadiusOrLength) must be positive");
    }
    if (config.rigidBodyRadiusOrLength > max_size) {
        throw std::runtime_error("Rigid body size (" + std::to_string(config.rigidBodyRadiusOrLength) + 
                                ") is too large for domain. Maximum recommended: " + std::to_string(max_size));
    }
    
    // Validate physical states
    auto validateState = [](const std::array<double, 4>& state, const std::string& name) {
        if (state[0] <= 0.0) {  // density
            throw std::runtime_error(name + " density (" + std::to_string(state[0]) + ") must be positive");
        }
        if (state[3] <= 0.0) {  // pressure
            throw std::runtime_error(name + " pressure (" + std::to_string(state[3]) + ") must be positive");
        }
        // Check for extreme values that might cause numerical issues
        if (state[0] > 100.0 || state[3] > 100.0) {
            std::cout << "Warning: " << name << " has large values (rho=" << state[0] 
                      << ", p=" << state[3] << ") that might cause numerical issues" << std::endl;
        }
    };
    
    validateState(config.leftState, "leftState");
    validateState(config.rightState, "rightState");
    
    // Validate timing parameters
    if (config.tStop <= 0.0) {
        throw std::runtime_error("Simulation time tStop (" + std::to_string(config.tStop) + ") must be positive");
    }
    if (config.tStop > 100.0) {
        std::cout << "Warning: Long simulation time (" << config.tStop 
                  << ") may take significant computational resources" << std::endl;
    }
    
    // Validate output parameters
    if (config.loggingFactor <= 0) {
        throw std::runtime_error("loggingFactor must be positive");
    }
    if (config.loggingFactor == 1) {
        std::cout << "Warning: loggingFactor=1 will generate many output files" << std::endl;
    }
}

void ConfigManager::parseYamlConfig(const std::string& configPath, SimulationConfig& config) {
    SimpleYamlParser parser;
    auto data = parser.parseFile(configPath);
    
    std::unordered_set<std::string> required_params = {
        "nCells_x", "nCells_y", "x0", "x1", "y0", "y1", 
        "c", "tStop", "loggingFactor",
        "rigidBodyType", "initRigidBodyVel", "initRigidBodyLoc",
        "rigidBodyRadiusOrLength", "leftState", "rightState",
        "leftRightStateBoundary"
    };
    
    // Parse required parameters
    if (parser.hasKey(data, "nCells_x")) {
        config.nCells_x = parser.getInt(data, "nCells_x");
        required_params.erase("nCells_x");
    }
    if (parser.hasKey(data, "nCells_y")) {
        config.nCells_y = parser.getInt(data, "nCells_y");
        required_params.erase("nCells_y");
    }
    if (parser.hasKey(data, "x0")) {
        config.x0 = parser.getDouble(data, "x0");
        required_params.erase("x0");
    }
    if (parser.hasKey(data, "x1")) {
        config.x1 = parser.getDouble(data, "x1");
        required_params.erase("x1");
    }
    if (parser.hasKey(data, "y0")) {
        config.y0 = parser.getDouble(data, "y0");
        required_params.erase("y0");
    }
    if (parser.hasKey(data, "y1")) {
        config.y1 = parser.getDouble(data, "y1");
        required_params.erase("y1");
    }
    if (parser.hasKey(data, "c")) {
        config.c = parser.getDouble(data, "c");
        required_params.erase("c");
    }
    if (parser.hasKey(data, "tStop")) {
        config.tStop = parser.getDouble(data, "tStop");
        required_params.erase("tStop");
    }
    if (parser.hasKey(data, "loggingFactor")) {
        config.loggingFactor = parser.getInt(data, "loggingFactor");
        required_params.erase("loggingFactor");
    }
    if (parser.hasKey(data, "rigidBodyType")) {
        config.rigidBodyType = parser.getInt(data, "rigidBodyType");
        required_params.erase("rigidBodyType");
    }
    if (parser.hasKey(data, "rigidBodyRadiusOrLength")) {
        config.rigidBodyRadiusOrLength = parser.getDouble(data, "rigidBodyRadiusOrLength");
        required_params.erase("rigidBodyRadiusOrLength");
    }
    if (parser.hasKey(data, "leftRightStateBoundary")) {
        config.leftRightStateBoundary = parser.getDouble(data, "leftRightStateBoundary");
        required_params.erase("leftRightStateBoundary");
    }
    
    // Parse arrays
    if (parser.hasKey(data, "initRigidBodyVel")) {
        auto vel = parser.getDoubleArray(data, "initRigidBodyVel");
        if (vel.size() != 2) {
            throw std::runtime_error("initRigidBodyVel must have exactly 2 elements");
        }
        config.initRigidBodyVel[0] = vel[0];
        config.initRigidBodyVel[1] = vel[1];
        required_params.erase("initRigidBodyVel");
    }
    if (parser.hasKey(data, "initRigidBodyLoc")) {
        auto loc = parser.getDoubleArray(data, "initRigidBodyLoc");
        if (loc.size() != 2) {
            throw std::runtime_error("initRigidBodyLoc must have exactly 2 elements");
        }
        config.initRigidBodyLoc[0] = loc[0];
        config.initRigidBodyLoc[1] = loc[1];
        required_params.erase("initRigidBodyLoc");
    }
    if (parser.hasKey(data, "leftState")) {
        auto state = parser.getDoubleArray(data, "leftState");
        if (state.size() != 4) {
            throw std::runtime_error("leftState must have exactly 4 elements [rho, u, v, p]");
        }
        config.leftState[0] = state[0];
        config.leftState[1] = state[1];
        config.leftState[2] = state[2];
        config.leftState[3] = state[3];
        required_params.erase("leftState");
    }
    if (parser.hasKey(data, "rightState")) {
        auto state = parser.getDoubleArray(data, "rightState");
        if (state.size() != 4) {
            throw std::runtime_error("rightState must have exactly 4 elements [rho, u, v, p]");
        }
        config.rightState[0] = state[0];
        config.rightState[1] = state[1];
        config.rightState[2] = state[2];
        config.rightState[3] = state[3];
        required_params.erase("rightState");
    }
    
    // Parse optional parameters
    if (parser.hasKey(data, "reinitFactor")) {
        config.reinitFactor = parser.getInt(data, "reinitFactor");
    }
    if (parser.hasKey(data, "acceleratingRigidBody")) {
        config.acceleratingRigidBody = parser.getInt(data, "acceleratingRigidBody");
    }
    if (parser.hasKey(data, "rigidBodyAdditionalFactor")) {
        config.rigidBodyAdditionalFactor = parser.getDouble(data, "rigidBodyAdditionalFactor");
    }
    if (parser.hasKey(data, "runName")) {
        config.runName = parser.getString(data, "runName");
    }
    
    // Check for missing required parameters
    if (!required_params.empty()) {
        std::string missing_params;
        for (const auto& param : required_params) {
            missing_params += param + ", ";
        }
        missing_params.erase(missing_params.length() - 2);
        throw std::runtime_error("Missing required parameters in YAML config: " + missing_params);
    }
}
