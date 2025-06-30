#ifndef INITIAL_CONDITIONS_HH
#define INITIAL_CONDITIONS_HH

#include "config_manager.hh"
#include "constants.hh"
#include <vector>
#include <array>

class InitialConditions {
public:
    // Set initial conditions based on configuration
    static void setInitialConditions(
        std::vector<std::vector<std::array<double, 4>>>& domain,
        std::vector<std::vector<double>>& levelSetDomain,
        const SimulationConfig& config
    );
    
    // Transform primitive variables to conservative form
    static void transformToConservativeForm(
        std::vector<std::vector<std::array<double, 4>>>& domain
    );

private:
    // Initialize level set function based on rigid body type
    static void initializeLevelSet(
        std::vector<std::vector<double>>& levelSetDomain,
        const SimulationConfig& config
    );
    
    // Initialize fluid state variables
    static void initializeFluidState(
        std::vector<std::vector<std::array<double, 4>>>& domain,
        const SimulationConfig& config
    );
    
    // Calculate grid coordinates
    static std::pair<double, double> getGridCoordinates(
        int i, int j, const SimulationConfig& config
    );
};

#endif