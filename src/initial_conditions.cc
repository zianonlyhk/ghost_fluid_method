#include "initial_conditions.hh"
#include "level_set_functions.hh"
#include <stdexcept>

void InitialConditions::setInitialConditions(
    std::vector<std::vector<std::array<double, 4>>>& domain,
    std::vector<std::vector<double>>& levelSetDomain,
    const SimulationConfig& config) {
    
    initializeLevelSet(levelSetDomain, config);
    initializeFluidState(domain, config);
}

void InitialConditions::transformToConservativeForm(
    std::vector<std::vector<std::array<double, 4>>>& domain) {
    
    const int nCell_y = domain.size() - Constants::TOTAL_GHOST_CELLS;
    const int nCell_x = domain[0].size() - Constants::TOTAL_GHOST_CELLS;

    for (int j = 0; j < nCell_y + Constants::TOTAL_GHOST_CELLS; ++j) {
        for (int i = 0; i < nCell_x + Constants::TOTAL_GHOST_CELLS; ++i) {
            const double rho = domain[j][i][0];
            const double u = domain[j][i][1];
            const double v = domain[j][i][2];
            const double p = domain[j][i][3];

            // Transform to conservative variables
            domain[j][i][1] = rho * u;  // momentum x
            domain[j][i][2] = rho * v;  // momentum y
            domain[j][i][3] = p / (Constants::GAMMA - 1.0) + 0.5 * rho * (u * u + v * v);  // total energy
        }
    }
}

void InitialConditions::initializeLevelSet(
    std::vector<std::vector<double>>& levelSetDomain,
    const SimulationConfig& config) {
    
    const int nCell_y = levelSetDomain.size() - Constants::TOTAL_GHOST_CELLS;
    const int nCell_x = levelSetDomain[0].size() - Constants::TOTAL_GHOST_CELLS;

    // Initialize all cells (including ghost cells)
    for (int j = 0; j < nCell_y + Constants::TOTAL_GHOST_CELLS; ++j) {
        for (int i = 0; i < nCell_x + Constants::TOTAL_GHOST_CELLS; ++i) {
            auto [x, y] = getGridCoordinates(i, j, config);
            
            switch (config.rigidBodyType) {
                case Constants::RigidBodyType::CIRCLE:
                    levelSetDomain[j][i] = LevelSetFunctions::singleCircle(
                        config.rigidBodyRadiusOrLength,
                        config.initRigidBodyLoc[0],
                        config.initRigidBodyLoc[1],
                        x, y
                    );
                    break;
                    
                case Constants::RigidBodyType::SQUARE:
                    levelSetDomain[j][i] = LevelSetFunctions::singleSquare(
                        config.rigidBodyRadiusOrLength,
                        config.initRigidBodyLoc[0],
                        config.initRigidBodyLoc[1],
                        x, y
                    );
                    break;
                    
                case Constants::RigidBodyType::DOUBLE_CIRCLE:
                    levelSetDomain[j][i] = LevelSetFunctions::doubleCircle(
                        config.rigidBodyRadiusOrLength,
                        config.rigidBodyRadiusOrLength,
                        config.initRigidBodyLoc[0],
                        config.initRigidBodyLoc[0],
                        config.initRigidBodyLoc[1] - config.rigidBodyAdditionalFactor,
                        config.initRigidBodyLoc[1] + config.rigidBodyAdditionalFactor,
                        x, y
                    );
                    break;
                    
                default:
                    throw std::invalid_argument("Invalid rigidBodyType: " + std::to_string(config.rigidBodyType));
            }
        }
    }
}

void InitialConditions::initializeFluidState(
    std::vector<std::vector<std::array<double, 4>>>& domain,
    const SimulationConfig& config) {
    
    const int nCell_y = domain.size() - Constants::TOTAL_GHOST_CELLS;
    const int nCell_x = domain[0].size() - Constants::TOTAL_GHOST_CELLS;

    // Initialize only interior cells (ghost cells will be handled by boundary conditions)
    for (int j = Constants::GHOST_CELLS; j < nCell_y + Constants::GHOST_CELLS; ++j) {
        for (int i = Constants::GHOST_CELLS; i < nCell_x + Constants::GHOST_CELLS; ++i) {
            auto [x, y] = getGridCoordinates(i, j, config);
            
            if (x <= config.leftRightStateBoundary) {
                domain[j][i] = config.leftState;
            } else {
                domain[j][i] = config.rightState;
            }
        }
    }
}

std::pair<double, double> InitialConditions::getGridCoordinates(
    int i, int j, const SimulationConfig& config) {
    
    const int nCell_x = config.nCells_x;
    const int nCell_y = config.nCells_y;
    
    const double dx = (config.x1 - config.x0) / nCell_x;
    const double dy = (config.y1 - config.y0) / nCell_y;
    
    const double x = config.x0 + (i - Constants::GHOST_CELLS) * dx;
    const double y = config.y0 + (j - Constants::GHOST_CELLS) * dy;
    
    return {x, y};
}