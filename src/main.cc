#include "inline/debug_tools.hh"
#include "gfm_2d_euler_solver.hh"
#include "config_manager.hh"
#include "initial_conditions.hh"
#include "constants.hh"
#include <iostream>

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [config_file.yaml]" << std::endl;
    std::cout << std::endl;
    std::cout << "Ghost Fluid Method 2D Euler Solver" << std::endl;
    std::cout << "Solves 2D Euler equations with immersed solid boundaries" << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  config_file.yaml    YAML configuration file (default: config.yaml)" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << programName << "                                    # Use default config.yaml" << std::endl;
    std::cout << "  " << programName << " examples/shock_circle_interaction.yaml  # Use specific example" << std::endl;
    std::cout << "  " << programName << " my_custom_config.yaml              # Use custom configuration" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -h, --help          Show this help message" << std::endl;
}

int main(int argc, char* argv[]) {
    try {
        // Parse command line arguments
        std::string configPath = "./config.yaml";  // Default configuration file
        
        if (argc > 1) {
            std::string arg1 = argv[1];
            if (arg1 == "-h" || arg1 == "--help") {
                printUsage(argv[0]);
                return 0;
            } else if (arg1.substr(0, 1) == "-") {
                std::cerr << "Error: Unknown option '" << arg1 << "'" << std::endl;
                std::cerr << "Use -h or --help for usage information." << std::endl;
                return 1;
            } else {
                configPath = arg1;
            }
        }
        
        if (argc > 2) {
            std::cerr << "Error: Too many arguments provided." << std::endl;
            std::cerr << "Use -h or --help for usage information." << std::endl;
            return 1;
        }
        
        // Load configuration
        ConfigManager configManager;
        SimulationConfig config = configManager.loadConfig(configPath);

        // Initialize computational domain
        std::vector<std::vector<std::array<double, 4>>> compDomain;
        std::vector<std::vector<double>> levelSetCompDomain;
        
        const int totalCells_x = config.nCells_x + Constants::TOTAL_GHOST_CELLS;
        const int totalCells_y = config.nCells_y + Constants::TOTAL_GHOST_CELLS;
        
        compDomain.resize(totalCells_y);
        levelSetCompDomain.resize(totalCells_y);
        for (int j = 0; j < totalCells_y; ++j) {
            compDomain[j].resize(totalCells_x);
            levelSetCompDomain[j].resize(totalCells_x);
        }

        // Set initial conditions
        InitialConditions::setInitialConditions(compDomain, levelSetCompDomain, config);
        InitialConditions::transformToConservativeForm(compDomain);

        // Initialize solver
        GFM_2D_EulerSolver gfmSolver(compDomain, config.nCells_x, config.nCells_y);
        gfmSolver.setBound(config.x0, config.x1, config.y0, config.y1, config.tStop);
        gfmSolver.setCFL(config.c);
        gfmSolver.setName(config.runName);
        gfmSolver.setRepoDir(config.repoDir);
        gfmSolver.setLevelSet(levelSetCompDomain);
        gfmSolver.setRigidBodyVel(config.initRigidBodyVel);
        gfmSolver.setRigidBodyCentreCoor(config.initRigidBodyLoc);

        // Initial setup
        gfmSolver.updateBoundaryTrans();
        gfmSolver.calculateMockSchliren();
        gfmSolver.initiateOutputLogging();
        
        std::cout << "1: 0 / " << gfmSolver.tStop() << std::endl;

        // Main simulation loop
        double t = 0.0;
        int numIter = 0;
        
        while (t < gfmSolver.tStop()) {
            ++numIter;

            // Update time step
            gfmSolver.updateMaxA(numIter);
            gfmSolver.updateDt();
            t += gfmSolver.dt();

            // Ghost cell operations
            gfmSolver.updateGhostCellBoundary();
            gfmSolver.propagateGhostCell();

            // Fluid dynamics sweeps
            gfmSolver.mhHllcSweepX();
            gfmSolver.updateBoundaryTrans();
            gfmSolver.mhHllcSweepY();
            gfmSolver.updateBoundaryTrans();

            gfmSolver.calculateMockSchliren();

            // Rigid body operations
            if (config.acceleratingRigidBody == Constants::AccelerationFlag::ENABLED) {
                gfmSolver.accelerateRigidBody_circ();
            }
            
            gfmSolver.advectLevelSet();
            gfmSolver.updateLevelSetBoundaryTrans();

            // Level set reinitialization
            if (numIter % config.reinitFactor == 0) {
                gfmSolver.reinitLevelSet();
            }

            // Output logging
            if (numIter % config.loggingFactor == 0) {
                gfmSolver.writeToFiles(t);
                std::cout << (numIter / config.loggingFactor + 1) << ": " << t << " / " << gfmSolver.tStop() << std::endl;
            }
        }

        gfmSolver.cleanUp();
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
