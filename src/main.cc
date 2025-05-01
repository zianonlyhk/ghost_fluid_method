#include "inline/debug_tools.hh"
#include "gfm_2d_euler_solver.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <math.h>
#include <filesystem>
#include <unordered_set>
#include <functional>

#define PI 3.14159265

void loadConfig(int &i_nCells_x, int &i_nCells_y, double &i_x0, double &i_x1, double &i_y0, double &i_y1, double &i_c, double &i_tStop, int &i_loggingFactor, int &i_reinitFactor,
                int &i_rigidBodyType, int &i_acceleratingRigidBody, std::array<double, 2> &i_initRigidBodyVel, std::array<double, 2> &i_initRigidBodyLoc,
                double &i_rigidBodyRadiusOrLength, double &i_rigidBodyAdditionalFactor,
                std::array<double, 4> &i_leftState, std::array<double, 4> &i_rightState, double &i_leftRightStateBoundary,
                std::string &i_runName, std::string &i_repoDir)
{
    // Set default values for optional parameters
    i_runName = "test";
    i_repoDir = std::filesystem::current_path().string();
    i_reinitFactor = 1;

    // Open config file
    std::ifstream in("./config.txt");
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open config.txt");
    }

    std::unordered_set<std::string> required_params = {
        "nCells_x", "nCells_y", "x0", "x1", "y0", "y1", 
        "c", "tStop", "loggingFactor",
        "rigidBodyType", "initRigidBodyVel", "initRigidBodyLoc",
        "rigidBodyRadiusOrLength", "leftState", "rightState",
        "leftRightStateBoundary"
    };

    auto readInt = [&](int &val, bool positive = false, int min = 0, int max = 0) {
        if (!(in >> val)) throw std::runtime_error("Invalid integer value");
        if (positive && val <= 0) throw std::runtime_error("Value must be positive");
        if (min != max && (val < min || val > max)) throw std::runtime_error("Value must be between " + std::to_string(min) + " and " + std::to_string(max));
    };

    auto readDouble = [&](double &val, bool positive = false, double min = 0.0, double max = 0.0) {
        if (!(in >> val)) throw std::runtime_error("Invalid double value");
        if (positive && val <= 0) throw std::runtime_error("Value must be positive");
        if (min != max && (val < min || val > max)) throw std::runtime_error("Value must be between " + std::to_string(min) + " and " + std::to_string(max));
    };

    auto readString = [&](std::string &val) {
        if (!(in >> val)) throw std::runtime_error("Invalid string value");
    };

    auto readArray = [&](auto &arr, size_t size) {
        for (size_t i = 0; i < size; ++i) {
            if (!(in >> arr[i])) throw std::runtime_error("Invalid array element at index " + std::to_string(i));
        }
    };

    // Define parameter handlers
    std::unordered_map<std::string, std::function<void()>> param_handlers = {
        {"nCells_x", [&](){ readInt(i_nCells_x, true); required_params.erase("nCells_x"); }},
        {"nCells_y", [&](){ readInt(i_nCells_y, true); required_params.erase("nCells_y"); }},
        {"x0", [&](){ readDouble(i_x0); required_params.erase("x0"); }},
        {"x1", [&](){ readDouble(i_x1); required_params.erase("x1"); }},
        {"y0", [&](){ readDouble(i_y0); required_params.erase("y0"); }},
        {"y1", [&](){ readDouble(i_y1); required_params.erase("y1"); }},
        {"c", [&](){ readDouble(i_c); required_params.erase("c"); }},
        {"tStop", [&](){ readDouble(i_tStop, true); required_params.erase("tStop"); }},
        {"loggingFactor", [&](){ readInt(i_loggingFactor, true); required_params.erase("loggingFactor"); }},
        {"reinitFactor", [&](){ readInt(i_reinitFactor); }},
        {"rigidBodyType", [&](){ readInt(i_rigidBodyType, false, 1, 3); required_params.erase("rigidBodyType"); }},
        {"acceleratingRigidBody", [&](){ readInt(i_acceleratingRigidBody); }},
        {"initRigidBodyVel", [&](){ readArray(i_initRigidBodyVel, 2); required_params.erase("initRigidBodyVel"); }},
        {"initRigidBodyLoc", [&](){ readArray(i_initRigidBodyLoc, 2); required_params.erase("initRigidBodyLoc"); }},
        {"rigidBodyRadiusOrLength", [&](){ readDouble(i_rigidBodyRadiusOrLength); required_params.erase("rigidBodyRadiusOrLength"); }},
        {"rigidBodyAdditionalFactor", [&](){ readDouble(i_rigidBodyAdditionalFactor); }},
        {"leftState", [&](){ readArray(i_leftState, 4); required_params.erase("leftState"); }},
        {"rightState", [&](){ readArray(i_rightState, 4); required_params.erase("rightState"); }},
        {"leftRightStateBoundary", [&](){ 
            readDouble(i_leftRightStateBoundary); 
            if (i_leftRightStateBoundary <= i_x0 || i_leftRightStateBoundary >= i_x1) {
                throw std::runtime_error("leftRightStateBoundary must be between x0 and x1");
            }
            required_params.erase("leftRightStateBoundary"); 
        }},
        {"runName", [&](){ readString(i_runName); }}
    };

    std::string parameter;
    while (in >> parameter) {
        try {
            auto handler = param_handlers.find(parameter);
            if (handler != param_handlers.end()) {
                handler->second();
            }
        } catch (const std::exception &e) {
            throw std::runtime_error("Error reading parameter '" + parameter + "': " + e.what());
        }
    }
    in.close();

    if (!required_params.empty()) {
        std::string missing_params;
        for (const auto& param : required_params) {
            missing_params += param + ", ";
        }
        missing_params.erase(missing_params.length() - 2);
        throw std::runtime_error("Missing required parameters in config.txt: " + missing_params);
    }
}

double singleCircleLevelSetFunc(double i_r, double i_centre_x, double i_centre_y, double i_x, double i_y)
{
    return -(i_r - pow((pow(i_x - i_centre_x, 2) + pow(i_y - i_centre_y, 2)), 0.5));
}

double doubleCircleLevelSetFunc(double i_r1, double i_r2, double i_centre_x1, double i_centre_x2, double i_centre_y1, double i_centre_y2, double i_x, double i_y)
{
    double phi1 = -(i_r1 - pow((pow(i_x - i_centre_x1, 2) + pow(i_y - i_centre_y1, 2)), 0.5));
    double phi2 = -(i_r2 - pow((pow(i_x - i_centre_x2, 2) + pow(i_y - i_centre_y2, 2)), 0.5));

    return std::min(phi1, phi2);
}

double singleSqaureLevelSetFunc(double i_l, double i_centre_x, double i_centre_y, double i_x, double i_y)
{
    double relative_x_coor = i_x - i_centre_x;
    double relative_y_coor = i_y - i_centre_y;

    double phi;

    if (relative_x_coor != 0.0)
    {
        if (relative_x_coor > 0)
        {
            phi = atan(relative_y_coor / relative_x_coor);
            if (phi < PI / 4 && phi > -PI / 4) // east
            {
                return -(i_l / 2 - relative_x_coor);
            }
            else if (phi <= -PI / 4) // south
            {
                return -(i_l / 2 + relative_y_coor);
            }
            else // north
            {
                return -(i_l / 2 - relative_y_coor);
            }
        }
        else
        {
            phi = atan(-relative_y_coor / relative_x_coor);
            if (phi < PI / 4 && phi > -PI / 4) // west
            {
                return -(i_l / 2 + relative_x_coor);
            }
            else if (phi <= -PI / 4) // south
            {
                return -(i_l / 2 + relative_y_coor);
            }
            else // north
            {
                return -(i_l / 2 - relative_y_coor);
            }
        }
    }
    else
    {
        if (relative_y_coor > 0) // north
        {
            return -(i_l / 2 - relative_y_coor);
        }
        else // south
        {
            return -(i_l / 2 + relative_y_coor);
        }
    }
}

double oneD_WallLevelSetFunc(double i_centre_x, double i_x, double i_y)
{
    return -i_x + i_centre_x;
}

void conservativeFormTransform(std::vector<std::vector<std::array<double, 4>>> &i_inputVec)
{
    double localGamma = 1.4;
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    double velX_Copy;
    double velY_Copy;
    double pressureCopy;

    for (int j = 0; j < nCell_y + 4; ++j)
    {
        for (int i = 0; i < nCell_x + 4; ++i)
        {
            velX_Copy = i_inputVec[j][i][1];
            velY_Copy = i_inputVec[j][i][2];
            pressureCopy = i_inputVec[j][i][3];

            i_inputVec[j][i][1] = velX_Copy * i_inputVec[j][i][0];
            i_inputVec[j][i][2] = velY_Copy * i_inputVec[j][i][0];
            i_inputVec[j][i][3] = pressureCopy / (localGamma - 1) + 0.5 * i_inputVec[j][i][0] * (velX_Copy * velX_Copy + velY_Copy * velY_Copy);
        }
    }
}

// #########################################################################################################################################################################################
// #########################################################################################################################################################################################

// manual initial state setting
void setInitialConditions(std::vector<std::vector<std::array<double, 4>>> &i_inputVec, std::vector<std::vector<double>> &i_levelSetFunc, double i_x0, double i_x1, double i_y0, double i_y1)
{
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    double dx = (i_x1 - i_x0) / nCell_x;
    double dy = (i_y1 - i_y0) / nCell_y;

    double currX;
    double currY;

    for (int j = 0; j < nCell_y + 4; ++j)
    {
        for (int i = 0; i < nCell_x + 4; ++i)
        {
            currX = i_x0 + (i - 2) * dx;
            currY = i_y0 + (j - 2) * dy;
            i_levelSetFunc[j][i] = 999;
        }
    }

    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            currX = i_x0 + (i - 2) * dx;
            currY = i_y0 + (j - 2) * dy;

            // if (currX <= 0.5)
            if (0.4 - sqrt((currX - 1) * (currX - 1) + (currY - 1) * (currY - 1)) > 0)
            {
                i_inputVec[j][i] = (std::array<double, 4>){1, 0.0, 0.0, 1};
            }
            else
            {
                i_inputVec[j][i] = (std::array<double, 4>){0.125, 0.0, 0.0, 0.1};
            }
        }
    }
}

// initial conditions read from config.txt
void setInitialConditions(std::vector<std::vector<std::array<double, 4>>> &i_inputVec, std::vector<std::vector<double>> &i_levelSetFunc, double i_x0, double i_x1, double i_y0, double i_y1,
                          int i_rigidBodyType, double i_rbRadiusOrLen, double i_rbAdditionalFactor, std::array<double, 2> i_initRigidBodyLoc,
                          std::array<double, 4> i_leftState, std::array<double, 4> i_righState, double i_lrStateBoundary)
{
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;
    double dx = (i_x1 - i_x0) / nCell_x;
    double dy = (i_y1 - i_y0) / nCell_y;
    double currX;
    double currY;
    for (int j = 0; j < nCell_y + 4; ++j)
    {
        for (int i = 0; i < nCell_x + 4; ++i)
        {
            currX = i_x0 + (i - 2) * dx;
            currY = i_y0 + (j - 2) * dy;
            if (i_rigidBodyType == 1)
            {
                i_levelSetFunc[j][i] = singleCircleLevelSetFunc(i_rbRadiusOrLen, i_initRigidBodyLoc[0], i_initRigidBodyLoc[1], currX, currY);
            }
            else if (i_rigidBodyType == 2)
            {
                i_levelSetFunc[j][i] = singleSqaureLevelSetFunc(i_rbRadiusOrLen, i_initRigidBodyLoc[0], i_initRigidBodyLoc[1], currX, currY);
            }
            else if (i_rigidBodyType == 3)
            {
                i_levelSetFunc[j][i] = doubleCircleLevelSetFunc(i_rbRadiusOrLen, i_rbRadiusOrLen,
                                                                i_initRigidBodyLoc[0], i_initRigidBodyLoc[0],
                                                                i_initRigidBodyLoc[1] - i_rbAdditionalFactor,
                                                                i_initRigidBodyLoc[1] + i_rbAdditionalFactor,
                                                                currX, currY);
            }
            else
            {
                throw std::invalid_argument("rigidBodyType can only be 1, 2 or 3");
            }
        }
    }

    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            currX = i_x0 + (i - 2) * dx;
            currY = i_y0 + (j - 2) * dy;
            if (currX <= i_lrStateBoundary)
            {
                i_inputVec[j][i] = i_leftState;
            }
            else
            {
                i_inputVec[j][i] = i_righState;
            }
        }
    }
}

int main()
{
    // these simulation variables are to be read from the ./cofig.txt file
    int nCells_x;
    int nCells_y;
    double x0;
    double x1;
    double y0;
    double y1;
    double c;
    double tStop;
    int loggingFactor;
    int reinitFactor;
    int rigidBodyType;
    int acceleratingRigidBody;
    std::array<double, 2> initRigidBodyVel;
    std::array<double, 2> initRigidBodyLoc;
    double rigidBodyRadiusOrLength;
    double rigidBodyAdditionalFactor;
    std::array<double, 4> leftState;
    std::array<double, 4> rightState;
    double leftRightStateBoundary;
    std::string runName;
    std::string repoDir;

    loadConfig(nCells_x, nCells_y, x0, x1, y0, y1, c, tStop, loggingFactor,
               reinitFactor, rigidBodyType, acceleratingRigidBody, initRigidBodyVel,
               initRigidBodyLoc, rigidBodyRadiusOrLength, rigidBodyAdditionalFactor,
               leftState, rightState, leftRightStateBoundary, runName, repoDir);

    // size the computational domain and the level set grid by {nCells_x, nCells_y}
    std::vector<std::vector<std::array<double, 4>>> compDomain;
    std::vector<std::vector<double>> levelSetCompDomain;
    compDomain.resize(nCells_y + 4);
    levelSetCompDomain.resize(nCells_y + 4);
    for (int j = 0; j < nCells_y + 4; ++j)
    {
        compDomain[j].resize(nCells_x + 4);
        levelSetCompDomain[j].resize(nCells_x + 4);
    }

    // setInitialConditions(compDomain, levelSetCompDomain, x0, x1, y0, y1);
    setInitialConditions(compDomain, levelSetCompDomain, x0, x1, y0, y1, rigidBodyType,
                         rigidBodyRadiusOrLength, rigidBodyAdditionalFactor, initRigidBodyLoc,
                         leftState, rightState, leftRightStateBoundary);
    conservativeFormTransform(compDomain);

    // constructing solver class and initiate with read states
    GFM_2D_EulerSolver gfmSolverClass(compDomain, nCells_x, nCells_y);
    gfmSolverClass.setBound(x0, x1, y0, y1, tStop);
    gfmSolverClass.setCFL(c);
    gfmSolverClass.setName(runName);
    gfmSolverClass.setRepoDir(repoDir);
    gfmSolverClass.setLevelSet(levelSetCompDomain);
    gfmSolverClass.setRigidBodyVel(initRigidBodyVel);
    gfmSolverClass.setRigidBodyCentreCoor(initRigidBodyLoc);

    gfmSolverClass.updateBoundaryTrans();
    // gfmSolverClass.updateBoundaryReflect();
    gfmSolverClass.calculateMockSchliren();

    gfmSolverClass.initiateOutputLogging();
    std::cout << 1 << ": " << 0 << " / " << gfmSolverClass.tStop() << std::endl;

    double t = 0.0;
    int numIter = 0;
    do
    {
        ++numIter;

        gfmSolverClass.updateMaxA(numIter);
        gfmSolverClass.updateDt();

        t += gfmSolverClass.dt();

        gfmSolverClass.updateGhostCellBoundary();
        gfmSolverClass.propagateGhostCell();

        gfmSolverClass.mhHllcSweepX();
        gfmSolverClass.updateBoundaryTrans();
        // gfmSolverClass.updateBoundaryReflect();
        gfmSolverClass.mhHllcSweepY();
        gfmSolverClass.updateBoundaryTrans();
        // gfmSolverClass.updateBoundaryReflect();

        gfmSolverClass.calculateMockSchliren();

        if (acceleratingRigidBody == 1)
        {
            gfmSolverClass.accelerateRigidBody_circ();
        }
        gfmSolverClass.advectLevelSet();
        gfmSolverClass.updateLevelSetBoundaryTrans();

        if (numIter % reinitFactor == 0)
        {
            gfmSolverClass.reinitLevelSet();
        }

        if (numIter % loggingFactor == 0)
        {
            gfmSolverClass.writeToFiles(t);
            std::cout << numIter / loggingFactor + 1 << ": " << t << " / " << gfmSolverClass.tStop() << std::endl;
        }

    } while (t < gfmSolverClass.tStop());

    gfmSolverClass.cleanUp();

    return 0;
}
