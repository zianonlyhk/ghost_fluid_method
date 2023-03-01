/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   exec_interface.cc                                 Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:54:13 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "inline/debug_tools.hh"
#include "gfm_2d_euler_solver.hh"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <math.h>

#define PI 3.14159265

void loadConfig(int &i_nCells_x, int &i_nCells_y, double &i_x0, double &i_x1, double &i_y0, double &i_y1, double &i_c, double &i_tStop, int &i_loggingFactor, int &i_reinitFactor,
                int &i_rigidBodyType, int &i_acceleratingRigidBody, std::array<double, 2> &i_initRigidBodyVel, std::array<double, 2> &i_initRigidBodyLoc,
                double &i_rigidBodyRadiusOrLength, double &i_rigidBodyAdditionalFactor,
                std::array<double, 4> &i_leftState, std::array<double, 4> &i_rightState, double &i_leftRightStateBoundary,
                std::string &i_runName, std::string &i_repoDir)
{
    std::ifstream in("./config.txt");

    std::string parameter;

    // dummie disposable variables
    double double_value;
    int int_value;
    std::string string_value;
    std::array<double, 2> initRigidBody_value;
    std::array<double, 4> leftRightState_value;

    while (!in.eof())
    {
        // if no read string fits the required parameter, do nothing
        in >> parameter;

        if (parameter == "nCells_x")
        {
            in >> int_value;
            i_nCells_x = int_value;
        }
        else if (parameter == "nCells_y")
        {
            in >> int_value;
            i_nCells_y = int_value;
        }
        else if (parameter == "x0")
        {
            in >> double_value;
            i_x0 = double_value;
        }
        else if (parameter == "x1")
        {
            in >> double_value;
            i_x1 = double_value;
        }
        else if (parameter == "y0")
        {
            in >> double_value;
            i_y0 = double_value;
        }
        else if (parameter == "y1")
        {
            in >> double_value;
            i_y1 = double_value;
        }
        else if (parameter == "c")
        {
            in >> double_value;
            i_c = double_value;
        }
        else if (parameter == "tStop")
        {
            in >> double_value;
            i_tStop = double_value;
        }
        else if (parameter == "loggingFactor")
        {
            in >> int_value;
            i_loggingFactor = int_value;
        }
        else if (parameter == "reinitFactor")
        {
            in >> int_value;
            i_reinitFactor = int_value;
        }
        else if (parameter == "rigidBodyType")
        {
            in >> int_value;
            i_rigidBodyType = int_value;
        }
        else if (parameter == "acceleratingRigidBody")
        {
            in >> int_value;
            i_acceleratingRigidBody = int_value;
        }
        else if (parameter == "initRigidBodyVel")
        {
            in >> double_value;
            i_initRigidBodyVel[0] = double_value;
            in >> double_value;
            i_initRigidBodyVel[1] = double_value;
        }
        else if (parameter == "initRigidBodyLoc")
        {
            in >> double_value;
            i_initRigidBodyLoc[0] = double_value;
            in >> double_value;
            i_initRigidBodyLoc[1] = double_value;
        }
        else if (parameter == "rigidBodyRadiusOrLength")
        {
            in >> double_value;
            i_rigidBodyRadiusOrLength = double_value;
        }
        else if (parameter == "rigidBodyAdditionalFactor")
        {
            in >> double_value;
            i_rigidBodyAdditionalFactor = double_value;
        }
        else if (parameter == "leftState")
        {
            in >> double_value;
            i_leftState[0] = double_value;
            in >> double_value;
            i_leftState[1] = double_value;
            in >> double_value;
            i_leftState[2] = double_value;
            in >> double_value;
            i_leftState[3] = double_value;
        }
        else if (parameter == "rightState")
        {
            in >> double_value;
            i_rightState[0] = double_value;
            in >> double_value;
            i_rightState[1] = double_value;
            in >> double_value;
            i_rightState[2] = double_value;
            in >> double_value;
            i_rightState[3] = double_value;
        }
        else if (parameter == "leftRightStateBoundary")
        {
            in >> double_value;
            i_leftRightStateBoundary = double_value;
        }
        else if (parameter == "runName")
        {
            in >> string_value;
            i_runName = string_value;
        }
        else if (parameter == "repoDir")
        {
            in >> string_value;
            i_repoDir = string_value;
        }
    }
    in.close();
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
    GFM_2D_EulerSolver testSolverClass(compDomain, nCells_x, nCells_y);
    testSolverClass.setBound(x0, x1, y0, y1, tStop);
    testSolverClass.setCFL(c);
    testSolverClass.setName(runName);
    testSolverClass.setRepoDir(repoDir);
    testSolverClass.setLevelSet(levelSetCompDomain);
    testSolverClass.setRigidBodyVel(initRigidBodyVel);
    testSolverClass.setRigidBodyCentreCoor(initRigidBodyLoc);

    testSolverClass.updateBoundaryTrans();
    // testSolverClass.updateBoundaryReflect();
    testSolverClass.calculateMockSchliren();

    testSolverClass.initiateDataLogging();
    std::cout << 1 << ": " << 0 << " / " << testSolverClass.tStop() << std::endl;

    double t = 0.0;
    int numIter = 0;
    do
    {
        ++numIter;

        testSolverClass.updateMaxA(numIter);
        testSolverClass.updateDt();

        t += testSolverClass.dt();

        testSolverClass.updateGhostCellBoundary();
        testSolverClass.propagateGhostCell();

        testSolverClass.mhHllcSweepX();
        testSolverClass.updateBoundaryTrans();
        // testSolverClass.updateBoundaryReflect();
        testSolverClass.mhHllcSweepY();
        testSolverClass.updateBoundaryTrans();
        // testSolverClass.updateBoundaryReflect();

        testSolverClass.calculateMockSchliren();

        if (acceleratingRigidBody == 1)
        {
            testSolverClass.accelerateRigidBody_circ();
        }
        testSolverClass.advectLevelSet();
        testSolverClass.updateLevelSetBoundaryTrans();

        if (numIter % reinitFactor == 0)
        {
            testSolverClass.reinitLevelSet();
        }

        if (numIter % loggingFactor == 0)
        {
            testSolverClass.writeToFiles(t);
            std::cout << numIter / loggingFactor + 1 << ": " << t << " / " << testSolverClass.tStop() << std::endl;
        }

    } while (t < testSolverClass.tStop());

    testSolverClass.cleanUp();

    return 0;
}
