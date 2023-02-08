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
#include <math.h>

#define PI 3.14159265

void loadConfig(int &i_nCells_x, int &i_nCells_y, double &i_x0, double &i_x1, double &i_y0, double &i_y1, double &i_c, double &i_tStop, int &i_loggingFactor)
{
    std::ifstream in("./config.txt");

    std::string parameter;
    double double_value;
    int int_value;

    while (!in.eof())
    {
        in >> parameter;
        if (parameter == "nCells_x")
        {
            in >> int_value;
            i_nCells_x = int_value;
        }
        if (parameter == "nCells_y")
        {
            in >> int_value;
            i_nCells_y = int_value;
        }
        if (parameter == "x0")
        {
            in >> double_value;
            i_x0 = double_value;
        }
        if (parameter == "x1")
        {
            in >> double_value;
            i_x1 = double_value;
        }
        if (parameter == "y0")
        {
            in >> double_value;
            i_y0 = double_value;
        }
        if (parameter == "y1")
        {
            in >> double_value;
            i_y1 = double_value;
        }
        if (parameter == "c")
        {
            in >> double_value;
            i_c = double_value;
        }
        if (parameter == "tStop")
        {
            in >> double_value;
            i_tStop = double_value;
        }
        if (parameter == "loggingFactor")
        {
            in >> int_value;
            i_loggingFactor = int_value;
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

    if (abs(phi1) < abs(phi2))
    {
        return phi1;
    }
    else
    {
        return phi2;
    }
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

// #############################################################################################################################################################################################
// things up there are probably fine -----------------------------------------------------------------------------------------------------------------------------------------------------------
// #############################################################################################################################################################################################

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

            // i_levelSetFunc[j][i] = singleCircleLevelSetFunc(0.2, 0.6, 0.5, currX, currY);
            // i_levelSetFunc[j][i] = singleSqaureLevelSetFunc(0.4, 0.6, 0.5, currX, currY);
            // i_levelSetFunc[j][i] = doubleCircleLevelSetFunc(0.2, 0.2, 0.6, 0.6, 0.25, 0.75, currX, currY);
            // i_levelSetFunc[j][i] = doubleCircleLevelSetFunc(0.2, 0.2, 0.6, 0.6, 0.35, 0.65, currX, currY);
            // i_levelSetFunc[j][i] = oneD_WallLevelSetFunc(0.6, currX, currY);

            i_levelSetFunc[j][i] = singleCircleLevelSetFunc(0.2, 2.5, 0.5, currX, currY);
        }
    }

    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            currX = i_x0 + (i - 2) * dx;
            currY = i_y0 + (j - 2) * dy;

            if (currX <= 0.2)
            {
                // i_inputVec[j][i] = (std::array<double, 4>){1.3764, 0.394, 0.0, 1.5698};
                i_inputVec[j][i] = (std::array<double, 4>){1, 0.0, 0.0, 1};
            }
            else
            {
                i_inputVec[j][i] = (std::array<double, 4>){1, 0.0, 0.0, 1};
            }
        }
    }
}

int main()
{
    int nCells_x;
    int nCells_y;
    double x0;
    double x1;
    double y0;
    double y1;
    double c;
    double tStop;

    int loggingFactor;

    loadConfig(nCells_x, nCells_y, x0, x1, y0, y1, c, tStop, loggingFactor);
    std::vector<std::vector<std::array<double, 4>>> compDomain;
    std::vector<std::vector<double>> levelSetCompDomain;
    compDomain.resize(nCells_y + 4);
    levelSetCompDomain.resize(nCells_y + 4);
    for (int j = 0; j < nCells_y + 4; ++j)
    {
        compDomain[j].resize(nCells_x + 4);
        levelSetCompDomain[j].resize(nCells_x + 4);
    }

    setInitialConditions(compDomain, levelSetCompDomain, x0, x1, y0, y1);
    conservativeFormTransform(compDomain);

    GFM_2D_EulerSolver testSolverClass(compDomain, nCells_x, nCells_y);
    testSolverClass.setBound(x0, x1, y0, y1, tStop);
    testSolverClass.setCFL(c);
    testSolverClass.setName((std::string) "test");
    testSolverClass.setRepoDir((std::string) "/Users/zianhuang/Room214N/dev/mphil/MPhil_writtenAssignment_GFM/");
    testSolverClass.setLevelSet(levelSetCompDomain);

    testSolverClass.setRigidBodyVel(std::array<double, 2>{-1.0, 0.0});

    testSolverClass.updateBoundaryTrans();

    testSolverClass.calculateMockSchliren();

    testSolverClass.initiateDataLogging();

    // DEBUG
    // testSolverClass.printBoundaryCoor();
    // printLevelSet(testSolverClass.levelSet());

    double t = 0.0;
    int numIter = 0;
    do
    {
        // DEBUG
        // printDomainDensity(testSolverClass.uVec());
        // printLevelSet(testSolverClass.levelSet());

        ++numIter;

        testSolverClass.updateMaxA(numIter);
        testSolverClass.updateDt();

        t += testSolverClass.dt();

        testSolverClass.updateGhostCellBoundary(true);
        // testSolverClass.updateGhostCellBoundary(false);
        testSolverClass.propagateGhostCell();

        testSolverClass.mhHllcSweepX();
        testSolverClass.updateBoundaryTrans();

        testSolverClass.mhHllcSweepY();
        testSolverClass.updateBoundaryTrans();

        testSolverClass.calculateMockSchliren();

        testSolverClass.advectLevelSet();
        testSolverClass.updateLevelSetBoundaryTrans();

        if (numIter % loggingFactor == 0)
        {
            testSolverClass.writeToFiles(t);
            std::cout << t << " / " << testSolverClass.tStop() << std::endl;
        }

    } while (t < testSolverClass.tStop());

    testSolverClass.cleanUp();

    return 0;
}
