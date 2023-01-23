/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   exec_interface.cc                                 Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/21 10:44:41 by Zian Huang                               */
/*   Updated: 2023/01/23 18:07:12 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "gfm_2d_euler_solver.hh"
#include <fstream>
#include <iostream>
#include <math.h>

void setInitialConditions(std::vector<std::vector<std::array<double, 4>>> &i_inputVec, std::vector<std::vector<double>> &i_levelSetFunc, double i_x0, double i_x1, double i_y0, double i_y1)
{
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    double dx = (i_x1 - i_x0) / nCell_x;
    double dy = (i_x1 - i_x0) / nCell_y;

    double currX;
    double currY;
    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            currX = i_x0 + (i - 2) * dx;
            currY = i_y0 + (j - 2) * dy;

            i_levelSetFunc[j][i] = 0.4 - (pow(currX - 1, 2) + pow(currY - 1, 2));

            if (pow(currX - 1, 2) + pow(currY - 1, 2) < 0.4)
            {
                i_inputVec[j][i] = (std::array<double, 4>){1, 0, 0, 1};
            }
            else
            {
                i_inputVec[j][i] = (std::array<double, 4>){0.125, 0, 0, 0.1};
            }
        }
    }
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
    nCells_x = 100;
    nCells_y = 100;
    x0 = 0.0;
    x1 = 2.0;
    y0 = 0.0;
    y1 = 2.0;
    c = 0.95;
    tStop = 0.25;

    std::vector<std::vector<std::array<double, 4>>> compDomain;
    compDomain.resize(nCells_y + 4);
    for (int i = 0; i < nCells_y + 4; ++i)
    {
        compDomain[i].resize(nCells_x + 4);
    }

    std::vector<std::vector<double>> levelSetcompDomain;
    levelSetcompDomain.resize(nCells_y + 4);
    for (int i = 0; i < nCells_y + 4; ++i)
    {
        levelSetcompDomain[i].resize(nCells_x + 4);
    }

    setInitialConditions(compDomain, levelSetcompDomain, x0, x1, y0, y1);
    conservativeFormTransform(compDomain);

    GFM_2D_EulerSolver testSolverClass(compDomain, nCells_x, nCells_y);
    testSolverClass.setBound(x0, x1, y0, y1, tStop);
    testSolverClass.setCFL(c);
    testSolverClass.setName((std::string) "test");
    testSolverClass.setRepoDir((std::string) "/Users/zianhuang/Room214N/dev/MPhil_writtenAssignment_GFM/");

    testSolverClass.updateBoundary_transmissive();

    testSolverClass.initiateDataLogging();

    double t = 0.0;
    int numIter = 0;
    do
    {
        ++numIter;

        testSolverClass.updateMaxA(numIter);
        testSolverClass.updateDt();

        t += testSolverClass.dt();

        testSolverClass.musclHancock_sweepX();
        testSolverClass.updateBoundary_transmissive();
        testSolverClass.musclHancock_sweepY();
        testSolverClass.updateBoundary_transmissive();

        if (numIter % 1 == 0)
        {
            testSolverClass.writeToFiles(t);
            std::cout << t << " / " << testSolverClass.tStop() << std::endl;
        }

    } while (t < testSolverClass.tStop());

    testSolverClass.cleanUp();

    return 0;
}
