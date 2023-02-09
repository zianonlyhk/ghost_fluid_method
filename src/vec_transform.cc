/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.cc                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:53:40 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "vec_transform.hh"
#include "inline/cell_operation.hh"
#include "inline/debug_tools.hh"
#include "inline/primitive_tran.hh"
#include <iostream>
#include <algorithm>
#include <math.h>

// Definitions #####################################################################################

VecTran::VecTran() {}

std::vector<std::vector<std::array<double, 4>>> VecTran::slicVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dx, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::array<double, 4> fluxNext;
    std::array<double, 4> fluxBefore;

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            fluxNext = fluxFunc.slicFlux_x(i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_inputU_Vec[iter_y][iter_x + 2], i_dx, i_dt);
            fluxBefore = fluxFunc.slicFlux_x(i_inputU_Vec[iter_y][iter_x - 2], i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_dx, i_dt);
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dx, diffCell(fluxNext, fluxBefore)));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::slicVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dy, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::array<double, 4> fluxNext;
    std::array<double, 4> fluxBefore;

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            fluxNext = fluxFunc.slicFlux_y(i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_inputU_Vec[iter_y + 2][iter_x], i_dy, i_dt);
            fluxBefore = fluxFunc.slicFlux_y(i_inputU_Vec[iter_y - 2][iter_x], i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_dy, i_dt);
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dy, diffCell(fluxNext, fluxBefore)));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::musclHancockVecTranHLLC_x(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dx, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::array<double, 4> fluxNext;
    std::array<double, 4> fluxBefore;

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            fluxNext = fluxFunc.musclHancockHllcFlux_x(i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_inputU_Vec[iter_y][iter_x + 2], i_dx, i_dt);
            fluxBefore = fluxFunc.musclHancockHllcFlux_x(i_inputU_Vec[iter_y][iter_x - 2], i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_dx, i_dt);
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dx, diffCell(fluxNext, fluxBefore)));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::musclHancockVecTranHLLC_y(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dy, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::array<double, 4> fluxNext;
    std::array<double, 4> fluxBefore;

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            fluxNext = fluxFunc.musclHancockHllcFlux_y(i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_inputU_Vec[iter_y + 2][iter_x], i_dy, i_dt);
            fluxBefore = fluxFunc.musclHancockHllcFlux_y(i_inputU_Vec[iter_y - 2][iter_x], i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_dy, i_dt);
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dy, diffCell(fluxNext, fluxBefore)));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::ghostCellBoundary(const std::vector<std::vector<std::array<double, 4>>> &i_uVec, const std::vector<std::vector<double>> &i_levelSet, double i_dx, double i_dy)
{
    int xVecLen = i_uVec[0].size();
    int yVecLen = i_uVec.size();
    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    toBeReturnVec = i_uVec;

    std::vector<std::array<int, 2>> boundaryCoorArr = getBoundaryCellCoor(i_levelSet);

    for (int i = 0; i < boundaryCoorArr.size(); ++i)
    {
        toBeReturnVec[boundaryCoorArr[i][1]][boundaryCoorArr[i][0]] = ghostFluidUtilities.ghostCellValues(i_levelSet, i_uVec, std::array<int, 2>{boundaryCoorArr[i][0], boundaryCoorArr[i][1]}, i_dx, i_dy);
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::ghostCellBoundary(const std::vector<std::vector<std::array<double, 4>>> &i_uVec, const std::vector<std::vector<double>> &i_levelSet, double i_dx, double i_dy, std::array<double, 2> i_rigidBodyVel)
{
    int xVecLen = i_uVec[0].size();
    int yVecLen = i_uVec.size();
    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    toBeReturnVec = i_uVec;

    std::vector<std::array<int, 2>> boundaryCoorArr = getBoundaryCellCoor(i_levelSet);

    for (int i = 0; i < boundaryCoorArr.size(); ++i)
    {
        toBeReturnVec[boundaryCoorArr[i][1]][boundaryCoorArr[i][0]] = ghostFluidUtilities.ghostCellValues(i_levelSet, i_uVec, std::array<int, 2>{boundaryCoorArr[i][0], boundaryCoorArr[i][1]}, i_dx, i_dy, i_rigidBodyVel);
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::propagateGhostInterface(const std::vector<std::vector<std::array<double, 4>>> &i_uVec, const std::vector<std::vector<double>> &i_levelSet, double i_dx, double i_dy)
{
    // DEBUG
    // printDomainDensity(i_uVec);

    int xVecLen = i_uVec[0].size();
    int yVecLen = i_uVec.size();
    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }
    // potential bug here
    toBeReturnVec = i_uVec;

    fastSweepingConstantPropagation(toBeReturnVec, i_uVec, i_levelSet, i_dx, i_dy);
    // fastSweepingConstantPropagation(toBeReturnVec, i_uVec, i_levelSet, i_dx, i_dy);
    // fastSweepingConstantPropagation(toBeReturnVec, i_uVec, i_levelSet, i_dx, i_dy);

    return toBeReturnVec;
}

std::vector<std::vector<double>> VecTran::mockSchlierenTrans(const std::vector<std::vector<std::array<double, 4>>> &i_uVec, double i_dx, double i_dy)
{
    int xVecLen = i_uVec[0].size();
    int yVecLen = i_uVec.size();
    std::vector<std::vector<double>> toBeReturnMsVec;
    toBeReturnMsVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnMsVec[i].resize(xVecLen);
    }

    double gradRho_x;
    double gradRho_y;
    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            gradRho_x = (i_uVec[iter_y][iter_x + 1][0] - i_uVec[iter_y][iter_x - 1][0]) / 2 / i_dx;
            gradRho_y = (i_uVec[iter_y + 1][iter_x][0] - i_uVec[iter_y - 1][iter_x][0]) / 2 / i_dy;

            toBeReturnMsVec[iter_y][iter_x] = exp((-20 * pow(gradRho_x * gradRho_x + gradRho_y * gradRho_y, 0.5)) / 1000 / i_uVec[iter_y][iter_x][0]);
        }
    }

    return toBeReturnMsVec;
}

std::vector<std::vector<double>> VecTran::levelSetAdvectionTransform(const std::vector<std::vector<double>> &i_levelSet, std::array<double, 2> i_rigidBodyVel, double i_dx, double i_dy, double i_dt)
{
    int xVecLen = i_levelSet[0].size();
    int yVecLen = i_levelSet.size();

    std::vector<std::vector<double>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    toBeReturnVec = i_levelSet;

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            if (i_rigidBodyVel[0] < 0)
            {
                toBeReturnVec[iter_y][iter_x] = i_levelSet[iter_y][iter_x] - i_rigidBodyVel[0] * i_dt / i_dx * (i_levelSet[iter_y][iter_x + 1] - i_levelSet[iter_y][iter_x]);
            }
            else
            {
                toBeReturnVec[iter_y][iter_x] = i_levelSet[iter_y][iter_x] + i_rigidBodyVel[0] * i_dt / i_dx * (i_levelSet[iter_y][iter_x - 1] - i_levelSet[iter_y][iter_x]);
            }
        }
    }

    std::vector<std::vector<double>> toBeReturnVecAgain;
    toBeReturnVecAgain.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVecAgain[i].resize(xVecLen);
    }

    toBeReturnVecAgain = toBeReturnVec;

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            if (i_rigidBodyVel[1] < 0)
            {
                toBeReturnVecAgain[iter_y][iter_x] = toBeReturnVec[iter_y][iter_x] - i_rigidBodyVel[1] * i_dt / i_dy * (toBeReturnVec[iter_y + 1][iter_x] - toBeReturnVec[iter_y][iter_x]);
            }
            else
            {
                toBeReturnVecAgain[iter_y][iter_x] = toBeReturnVec[iter_y][iter_x] + i_rigidBodyVel[1] * i_dt / i_dy * (toBeReturnVec[iter_y - 1][iter_x] - toBeReturnVec[iter_y][iter_x]);
            }
        }
    }

    return toBeReturnVecAgain;
}

std::vector<std::vector<double>> VecTran::reinitLevelSetInsideRB(const std::vector<std::vector<double>> &i_levelSet, double i_dx, double i_dy)
{
    int xVecLen = i_levelSet[0].size();
    int yVecLen = i_levelSet.size();

    std::vector<std::vector<double>> toBeReturnLevelSet;
    toBeReturnLevelSet.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnLevelSet[i].resize(xVecLen);
    }

    toBeReturnLevelSet = i_levelSet;

    bool insideRigidBody;

    std::vector<std::array<int, 2>> boundaryCoorArr = getBoundaryCellCoor(i_levelSet);

    bool updateLocalLevelSet;

    // +ve x sweep
    for (int j = 2; j < yVecLen - 2; ++j)
    {
        insideRigidBody = false;
        for (int i = 2; i < xVecLen - 2; ++i)
        {
            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                insideRigidBody = true;
                continue;
            }

            if (insideRigidBody)
            {
                if (i_levelSet[j][i + 1] > 0)
                {
                    break;
                }
                updateLocalLevelSet = true;

                for (int dummie = 0; dummie < boundaryCoorArr.size(); ++dummie)
                {
                    if (boundaryCoorArr[dummie][0] == i && boundaryCoorArr[dummie][1] == j)
                    {
                        updateLocalLevelSet = false;
                        break;
                    }
                }
                if (updateLocalLevelSet)
                {
                    toBeReturnLevelSet[j][i] = ghostFluidUtilities.solveForLevelSetReinit(i_levelSet, std::array<int, 2>{i, j}, i_dx, i_dy);
                }
            }
        }
    }

    // -ve x sweep
    for (int j = 2; j < yVecLen - 2; ++j)
    {
        insideRigidBody = false;
        for (int i = xVecLen - 3; i > 1; --i)
        {
            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                insideRigidBody = true;
                continue;
            }

            if (insideRigidBody)
            {
                if (i_levelSet[j][i + 1] > 0)
                {
                    break;
                }
                updateLocalLevelSet = true;

                for (int dummie = 0; dummie < boundaryCoorArr.size(); ++dummie)
                {
                    if (boundaryCoorArr[dummie][0] == i && boundaryCoorArr[dummie][1] == j)
                    {
                        updateLocalLevelSet = false;
                        break;
                    }
                }
                if (updateLocalLevelSet)
                {
                    toBeReturnLevelSet[j][i] = ghostFluidUtilities.solveForLevelSetReinit(i_levelSet, std::array<int, 2>{i, j}, i_dx, i_dy);
                }
            }
        }
    }

    for (int i = 2; i < xVecLen - 2; ++i)
    {
        insideRigidBody = false;
        for (int j = 2; j < yVecLen - 2; ++j)
        {
            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                insideRigidBody = true;
                continue;
            }

            if (insideRigidBody)
            {
                if (i_levelSet[j][i + 1] > 0)
                {
                    break;
                }
                updateLocalLevelSet = true;

                for (int dummie = 0; dummie < boundaryCoorArr.size(); ++dummie)
                {
                    if (boundaryCoorArr[dummie][0] == i && boundaryCoorArr[dummie][1] == j)
                    {
                        updateLocalLevelSet = false;
                        break;
                    }
                }
                if (updateLocalLevelSet)
                {
                    toBeReturnLevelSet[j][i] = ghostFluidUtilities.solveForLevelSetReinit(i_levelSet, std::array<int, 2>{i, j}, i_dx, i_dy);
                }
            }
        }
    }

    for (int i = 2; i < xVecLen - 2; ++i)
    {
        insideRigidBody = false;
        for (int j = yVecLen - 3; j > 1; --j)
        {
            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                insideRigidBody = true;
                continue;
            }

            if (insideRigidBody)
            {
                if (i_levelSet[j][i + 1] > 0)
                {
                    break;
                }
                updateLocalLevelSet = true;

                for (int dummie = 0; dummie < boundaryCoorArr.size(); ++dummie)
                {
                    if (boundaryCoorArr[dummie][0] == i && boundaryCoorArr[dummie][1] == j)
                    {
                        updateLocalLevelSet = false;
                        break;
                    }
                }
                if (updateLocalLevelSet)
                {
                    toBeReturnLevelSet[j][i] = ghostFluidUtilities.solveForLevelSetReinit(i_levelSet, std::array<int, 2>{i, j}, i_dx, i_dy);
                }
            }
        }
    }

    return toBeReturnLevelSet;
}

std::vector<std::array<int, 2>> VecTran::getBoundaryCellCoor(const std::vector<std::vector<double>> &i_levelSet)
{
    return ghostFluidUtilities.ghostBoundaryCellCoor(i_levelSet);
}

void VecTran::fastSweepingConstantPropagation(std::vector<std::vector<std::array<double, 4>>> &i_toBeReturned, const std::vector<std::vector<std::array<double, 4>>> &i_uVec, const std::vector<std::vector<double>> &i_levelSet, double i_dx, double i_dy)
{
    int xVecLen = i_uVec[0].size();
    int yVecLen = i_uVec.size();
    bool insideRigidBody;
    double maxPhi;
    std::array<double, 4> tempArr;

    // +ve x axis sweep

    // DEBUG
    // std::cout << "starting to sweeping the +ve x direction" << std::endl;

    for (int j = 2; j < yVecLen - 2; ++j)
    {
        insideRigidBody = false;

        for (int i = 2; i < xVecLen - 2; ++i)
        {
            // DEBUG
            // std::cout << "at coor " << i << ", " << j << std::endl;

            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                // DEBUG
                // std::cout << "hitting the interface at coor " << i << ", " << j << std::endl;

                insideRigidBody = true;
                maxPhi = abs(i_levelSet[j][i]);
                continue;
            }

            if (insideRigidBody) // starting to propagate each cell
            {
                // only solve if the next phi is bigger so not maxed yet
                if (abs(i_levelSet[j][i]) > maxPhi)
                {
                    // DEBUG
                    // std::cout << "solving at coor " << i << ", " << j << std::endl;

                    maxPhi = abs(i_levelSet[j][i]);
                    tempArr = ghostFluidUtilities.solveForConstantExtrapolation(i_levelSet, i_uVec, std::array<int, 2>{i, j}, i_dx, i_dy);

                    // if (abs(tempArr[0]) > abs(i_uVec[j][i][0]))
                    // {
                    //     tempArr[0] = i_uVec[j][i][0];
                    // }
                    // if (abs(tempArr[1]) > abs(i_uVec[j][i][1]))
                    // {
                    //     tempArr[1] = i_uVec[j][i][1];
                    // }
                    // if (abs(tempArr[2]) > abs(i_uVec[j][i][2]))
                    // {
                    //     tempArr[2] = i_uVec[j][i][2];
                    // }
                    // if (abs(tempArr[3]) > abs(i_uVec[j][i][3]))
                    // {
                    //     tempArr[3] = i_uVec[j][i][3];
                    // }

                    std::copy(std::begin(tempArr), std::end(tempArr), std::begin(i_toBeReturned[j][i]));
                }
                else
                {
                    // DEBUG
                    // std::cout << "just crossed the max, breaking" << std::endl;
                    break;
                }
            }
        }
    }

    // -ve x axis sweep

    // DEBUG
    // std::cout << "starting to sweeping the -ve x direction" << std::endl;

    for (int j = 2; j < yVecLen - 2; ++j)
    {
        insideRigidBody = false;

        for (int i = xVecLen - 3; i > 1; --i)
        {
            // DEBUG
            // std::cout << "at coor " << i << ", " << j << std::endl;

            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                // DEBUG
                // std::cout << "hitting the interface at coor " << i << ", " << j << std::endl;

                insideRigidBody = true;
                maxPhi = abs(i_levelSet[j][i]);
                continue;
            }

            if (insideRigidBody) // starting to propagate each cell
            {
                // only solve if the next phi is bigger so not maxed yet
                if (abs(i_levelSet[j][i]) > maxPhi)
                {
                    // DEBUG
                    // std::cout << "solving at coor " << i << ", " << j << std::endl;

                    maxPhi = abs(i_levelSet[j][i]);
                    tempArr = ghostFluidUtilities.solveForConstantExtrapolation(i_levelSet, i_uVec, std::array<int, 2>{i, j}, i_dx, i_dy);

                    std::copy(std::begin(tempArr), std::end(tempArr), std::begin(i_toBeReturned[j][i]));
                }
                else
                {
                    // DEBUG
                    // std::cout << "just crossed the max, breaking" << std::endl;
                    break;
                }
            }
        }
    }

    // +ve y axis sweep

    // DEBUG
    // std::cout << "starting to sweeping the +ve y direction" << std::endl;

    for (int i = 2; i < xVecLen - 2; ++i)
    {
        insideRigidBody = false;

        for (int j = 2; j < yVecLen - 2; ++j)
        {
            // DEBUG
            // std::cout << "at coor " << i << ", " << j << std::endl;

            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                // DEBUG
                // std::cout << "hitting the interface at coor " << i << ", " << j << std::endl;

                insideRigidBody = true;
                maxPhi = abs(i_levelSet[j][i]);
                continue;
            }

            if (insideRigidBody) // starting to propagate each cell
            {
                // only solve if the next phi is bigger so not maxed yet
                if (abs(i_levelSet[j][i]) > maxPhi)
                {
                    // DEBUG
                    // std::cout << "solving at coor " << i << ", " << j << std::endl;

                    maxPhi = abs(i_levelSet[j][i]);
                    tempArr = ghostFluidUtilities.solveForConstantExtrapolation(i_levelSet, i_uVec, std::array<int, 2>{i, j}, i_dx, i_dy);

                    std::copy(std::begin(tempArr), std::end(tempArr), std::begin(i_toBeReturned[j][i]));
                }
                else
                {
                    // DEBUG
                    // std::cout << "just crossed the max, breaking" << std::endl;
                    break;
                }
            }
        }
    }

    // -ve y axis sweep

    // DEBUG
    // std::cout << "starting to sweeping the -ve y direction" << std::endl;

    for (int i = 2; i < xVecLen - 2; ++i)
    {
        insideRigidBody = false;

        for (int j = yVecLen - 3; j > 1; --j)
        {
            // DEBUG
            // std::cout << "at coor " << i << ", " << j << std::endl;

            if (i_levelSet[j][i] <= 0 && !insideRigidBody)
            {
                // DEBUG
                // std::cout << "hitting the interface at coor " << i << ", " << j << std::endl;

                insideRigidBody = true;
                maxPhi = abs(i_levelSet[j][i]);
                continue;
            }

            if (insideRigidBody) // starting to propagate each cell
            {
                // only solve if the next phi is bigger so not maxed yet
                if (abs(i_levelSet[j][i]) > maxPhi)
                {
                    // DEBUG
                    // std::cout << "solving at coor " << i << ", " << j << std::endl;

                    maxPhi = abs(i_levelSet[j][i]);
                    tempArr = ghostFluidUtilities.solveForConstantExtrapolation(i_levelSet, i_uVec, std::array<int, 2>{i, j}, i_dx, i_dy);

                    std::copy(std::begin(tempArr), std::end(tempArr), std::begin(i_toBeReturned[j][i]));
                }
                else
                {
                    // DEBUG
                    // std::cout << "just crossed the max, breaking" << std::endl;
                    break;
                }
            }
        }
    }
}