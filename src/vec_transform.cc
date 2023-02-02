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

std::vector<std::vector<std::array<double, 4>>> VecTran::propagateGhostInterface(const std::vector<std::vector<std::array<double, 4>>> &i_uVec, const std::vector<std::vector<double>> &i_levelSet, double i_dx, double i_dy)
{
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

    // DEBUG
    std::cout << "before propagation:" << std::endl;
    printDomainDensity(i_uVec);

    fastSweepingConstantPropagation(toBeReturnVec, i_uVec, i_levelSet, i_dx, i_dy);

    // DEBUG
    std::cout << "after propagation:" << std::endl;
    printDomainDensity(toBeReturnVec);
    std::cout << "##################################################################################################################" << std::endl;

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::fillGhostRegionWithConstant(const std::vector<std::vector<std::array<double, 4>>> &i_uVec, const std::vector<std::vector<double>> &i_levelSet)
{
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

    std::vector<std::array<int, 2>> boundaryCoorArr = getBoundaryCellCoor(i_levelSet);

    bool dummie_bool;
    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            if (i_levelSet[iter_y][iter_x] < 0)
            {
                dummie_bool = false;
                for (int pt = 0; pt < boundaryCoorArr.size(); ++pt)
                {
                    if (iter_x == boundaryCoorArr[pt][0] && iter_y == boundaryCoorArr[pt][1])
                    {
                        dummie_bool = true;
                    }
                }
                if (!dummie_bool)
                {
                    toBeReturnVec[iter_y][iter_x][0] = 99999.9;
                    toBeReturnVec[iter_y][iter_x][1] = 99999.9;
                    toBeReturnVec[iter_y][iter_x][2] = 99999.9;
                    toBeReturnVec[iter_y][iter_x][3] = 99999.9;
                }
            }
        }
    }

    return toBeReturnVec;
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

                    // DEBUG
                    // std::cout << "before comparing min" << std::endl;
                    // std::cout << "tempArr = (" << tempArr[0] << ' ' << tempArr[1] << ' ' << tempArr[2] << ' ' << tempArr[3] << ')' << std::endl;

                    tempArr[0] = std::min(abs(tempArr[0]), abs(i_uVec[j][i][0]));
                    tempArr[1] = std::min(abs(tempArr[1]), abs(i_uVec[j][i][1]));
                    tempArr[2] = std::min(abs(tempArr[2]), abs(i_uVec[j][i][2]));
                    tempArr[3] = std::min(abs(tempArr[3]), abs(i_uVec[j][i][3]));

                    // DEBUG
                    // std::cout << "after comparing min" << std::endl;
                    // std::cout << "tempArr = (" << tempArr[0] << ' ' << tempArr[1] << ' ' << tempArr[2] << ' ' << tempArr[3] << ')' << std::endl;

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

                    tempArr[0] = std::min(abs(tempArr[0]), abs(i_uVec[j][i][0]));
                    tempArr[1] = std::min(abs(tempArr[1]), abs(i_uVec[j][i][1]));
                    tempArr[2] = std::min(abs(tempArr[2]), abs(i_uVec[j][i][2]));
                    tempArr[3] = std::min(abs(tempArr[3]), abs(i_uVec[j][i][3]));

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

                    tempArr[0] = std::min(abs(tempArr[0]), abs(i_uVec[j][i][0]));
                    tempArr[1] = std::min(abs(tempArr[1]), abs(i_uVec[j][i][1]));
                    tempArr[2] = std::min(abs(tempArr[2]), abs(i_uVec[j][i][2]));
                    tempArr[3] = std::min(abs(tempArr[3]), abs(i_uVec[j][i][3]));

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

                    tempArr[0] = std::min(abs(tempArr[0]), abs(i_uVec[j][i][0]));
                    tempArr[1] = std::min(abs(tempArr[1]), abs(i_uVec[j][i][1]));
                    tempArr[2] = std::min(abs(tempArr[2]), abs(i_uVec[j][i][2]));
                    tempArr[3] = std::min(abs(tempArr[3]), abs(i_uVec[j][i][3]));

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