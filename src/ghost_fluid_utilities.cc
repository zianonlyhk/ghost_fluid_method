/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   ghost_fluid_utilities.cc                          Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:53:49 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "ghost_fluid_utilities.hh"
#include "inline/cell_operation.hh"
#include "inline/primitive_tran.hh"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <assert.h>

// public:
// #########################################################################################################################################################################
// 楚河 =====================================================================================================================================================================
// #########################################################################################################################################################################

GhostFluidUtilities::GhostFluidUtilities() {}

std::vector<std::array<int, 2>> GhostFluidUtilities::ghostBoundaryCellCoor(const std::vector<std::vector<double>> &i_levelSet)
{
    std::vector<std::array<int, 2>> toBeReturn;

    int nCell_y = i_levelSet.size() - 4;
    int nCell_x = i_levelSet[0].size() - 4;

    double currPhi;

    std::vector<std::vector<int>> boundaryCellCheck;
    boundaryCellCheck.resize(nCell_y + 4);
    for (int j = 0; j < nCell_y + 4; ++j)
    {
        boundaryCellCheck[j].resize(nCell_x + 4);
        for (int i = 0; i < nCell_x + 4; ++i)
        {
            boundaryCellCheck[j][i] = 0;
        }
    }

    // positive x-sweep for +ve level set boundary
    for (int j = 2; j < nCell_y + 2; ++j)
    {
        currPhi = i_levelSet[j][2];

        for (int i = 2; i < nCell_x + 2; ++i)
        {
            if (currPhi > 0 && i_levelSet[j][i] <= 0)
            {
                if (boundaryCellCheck[j][i] == 0)
                {
                    toBeReturn.push_back(std::array<int, 2>{i, j});
                    boundaryCellCheck[j][i] = 1;
                }
            }
            currPhi = i_levelSet[j][i];
        }
    }

    // negative x-sweep for +ve level set boundary
    for (int j = 2; j < nCell_y + 2; ++j)
    {
        currPhi = i_levelSet[j][nCell_x + 1];

        for (int i = nCell_x + 1; i > 1; --i)
        {
            if (currPhi > 0 && i_levelSet[j][i] < 0)
            {
                if (boundaryCellCheck[j][i] == 0)
                {
                    toBeReturn.push_back(std::array<int, 2>{i, j});
                    boundaryCellCheck[j][i] = 1;
                }
            }
            currPhi = i_levelSet[j][i];
        }
    }

    // positive y-sweep for +ve level set boundary
    for (int i = 2; i < nCell_x + 2; ++i)
    {
        currPhi = i_levelSet[2][i];

        for (int j = 2; j < nCell_y + 2; ++j)
        {
            if (currPhi > 0 && i_levelSet[j][i] < 0)
            {
                if (boundaryCellCheck[j][i] == 0)
                {
                    toBeReturn.push_back(std::array<int, 2>{i, j});
                    boundaryCellCheck[j][i] = 1;
                }
            }
            currPhi = i_levelSet[j][i];
        }
    }

    // negative y-sweep for +ve level set boundary
    for (int i = 2; i < nCell_x + 2; ++i)
    {
        currPhi = i_levelSet[nCell_y + 1][i];

        for (int j = nCell_y + 1; j > 1; --j)
        {
            if (currPhi > 0 && i_levelSet[j][i] < 0)
            {
                if (boundaryCellCheck[j][i] == 0)
                {
                    toBeReturn.push_back(std::array<int, 2>{i, j});
                    boundaryCellCheck[j][i] = 1;
                }
            }
            currPhi = i_levelSet[j][i];
        }
    }

    return toBeReturn;
}

std::array<double, 4> GhostFluidUtilities::ghostCellValues(const std::vector<std::vector<double>> &i_levelSet, const std::vector<std::vector<std::array<double, 4>>> &i_compDomain, std::array<int, 2> i_coor, double i_dx, double i_dy)
{
    double local_gamma = 1.4;

    std::array<double, 2> normalVec_ghostCell = normalUnitVector(i_levelSet, i_coor[0], i_coor[1], i_dx, i_dy);
    std::array<double, 2> normalVec_realFluid = {normalVec_ghostCell[0], normalVec_ghostCell[1]};
    std::array<double, 4> mirrorState = getBilinearlyProbedCell(i_compDomain, i_levelSet, i_coor, i_dx, i_dy);

    // turn into primitive form here
    double mirrorVelX = primitiveX_Vel(mirrorState);
    double mirrorVelY = primitiveY_Vel(mirrorState);
    double mirrorP = primitivePressure(mirrorState);

    double dotProduct = normalVec_realFluid[0] * mirrorVelX + normalVec_realFluid[1] * mirrorVelY;

    std::array<double, 2> normalComponent = {dotProduct * normalVec_realFluid[0], dotProduct * normalVec_realFluid[1]};
    std::array<double, 2> tangentialComponent = {mirrorVelX - normalComponent[0], mirrorVelY - normalComponent[1]};

    std::array<double, 3> riemannLeftState = {mirrorState[0], dotProduct, mirrorP};
    std::array<double, 3> riemannRightState = {mirrorState[0], -dotProduct, mirrorP};

    std::array<double, 3> starredState = HLLC_1D(riemannLeftState, riemannRightState);
    double finalRho = starredState[0];
    double finalMomentumX = starredState[1] * normalVec_realFluid[0] + tangentialComponent[0] * finalRho;
    double finalMomentumY = starredState[1] * normalVec_realFluid[1] + tangentialComponent[1] * finalRho;
    double finalEnergy = starredState[2] + 0.5 * finalRho * (tangentialComponent[0] * tangentialComponent[0] + tangentialComponent[1] * tangentialComponent[1]);

    std::array<double, 4> toBeReturned = {finalRho, finalMomentumX, finalMomentumY, finalEnergy};

    // double finalRho = mirrorState[0];
    // double finalVelX = -normalComponent[0] + tangentialComponent[0];
    // double finalVelY = -normalComponent[1] + tangentialComponent[1];
    // double finalMomentumX = finalVelX * finalRho;
    // double finalMomentumY = finalVelY * finalRho;
    // double finalEnergy = mirrorState[3] / (local_gamma - 1) + 0.5 * finalRho * (finalVelX * finalVelX + finalVelY * finalVelY);

    // DEBUG
    // std::cout << "mirrorState is: (" << mirrorState[0] << ", " << mirrorState[1] << ", " << mirrorState[2] << ", " << mirrorState[3] << ')' << std::endl;
    // std::cout << "toBeReturned is: (" << toBeReturned[0] << ", " << toBeReturned[1] << ", " << toBeReturned[2] << ", " << toBeReturned[3] << ')' << std::endl;
    // if (i_coor[0] == 34 && i_coor[1] == 28)
    // {
    //     std::cout << "mirrorState is: (" << mirrorState[0] << ", " << mirrorState[1] << ", " << mirrorState[2] << ", " << mirrorState[3] << ')' << std::endl;
    //     std::cout << "toBeReturned is: (" << toBeReturned[0] << ", " << toBeReturned[1] << ", " << toBeReturned[2] << ", " << toBeReturned[3] << ')' << std::endl;
    // }

    return toBeReturned;
}

std::array<double, 4> GhostFluidUtilities::ghostCellValues(const std::vector<std::vector<double>> &i_levelSet, const std::vector<std::vector<std::array<double, 4>>> &i_compDomain, std::array<int, 2> i_coor, double i_dx, double i_dy, std::array<double, 2> i_rigidBodyVel)
{
    // double local_gamma = 1.4;

    // std::array<double, 2> normalVec_ghostCell = normalUnitVector(i_levelSet, i_coor[0], i_coor[1], i_dx, i_dy);
    // std::array<double, 2> normalVec_realFluid = {normalVec_ghostCell[0], normalVec_ghostCell[1]};
    // std::array<double, 4> mirrorState = getBilinearlyProbedCell(i_compDomain, i_levelSet, i_coor, i_dx, i_dy);

    // // turn into primitive form here
    // double mirrorVelX = primitiveX_Vel(mirrorState);
    // double mirrorVelY = primitiveY_Vel(mirrorState);
    // double mirrorP = primitivePressure(mirrorState);

    // double dotProduct = normalVec_realFluid[0] * mirrorVelX + normalVec_realFluid[1] * mirrorVelY;

    // std::array<double, 2> normalComponent = {dotProduct * normalVec_realFluid[0], dotProduct * normalVec_realFluid[1]};
    // std::array<double, 2> tangentialComponent = {mirrorVelX - normalComponent[0], mirrorVelY - normalComponent[1]};

    // std::array<double, 3> riemannLeftState = {mirrorState[0], dotProduct, mirrorP};
    // std::array<double, 3> riemannRightState = {mirrorState[0], -dotProduct, mirrorP};

    // std::array<double, 3> starredState = HLLC_1D(riemannLeftState, riemannRightState);
    // double finalRho = starredState[0];
    // double finalMomentumX = starredState[1] * normalVec_realFluid[0] + tangentialComponent[0] * finalRho;
    // double finalMomentumY = starredState[1] * normalVec_realFluid[1] + tangentialComponent[1] * finalRho;
    // double finalEnergy = starredState[2] + 0.5 * finalRho * (tangentialComponent[0] * tangentialComponent[0] + tangentialComponent[1] * tangentialComponent[1]);

    // std::array<double, 4> toBeReturned = {finalRho, finalMomentumX, finalMomentumY, finalEnergy};

    double local_gamma = 1.4;

    std::array<double, 2> normalVec_ghostCell = normalUnitVector(i_levelSet, i_coor[0], i_coor[1], i_dx, i_dy);
    std::array<double, 2> normalVec_realFluid = {normalVec_ghostCell[0], normalVec_ghostCell[1]};
    std::array<double, 4> mirrorState = getBilinearlyProbedCell(i_compDomain, i_levelSet, i_coor, i_dx, i_dy);

    // turn into primitive form here
    double mirrorVelX = primitiveX_Vel(mirrorState);
    double mirrorVelY = primitiveY_Vel(mirrorState);
    double mirrorP = primitivePressure(mirrorState);

    double dotProduct = normalVec_realFluid[0] * mirrorVelX + normalVec_realFluid[1] * mirrorVelY;

    double rigidBodyVelDotProd = normalVec_ghostCell[0] * i_rigidBodyVel[0] + normalVec_ghostCell[1] * i_rigidBodyVel[1];
    std::array<double, 2> normalComponent = {dotProduct * mirrorVelX + rigidBodyVelDotProd * i_rigidBodyVel[0], dotProduct * mirrorVelX + rigidBodyVelDotProd * i_rigidBodyVel[1]};

    // std::array<double, 2> normalComponent = {dotProduct * normalVec_realFluid[0], dotProduct * normalVec_realFluid[1]};
    std::array<double, 2> tangentialComponent = {mirrorVelX + i_rigidBodyVel[0] - normalComponent[0], mirrorVelY + i_rigidBodyVel[1] - normalComponent[1]};

    std::array<double, 3> riemannLeftState = {mirrorState[0], dotProduct, mirrorP};
    std::array<double, 3> riemannRightState = {mirrorState[0], -dotProduct - rigidBodyVelDotProd, mirrorP};

    std::array<double, 3> starredState = HLLC_1D(riemannLeftState, riemannRightState);
    double finalRho = starredState[0];
    double finalMomentumX = starredState[1] * normalVec_realFluid[0] - tangentialComponent[0] * normalVec_realFluid[0];
    double finalMomentumY = starredState[1] * normalVec_realFluid[1] - tangentialComponent[1] * normalVec_realFluid[0];
    double finalEnergy = starredState[2] + 0.5 * finalRho * (tangentialComponent[0] * tangentialComponent[0] + tangentialComponent[1] * tangentialComponent[1]);

    std::array<double, 4> toBeReturned = {finalRho, finalMomentumX, finalMomentumY, finalEnergy};

    return toBeReturned;
}

std::array<double, 4> GhostFluidUtilities::solveForConstantExtrapolation(const std::vector<std::vector<double>> &i_levelSet, const std::vector<std::vector<std::array<double, 4>>> &i_compDomain, std::array<int, 2> i_coor, double i_dx, double i_dy)
{
    std::array<double, 2> normalVec = normalUnitVector(i_levelSet, i_coor[0], i_coor[1], i_dx, i_dy);

    std::array<double, 4> referenceCell_x;
    std::array<double, 4> referenceCell_y;

    if (normalVec[0] > 0)
    {
        referenceCell_x = i_compDomain[i_coor[1]][i_coor[0] + 1];
    }
    else
    {
        referenceCell_x = i_compDomain[i_coor[1]][i_coor[0] - 1];
    }
    if (normalVec[1] > 0)
    {
        referenceCell_y = i_compDomain[i_coor[1] + 1][i_coor[0]];
    }
    else
    {
        referenceCell_y = i_compDomain[i_coor[1] - 1][i_coor[0]];
    }

    std::array<double, 4> toBeReturned;
    // double scalingConstant = normalVec[0] / i_dx + normalVec[1] / i_dy;
    // toBeReturned = scalingCell(1 / scalingConstant, sumCell(scalingCell(normalVec[0] / i_dx, referenceCell_x), scalingCell(normalVec[1] / i_dy, referenceCell_y)));
    toBeReturned = sumCell(scalingCell(normalVec[0] * normalVec[0], referenceCell_x), scalingCell(normalVec[1] * normalVec[1], referenceCell_y));

    return toBeReturned;
}

// private:
// #########################################################################################################################################################################
// 漢界 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
// #########################################################################################################################################################################

std::array<double, 2> GhostFluidUtilities::normalUnitVector(const std::vector<std::vector<double>> &i_levelSet, int i_x, int i_y, double i_dx, double i_dy)
{
    double dphidx = (i_levelSet[i_y][i_x + 1] - i_levelSet[i_y][i_x - 1]) / (2 * i_dx);
    double dphidy = (i_levelSet[i_y + 1][i_x] - i_levelSet[i_y - 1][i_x]) / (2 * i_dy);

    double scalingFactor = sqrt(dphidx * dphidx + dphidy * dphidy);

    std::array<double, 2> toBeReturn;
    toBeReturn = std::array<double, 2>{dphidx / scalingFactor, dphidy / scalingFactor};

    return toBeReturn;
}

std::array<double, 2> GhostFluidUtilities::probeCoor(const std::vector<std::vector<double>> &i_levelSet, std::array<int, 2> i_currCoor, double i_dx, double i_dy)
{
    double interfaceX;
    double interfaceY;
    double currPhi = abs(i_levelSet[i_currCoor[1]][i_currCoor[0]]);
    std::array<double, 2> i_normalVector = normalUnitVector(i_levelSet, i_currCoor[0], i_currCoor[1], i_dx, i_dy);

    interfaceX = i_currCoor[0] + currPhi * i_normalVector[0];
    interfaceY = i_currCoor[1] + currPhi * i_normalVector[1];

    assert(interfaceX + 1.5 * i_normalVector[0] > 0 && interfaceY + 1.5 * i_normalVector[1] > 0);

    // DEBUG
    // std::cout << "at the coor (" << i_currCoor[0] << ", " << i_currCoor[1] << ')' << std::endl;
    // std::cout << "it is pointed to (" << interfaceX + 1.5 * i_normalVector[0] << ", " << interfaceY + 1.5 * i_normalVector[1] << ')' << std::endl;

    return std::array<double, 2>{interfaceX + 1.5 * i_normalVector[0], interfaceY + 1.5 * i_normalVector[1]};
}

std::array<double, 4> GhostFluidUtilities::getBilinearlyProbedCell(const std::vector<std::vector<std::array<double, 4>>> &i_compDomain, const std::vector<std::vector<double>> &i_levelSet, std::array<int, 2> i_coor, double i_dx, double i_dy)
{
    std::array<double, 2> exactProbeCoor = probeCoor(i_levelSet, i_coor, i_dx, i_dy);

    int lowerBound_x = std::floor(exactProbeCoor[0]);
    int upperBound_x = lowerBound_x + 1;
    int lowerBound_y = std::floor(exactProbeCoor[1]);
    int upperBound_y = lowerBound_y + 1;

    std::array<int, 2> nokia_1 = {lowerBound_x, upperBound_y};
    std::array<int, 2> nokia_3 = {upperBound_x, upperBound_y};
    std::array<int, 2> nokia_7 = {lowerBound_x, lowerBound_y};
    std::array<int, 2> nokia_9 = {upperBound_x, upperBound_y};

    std::array<double, 4> nokia_1_cell = i_compDomain[nokia_1[1]][nokia_1[0]];
    std::array<double, 4> nokia_3_cell = i_compDomain[nokia_3[1]][nokia_3[0]];
    std::array<double, 4> nokia_7_cell = i_compDomain[nokia_7[1]][nokia_7[0]];
    std::array<double, 4> nokia_9_cell = i_compDomain[nokia_9[1]][nokia_9[0]];

    double del_x = exactProbeCoor[0] - lowerBound_x;
    double del_y = exactProbeCoor[1] - lowerBound_y;

    std::array<double, 4> upper = sumCell(scalingCell(1 - del_x, nokia_1_cell), scalingCell(del_x, nokia_3_cell));
    std::array<double, 4> lower = sumCell(scalingCell(1 - del_x, nokia_7_cell), scalingCell(del_x, nokia_9_cell));

    return sumCell(scalingCell(1 - del_y, lower), scalingCell(del_y, upper));
}

std::array<double, 3> GhostFluidUtilities::HLLC_1D(std::array<double, 3> i_W_L, std::array<double, 3> i_W_R)
{
    double local_gamma = 1.4;

    // transforming into primitive
    double rho_L = i_W_L[0];
    double rho_R = i_W_R[0];
    double u_L = i_W_L[1];
    double u_R = i_W_R[1];
    double p_L = i_W_L[2];
    double p_R = i_W_R[2];

    // Toro Step 1
    double rho_bar = 0.5 * (rho_L + rho_R);
    double a_L = sqrt(local_gamma * p_L / rho_L);
    double a_R = sqrt(local_gamma * p_R / rho_R);
    double a_bar = 0.5 * (a_L + a_R);
    double p_pvrs = 0.5 * (p_L + p_R) - 0.5 * (u_R - u_L) * rho_bar * a_bar;
    double p_star = std::max(0.0, p_pvrs);

    // Toro Step 2
    double q_L;
    double q_R;

    if (p_star > p_L)
    {
        q_L = sqrt(1 + (local_gamma + 1) / (2 * local_gamma) * (p_star / p_L - 1));
    }
    else
    {
        q_L = 1;
    }

    if (p_star > p_R)
    {
        q_R = sqrt(1 + (local_gamma + 1) / (2 * local_gamma) * (p_star / p_R - 1));
    }
    else
    {
        q_R = 1;
    }

    double S_L = std::min(u_L - a_L * q_L, u_R - a_R * q_R);
    double S_R = std::min(u_L + a_L * q_L, u_R + a_R * q_R);

    double S_star = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) / (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));

    // Toro Step 3
    double toroScaler_L = rho_L * (S_L - u_L) / (S_L - S_star);
    std::array<double, 3> U_star_L;

    double energy_L = p_L / (local_gamma - 1) + 0.5 * rho_L * (u_L * u_L);

    U_star_L[0] = 1;
    U_star_L[1] = S_star;
    U_star_L[2] = energy_L / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L)));

    U_star_L = std::array<double, 3>{toroScaler_L * U_star_L[0], toroScaler_L * U_star_L[1], toroScaler_L * U_star_L[2]};

    return U_star_L;
}
