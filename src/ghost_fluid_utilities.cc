/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   ghost_fluid_utilities.cc                          Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/24 12:32:28 by Zian Huang                               */
/*   Updated: 2023/01/27 11:00:41 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "ghost_fluid_utilities.hh"
#include <math.h>
#include <algorithm>
#include <iostream>

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
            if (currPhi < 0 && i_levelSet[j][i] > 0)
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
            if (currPhi < 0 && i_levelSet[j][i] > 0)
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

        for (int j = 2; j < nCell_y + 2; ++i)
        {
            if (currPhi < 0 && i_levelSet[j][i] > 0)
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

        for (int j = nCell_x + 1; j > 1; --i)
        {
            if (currPhi < 0 && i_levelSet[j][i] > 0)
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

// private:
// #########################################################################################################################################################################
// 漢界 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
// #########################################################################################################################################################################

std::array<double, 2> GhostFluidUtilities::normalUnitVector(const std::vector<std::vector<double>> &i_levelSet, int i_x, int i_y, double i_dx, double i_dy)
{
    double dphidx = (i_levelSet[i_y][i_x + 1] - i_levelSet[i_y][i_x - 1]) / (2 * i_dx);
    double dphidy = (i_levelSet[i_y + 1][i_x] - i_levelSet[i_y - 1][i_x]) / (2 * i_dy);

    std::array<double, 2> toBeReturn;
    toBeReturn = std::array<double, 2>{dphidx, dphidy};

    return toBeReturn;
}