/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   debug_tools.hh                                    Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/24 10:58:05 by Zian Huang                               */
/*   Updated: 2023/02/01 15:14:01 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef DEBUG_TOOLS_HH
#define DEBUG_TOOLS_HH

#include <iostream>
#include <vector>
#include <array>

inline void printDomainDensity(std::vector<std::vector<std::array<double, 4>>> i_compDomain)
{

    std::cout.precision(2);
    std::cout << std::fixed;

    int xVecLen = i_compDomain[0].size();
    int yVecLen = i_compDomain.size();

    for (int iter_y = 0; iter_y < yVecLen; ++iter_y)
    {
        for (int iter_x = 0; iter_x < xVecLen; ++iter_x)
        {
            std::cout << i_compDomain[iter_y][iter_x][0] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

inline void printDomainMomentumY(std::vector<std::vector<std::array<double, 4>>> i_compDomain)
{

    std::cout.precision(2);
    std::cout << std::fixed;

    int xVecLen = i_compDomain[0].size();
    int yVecLen = i_compDomain.size();

    for (int iter_y = 0; iter_y < yVecLen; ++iter_y)
    {
        for (int iter_x = 0; iter_x < xVecLen; ++iter_x)
        {

            if (i_compDomain[iter_y][iter_x][2] < 0)
            {
                std::cout << i_compDomain[iter_y][iter_x][2] << ' ';
            }
            else
            {
                std::cout << ' ' << i_compDomain[iter_y][iter_x][2] << ' ';
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

inline void printLevelSet(const std::vector<std::vector<double>> &i_levelSet)
{

    std::cout.precision(2);
    std::cout << std::fixed;

    int xVecLen = i_levelSet[0].size();
    int yVecLen = i_levelSet.size();

    for (int iter_y = 0; iter_y < yVecLen; ++iter_y)
    {
        for (int iter_x = 0; iter_x < xVecLen; ++iter_x)
        {
            if (i_levelSet[iter_y][iter_x] < 0)
            {
                std::cout << i_levelSet[iter_y][iter_x] << ' ';
            }
            else
            {
                std::cout << ' ' << i_levelSet[iter_y][iter_x] << ' ';
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

inline void printBoundaryCellCoor(const std::vector<std::array<int, 2>> &coorList)
{
    int listLen = coorList.size();

    for (int i = 0; i < listLen; ++i)
    {

        std::cout << coorList[i][0] << ' ' << coorList[i][1];
        std::cout << std::endl;
    }
}

#endif
