/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   debug_tools.hh                                    Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/24 10:58:05 by Zian Huang                               */
/*   Updated: 2023/01/24 11:10:09 by Zian Huang                               */
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

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            std::cout << i_compDomain[iter_y][iter_x][0] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

#endif