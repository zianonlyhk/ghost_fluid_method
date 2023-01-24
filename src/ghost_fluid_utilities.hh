/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   ghost_fluid_utilities.hh                          Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/24 12:32:25 by Zian Huang                               */
/*   Updated: 2023/01/24 13:14:07 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef GHOST_FLUID_UTILITIES_HH
#define GHOST_FLUID_UTILITIES_HH

#include <array>
#include <vector>

class GhostFluidUtilities
{
public:
    GhostFluidUtilities();

    std::vector<std::array<int, 2>> ghostBoundaryCellCoor(std::vector<std::vector<double>> levelSet);

private:
};

#endif