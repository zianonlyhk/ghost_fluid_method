/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   ghost_fluid_utilities.hh                          Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/24 12:32:25 by Zian Huang                               */
/*   Updated: 2023/01/30 15:14:46 by Zian Huang                               */
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

    std::vector<std::array<int, 2>> ghostBoundaryCellCoor(const std::vector<std::vector<double>> &levelSet);

    std::array<double, 4> ghostCellValues(const std::vector<std::vector<double>> &compDomain, const std::vector<std::vector<std::array<double, 4>>> &levelSet, std::array<int, 2> coor);

private:
    std::array<double, 2> normalUnitVector(const std::vector<std::vector<double>> &levelSet, int x, int y, double dx, double dy);

    std::array<double, 4> binearInterpolation(std::array<int, 4> p1, std::array<int, 4> p2, std::array<int, 4> p3, std::array<int, 4> p4);

    std::array<double, 4> ghostCellReflectiveBoundary();
    std::array<double, 2> tangentialVel();
    std::array<std::array<double, 3>, 2> normalRiemannStates();
};

#endif