/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   ghost_fluid_utilities.hh                          Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/24 12:32:25 by Zian Huang                               */
/*   Updated: 2023/01/30 22:30:11 by Zian Huang                               */
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

    std::array<double, 4> ghostCellValues(const std::vector<std::vector<double>> &levelSet, const std::vector<std::vector<std::array<double, 4>>> &compDomain, std::array<int, 2> coor);
    std::array<double, 4> ghostCellValues(const std::vector<std::vector<double>> &levelSet, const std::vector<std::vector<std::array<double, 4>>> &compDomain, std::array<int, 2> coor, std::array<double, 2> rigidBodyVel);

private:
    std::array<double, 2> normalUnitVector(const std::vector<std::vector<double>> &levelSet, int x, int y);
    std::array<double, 2> probeCoor(const std::vector<std::vector<double>> &levelSet, std::array<int, 2> currCoor);
    std::array<double, 4> getBilinearlyProbedCell(const std::vector<std::vector<std::array<double, 4>>> &compDomain, const std::vector<std::vector<double>> &levelSet, std::array<int, 2> coor);

    std::array<double, 3> HLLC_1D(std::array<double, 3> W_L, std::array<double, 3> W_R);
};

#endif