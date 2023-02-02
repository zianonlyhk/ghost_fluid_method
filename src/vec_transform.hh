/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.hh                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/21 10:45:20 by Zian Huang                               */
/*   Updated: 2023/02/01 14:15:04 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef VEC_TRANSFORM_HH
#define VEC_TRANSFORM_HH

#include "flux_func.hh"
#include "ghost_fluid_utilities.hh"
#include <vector>
#include <array>

class VecTran
{
public:
    VecTran();

    // FVM conservation time update
    std::vector<std::vector<std::array<double, 4>>> slicVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dx, double dt);
    std::vector<std::vector<std::array<double, 4>>> slicVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dy, double dt);
    std::vector<std::vector<std::array<double, 4>>> musclHancockVecTranHLLC_x(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dx, double dt);
    std::vector<std::vector<std::array<double, 4>>> musclHancockVecTranHLLC_y(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dy, double dt);

    // GFM transformation functions
    std::vector<std::vector<std::array<double, 4>>> ghostCellBoundary(const std::vector<std::vector<std::array<double, 4>>> &uVec, const std::vector<std::vector<double>> &levelSet, double dx, double dy);
    std::vector<std::vector<std::array<double, 4>>> propagateGhostInterface(const std::vector<std::vector<std::array<double, 4>>> &uVec, const std::vector<std::vector<double>> &levelSet, double dx, double dy);

    std::vector<std::array<int, 2>> getBoundaryCellCoor(const std::vector<std::vector<double>> &levelSet);

private:
    FluxFunc fluxFunc;
    GhostFluidUtilities ghostFluidUtilities;

    void fastSweepingConstantPropagation(std::vector<std::vector<std::array<double, 4>>> &toBeReturned, const std::vector<std::vector<std::array<double, 4>>> &uVec, const std::vector<std::vector<double>> &levelSet, double dx, double dy);
};

#endif