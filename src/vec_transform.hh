/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.hh                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/21 10:45:20 by Zian Huang                               */
/*   Updated: 2023/01/24 13:10:44 by Zian Huang                               */
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

private:
    FluxFunc fluxFunc;
    GhostFluidUtilities ghostFluidUtilities;
};

#endif