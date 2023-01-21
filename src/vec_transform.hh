/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.hh                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/21 10:45:20 by Zian Huang                               */
/*   Updated: 2023/01/21 18:00:34 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef VEC_TRANSFORM_HH
#define VEC_TRANSFORM_HH

#include <vector>
#include <array>
#include "flux_func.hh"

class VecTran
{
public:
    VecTran();
    FluxFVM fluxFunc;

    std::vector<std::vector<std::array<double, 4>>> slicVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dx, double dt);
    std::vector<std::vector<std::array<double, 4>>> slicVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dy, double dt);

    std::vector<std::vector<std::array<double, 4>>> musclHancockVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dx, double dt);
    std::vector<std::vector<std::array<double, 4>>> musclHancockVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dy, double dt);
};

#endif