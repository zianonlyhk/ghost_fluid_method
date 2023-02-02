/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   flux_func.hh                                      Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:54:01 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef FLUX_FUNC_HH
#define FLUX_FUNC_HH

#include <array>

class FluxFunc
{
public:
    FluxFunc();

    // FVM conservation flux
    std::array<double, 4> conservationFlux_x(std::array<double, 4>);
    std::array<double, 4> conservationFlux_y(std::array<double, 4>);

    // FORCE flux
    std::array<double, 4> forceFlux_x(std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<double, 4> forceFlux_y(std::array<double, 4>, std::array<double, 4>, double, double);

    // SLIC flux
    std::array<double, 4> slicFlux_x(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<double, 4> slicFlux_y(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);

    // MUSCL-Hancock HLLC flux
    std::array<double, 4> musclHancockHllcFlux_x(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<double, 4> musclHancockHllcFlux_y(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);

private:
    // Slope limiting for SLIC and MUSCL-Hancock
    std::array<std::array<double, 4>, 2> slopeLimitedLR_U_x(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<std::array<double, 4>, 2> slopeLimitedLR_U_y(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);

    // Approx Riemann Solvers
    std::array<double, 4> HLL_Riemannflux_x(std::array<double, 4>, std::array<double, 4>);
    std::array<double, 4> HLL_Riemannflux_y(std::array<double, 4>, std::array<double, 4>);
    std::array<double, 4> HLLC_Riemannflux_x(std::array<double, 4>, std::array<double, 4>);
    std::array<double, 4> HLLC_Riemannflux_y(std::array<double, 4>, std::array<double, 4>);
};

#endif