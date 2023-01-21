/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   flux_func.hh                                      Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/21 16:12:07 by Zian Huang                               */
/*   Updated: 2023/01/21 16:18:02 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef FLUX_FUNC_HH
#define FLUX_FUNC_HH

#include <array>

class FluxFVM
{
public:
    FluxFVM();

    // instantaneous flux based on the state of cell
    std::array<double, 4> conservationFlux_x(std::array<double, 4>);
    std::array<double, 4> conservationFlux_y(std::array<double, 4>);

    // centred scheme flux that requires the 2 neighbours and the space time distances between them
    std::array<double, 4> forceFlux_x(std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<double, 4> forceFlux_y(std::array<double, 4>, std::array<double, 4>, double, double);

    // SLIC scheme flux requires the 5 neighbours and the space time distances between them
    std::array<double, 4> slicFlux_x(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<double, 4> slicFlux_y(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);
};

#endif