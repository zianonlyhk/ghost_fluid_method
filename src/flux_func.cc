/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   flux_func.cc                                      Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/21 16:12:23 by Zian Huang                               */
/*   Updated: 2023/01/21 16:20:28 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "flux_func.hh"
#include "inline/primitive_tran.hh"
#include "inline/cell_operation.hh"

FluxFVM::FluxFVM() {}

std::array<double, 4> FluxFVM::conservationFlux_x(std::array<double, 4> i_inputVec)
{
    std::array<double, 4> arrayToBeReturned;

    arrayToBeReturned[0] = i_inputVec[1];
    arrayToBeReturned[1] = i_inputVec[1] * primitiveX_Vel(i_inputVec) + primitivePressure(i_inputVec);
    arrayToBeReturned[2] = i_inputVec[1] * primitiveY_Vel(i_inputVec);
    arrayToBeReturned[3] = primitiveX_Vel(i_inputVec) * (i_inputVec[3] + primitivePressure(i_inputVec));

    return arrayToBeReturned;
}

std::array<double, 4> FluxFVM::conservationFlux_y(std::array<double, 4> i_inputVec)
{
    std::array<double, 4> arrayToBeReturned;

    arrayToBeReturned[0] = i_inputVec[2];
    arrayToBeReturned[1] = i_inputVec[2] * primitiveX_Vel(i_inputVec);
    arrayToBeReturned[2] = i_inputVec[2] * primitiveY_Vel(i_inputVec) + primitivePressure(i_inputVec);
    arrayToBeReturned[3] = primitiveY_Vel(i_inputVec) * (i_inputVec[3] + primitivePressure(i_inputVec));

    return arrayToBeReturned;
}

std::array<double, 4> FluxFVM::forceFlux_x(std::array<double, 4> i_uVector_i, std::array<double, 4> i_uVector_i_next, double i_dx, double i_dt)
{
    std::array<double, 4> lfFluxArr;
    std::array<double, 4> riTempArr;

    lfFluxArr = sumCell(scalingCell(0.5 * i_dx / i_dt, diffCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5, sumCell(conservationFlux_x(i_uVector_i_next), conservationFlux_x(i_uVector_i))));
    riTempArr = diffCell(scalingCell(0.5, sumCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5 * i_dt / i_dx, diffCell(conservationFlux_x(i_uVector_i_next), conservationFlux_x(i_uVector_i))));

    return scalingCell(0.5, sumCell(lfFluxArr, conservationFlux_x(riTempArr)));
}

std::array<double, 4> FluxFVM::forceFlux_y(std::array<double, 4> i_uVector_i, std::array<double, 4> i_uVector_i_next, double i_dy, double i_dt)
{
    std::array<double, 4> lfFluxArr;
    std::array<double, 4> riTempArr;

    lfFluxArr = sumCell(scalingCell(0.5 * i_dy / i_dt, diffCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5, sumCell(conservationFlux_y(i_uVector_i_next), conservationFlux_y(i_uVector_i))));
    riTempArr = diffCell(scalingCell(0.5, sumCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5 * i_dt / i_dy, diffCell(conservationFlux_y(i_uVector_i_next), conservationFlux_y(i_uVector_i))));

    return scalingCell(0.5, sumCell(lfFluxArr, conservationFlux_y(riTempArr)));
}

std::array<double, 4> FluxFVM::slicFlux_x(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dx, double i_dt)
{
    const double local_omega = 0.0;
    std::array<double, 4> uVec_1_L;
    std::array<double, 4> uVec_2_L;
    std::array<double, 4> uVec_3_L;
    std::array<double, 4> uVec_1_R;
    std::array<double, 4> uVec_2_R;
    std::array<double, 4> uVec_3_R;

    // interpret left and right flux
    uVec_1_L = diffCell(i_uVector_1, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_1, i_uVector_0)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_2, i_uVector_1)))));
    uVec_1_R = sumCell(i_uVector_1, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_1, i_uVector_0)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_2, i_uVector_1)))));
    uVec_2_L = diffCell(i_uVector_2, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_2, i_uVector_1)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_3, i_uVector_2)))));
    uVec_2_R = sumCell(i_uVector_2, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_2, i_uVector_1)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_3, i_uVector_2)))));

    // locally half time evolve
    uVec_1_R = diffCell(uVec_1_R, scalingCell(0.5 * i_dt / i_dx, diffCell(conservationFlux_x(uVec_1_R), conservationFlux_x(uVec_1_L))));
    uVec_2_L = diffCell(uVec_2_L, scalingCell(0.5 * i_dt / i_dx, diffCell(conservationFlux_x(uVec_2_R), conservationFlux_x(uVec_2_L))));

    return forceFlux_x(uVec_1_R, uVec_2_L, i_dx, i_dt);
}

std::array<double, 4> FluxFVM::slicFlux_y(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dy, double i_dt)
{
    // const parameter
    const double local_omega = 0.0;
    std::array<double, 4> uVec_1_L;
    std::array<double, 4> uVec_2_L;
    std::array<double, 4> uVec_3_L;
    std::array<double, 4> uVec_1_R;
    std::array<double, 4> uVec_2_R;
    std::array<double, 4> uVec_3_R;

    // interpret left and right flux
    uVec_1_L = diffCell(i_uVector_1, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_1, i_uVector_0)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_2, i_uVector_1)))));
    uVec_1_R = sumCell(i_uVector_1, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_1, i_uVector_0)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_2, i_uVector_1)))));
    uVec_2_L = diffCell(i_uVector_2, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_2, i_uVector_1)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_3, i_uVector_2)))));
    uVec_2_R = sumCell(i_uVector_2, scalingCell(0.5, sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_2, i_uVector_1)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_3, i_uVector_2)))));

    // locally half time evolve
    uVec_1_R = diffCell(uVec_1_R, scalingCell(0.5 * i_dt / i_dy, diffCell(conservationFlux_y(uVec_1_R), conservationFlux_y(uVec_1_L))));
    uVec_2_L = diffCell(uVec_2_L, scalingCell(0.5 * i_dt / i_dy, diffCell(conservationFlux_y(uVec_2_R), conservationFlux_y(uVec_2_L))));

    return forceFlux_y(uVec_1_R, uVec_2_L, i_dy, i_dt);
}
