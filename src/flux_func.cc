/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   flux_func.cc                                      Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:54:06 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "flux_func.hh"
#include "inline/primitive_tran.hh"
#include "inline/cell_operation.hh"
#include <math.h>
#include <algorithm>
#include <iostream>

std::array<double, 4> vanLeerLimiter(std::array<double, 4> i_arr0, std::array<double, 4> i_arr1, std::array<double, 4> i_arr2)
{
    std::array<double, 4> delta_formerHalf = diffCell(i_arr1, i_arr0);
    std::array<double, 4> delta_laterHalf = diffCell(i_arr2, i_arr1);

    std::array<double, 4> r = divisionCell(delta_formerHalf, delta_laterHalf);

    std::array<double, 4> e_L = divisionCell(scalingCell(2, r), scalarAdditionCell(1, r));
    std::array<double, 4> e_R = scalarDivisionCell(2, scalarAdditionCell(1, r));

    std::array<double, 4> cellToBeReturned;
    for (int i = 0; i < 4; ++i)
    {
        if (r[i] < 0.0 || r[i] != r[i])
        {
            cellToBeReturned[i] = 0;
        }
        else
        {
            cellToBeReturned[i] = std::min(e_L[i], e_R[i]);
        }
    }

    return cellToBeReturned;
};

std::array<double, 4> minbeeLimiter(std::array<double, 4> i_arr0, std::array<double, 4> i_arr1, std::array<double, 4> i_arr2)
{
    std::array<double, 4> delta_formerHalf = diffCell(i_arr1, i_arr0);
    std::array<double, 4> delta_laterHalf = diffCell(i_arr2, i_arr1);

    std::array<double, 4> r = divisionCell(delta_formerHalf, delta_laterHalf);

    std::array<double, 4> e_L = divisionCell(scalingCell(2, r), scalarAdditionCell(1, r));
    std::array<double, 4> e_R = scalarDivisionCell(2, scalarAdditionCell(1, r));

    std::array<double, 4> cellToBeReturned;
    for (int i = 0; i < 4; ++i)
    {
        if (r[i] < 0.0 || r[i] != r[i])
        {
            cellToBeReturned[i] = 0;
        }
        else if (r[i] > 0 && r[i] <= 1)
        {
            cellToBeReturned[i] = r[i];
        }
        else
        {
            cellToBeReturned[i] = std::min(1.0, e_R[i]);
        }
    }

    return cellToBeReturned;
};

std::array<double, 4> superbeeLimiter(std::array<double, 4> i_arr0, std::array<double, 4> i_arr1, std::array<double, 4> i_arr2)
{
    std::array<double, 4> delta_formerHalf = diffCell(i_arr1, i_arr0);
    std::array<double, 4> delta_laterHalf = diffCell(i_arr2, i_arr1);

    std::array<double, 4> r = divisionCell(delta_formerHalf, delta_laterHalf);

    std::array<double, 4> e_L = divisionCell(scalingCell(2, r), scalarAdditionCell(1, r));
    std::array<double, 4> e_R = scalarDivisionCell(2, scalarAdditionCell(1, r));

    std::array<double, 4> cellToBeReturned;
    for (int i = 0; i < 4; ++i)
    {
        if (r[i] < 0.0 || r[i] != r[i])
        {
            cellToBeReturned[i] = 0;
        }
        else if (r[i] > 0 && r[i] <= 0.5)
        {
            cellToBeReturned[i] = 2 * r[i];
        }
        else if (r[i] > 0.5 && r[i] <= 1)
        {
            cellToBeReturned[i] = 1;
        }
        else
        {
            cellToBeReturned[i] = std::min(2.0, std::min(r[i], e_R[i]));
        }
    }

    return cellToBeReturned;
}

// public:
// #########################################################################################################################################################################
// #########################################################################################################################################################################

FluxFunc::FluxFunc() {}

std::array<double, 4> FluxFunc::conservationFlux_x(std::array<double, 4> i_inputVec)
{
    std::array<double, 4> arrayToBeReturned;

    arrayToBeReturned[0] = i_inputVec[1];
    arrayToBeReturned[1] = i_inputVec[1] * primitiveX_Vel(i_inputVec) + primitivePressure(i_inputVec);
    arrayToBeReturned[2] = i_inputVec[1] * primitiveY_Vel(i_inputVec);
    arrayToBeReturned[3] = primitiveX_Vel(i_inputVec) * (i_inputVec[3] + primitivePressure(i_inputVec));

    return arrayToBeReturned;
}

std::array<double, 4> FluxFunc::conservationFlux_y(std::array<double, 4> i_inputVec)
{
    std::array<double, 4> arrayToBeReturned;

    arrayToBeReturned[0] = i_inputVec[2];
    arrayToBeReturned[1] = i_inputVec[2] * primitiveX_Vel(i_inputVec);
    arrayToBeReturned[2] = i_inputVec[2] * primitiveY_Vel(i_inputVec) + primitivePressure(i_inputVec);
    arrayToBeReturned[3] = primitiveY_Vel(i_inputVec) * (i_inputVec[3] + primitivePressure(i_inputVec));

    return arrayToBeReturned;
}

std::array<double, 4> FluxFunc::forceFlux_x(std::array<double, 4> i_uVector_i, std::array<double, 4> i_uVector_i_next, double i_dx, double i_dt)
{
    std::array<double, 4> lfFluxArr;
    std::array<double, 4> riTempArr;

    lfFluxArr = sumCell(scalingCell(0.5 * i_dx / i_dt, diffCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5, sumCell(conservationFlux_x(i_uVector_i_next), conservationFlux_x(i_uVector_i))));
    riTempArr = diffCell(scalingCell(0.5, sumCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5 * i_dt / i_dx, diffCell(conservationFlux_x(i_uVector_i_next), conservationFlux_x(i_uVector_i))));

    return scalingCell(0.5, sumCell(lfFluxArr, conservationFlux_x(riTempArr)));
}

std::array<double, 4> FluxFunc::forceFlux_y(std::array<double, 4> i_uVector_i, std::array<double, 4> i_uVector_i_next, double i_dy, double i_dt)
{
    std::array<double, 4> lfFluxArr;
    std::array<double, 4> riTempArr;

    lfFluxArr = sumCell(scalingCell(0.5 * i_dy / i_dt, diffCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5, sumCell(conservationFlux_y(i_uVector_i_next), conservationFlux_y(i_uVector_i))));
    riTempArr = diffCell(scalingCell(0.5, sumCell(i_uVector_i, i_uVector_i_next)), scalingCell(0.5 * i_dt / i_dy, diffCell(conservationFlux_y(i_uVector_i_next), conservationFlux_y(i_uVector_i))));

    return scalingCell(0.5, sumCell(lfFluxArr, conservationFlux_y(riTempArr)));
}

std::array<double, 4> FluxFunc::slicFlux_x(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dx, double i_dt)
{
    std::array<std::array<double, 4>, 2> leftRightHalftimeLimitedState = slopeLimitedLR_U_x(i_uVector_0, i_uVector_1, i_uVector_2, i_uVector_3, i_dx, i_dt);

    return forceFlux_x(leftRightHalftimeLimitedState[0], leftRightHalftimeLimitedState[1], i_dx, i_dt);
}

std::array<double, 4> FluxFunc::slicFlux_y(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dy, double i_dt)
{
    std::array<std::array<double, 4>, 2> leftRightHalftimeLimitedState = slopeLimitedLR_U_y(i_uVector_0, i_uVector_1, i_uVector_2, i_uVector_3, i_dy, i_dt);

    return forceFlux_y(leftRightHalftimeLimitedState[0], leftRightHalftimeLimitedState[1], i_dy, i_dt);
}

std::array<double, 4> FluxFunc::musclHancockHllcFlux_x(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dx, double i_dt)
{
    std::array<std::array<double, 4>, 2> leftRightHalftimeLimitedState = slopeLimitedLR_U_x(i_uVector_0, i_uVector_1, i_uVector_2, i_uVector_3, i_dx, i_dt);

    return HLLC_Riemannflux_x(leftRightHalftimeLimitedState[0], leftRightHalftimeLimitedState[1]);
}

std::array<double, 4> FluxFunc::musclHancockHllcFlux_y(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dy, double i_dt)
{
    std::array<std::array<double, 4>, 2> leftRightHalftimeLimitedState = slopeLimitedLR_U_y(i_uVector_0, i_uVector_1, i_uVector_2, i_uVector_3, i_dy, i_dt);

    return HLLC_Riemannflux_y(leftRightHalftimeLimitedState[0], leftRightHalftimeLimitedState[1]);
}

// private:
// #########################################################################################################################################################################
// #########################################################################################################################################################################

std::array<std::array<double, 4>, 2> FluxFunc::slopeLimitedLR_U_x(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dx, double i_dt)
{
    // const parameter
    const double local_omega = 0.0;
    std::array<double, 4> uVec_1_L;
    std::array<double, 4> uVec_2_L;
    std::array<double, 4> uVec_3_L;
    std::array<double, 4> uVec_1_R;
    std::array<double, 4> uVec_2_R;
    std::array<double, 4> uVec_3_R;

    std::array<double, 4> Delta_i_1 = sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_1, i_uVector_0)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_2, i_uVector_1)));
    std::array<double, 4> Delta_i_2 = sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_2, i_uVector_1)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_3, i_uVector_2)));

    std::array<double, 4> limiter_1 = superbeeLimiter(i_uVector_0, i_uVector_1, i_uVector_2);
    std::array<double, 4> limiter_2 = superbeeLimiter(i_uVector_1, i_uVector_2, i_uVector_3);

    // interpret left and right flux
    uVec_1_L = diffCell(i_uVector_1, scalingCell(0.5, productCell(limiter_1, Delta_i_1)));
    uVec_1_R = sumCell(i_uVector_1, scalingCell(0.5, productCell(limiter_1, Delta_i_1)));
    uVec_2_L = diffCell(i_uVector_2, scalingCell(0.5, productCell(limiter_2, Delta_i_2)));
    uVec_2_R = sumCell(i_uVector_2, scalingCell(0.5, productCell(limiter_2, Delta_i_2)));

    // locally half time evolve
    uVec_1_R = diffCell(uVec_1_R, scalingCell(0.5 * i_dt / i_dx, diffCell(conservationFlux_x(uVec_1_R), conservationFlux_x(uVec_1_L))));
    uVec_2_L = diffCell(uVec_2_L, scalingCell(0.5 * i_dt / i_dx, diffCell(conservationFlux_x(uVec_2_R), conservationFlux_x(uVec_2_L))));

    return std::array<std::array<double, 4>, 2>{uVec_1_R, uVec_2_L};
}

std::array<std::array<double, 4>, 2> FluxFunc::slopeLimitedLR_U_y(std::array<double, 4> i_uVector_0, std::array<double, 4> i_uVector_1, std::array<double, 4> i_uVector_2, std::array<double, 4> i_uVector_3, double i_dy, double i_dt)
{
    // const parameter
    const double local_omega = 0.0;
    std::array<double, 4> uVec_1_L;
    std::array<double, 4> uVec_2_L;
    std::array<double, 4> uVec_3_L;
    std::array<double, 4> uVec_1_R;
    std::array<double, 4> uVec_2_R;
    std::array<double, 4> uVec_3_R;

    std::array<double, 4> Delta_i_1 = sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_1, i_uVector_0)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_2, i_uVector_1)));
    std::array<double, 4> Delta_i_2 = sumCell(scalingCell(0.5 * (1 + local_omega), diffCell(i_uVector_2, i_uVector_1)), scalingCell(0.5 * (1 - local_omega), diffCell(i_uVector_3, i_uVector_2)));

    std::array<double, 4> limiter_1 = superbeeLimiter(i_uVector_0, i_uVector_1, i_uVector_2);
    std::array<double, 4> limiter_2 = superbeeLimiter(i_uVector_1, i_uVector_2, i_uVector_3);

    // interpret left and right flux
    uVec_1_L = diffCell(i_uVector_1, scalingCell(0.5, productCell(limiter_1, Delta_i_1)));
    uVec_1_R = sumCell(i_uVector_1, scalingCell(0.5, productCell(limiter_1, Delta_i_1)));
    uVec_2_L = diffCell(i_uVector_2, scalingCell(0.5, productCell(limiter_2, Delta_i_2)));
    uVec_2_R = sumCell(i_uVector_2, scalingCell(0.5, productCell(limiter_2, Delta_i_2)));

    // locally half time evolve
    uVec_1_R = diffCell(uVec_1_R, scalingCell(0.5 * i_dt / i_dy, diffCell(conservationFlux_y(uVec_1_R), conservationFlux_y(uVec_1_L))));
    uVec_2_L = diffCell(uVec_2_L, scalingCell(0.5 * i_dt / i_dy, diffCell(conservationFlux_y(uVec_2_R), conservationFlux_y(uVec_2_L))));

    return std::array<std::array<double, 4>, 2>{uVec_1_R, uVec_2_L};
}

std::array<double, 4> FluxFunc::HLL_Riemannflux_x(std::array<double, 4> i_UL, std::array<double, 4> i_UR)
{
    double local_gamma = 1.4;

    double rho_L = i_UL[0];
    double rho_R = i_UR[0];
    double u_L = primitiveX_Vel(i_UL);
    double u_R = primitiveX_Vel(i_UR);
    double v_L = primitiveY_Vel(i_UL);
    double v_R = primitiveY_Vel(i_UR);
    double p_L = primitivePressure(i_UL);
    double p_R = primitivePressure(i_UR);
    double a_L = sqrt(local_gamma * p_L / rho_L);
    double a_R = sqrt(local_gamma * p_R / rho_R);

    double S_L = std::min(u_L - a_L, u_R - a_R);
    double S_R = std::max(u_L + a_L, u_R + a_R);

    if (S_L >= 0)
    {
        return conservationFlux_x(i_UL);
    }
    else if (S_R <= 0)
    {
        return conservationFlux_x(i_UR);
    }
    else
    {
        return scalingCell(1 / (S_R - S_L), sumCell(diffCell(scalingCell(S_R, conservationFlux_x(i_UL)), scalingCell(S_L, conservationFlux_x(i_UR))), scalingCell(S_L * S_R, diffCell(i_UR, i_UL))));
    }
}

std::array<double, 4> FluxFunc::HLL_Riemannflux_y(std::array<double, 4> i_UL, std::array<double, 4> i_UR)
{
    double local_gamma = 1.4;

    double rho_L = i_UL[0];
    double rho_R = i_UR[0];
    double u_L = primitiveX_Vel(i_UL);
    double u_R = primitiveX_Vel(i_UR);
    double v_L = primitiveY_Vel(i_UL);
    double v_R = primitiveY_Vel(i_UR);
    double p_L = primitivePressure(i_UL);
    double p_R = primitivePressure(i_UR);
    double a_L = sqrt(local_gamma * p_L / rho_L);
    double a_R = sqrt(local_gamma * p_R / rho_R);

    double S_L = std::min(v_L - a_L, v_R - a_R);
    double S_R = std::max(v_L + a_L, v_R + a_R);

    if (S_L >= 0)
    {
        return conservationFlux_y(i_UL);
    }
    else if (S_R <= 0)
    {
        return conservationFlux_y(i_UR);
    }
    else
    {
        return scalingCell(1 / (S_R - S_L), sumCell(diffCell(scalingCell(S_R, conservationFlux_y(i_UL)), scalingCell(S_L, conservationFlux_y(i_UR))), scalingCell(S_L * S_R, diffCell(i_UR, i_UL))));
    }
}

std::array<double, 4> FluxFunc::HLLC_Riemannflux_x(std::array<double, 4> i_UL, std::array<double, 4> i_UR)
{
    double local_gamma = 1.4;

    // transforming into primitive
    double rho_L = i_UL[0];
    double rho_R = i_UR[0];
    double u_L = primitiveX_Vel(i_UL);
    double u_R = primitiveX_Vel(i_UR);
    double v_L = primitiveY_Vel(i_UL);
    double v_R = primitiveY_Vel(i_UR);
    double p_L = primitivePressure(i_UL);
    double p_R = primitivePressure(i_UR);

    // Toro Step 1
    double rho_bar = 0.5 * (rho_L + rho_R);
    double a_L = sqrt(local_gamma * p_L / rho_L);
    double a_R = sqrt(local_gamma * p_R / rho_R);
    double a_bar = 0.5 * (a_L + a_R);
    double p_pvrs = 0.5 * (p_L + p_R) - 0.5 * (u_R - u_L) * rho_bar * a_bar;
    double p_star = std::max(0.0, p_pvrs);

    // Toro Step 2
    double q_L;
    double q_R;

    if (p_star > p_L)
    {
        q_L = sqrt(1 + (local_gamma + 1) / (2 * local_gamma) * (p_star / p_L - 1));
    }
    else
    {
        q_L = 1;
    }

    if (p_star > p_R)
    {
        q_R = sqrt(1 + (local_gamma + 1) / (2 * local_gamma) * (p_star / p_R - 1));
    }
    else
    {
        q_R = 1;
    }

    double S_L = std::min(u_L - a_L * q_L, u_R - a_R * q_R);
    double S_R = std::min(u_L + a_L * q_L, u_R + a_R * q_R);

    double S_star = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) / (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));

    // Toro Step 3
    if (S_L >= 0)
    {
        return conservationFlux_x(i_UL);
    }
    else if (S_R <= 0)
    {
        return conservationFlux_x(i_UR);
    }
    else
    {
        double toroScaler_L = rho_L * (S_L - u_L) / (S_L - S_star);
        double toroScaler_R = rho_R * (S_R - u_R) / (S_R - S_star);

        std::array<double, 4> F_star_L;
        std::array<double, 4> F_star_R;

        std::array<double, 4> U_star_L;
        std::array<double, 4> U_star_R;

        U_star_L[0] = 1;
        U_star_L[1] = S_star;
        U_star_L[2] = v_L;
        U_star_L[3] = i_UL[3] / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L)));
        U_star_L = scalingCell(toroScaler_L, U_star_L);

        U_star_R[0] = 1;
        U_star_R[1] = S_star;
        U_star_R[2] = v_R;
        U_star_R[3] = i_UR[3] / rho_R + (S_star - u_R) * (S_star + p_R / (rho_R * (S_R - u_R)));
        U_star_R = scalingCell(toroScaler_R, U_star_R);

        F_star_L = sumCell(conservationFlux_x(i_UL), scalingCell(S_L, diffCell(U_star_L, i_UL)));
        F_star_R = sumCell(conservationFlux_x(i_UR), scalingCell(S_R, diffCell(U_star_R, i_UR)));

        if (S_L <= 0 && S_star >= 0)
        {
            return F_star_L;
        }
        else
        {
            return F_star_R;
        }
    }
}

std::array<double, 4> FluxFunc::HLLC_Riemannflux_y(std::array<double, 4> i_UL, std::array<double, 4> i_UR)
{
    double local_gamma = 1.4;

    // transforming into primitive
    double rho_L = i_UL[0];
    double rho_R = i_UR[0];
    double u_L = primitiveX_Vel(i_UL);
    double u_R = primitiveX_Vel(i_UR);
    double v_L = primitiveY_Vel(i_UL);
    double v_R = primitiveY_Vel(i_UR);
    double p_L = primitivePressure(i_UL);
    double p_R = primitivePressure(i_UR);

    // Toro Step 1
    double rho_bar = 0.5 * (rho_L + rho_R);
    double a_L = sqrt(local_gamma * p_L / rho_L);
    double a_R = sqrt(local_gamma * p_R / rho_R);
    double a_bar = 0.5 * (a_L + a_R);
    double p_pvrs = 0.5 * (p_L + p_R) - 0.5 * (v_R - v_L) * rho_bar * a_bar;
    double p_star = std::max(0.0, p_pvrs);

    // Toro Step 2
    double q_L;
    double q_R;

    if (p_star > p_L)
    {
        q_L = sqrt(1 + (local_gamma + 1) / (2 * local_gamma) * (p_star / p_L - 1));
    }
    else
    {
        q_L = 1;
    }

    if (p_star > p_R)
    {
        q_R = sqrt(1 + (local_gamma + 1) / (2 * local_gamma) * (p_star / p_R - 1));
    }
    else
    {
        q_R = 1;
    }

    double S_L = v_L - a_L * q_L;
    double S_R = v_R + a_R * q_R;

    double S_star = (p_R - p_L + rho_L * v_L * (S_L - v_L) - rho_R * v_R * (S_R - v_R)) / (rho_L * (S_L - v_L) - rho_R * (S_R - v_R));

    // Toro Step 3
    if (S_L >= 0)
    {
        return conservationFlux_y(i_UL);
    }
    else if (S_R <= 0)
    {
        return conservationFlux_y(i_UR);
    }
    else
    {
        double toroScaler_L = rho_L * (S_L - v_L) / (S_L - S_star);
        double toroScaler_R = rho_R * (S_R - v_R) / (S_R - S_star);

        std::array<double, 4> F_star_L;
        std::array<double, 4> F_star_R;

        std::array<double, 4> U_star_L;
        std::array<double, 4> U_star_R;

        U_star_L[0] = 1;
        U_star_L[1] = u_L;
        U_star_L[2] = S_star;
        U_star_L[3] = i_UL[3] / rho_L + (S_star - v_L) * (S_star + p_L / (rho_L * (S_L - v_L)));
        U_star_L = scalingCell(toroScaler_L, U_star_L);

        U_star_R[0] = 1;
        U_star_R[1] = u_R;
        U_star_R[2] = S_star;
        U_star_R[3] = i_UR[3] / rho_R + (S_star - v_R) * (S_star + p_R / (rho_R * (S_R - v_R)));
        U_star_R = scalingCell(toroScaler_R, U_star_R);

        F_star_L = sumCell(conservationFlux_y(i_UL), scalingCell(S_L, diffCell(U_star_L, i_UL)));
        F_star_R = sumCell(conservationFlux_y(i_UR), scalingCell(S_R, diffCell(U_star_R, i_UR)));

        if (S_L <= 0 && S_star >= 0)
        {
            return F_star_L;
        }
        else
        {
            return F_star_R;
        }
    }
}
