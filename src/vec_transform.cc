/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.cc                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/01/21 10:45:26 by Zian Huang                               */
/*   Updated: 2023/01/21 18:00:52 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "vec_transform.hh"
#include "inline/cell_operation.hh"
#include "inline/primitive_tran.hh"

// Definitions #####################################################################################

VecTran::VecTran() {}

std::vector<std::vector<std::array<double, 4>>> VecTran::slicVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dx, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dx, diffCell(fluxFunc.slicFlux_x(i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_inputU_Vec[iter_y][iter_x + 2], i_dx, i_dt), fluxFunc.slicFlux_x(i_inputU_Vec[iter_y][iter_x - 2], i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_dx, i_dt))));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::slicVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dy, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dy, diffCell(fluxFunc.slicFlux_y(i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_inputU_Vec[iter_y + 2][iter_x], i_dy, i_dt), fluxFunc.slicFlux_y(i_inputU_Vec[iter_y - 2][iter_x], i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_dy, i_dt))));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::musclHancockVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dx, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dx, diffCell(fluxFunc.musclHancockFlux_x(i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_inputU_Vec[iter_y][iter_x + 2], i_dx, i_dt), fluxFunc.musclHancockFlux_x(i_inputU_Vec[iter_y][iter_x - 2], i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_dx, i_dt))));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::musclHancockVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dy, double i_dt)
{
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dy, diffCell(fluxFunc.musclHancockFlux_y(i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_inputU_Vec[iter_y + 2][iter_x], i_dy, i_dt), fluxFunc.musclHancockFlux_y(i_inputU_Vec[iter_y - 2][iter_x], i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_dy, i_dt))));
        }
    }

    return toBeReturnVec;
}
