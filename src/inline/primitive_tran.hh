/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   primitive_tran.hh                                 Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/10 10:02:28 by Zian Huang                               */
/*   Updated: 2023/01/21 12:22:19 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef INLINE_PRIMITIVE_TRAN_HH
#define INLINE_PRIMITIVE_TRAN_HH

#include <array>

inline double primitiveX_Vel(std::array<double, 4> i_uVector_i)
{
    return i_uVector_i[1] / i_uVector_i[0];
}

inline double primitiveY_Vel(std::array<double, 4> i_uVector_i)
{
    return i_uVector_i[2] / i_uVector_i[0];
}

inline double primitivePressure(std::array<double, 4> i_uVector_i)
{
    // only applicable to air modelling
    double localGamma = 1.4;

    return (localGamma - 1) * (i_uVector_i[3] - (i_uVector_i[1] * i_uVector_i[1] + i_uVector_i[2] * i_uVector_i[2]) / 2 / i_uVector_i[0]);
}

#endif