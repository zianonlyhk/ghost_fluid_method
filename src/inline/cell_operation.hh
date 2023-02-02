/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   cell_operation.hh                                 Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:54:22 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef INLINE_CELL_OPERATION_HH
#define INLINE_CELL_OPERATION_HH

#include <array>

inline std::array<double, 4> sumCell(std::array<double, 4> i_firstCell, std::array<double, 4> i_secondCell)
{
    std::array<double, 4> cellToBeReturned;

    for (int i = 0; i < 4; ++i)
    {
        cellToBeReturned[i] = i_firstCell[i] + i_secondCell[i];
    }

    return cellToBeReturned;
}

inline std::array<double, 4> diffCell(std::array<double, 4> i_firstCell, std::array<double, 4> i_secondCell)
{
    std::array<double, 4> cellToBeReturned;

    for (int i = 0; i < 4; ++i)
    {
        cellToBeReturned[i] = i_firstCell[i] - i_secondCell[i];
    }

    return cellToBeReturned;
}

inline std::array<double, 4> productCell(std::array<double, 4> i_firstCell, std::array<double, 4> i_secondCell)
{
    std::array<double, 4> cellToBeReturned;

    for (int i = 0; i < 4; ++i)
    {
        cellToBeReturned[i] = i_firstCell[i] * i_secondCell[i];
    }

    return cellToBeReturned;
}

inline std::array<double, 4> divisionCell(std::array<double, 4> i_firstCell, std::array<double, 4> i_secondCell)
{
    std::array<double, 4> cellToBeReturned;

    for (int i = 0; i < 4; ++i)
    {
        cellToBeReturned[i] = i_firstCell[i] / i_secondCell[i];
    }

    return cellToBeReturned;
}

inline std::array<double, 4> scalingCell(double i_factor, std::array<double, 4> i_inputCell)
{
    std::array<double, 4> cellToBeReturned;

    for (int i = 0; i < 4; ++i)
    {
        cellToBeReturned[i] = i_factor * i_inputCell[i];
    }

    return cellToBeReturned;
}

inline std::array<double, 4> scalarAdditionCell(double i_added, std::array<double, 4> i_inputCell)
{
    std::array<double, 4> cellToBeReturned;

    for (int i = 0; i < 4; ++i)
    {
        cellToBeReturned[i] = i_added + i_inputCell[i];
    }

    return cellToBeReturned;
}

inline std::array<double, 4> scalarDivisionCell(double i_nominator, std::array<double, 4> i_inputCell)
{
    std::array<double, 4> cellToBeReturned;

    for (int i = 0; i < 4; ++i)
    {
        cellToBeReturned[i] = i_nominator / i_inputCell[i];
    }

    return cellToBeReturned;
}

#endif