/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   gfm_2d_euler_solver.cc                            Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:53:57 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "gfm_2d_euler_solver.hh"
#include "inline/primitive_tran.hh"
#include "inline/debug_tools.hh"
#include <math.h>

double calcMaxA_eachArr(std::array<double, 4> i_uVec, double i_gamma)
{
    double maxA = sqrt(primitiveX_Vel(i_uVec) * primitiveX_Vel(i_uVec) + primitiveY_Vel(i_uVec) * primitiveY_Vel(i_uVec)) + sqrt(i_gamma * primitivePressure(i_uVec) / i_uVec[0]);

    return maxA;
}

void writeToFileStream(std::ofstream &i_fstream, std::vector<std::vector<std::array<double, 4>>> const &i_inputVec, double i_x0, double i_dx, double i_y0, double i_dy, double i_t, int i_idx)
{
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            i_fstream << i_t << ' ' << i_x0 + (i - 2) * i_dx << ' ' << i_y0 + (j - 2) * i_dy << ' ' << i_inputVec[j][i][i_idx] << std::endl;
        }

        i_fstream << std::endl;
    }
    i_fstream << std::endl;
}

void writeToFileStream(std::ofstream &i_fstream, std::vector<std::vector<double>> const &i_mockschlieren, double i_x0, double i_dx, double i_y0, double i_dy, double i_t)
{
    int nCell_y = i_mockschlieren.size() - 4;
    int nCell_x = i_mockschlieren[0].size() - 4;

    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            i_fstream << i_t << ' ' << i_x0 + (i - 2) * i_dx << ' ' << i_y0 + (j - 2) * i_dy << ' ' << i_mockschlieren[j][i] << std::endl;
        }

        i_fstream << std::endl;
    }
    i_fstream << std::endl;
}

// Definitions #####################################################################################

GFM_2D_EulerSolver::GFM_2D_EulerSolver() {}

GFM_2D_EulerSolver::GFM_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> i_uVec, int i_nCell_x, int i_nCell_y)
{
    m_nCell_x = i_nCell_x;
    m_nCell_y = i_nCell_y;
    m_uVec = i_uVec;
}

GFM_2D_EulerSolver::GFM_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> i_uVec, int i_nCell_x, int i_nCell_y, double i_x0, double i_x1, double i_y0, double i_y1, double i_tStop, double i_c)
{
    m_nCell_x = i_nCell_x;
    m_nCell_y = i_nCell_y;
    m_x0 = i_x0;
    m_x1 = i_x1;
    m_y0 = i_y0;
    m_y1 = i_y1;
    m_c = i_c;
    m_tStop = i_tStop;
    m_dx = (i_x1 - i_x0) / i_nCell_x;
    m_dy = (i_y1 - i_y0) / i_nCell_y;
    m_uVec = i_uVec;
}

void GFM_2D_EulerSolver::setBound(double i_x0, double i_x1, double i_y0, double i_y1, double i_tStop)
{
    m_x0 = i_x0;
    m_x1 = i_x1;
    m_y0 = i_y0;
    m_y1 = i_y1;
    m_tStop = i_tStop;

    // implied values
    m_dx = (i_x1 - i_x0) / m_nCell_x;
    m_dy = (i_y1 - i_y0) / m_nCell_y;
}

void GFM_2D_EulerSolver::setCFL(double i_c)
{
    m_c = i_c;
}

void GFM_2D_EulerSolver::setLevelSet(std::vector<std::vector<double>> i_levelSet)
{
    m_levelSet = i_levelSet;
}

void GFM_2D_EulerSolver::setName(std::string i_name)
{
    m_name = i_name;
}

void GFM_2D_EulerSolver::setRepoDir(std::string repoDir)
{
    m_repoDir = repoDir;
}

void GFM_2D_EulerSolver::setRigidBodyVel(std::array<double, 2> velArr)
{
    m_rigidBodyVel = velArr;
}

void GFM_2D_EulerSolver::setRigidBodyCentreCoor(std::array<double, 2> coorArr)
{
    m_rigidBodyCentreCoor = coorArr;
}

void GFM_2D_EulerSolver::updateMaxA(int i_numIter)
{
    double localMaxA = 0.0;
    double dummieMaxA;

    for (int j = 2; j < m_nCell_y + 2; ++j)
    {
        for (int i = 2; i < m_nCell_x + 2; ++i)
        {
            dummieMaxA = calcMaxA_eachArr(m_uVec[j][i], m_gamma);

            if (dummieMaxA > localMaxA)
            {
                localMaxA = dummieMaxA;
            }
        }
    }

    // in the first 10 iterations of the numerical scheme, bigger aMax
    if (i_numIter < 10)
    {
        m_aMax = localMaxA * 3;
    }
    else
    {
        m_aMax = localMaxA;
    }
}

void GFM_2D_EulerSolver::updateDt()
{
    double smallerD;
    if (m_dx < m_dy)
    {
        smallerD = m_dx;
    }
    else
    {
        smallerD = m_dy;
    }

    m_dt = m_c * smallerD / m_aMax;
}

void GFM_2D_EulerSolver::advectLevelSet()
{
    m_levelSet = vecTran.levelSetAdvectionTransform(m_levelSet, m_rigidBodyVel, m_dx, m_dy, m_dt);
    m_rigidBodyCentreCoor[0] += m_rigidBodyVel[0] * m_dt;
    m_rigidBodyCentreCoor[1] += m_rigidBodyVel[1] * m_dt;
}

void GFM_2D_EulerSolver::updateLevelSetBoundaryTrans()
{
    for (int i = 2; i < m_nCell_y + 2; ++i)
    {
        m_levelSet[i][0] = m_levelSet[i][3];
        m_levelSet[i][1] = m_levelSet[i][2];
        m_levelSet[i][m_nCell_x + 3] = m_levelSet[i][m_nCell_x];
        m_levelSet[i][m_nCell_x + 2] = m_levelSet[i][m_nCell_x + 1];
    }

    for (int i = 2; i < m_nCell_x + 2; ++i)
    {
        m_levelSet[0][i] = m_levelSet[3][i];
        m_levelSet[1][i] = m_levelSet[2][i];
        m_levelSet[m_nCell_y + 3][i] = m_levelSet[m_nCell_y][i];
        m_levelSet[m_nCell_y + 2][i] = m_levelSet[m_nCell_y + 1][i];
    }
}

void GFM_2D_EulerSolver::accelerateRigidBody_circ()
{
    std::array<double, 2> rVec = {(m_x1 + m_x0) / 2 - m_rigidBodyCentreCoor[0], (m_y1 + m_y0) / 2 - m_rigidBodyCentreCoor[1]};
    double rNorm = sqrt(rVec[0] * rVec[0] + rVec[1] * rVec[1]);
    std::array<double, 2> rVecNormalised = {rVec[0] / rNorm, rVec[1] / rNorm};

    double aNorm = (m_rigidBodyVel[0] * m_rigidBodyVel[0] + m_rigidBodyVel[1] * m_rigidBodyVel[1]) / rNorm;
    std::array<double, 2> acceleration = {aNorm * rVecNormalised[0], aNorm * rVecNormalised[1]};

    m_rigidBodyVel[0] += acceleration[0] * m_dt;
    m_rigidBodyVel[1] += acceleration[1] * m_dt;
}

void GFM_2D_EulerSolver::reinitLevelSet()
{
    m_levelSet = vecTran.reinitLevelSetInsideRB(m_levelSet, m_dx, m_dy);
}

void GFM_2D_EulerSolver::updateGhostCellBoundary()
{
    m_uVec = vecTran.ghostCellBoundary(m_uVec, m_levelSet, m_dx, m_dy, m_rigidBodyVel);
}

void GFM_2D_EulerSolver::propagateGhostCell()
{
    m_uVec = vecTran.propagateGhostInterface(m_uVec, m_levelSet, m_dx, m_dy);
}

void GFM_2D_EulerSolver::updateBoundaryTrans()
{
    for (int i = 2; i < m_nCell_y + 2; ++i)
    {
        m_uVec[i][0] = m_uVec[i][3];
        m_uVec[i][1] = m_uVec[i][2];
        m_uVec[i][m_nCell_x + 3] = m_uVec[i][m_nCell_x];
        m_uVec[i][m_nCell_x + 2] = m_uVec[i][m_nCell_x + 1];
    }

    for (int i = 2; i < m_nCell_x + 2; ++i)
    {
        m_uVec[0][i] = m_uVec[3][i];
        m_uVec[1][i] = m_uVec[2][i];
        m_uVec[m_nCell_y + 3][i] = m_uVec[m_nCell_y][i];
        m_uVec[m_nCell_y + 2][i] = m_uVec[m_nCell_y + 1][i];
    }
}

void GFM_2D_EulerSolver::mhHllcSweepX()
{
    m_uVec = vecTran.musclHancockVecTranHLLC_x(m_uVec, m_dx, m_dt);
}

void GFM_2D_EulerSolver::mhHllcSweepY()
{
    m_uVec = vecTran.musclHancockVecTranHLLC_y(m_uVec, m_dy, m_dt);
}

void GFM_2D_EulerSolver::slicSweepX()
{
    m_uVec = vecTran.slicVecTran_x(m_uVec, m_dx, m_dt);
}

void GFM_2D_EulerSolver::slicSweepY()
{
    m_uVec = vecTran.slicVecTran_y(m_uVec, m_dy, m_dt);
}

void GFM_2D_EulerSolver::calculateMockSchliren()
{
    m_mockschlieren = vecTran.mockSchlierenTrans(m_uVec, m_dx, m_dy);
}

void GFM_2D_EulerSolver::initiateDataLogging()
{
    m_rhoResults.open(m_repoDir + (std::string) "data/" + m_name + (std::string) "_rhoResults.dat");
    m_momentumX_Results.open(m_repoDir + (std::string) "data/" + m_name + (std::string) "_momentumX_Results.dat");
    m_momentumY_Results.open(m_repoDir + (std::string) "data/" + m_name + (std::string) "_momentumY_Results.dat");
    m_energyResults.open(m_repoDir + (std::string) "data/" + m_name + (std::string) "_energyResults.dat");
    m_msResults.open(m_repoDir + (std::string) "data/" + m_name + (std::string) "_msResults.dat");
    m_levelSetResults.open(m_repoDir + (std::string) "data/" + m_name + (std::string) "_levelSetResults.dat");

    writeToFileStream(m_rhoResults, m_uVec, m_x0, m_dx, m_y0, m_dy, 0, 0);
    writeToFileStream(m_momentumX_Results, m_uVec, m_x0, m_dx, m_y0, m_dy, 0, 1);
    writeToFileStream(m_momentumY_Results, m_uVec, m_x0, m_dx, m_y0, m_dy, 0, 2);
    writeToFileStream(m_energyResults, m_uVec, m_x0, m_dx, m_y0, m_dy, 0, 3);
    writeToFileStream(m_msResults, m_mockschlieren, m_x0, m_dx, m_y0, m_dy, 0);
    writeToFileStream(m_levelSetResults, m_levelSet, m_x0, m_dx, m_y0, m_dy, 0);
}

void GFM_2D_EulerSolver::writeToFiles(double i_time)
{
    writeToFileStream(m_rhoResults, m_uVec, m_x0, m_dx, m_y0, m_dy, i_time, 0);
    writeToFileStream(m_momentumX_Results, m_uVec, m_x0, m_dx, m_y0, m_dy, i_time, 1);
    writeToFileStream(m_momentumY_Results, m_uVec, m_x0, m_dx, m_y0, m_dy, i_time, 2);
    writeToFileStream(m_energyResults, m_uVec, m_x0, m_dx, m_y0, m_dy, i_time, 3);
    writeToFileStream(m_msResults, m_mockschlieren, m_x0, m_dx, m_y0, m_dy, i_time);
    writeToFileStream(m_levelSetResults, m_levelSet, m_x0, m_dx, m_y0, m_dy, i_time);
}

void GFM_2D_EulerSolver::cleanUp()
{
    m_rhoResults.close();
    m_momentumX_Results.close();
    m_momentumY_Results.close();
    m_energyResults.close();
    m_msResults.close();
    m_levelSetResults.close();
}

const int GFM_2D_EulerSolver::nCell_x() { return m_nCell_x; }
const int GFM_2D_EulerSolver::nCell_y() { return m_nCell_y; }
const double GFM_2D_EulerSolver::x0() { return m_x0; }
const double GFM_2D_EulerSolver::x1() { return m_x1; }
const double GFM_2D_EulerSolver::dx() { return m_dx; }
const double GFM_2D_EulerSolver::y0() { return m_y0; }
const double GFM_2D_EulerSolver::y1() { return m_y1; }
const double GFM_2D_EulerSolver::dy() { return m_dy; }
const double GFM_2D_EulerSolver::tStart() { return m_tStart; }
const double GFM_2D_EulerSolver::tStop() { return m_tStop; }
const double GFM_2D_EulerSolver::c() { return m_c; }
const double GFM_2D_EulerSolver::gamma() { return m_gamma; }
const double GFM_2D_EulerSolver::aMax() { return m_aMax; }
const double GFM_2D_EulerSolver::dt() { return m_dt; }
const std::array<double, 2> GFM_2D_EulerSolver::rigidBodyVel() { return m_rigidBodyVel; }
const std::vector<std::vector<std::array<double, 4>>> GFM_2D_EulerSolver::uVec() { return m_uVec; }
const std::vector<std::vector<double>> GFM_2D_EulerSolver::levelSet() { return m_levelSet; }

void GFM_2D_EulerSolver::printBoundaryCoor()
{
    printBoundaryCellCoor(vecTran.getBoundaryCellCoor(m_levelSet));
}