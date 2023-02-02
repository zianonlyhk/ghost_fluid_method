/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   gfm_2d_euler_solver.hh                            Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2023/02/02 14:53:53 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef GFM_2D_EULER_SOLVER_HH
#define GFM_2D_EULER_SOLVER_HH

#include "vec_transform.hh"
#include <string>
#include <fstream>

class GFM_2D_EulerSolver
{
public:
    GFM_2D_EulerSolver();
    GFM_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> uVec, int nCell_x, int nCell_y);
    GFM_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> uVec, int nCell_x, int nCell_y, double x0, double x1, double y0, double y1, double tStop, double c);

    // ADDITIONAL CONSTRUCTION
    void setBound(double x0, double x1, double y0, double y1, double tStop);
    void setCFL(double c);
    void setLevelSet(std::vector<std::vector<double>> levelSet);
    void setName(std::string name);
    void setRepoDir(std::string repoDir);

    // NUMERICAL SCHEME TOOLBOX
    VecTran vecTran;

    void updateMaxA(int numIter);
    void updateDt();

    void updateGhostCellBoundary();
    void cleanupGhostRegion();
    void propagateGhostCell();

    void updateBoundaryTrans();

    void mhHllcSweepX();
    void mhHllcSweepY();
    void slicSweepX();
    void slicSweepY();

    void initiateDataLogging();
    void writeToFiles(double time);
    void cleanUp();

    void printBoundaryCoor();

    // ACCESSING PRIVATE MEMBERS
    const int nCell_x();
    const int nCell_y();
    const double x0();
    const double x1();
    const double dx();
    const double y0();
    const double y1();
    const double dy();
    const double tStart();
    const double tStop();
    const double c();
    const double gamma();
    const double aMax();
    const double dt();
    const std::vector<std::vector<std::array<double, 4>>> uVec();
    const std::vector<std::vector<double>> levelSet();

private:
    // CONSTANT ATTRIBUTES
    int m_nCell_x;
    int m_nCell_y;
    double m_x0;
    double m_x1;
    double m_dx;
    double m_y0;
    double m_y1;
    double m_dy;
    double m_tStart = 0.0;
    double m_tStop;
    double m_c;
    double m_gamma = 1.4;

    std::string m_name;
    std::string m_repoDir;
    std::ofstream m_rhoResults;
    std::ofstream m_momentumX_Results;
    std::ofstream m_momentumY_Results;
    std::ofstream m_energyResults;

    // DYNAMIC ATTRIBUTES
    double m_aMax;
    double m_dt;

    std::vector<std::vector<std::array<double, 4>>> m_uVec;
    std::vector<std::vector<double>> m_levelSet;
};

#endif