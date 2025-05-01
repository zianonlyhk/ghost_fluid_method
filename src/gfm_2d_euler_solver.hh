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
    void setRigidBodyVel(std::array<double, 2> velArr);
    void setRigidBodyCentreCoor(std::array<double, 2> coorArr);

    // NUMERICAL SCHEME TOOLBOX
    VecTran vecTran;

    void updateMaxA(int numIter);
    void updateDt();

    void advectLevelSet();
    void updateLevelSetBoundaryTrans();
    void accelerateRigidBody_circ();
    void reinitLevelSet();

    void updateGhostCellBoundary();
    void propagateGhostCell();

    void updateBoundaryTrans();
    void updateBoundaryReflect();

    void mhHllcSweepX();
    void mhHllcSweepY();
    void slicSweepX();
    void slicSweepY();

    void calculateMockSchliren();

    void initiateOutputLogging();
    void writeToFiles(double time);
    void cleanUp();

    // DEBUG
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
    const std::array<double, 2> rigidBodyVel();
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
    std::ofstream m_velMagResults;
    std::ofstream m_pressureResults;
    std::ofstream m_itnEnergyResults;
    std::ofstream m_msResults;
    std::ofstream m_levelSetResults;

    // DYNAMIC ATTRIBUTES
    double m_aMax;
    double m_dt;

    std::array<double, 2> m_rigidBodyCentreCoor;
    std::array<double, 2> m_rigidBodyVel;

    std::vector<std::vector<std::array<double, 4>>> m_uVec;
    std::vector<std::vector<double>> m_levelSet;
    std::vector<std::vector<double>> m_mockschlieren;
};

#endif