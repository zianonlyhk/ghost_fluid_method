#ifndef CONSTANTS_HH
#define CONSTANTS_HH

namespace Constants {
    // Mathematical constants
    constexpr double PI = 3.14159265358979323846;
    
    // Physical constants
    constexpr double GAMMA = 1.4;  // Heat capacity ratio for air
    
    // Numerical constants
    constexpr double LARGE_VALUE = 999.0;  // Used for level set initialization
    constexpr int GHOST_CELLS = 2;  // Number of ghost cells on each side
    constexpr int TOTAL_GHOST_CELLS = 4;  // Total ghost cells (both sides)
    
    // Rigid body types
    namespace RigidBodyType {
        constexpr int CIRCLE = 1;
        constexpr int SQUARE = 2;
        constexpr int DOUBLE_CIRCLE = 3;
    }
    
    // Acceleration flags
    namespace AccelerationFlag {
        constexpr int DISABLED = 0;
        constexpr int ENABLED = 1;
    }
}

#endif