#ifndef LEVEL_SET_FUNCTIONS_HH
#define LEVEL_SET_FUNCTIONS_HH

#include "constants.hh"

class LevelSetFunctions {
public:
    // Single circle level set function
    static double singleCircle(double radius, double center_x, double center_y, double x, double y);
    
    // Double circle level set function
    static double doubleCircle(double r1, double r2, 
                              double center_x1, double center_x2, 
                              double center_y1, double center_y2, 
                              double x, double y);
    
    // Single square level set function
    static double singleSquare(double side_length, double center_x, double center_y, double x, double y);
    
    // 1D wall level set function
    static double oneDWall(double center_x, double x, double y);
    
private:
    // Helper function for square level set calculation
    static double calculateSquareDistance(double relative_x, double relative_y, double half_length);
};

#endif