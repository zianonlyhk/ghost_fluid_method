#include "level_set_functions.hh"
#include <cmath>
#include <algorithm>

double LevelSetFunctions::singleCircle(double radius, double center_x, double center_y, double x, double y) {
    double dx = x - center_x;
    double dy = y - center_y;
    double distance = std::sqrt(dx * dx + dy * dy);
    return -(radius - distance);
}

double LevelSetFunctions::doubleCircle(double r1, double r2, 
                                      double center_x1, double center_x2, 
                                      double center_y1, double center_y2, 
                                      double x, double y) {
    double phi1 = singleCircle(r1, center_x1, center_y1, x, y);
    double phi2 = singleCircle(r2, center_x2, center_y2, x, y);
    return std::min(phi1, phi2);
}

double LevelSetFunctions::singleSquare(double side_length, double center_x, double center_y, double x, double y) {
    double relative_x = x - center_x;
    double relative_y = y - center_y;
    double half_length = side_length / 2.0;
    
    return calculateSquareDistance(relative_x, relative_y, half_length);
}

double LevelSetFunctions::oneDWall(double center_x, double x, double y) {
    return -x + center_x;
}

double LevelSetFunctions::calculateSquareDistance(double relative_x, double relative_y, double half_length) {
    if (std::abs(relative_x) < 1e-10) {  // relative_x == 0
        if (relative_y > 0) {  // north
            return -(half_length - relative_y);
        } else {  // south
            return -(half_length + relative_y);
        }
    }
    
    double phi;
    if (relative_x > 0) {
        phi = std::atan(relative_y / relative_x);
        if (phi < Constants::PI / 4 && phi > -Constants::PI / 4) {  // east
            return -(half_length - relative_x);
        } else if (phi <= -Constants::PI / 4) {  // south
            return -(half_length + relative_y);
        } else {  // north
            return -(half_length - relative_y);
        }
    } else {
        phi = std::atan(-relative_y / relative_x);
        if (phi < Constants::PI / 4 && phi > -Constants::PI / 4) {  // west
            return -(half_length + relative_x);
        } else if (phi <= -Constants::PI / 4) {  // south
            return -(half_length + relative_y);
        } else {  // north
            return -(half_length - relative_y);
        }
    }
}