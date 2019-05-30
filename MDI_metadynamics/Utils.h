#pragma once

#include <array>
#include <cmath>

typedef std::array<int, 4> array4dint;
typedef std::array<double, 3> array3d;
typedef std::array<double, 2> array2d;
const double bohr_to_ang = 0.529177;

double Dotp(array3d a, array3d b);
array3d Crossp (array3d a, array3d b);
array3d Minimum_Image(array3d xyz, array3d box_length);
double Norm(array3d array);
array3d Normalize(array3d array);
double Gaussian(double ds, double width, double height);
double Gaussian_derv(double ds, double width, double height);
