#include "Utils.h"
#include <iostream>

using namespace std;

// double Dotp(array3d& a, array3d& b) {
// 	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
// 
// }

array3d Crossp (array3d a, array3d b) {

	array3d ret;

	ret[0] = a[1]*b[2] - a[2]*b[1];
	ret[1] = - (a[0]*b[2] - a[2]*b[0]);
	ret[2] = a[0]*b[1] - a[1]*b[0];

	return ret;
}

array3d Minimum_Image(array3d xyz, array3d box_length) {

    xyz[0] = xyz[0] - box_length[0] * round(xyz[0] / box_length[0]);
    xyz[1] = xyz[1] - box_length[1] * round(xyz[1] / box_length[1]);
    xyz[2] = xyz[2] - box_length[2] * round(xyz[2] / box_length[2]);
	return xyz;
}


double Norm(array3d& array) {
	
	double size = array[0]*array[0] + array[1]*array[1] + array[2]*array[2];
	size = sqrt(size);

	return size;
}

// array3d Normalize(array3d array) {
// 	
// 	double size = Norm(array);
// 
// 	for (int i = 0; i<=3; i++) {
// 		array[i] = array[i] / size;
// 	}
// 
// 	return array;
// }

double Gaussian(double ds, double width, double height) {
	double arg = (ds * ds) / (2. * width * width);
	return height * exp(-arg);
}

double Gaussian_derv(double ds, double width, double height) {
	double arg = (ds * ds) / (2. * width * width);
	double pre = - ds / (width * width);
	return height * pre * exp(-arg);
//	double gauss = height * exp(-ds*ds/(2.*width*width));
//	return -gauss * (ds/2.0/(width*width));
}
