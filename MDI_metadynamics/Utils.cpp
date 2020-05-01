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

	array3d h_box_length;
	h_box_length[0] = box_length[0] / 2;
	h_box_length[1] = box_length[1] / 2;
	h_box_length[2] = box_length[2] / 2;

	if (xyz[0] > h_box_length[0]) {
		xyz[0] -= box_length[0];
	}
	else if (xyz[0] < - h_box_length[0]) {
		xyz[0] += box_length[0];
	}

	if (xyz[1] > h_box_length[1]) {
		xyz[1] -= box_length[1];
	}
	else if (xyz[1] < - h_box_length[1]) {
		xyz[1] += box_length[1];
	}

	if (xyz[2] > h_box_length[2]) {
		xyz[2] -= box_length[2];
	}
	else if (xyz[2] < - h_box_length[2]) {
		xyz[2] += box_length[2];
	}
	
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
}
