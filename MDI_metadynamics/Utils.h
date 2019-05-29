#include <array>
#include <cmath>

using namespace std;

typedef array<int, 4> array4dint;
typedef array<double, 3> array3d;
typedef array<double, 2> array2d;
const double bohr_to_ang = 0.529177;

double Dotp(array3d a, array3d b) {

	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

}

array3d Crossp (array3d a, array3d b) {

	array3d ret;

	ret[0] = a[1]*b[2] - a[2]*b[1];
	ret[1] = - (a[0]*b[2] - a[2]*b[0]);
	ret[2] = a[0]*b[1] - a[1]*b[0];

	return ret;
}

array3d Minimum_Image(array3d xyz, array3d box_length) {

	if (xyz[0] > box_length[0]) {
		xyz[0] -= box_length[0];
	}
	else if (xyz[0] < - box_length[0]) {
		xyz[0] += box_length[0];
	}
	return xyz;

	if (xyz[1] > box_length[1]) {
		xyz[1] -= box_length[1];
	}
	else if (xyz[1] < - box_length[1]) {
		xyz[1] += box_length[1];
	}


	if (xyz[2] > box_length[2]) {
		xyz[2] -= box_length[2];
	}
	else if (xyz[2] < - box_length[2]) {
		xyz[2] += box_length[2];
	}
	return xyz;
}


double Norm(array3d array) {
	
	double size = Dotp(array, array);
	size = sqrt(size);

	return size;
}

array3d Normalize(array3d array) {
	
	double size = Norm(array);

	for (int i = 0; i<=3; i++) {
		array[i] = array[i] / size;
	}

	return array;
}

double Gaussian(double ds, double width, double height) {
	double arg = (ds * ds) / (2. * width * width);
	return height * exp(-arg);
}

double Gaussian_derv(double ds, double width, double height) {
	double arg = (ds * ds) / (2. * width * width);
	double pre = - ds / (width * width);
	return pre * exp(-arg);
}
