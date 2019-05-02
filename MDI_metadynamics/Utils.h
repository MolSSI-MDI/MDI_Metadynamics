#include <array>
#include <cmath>

using namespace std;

typedef std::array<double, 3> array3d;
typedef std::array<double, 2> array2d;

double dotp(array3d a, array3d b) {

	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

}

array3d crossp (array3d a, array3d b) {

	array3d ret;

	ret[0] = a[1]*b[2] - a[2]*b[1];
	ret[1] = - (a[0]*b[2] - a[2]*b[0]);
	ret[2] = a[0]*b[1] - a[1]*b[0];

	return ret;
}

array3d Minimum_Image(array3d array, double box_length) {

	if (array[0] > box_length) {
		array[0] -= box_length;
	}
	else if (array[0] < - box_length) {
		array[0] += box_length;
	}
	return array;

	if (array[1] > box_length) {
		array[1] -= box_length;
	}
	else if (array[1] < - box_length) {
		array[1] += box_length;
	}


	if (array[2] > box_length) {
		array[2] -= box_length;
	}
	else if (array[2] < - box_length) {
		array[2] += box_length;
	}
	return array;
}

array3d Normalize(array3d array) {
	
	double norm = array[0] * array[0] + array[1] * array[1] + array[2] * array[2];
	norm = sqrt(norm);

	for (int i = 0; i<=3; i++) {
		array[i] = array[i] / norm;
	}

	return array;
}
