#include "Distance.h"
#include <iostream>

using namespace std;

Distance::Distance(int a1, int a2) {
	natoms_ = 2;
    atomi_ = a1-1;
    atomj_ = a2-1;
	gradi_ = {0.0, 0.0, 0.0};			
	gradj_ = {0.0, 0.0, 0.0};
	value_ = 0.0;
}

array2dint Distance::Get_Atoms()  {

	array2dint atoms;
	atoms[0] = atomi_;
	atoms[1] = atomj_;
	return atoms;
}

double Distance::Get_Value()  {
	return value_;
}

array<array3d, 2> Distance::Get_Gradient()  {
	return {gradi_,gradj_};
}


void Distance::Evaluate(const double xyz[], int natoms, array3d box_len)  {

	array3d atomi_xyz, atomj_xyz, atomk_xyz, atoml_xyz;

	atomi_xyz[0]=*(xyz+3*atomi_);
	atomi_xyz[1]=*(xyz+3*atomi_+1);
	atomi_xyz[2]=*(xyz+3*atomi_+2);

	atomj_xyz[0]=*(xyz+3*atomj_);
	atomj_xyz[1]=*(xyz+3*atomj_+1);
	atomj_xyz[2]=*(xyz+3*atomj_+2);

	array3d rij_, m_, n_; 

	rij_[0] = atomi_xyz[0] - atomj_xyz[0];
	rij_[1] = atomi_xyz[1] - atomj_xyz[1];
	rij_[2] = atomi_xyz[2] - atomj_xyz[2];

	rij_ = Minimum_Image(rij_, box_len);

	value_ = Norm(rij_);

	for (int i=0; i < 3; i++) {
		gradi_[i] = rij_[i] / value_;
		gradj_[i] = - rij_[i] / value_;
	}
}
