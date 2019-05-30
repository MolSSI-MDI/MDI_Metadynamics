#include "Dihedral.h"

using namespace std;

Dihedral::Dihedral(int a1, int a2, int a3, int a4) {
	natoms_ = 4;
        atomi_ = a1;
        atomj_ = a2;
        atomk_ = a3;
        atoml_ = a4;
	gradi_ = {0.0, 0.0, 0.0};			
	gradj_ = {0.0, 0.0, 0.0};			
	gradk_ = {0.0, 0.0, 0.0};			
	gradl_ = {0.0, 0.0, 0.0};			

}

array4dint Dihedral::Get_Atoms()  {

	array4dint atoms;
	atoms[0] = atomi_;
	atoms[1] = atomj_;
	atoms[2] = atomk_;
	atoms[3] = atoml_;
	return atoms;
}

double Dihedral::Get_Value()  {
	return value_;
}

array<array3d, 4> Dihedral::Get_Gradient()  {
	return {gradi_,gradj_,gradk_,gradl_};
}


void Dihedral::Evaluate(const double xyz[], int natoms, array3d box_len)  {


	array3d atomi_xyz, atomj_xyz, atomk_xyz, atoml_xyz;

	atomi_xyz[0]=*(xyz+3*atomi_)*bohr_to_ang;
	atomi_xyz[1]=*(xyz+3*atomi_+1)*bohr_to_ang;
	atomi_xyz[2]=*(xyz+3*atomi_+2)*bohr_to_ang;

	atomj_xyz[0]=*(xyz+3*atomj_)*bohr_to_ang;
	atomj_xyz[1]=*(xyz+3*atomj_+1)*bohr_to_ang;
	atomj_xyz[2]=*(xyz+3*atomj_+2)*bohr_to_ang;

	atomk_xyz[0]=*(xyz+3*atomk_)*bohr_to_ang;
	atomk_xyz[1]=*(xyz+3*atomk_+1)*bohr_to_ang;
	atomk_xyz[2]=*(xyz+3*atomk_+2)*bohr_to_ang;
	
	atoml_xyz[0]=*(xyz+3*atoml_)*bohr_to_ang;
	atoml_xyz[1]=*(xyz+3*atoml_+1)*bohr_to_ang;
	atoml_xyz[2]=*(xyz+3*atoml_+2)*bohr_to_ang;

	array3d rij_, rkj_, rkl_, m_, n_; 

	rij_[0] = atomi_xyz[0] - atomj_xyz[0];
	rij_[1] = atomi_xyz[1] - atomj_xyz[1];
	rij_[2] = atomi_xyz[2] - atomj_xyz[2];

	rkj_[0] = atomk_xyz[0] - atomj_xyz[0];
	rkj_[1] = atomk_xyz[1] - atomj_xyz[1];
	rkj_[2] = atomk_xyz[2] - atomj_xyz[2];

	rkl_[0] = atomk_xyz[0] - atoml_xyz[0];
	rkl_[1] = atomk_xyz[1] - atoml_xyz[1];
	rkl_[2] = atomk_xyz[2] - atoml_xyz[2];

	rij_ = Minimum_Image(rij_, box_len);
	rkj_ = Minimum_Image(rkj_, box_len);
	rkl_ = Minimum_Image(rkl_, box_len);

	m_=Crossp(rij_, rkj_);
	n_=Crossp(rkj_, rkl_);

	array3d rkj_hat_ = Normalize(rkj_);

	double X = Dotp(Crossp(m_, n_), rkj_hat_);
        double Y = Dotp(m_, n_);	

	value_ = atan2(X, Y);

	// Now, get the gradient
	
	double rkj_norm = Norm(rkj_);
	double mdotm = Dotp(m_, m_);
	double ndotn = Dotp(n_, n_);
	array3d rik_, rlj_;
	for (int i=0; i < 3; i++) {
		rik_[i] = rij_[i] - rkj_[i];
		rlj_[i] = rkl_[i] - rkj_[i];
	}
	array3d m_over_mbar2;
	array3d n_over_nbar2;
	double gradi_pref1;
	double gradj_pref1, gradj_pref2;
	double gradk_pref1, gradk_pref2;
	double gradl_pref1;

	gradj_pref1 = Dotp(rik_, rkj_) / rkj_norm;
	gradj_pref2 = Dotp(rkl_, rkj_) / rkj_norm;

	gradk_pref1 = - Dotp(rlj_, rkj_) / rkj_norm;
	gradk_pref2 = - Dotp(rij_, rkj_) / rkj_norm;

	for (int i=0; i < 3; i++) {
		m_over_mbar2[i] = m_[i] / mdotm;
		n_over_nbar2[i] = n_[i] / ndotn;

		gradi_[i] = rkj_norm * m_over_mbar2[i];
		gradl_[i] = - rkj_norm * n_over_nbar2[i];

		gradj_[i] = gradj_pref1 * m_over_mbar2[i];
		gradj_[i] += gradj_pref2 * n_over_nbar2[i];

		gradk_[i] = gradk_pref1 * n_over_nbar2[i];
		gradk_[i] += gradk_pref2 * m_over_mbar2[i];
	}
}