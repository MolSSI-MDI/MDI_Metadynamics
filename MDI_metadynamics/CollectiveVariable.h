#include <array>
#include <cmath>
#include "Utils.h"

using namespace std;

class CollectiveVariable
{
	protected:
		array2d bounds_;
		double value_;

	public:

		CollectiveVariable(): value_(0.0), bounds_{0,0} {}
		virtual void Evaluate(const double xyz[], int, array3d) = 0;
		virtual void Compute_gradient(const double xyz[], int, array3d) = 0;
		virtual double Get_Value() = 0;
		virtual array<array3d, 4> Get_Gradient() = 0;

};


class Dihedral: public CollectiveVariable
{
        private:
                int atomi_, atomj_, atomk_, atoml_;
		array3d rij_, rkj_, rkl_, m_, n_; 
		array3d gradi_, gradj_, gradk_, gradl_;


        public:

                Dihedral(int a1, int a2, int a3, int a4) {
                        atomi_ = a1;
                        atomj_ = a2;
                        atomk_ = a3;
                        atoml_ = a4;
			gradi_ = {0.0, 0.0, 0.0};			
			gradj_ = {0.0, 0.0, 0.0};			
			gradk_ = {0.0, 0.0, 0.0};			
			gradl_ = {0.0, 0.0, 0.0};			

                }

		double Get_Value() override {
			return value_;
		}

		array<array3d, 4> Get_Gradient() override {
			return {gradi_,gradj_,gradk_,gradl_};
		}


                void Evaluate(const double xyz[], int natoms, array3d box_len) override {
			const double bohr_to_ang = 0.529177;
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


			rij_[0] = atomi_xyz[0] - atomj_xyz[0];
			rij_[1] = atomi_xyz[1] - atomj_xyz[1];
			rij_[2] = atomi_xyz[2] - atomj_xyz[2];

			rkj_[0] = atomk_xyz[0] - atomj_xyz[0];
			rkj_[1] = atomk_xyz[1] - atomj_xyz[1];
			rkj_[2] = atomk_xyz[2] - atomj_xyz[2];

			rkl_[0] = atomk_xyz[0] - atoml_xyz[0];
			rkl_[1] = atomk_xyz[1] - atoml_xyz[1];
			rkl_[2] = atomk_xyz[2] - atoml_xyz[2];

			rij_ = Minimum_Image(rij_, 10.0);
			rkj_ = Minimum_Image(rkj_, 10.0);
			rkl_ = Minimum_Image(rkl_, 10.0);

			m_=Crossp(rij_, rkj_);
			n_=Crossp(rkj_, rkl_);

			array3d rkj_hat_ = Normalize(rkj_);

			double X = Dotp(Crossp(m_, n_), rkj_hat_);
		        double Y = Dotp(m_, n_);	

			value_ = atan2(X, Y);
		};


                void Compute_gradient(const double xyz[], int natoms, array3d box_len) override {

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
			
		};


};
