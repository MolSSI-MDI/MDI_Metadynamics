#include <array>
#include <cmath>
#include "Utils.h"

using namespace std;

class CollectiveVariable
{
	protected:
		array2d bounds_;
		// array3d grad_;
		double value_;
	public:

		CollectiveVariable(): value_(0.0), bounds_{{0,0}} {}
		virtual double Evaluate(const double xyz[], int) = 0;

};


class Dihedral: public CollectiveVariable
{
        private:
                int atomi_, atomj_, atomk_, atoml_;


        public:

                Dihedral(int a1, int a2, int a3, int a4) {
                        atomi_ = a1;
                        atomj_ = a2;
                        atomk_ = a3;
                        atoml_ = a4;
                }

                double Evaluate(const double xyz[], int natoms) override {
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

			std:array<double, 3> vij, vjk, vkl; 

			vij[0] = atomj_xyz[0] - atomi_xyz[0];
			vij[1] = atomj_xyz[1] - atomi_xyz[1];
			vij[2] = atomj_xyz[2] - atomi_xyz[2];

			vjk[0] = atomk_xyz[0] - atomj_xyz[0];
			vjk[1] = atomk_xyz[1] - atomj_xyz[1];
			vjk[2] = atomk_xyz[2] - atomj_xyz[2];

			vkl[0] = atoml_xyz[0] - atomk_xyz[0];
			vkl[1] = atoml_xyz[1] - atomk_xyz[1];
			vkl[2] = atoml_xyz[2] - atomk_xyz[2];

			vij = Minimum_Image(vij, 10.0);
			vjk = Minimum_Image(vjk, 10.0);
			vkl = Minimum_Image(vkl, 10.0);

			array3d A=crossp(vij, vjk);
			array3d B=crossp(vjk, vkl);

			array3d vjk_norm = Normalize(vjk);

			double X = dotp(crossp(A, B), vjk_norm);
		        double Y = dotp(A, B);	

			double phi = atan2(X, Y);
			return phi; 
		};
};
