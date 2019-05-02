#include <array>
#include <cmath>

using namespace std;

class CollectiveVariable
{
	protected:
		std::array<double, 2> bounds_;
		// std::array<double, 3> grad_;
		double value_;
	public:

		CollectiveVariable(): value_(0.0), bounds_{{0,0}} {}
		virtual double Evaluate(const double xyz[], int) = 0;

};


class Dihedral: public CollectiveVariable
{
        private:
                int atom1_, atom2_, atom3_, atom4_;


        public:

                Dihedral(int a1, int a2, int a3, int a4) {
                        atom1_ = a1;
                        atom2_ = a2;
                        atom3_ = a3;
                        atom4_ = a4;
                }

                double Evaluate(const double xyz[], int natoms) override {
			const double bohr_to_ang = 0.529177;
			std::array<double, 3> atom1_xyz, atom2_xyz, atom3_xyz, atom4_xyz;
			atom1_xyz[0]=*(xyz+3*atom1_)*bohr_to_ang;
			atom1_xyz[1]=*(xyz+3*atom1_+1)*bohr_to_ang;
			atom1_xyz[2]=*(xyz+3*atom1_+2)*bohr_to_ang;

			atom2_xyz[0]=*(xyz+3*atom2_)*bohr_to_ang;
			atom2_xyz[1]=*(xyz+3*atom2_+1)*bohr_to_ang;
			atom2_xyz[2]=*(xyz+3*atom2_+2)*bohr_to_ang;

			atom3_xyz[0]=*(xyz+3*atom3_)*bohr_to_ang;
			atom3_xyz[1]=*(xyz+3*atom3_+1)*bohr_to_ang;
			atom3_xyz[2]=*(xyz+3*atom3_+2)*bohr_to_ang;
			
			atom4_xyz[0]=*(xyz+3*atom4_)*bohr_to_ang;
			atom4_xyz[1]=*(xyz+3*atom4_+1)*bohr_to_ang;
			atom4_xyz[2]=*(xyz+3*atom4_+2)*bohr_to_ang;

			double rx12, rx32, rx34;
			double ry12, ry32, ry34;
			double rz12, rz32, rz34;

			rx12 = atom1_xyz[0] - atom2_xyz[0];
			ry12 = atom1_xyz[1] - atom2_xyz[1];
			rz12 = atom1_xyz[2] - atom2_xyz[2];

			rx32 = atom3_xyz[0] - atom2_xyz[0];
			ry32 = atom3_xyz[1] - atom2_xyz[1];
			rz32 = atom3_xyz[2] - atom2_xyz[2];

			rx34 = atom3_xyz[0] - atom4_xyz[0];
			ry34 = atom3_xyz[1] - atom4_xyz[1];
			rz34 = atom3_xyz[2] - atom4_xyz[2];

//			cout << rx12 << endl;
//			cout << ry12 << endl;
//			cout << rz12 << endl;
//
//			cout << rx32 << endl;
//			cout << ry32 << endl;
//			cout << rz32 << endl;
//
//			cout << rx34 << endl;
//			cout << ry34 << endl;
//			cout << rz34 << endl;
			double mx, my, mz; 

			mx =  ry12*rz32 - ry32*rz12;
			my = -rx12*rz32 + rz12*rx32;
			mz =  rx12*ry32 - ry12*rx32;

			double nx, ny, nz;

			nx =  ry32*rz34 - rz32*ry34;
			ny = -rx32*rz34 + rz32*rx34;
			nz =  rx32*ry34 - ry32*rx34;

			//cout << mx << endl;
			//cout << my << endl;
			//cout << mz << endl;
			//cout << nx << endl;
			//cout << ny << endl;
			//cout << nz << endl;

			double msq, abs_m, nsq, abs_n, mdn, r12dn, cosphi;

			msq = mx*mx + my*my + mz*mz;
			abs_m = std::sqrt(msq);
			nsq = nx*nx + ny*ny + nz*nz;
			abs_n = std::sqrt(nsq);
		
			// cout << abs_m << endl;
			// cout << abs_n << endl;

			mdn = mx*nx + my*ny + mz*nz;

			// cout << mdn << endl;

			cosphi = mdn/(abs_m * abs_n);

			// cout << cosphi << endl;

			cosphi = std::min(cosphi, 1.0);
			cosphi = std::max(cosphi, -1.0);
			r12dn = rx12*nx + ry12*ny + rz12*nz;

			//bool sign = r12dn > 0 ? true : false;

			double phi = std::acos(cosphi);
			//phi = sign > true ? -phi : phi;
			return phi;
		};
};
