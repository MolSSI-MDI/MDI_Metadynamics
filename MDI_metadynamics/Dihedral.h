#include "CollectiveVariable.h"

using namespace std;


class Dihedral: public CollectiveVariable
{
	private:
		int atom1, atom2, atom3, atom4;
	public:

		Dihedral(int a1, int a2, int a3, int a4) {
			a1 = a1;
			a2 = a2;
			a3 = a3;
			a4 = a4;
		}

		void Evaluate() override {
			value_ = 1.0;
		}
};

