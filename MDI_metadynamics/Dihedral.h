#include "CollectiveVariable.h"

class Dihedral: public CollectiveVariable
{
        private:
                int atomi_, atomj_, atomk_, atoml_;
		array3d gradi_, gradj_, gradk_, gradl_;

        public:

                Dihedral(int a1, int a2, int a3, int a4);
                array4dint Get_Atoms() override;
		double Get_Value() override;
		std::array<array3d, 4> Get_Gradient() override;
                void Evaluate(const double xyz[], int natoms, array3d box_len) override;
};
