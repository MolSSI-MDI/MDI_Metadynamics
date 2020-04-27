#include "CollectiveVariable.h"

class Distance: public CollectiveVariable
{
        private:
                int atomi_, atomj_;
				array3d gradi_, gradj_;

        public:

                Distance(int a1, int a2);
                array2dint Get_Atoms() override;
				double Get_Value() override;
				std::array<array3d, 2> Get_Gradient() override;
                void Evaluate(const double xyz[], int natoms, array3d box_len) override;
};
