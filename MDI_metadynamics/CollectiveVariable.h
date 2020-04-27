#pragma once

#include "Utils.h"

class CollectiveVariable
{
	protected:
		double value_;
		int natoms_;

	public:

		CollectiveVariable(): value_(0.0) {}
		virtual void Evaluate(const double xyz[], int, array3d) = 0;
		virtual double Get_Value() = 0;
		virtual std::array<array3d, 2> Get_Gradient() = 0;
		virtual array2dint Get_Atoms() = 0;
		~CollectiveVariable();

};
