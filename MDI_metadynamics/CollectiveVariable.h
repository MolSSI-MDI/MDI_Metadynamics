#include <array>

using namespace std;

class CollectiveVariable
{
	protected:
		std::array<double, 2> bounds_;
		// std::array<double, 3> grad_;
		double value_;
	public:

		CollectiveVariable(): value_(0.0), bounds_{{0,0}} {}

		virtual void Evaluate() = 0;

		double GetValue() const
		{
			return value_;
		}
};
