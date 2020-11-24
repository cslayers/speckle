#include <iostream>
#include <Eigen/Dense>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = NX,
		ValuesAtCompileTime = NY
	};
	typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

	int m_inputs, m_values;

	Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
	Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

	int inputs() const { return m_inputs; }
	int values() const { return m_values; }

};

struct my_functor : Functor<double>
{
	my_functor(void) : Functor<double>(1, 1) {}
	int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
	{
		auto x0 = x(0);

		fvec(0) = x0 + sin(x0) + sin(2 * x0) - 2.57;

		return 0;
	}
};


int eigen_examplex()
{
	Eigen::VectorXd x(1);
	x(0) = 2;
	std::cout << "x: " << x << std::endl;

	my_functor functor;
	Eigen::NumericalDiff<my_functor> numDiff(functor);
	Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor>, double> lm(numDiff);

	lm.parameters.maxfev = 2000;
	lm.parameters.xtol = 1.0e-10;
	int ret = lm.minimize(x);

	std::cout << "iter=";
	std::cout << lm.iter << std::endl;
	std::cout << "return=" << ret << std::endl;


	std::cout << "x that minimizes the function: " << x << std::endl;
	std::cout << "press [ENTER] to continue " << std::endl;
	std::cin.get();
	return 0;
}