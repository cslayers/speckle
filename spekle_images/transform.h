#ifndef _TRANSFORM_H_SPECKLE_
#define _TRANSFORM_H_SPECKLE_
#include "util.h"
#include <unsupported\Eigen\NumericalDiff>
#include <unsupported\Eigen\NonLinearOptimization>


namespace speckle {

	inline void inverse_matrix(float input[3][3], float result[3][3]) {
		int i, col, row;
		int temp;

		for (row = 0; row < 3; row++)
		{
			for (col = 0; col < 3; col++)
			{
				result[row][col] = col == row ? 1.0f : 0.0f;
			}
		}

		for (col = 0; col < 3; col++)
		{
			// Find pivot (maximum lth column element) in the rest (3-l) rows
			temp = col;
			for (row = col + 1; row < 3; row++)
			{
				if (input[row][col] > input[temp][col])
				{
					temp = row;
				}
			}
			if (std::abs(input[temp][col]) == 0.0f)
			{
				return;
			}
			// Swap the row which has maximum lth column element
			if (temp != col)
			{
				for (i = 0; i < 3; i++)
				{
					float buf = input[col][i];
					input[col][i] = input[temp][i];
					input[temp][i] = buf;

					buf = result[col][i];
					result[col][i] = result[temp][i];
					result[temp][i] = buf;
				}
			}
			// Perform row operation to form required identity matrix out of the Hessian matrix
			for (row = 0; row < 3; row++)
			{
				float curRowEle = input[row][col];
				if (row != col)
				{
					for (i = 0; i < 3; i++)
					{
						result[row][i] -= result[col][i] * curRowEle / input[col][col];
						input[row][i] -= input[col][i] * curRowEle / input[col][col];
					}
				}
				else
				{
					for (i = 0; i < 3; i++)
					{
						result[row][i] /= curRowEle;
						input[row][i] /= curRowEle;
					}
				}
			}
		}
	}
	//!==================================================================================
	// Eigen non-linear optimize util
	//!==================================================================================
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
	//!==================================================================================




	struct sin_functor : Functor<double> //usage
	{
		sin_functor(float_t xp, float_t a, float_t b, float_t w, float_t b2 = 0, float_t w2 = 0) :
			_xp(xp), _a(a), _b(b), _w(w), _b2(b2), _w2(w2),
			Functor<double>(1, 1) {}


		int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
		{
			auto x0 = x(0);

			fvec(0) = x0 + _a + _b * sin(_w * x0) + _b2 * sin(_w2 * x0) - _xp;

			return 0;
		}


		float_t _xp;
		float_t _a, _b, _w;    // xp = x + a + b sin (w * x)  + ...
		float_t _b2, _w2;

	};




	class mapbase {
	public:
		mapbase() {}
		virtual point2d map(const point2d& coor) const = 0;
		virtual point2d map_inv(const point2d& coor) const = 0;
	};

	class affine : public mapbase {
	public:
		virtual point2d map(const point2d& coor) const {


			float_t x = coor.x + u + ux * coor.x + uy * coor.y;
			float_t y = coor.y + v + vx * coor.x + vy * coor.y;

			return point2d{ x,y };
		}


		virtual point2d map_inv(const point2d& coor)const {

			float_t x = invm[0][0] * coor.x + invm[0][1] * coor.y + invm[0][2];
			float_t y = invm[1][0] * coor.x + invm[1][1] * coor.y + invm[1][2];
			return point2d{ x,y };
		}

		affine(float u, float ux, float uy, float v, float vx, float vy) {
			u = u;
			ux = ux;
			uy = uy;
			v = v;
			vx = vx;
			vy = vy;

			float m[3][3] = { 0.f };
			m[0][0] = 1 + ux;
			m[0][1] = uy;
			m[0][2] = u;
			m[1][0] = vx;
			m[1][1] = 1 + vy;
			m[1][2] = v;
			m[2][0] = 0.f;
			m[2][1] = 0.f;
			m[2][2] = 1.f;

			inverse_matrix(m, invm);
		}
	private:
		float invm[3][3];
		float u, ux, uy, v, vx, vy;
	};

	class sinosoidal_x : public mapbase {
	public:
		virtual point2d map(const point2d& coor) const {
			float_t x_prime = (*this)(coor.x);
			return { x_prime, coor.y };
		}


		virtual point2d map_inv(const point2d& coor)const {
			float_t _xp = coor.x;
			sin_functor functor(_xp, _a, _b, _w, _b2, _w2);

			Eigen::NumericalDiff<sin_functor> numDiff(functor);
			Eigen::LevenbergMarquardt<Eigen::NumericalDiff<sin_functor>, double> lm(numDiff);

			Eigen::VectorXd x(1);
			x(0) = _xp - _a;

			int ret = lm.minimize(x);
			//std::cout << "r: " << ret << std::endl;

			return { (float_t)x(0),coor.y };
		}


	public:
		float_t operator()(float_t x) const {
			return x + _a + _b * sin(_w * x) + _b2 * sin(_w2 * x);
		}

	public:
		sinosoidal_x(float_t a, float_t b, float_t w, float_t b2 = 0, float_t w2 = 0) :
			_a(a), _b(b), _w(w), _b2(b2), _w2(w2) {}

	private:
		float_t _a, _b, _w; /*x' = x + a + b sin (w * x)*/
		float_t _b2, _w2;


	};



	inline void test_transform() {
		/*
				sinosoidal_x m(10, 10, 2 * 2 * pi / 200);

				point2d p = { 207.7843f,100 };
				point2d mp = m.map(p);
				point2d pp = m.map_inv(mp);
				point2d ip = m.map_inv(p);

				cout << p << endl;
				cout << mp << endl;
				cout << pp - p << endl;
				cout << ip << endl;
		*/
		int L = 512;
		int N = 2;
		int N2 = N * 10;
		float_t A = (float_t)10;

		float_t B = (float_t)10;
		float_t B2 = (float_t)1;

		float_t W = 2 * pi * N / L;
		float_t W2 = 2 * pi * N2 / L;


		std::cout << "A " << std::endl; std::cin >> A;
		std::cout << "N " << std::endl; std::cin >> N;
		std::cout << "B " << std::endl; std::cin >> B;
		std::cout << "N2" << std::endl; std::cin >> N2;
		std::cout << "B2" << std::endl; std::cin >> B2;


		sinosoidal_x m(A, B, W, B2, W2);

		point2d p = { 101,101 };

		auto mp = m.map(p);
		auto pp = m.map_inv(mp);

		std::cout << p << std::endl;
		std::cout << mp << std::endl;
		std::cout << pp << std::endl;
		std::cout << p - pp << std::endl;
		std::cout << std::endl;

		auto ip = m.map_inv(p);
		auto mip = m.map(ip);

		std::cout << p << std::endl;
		std::cout << ip << std::endl;
		std::cout << mip << std::endl;
		std::cout << mip - p << std::endl;

	}

}
#endif
