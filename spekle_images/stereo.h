#ifndef _STEREO_H_SPECKLE_ 
#define _STEREO_H_SPECKLE_

#include <vector>
#include <iostream>
#include <iomanip>
#include "util.h"

namespace speckle {

	class instrinsic_parameter
	{
		matrix3d m;
	public:
		instrinsic_parameter(float_t fx, float_t fs, float_t cx, float_t fy, float_t cy) {
			m.data[0][0] = fx; m.data[0][1] = fs; m.data[0][2] = cx;
			m.data[1][0] = 0.f; m.data[1][1] = fy; m.data[1][2] = cy;
			m.data[2][0] = 0.f; m.data[2][1] = 0.f; m.data[2][2] = 1.f;
		}

		point2d projection(const point3d& p3) const {
			float_t lambda = p3.z;
			float_t x = m.data[0][2] * p3.z + m.data[0][1] * p3.y + m.data[0][0] * p3.x;
			float_t y = m.data[1][2] * p3.z + m.data[1][1] * p3.y;
			x /= lambda;
			y /= lambda;
			return { x,y };
		}

		point3d inverse_projection(const point2d& p2, float_t lambda) const {
			float_t Z = lambda;
			float_t Y = lambda * (p2.y - m.data[1][2]) / m.data[1][1];
			float_t X = lambda * (p2.x - m.data[0][2]) / m.data[0][0];

			if (m.data[0][1] != 0.f) {
				float_t part = m.data[0][1] / m.data[0][0] * Y;
				X -= part;
			}
			return { X,Y,Z };
		}

	}; // end class instrinsic_parameter


#include <iostream>
	using namespace std;

	void test_3d_to_2d(const instrinsic_parameter& k, const vector<point3d> p3s) {

		for (auto p3 : p3s)
		{
			point2d p2 = k.projection(p3);
			point3d p3_ = k.inverse_projection(p2, p3.z);

			cout << p3 << " > " << p2 << " > " << p3_ << endl;
		}

	}


	void test_2d_3d(const instrinsic_parameter& k, const vector<point2d> p2s, const float_t z) {
		cout << ios::fixed << setprecision(8) << endl;
		for (auto p2 : p2s)
		{
			point3d p3 = k.inverse_projection(p2, z);
			point2d p2_ = k.projection(p3);

			cout << p2 << " > " << p3 << " > " << p2_ << endl;
		}
	}




	static void test_intrinsic_parameters() {

		float_t fx, fy, fs, cx, cy;

		fx = 1024.f;
		fy = 1024.f;
		cx = 256.f;
		cy = 256.f;
		fs = 0.f;

		instrinsic_parameter k(fx, fs, cx, fy, cy);



		const float_t z = 800.f;  //Plane Z.

		vector<point3d> p3s =
		{
			{ -100,-100,z },
			{ 0, 0, z },
			{ 100,100,z }
		};

		test_3d_to_2d(k, p3s);




		vector<point2d> p2s = {
			{129,129},
			{255,257},
			{383,385},
		};

		test_2d_3d(k, p2s, z);


	}

}// end namespace spekle



#endif // !
