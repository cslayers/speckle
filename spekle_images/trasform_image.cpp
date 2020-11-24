#include <iostream>
using namespace std;



#include "image.h"
#include "builder.h"
#include "transform.h"
#include "interp.h"
using namespace speckle;


void transform_double_sin(image& src, image& dst) {
	int L = src.cols();
	int N = 2;
	int N2 = N * 10;
	float_t A = (float_t)10;

	float_t B = (float_t)10;
	float_t B2 = (float_t)1;



	std::cout << "u0 :"; std::cin >> A; cout << endl;
	std::cout << "N :"; std::cin >> N; cout << endl;
	std::cout << "um :"; std::cin >> B; cout << endl;
	std::cout << "N2:"; std::cin >> N2; cout << endl;
	std::cout << "uh:"; std::cin >> B2; cout << endl;


	float_t W = 2 * pi * N / L;
	float_t W2 = 2 * pi * N2 / L;


	string append_filename = "--A" + to_string(A);
	append_filename += " N-" + to_string(N);
	append_filename += " B-" + to_string(B);
	append_filename += " N2-" + to_string(N2);
	append_filename += " B2-" + to_string(B2);
	append_filename += ".bmp";

	dst.setpath(src.path() + append_filename);
	dst.resize(src.size());



	sinosoidal_x m(A, B, W, B2, W2);

	bspline b_interp(&src);

	for (int r = 0; r < dst.rows(); ++r)
		for (int c = 0; c < dst.cols(); ++c) {
			point2d ref_coor = m.map_inv({ (float_t)c,(float_t)r });
			dst[r][c] = (float_t)b_interp.interp(ref_coor);
		}


}


void transform_affine(image& src, image& dst) {



	string append_filename = "trans.bmp";
	dst.setpath(src.path() + append_filename);
	dst.resize(src.size());

	float u = 0, v = 0;
	float ux = -0.1f, vx = -0.1f;
	float uy = -0.1f, vy = -0.1f;
	affine m(u, ux, uy, v, vx, vy);


	bspline b_interp(&src);
	const int cx = src.cols() / 2;
	const int cy = src.rows() / 2;

	for (int r = 0; r < dst.rows(); ++r)
		for (int c = 0; c < dst.cols(); ++c) {
			point2d ref_coor = m.map_inv({ (float_t)c - cx,(float_t)r - cy });
			ref_coor.x += cx;
			ref_coor.y += cy;
			dst[r][c] = (float_t)b_interp.interp(ref_coor);
		}
}


int run_transform() {



	config c = speckle::default_config();
	c.size.cols = 512;

	image orgin;
	//image_builder::build_image(c, orgin);

	im_read("C:\\Users\\205\\Desktop\\work\\gen\\apple.tif", orgin);


	image transformed;

	transform_affine(orgin, transformed);

	im_write(transformed, transformed.path());


	return 0;
}