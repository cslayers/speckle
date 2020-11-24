#ifndef _INTERP_H_SPECKLE_
#define _INTERP_H_SPECKLE_

#include "image.h"
#include "util.h"

namespace speckle {


	class image_interpolator {
	public:
		virtual float_t interp(const point2d& subpixel) = 0;
		virtual ~image_interpolator() {}
	};



	class bilinear :public image_interpolator {
	public:
		float_t interp(const point2d& subpixel) {

			float_t x = subpixel.x;
			float_t y = subpixel.y;

			if (x < 0 || x >= _pimg->cols() - 1 || y < 0 || y >= _pimg->rows() - 1)
				return (float_t)255;

			int ix = (int)floor(x);
			int iy = (int)floor(y);

			float_t dx = x - ix;
			float_t dy = y - iy;


			auto v00 = (*_pimg)(ix, iy);
			auto v01 = (*_pimg)(ix, iy + 1);
			auto v10 = (*_pimg)(ix + 1, iy);
			auto v11 = (*_pimg)(ix + 1, iy + 1);

			float_t res = (float_t)0.f;
			res += v00 * (1 - dx) * (1 - dy);
			res += v01 * dx * (1 - dy);
			res += v10 * (1 - dx) * dy;
			res += v11 * dx * dy;
			return res;
		}

		explicit bilinear(image* pimg) :_pimg(pimg) {

		}
	private:
		image* _pimg;
	};


	class bspline :public image_interpolator {
	public:
		explicit bspline(image* pimg) :_pimg(pimg) {
			prepare();
		}

		virtual float_t interp(const point2d& subpixel) {
			int ir = (int)subpixel.y;
			int ic = (int)subpixel.x;

			if (ir <= 0 || ic <= 0) return 255;
			if (ir > _pimg->rows() - 3 || ic > _pimg->cols() - 3) return 255;


			float_t x1 = subpixel.x - ic;
			float_t y1 = subpixel.y - ir;

			float_t x2 = x1 * x1;
			float_t y2 = y1 * y1;

			float_t x3 = x2 * x1;
			float_t y3 = y2 * y1;


			const vector<vector<float_t>>& alpha = alphas[ir][ic];

			float_t res = 0.f;
			res += alpha[0][0];
			res += alpha[0][1] * x1;
			res += alpha[0][2] * x2;
			res += alpha[0][3] * x3;

			res += alpha[1][0] * y1;
			res += alpha[1][1] * y1 * x1;
			res += alpha[1][2] * y1 * x2;
			res += alpha[1][3] * y1 * x3;

			res += alpha[2][0] * y2;
			res += alpha[2][1] * y2 * x1;
			res += alpha[2][2] * y2 * x2;
			res += alpha[2][3] * y2 * x3;

			res += alpha[3][0] * y3;
			res += alpha[3][1] * y3 * x1;
			res += alpha[3][2] * y3 * x2;
			res += alpha[3][3] * y3 * x3;

			if (res > (float_t)image::intensity_max)
				res = (float_t)image::intensity_max;
			if (res < 0) res = 0;

			return res;
		}

		virtual ~bspline() {}


	private:
		image* _pimg;

	private:
		vector<vector<vector<vector<float_t >>>> alphas;

	private:
		const float controlMatrix[4][4] = {
			{ 71.0f / 56.0f, -19.0f / 56.0f, 5 / 56.0f, -1.0f / 56.0f },
			{ -19.0f / 56.0f, 95.0f / 56.0f, -25 / 56.0f, 5.0f / 56.0f },
			{ 5.0f / 56.0f, -25.0f / 56.0f, 95 / 56.0f, -19.0f / 56.0f },
			{ -1.0f / 56.0f, 5.0f / 56.0f, -19 / 56.0f, 71.0f / 56.0f }
		};

		const float functionMatrix[4][4] = {
			{ -1.0f / 6.0f, 3.0f / 6.0f, -3.0f / 6.0f, 1.0f / 6.0f },
			{ 3.0f / 6.0f, -6.0f / 6.0f, 3.0f / 6.0f, 0.0f },
			{ -3.0f / 6.0f, 0.0f, 3.0f / 6.0f, 0.0f },
			{ 1.0f / 6.0f, 4.0f / 6.0f, 1.0f / 6.0f, 0.0f }
		};

	private:
		void prepare() {

			auto rows = _pimg->rows();
			auto cols = _pimg->cols();

			vector<float_t> e(4);
			vector<vector<float_t>> ele(4, e);
			vector<vector<vector<float_t>>> row_data(cols, ele);
			vector<vector<vector<vector<float_t >>>> temp(rows, row_data);
			alphas.swap(temp);

			for (int r = 1; r < rows - 2; r++) {
				for (int c = 1; c < cols - 2; c++) {

					float omiga[4][4];
					float beta[4][4];

					for (int i = 0; i < 4; i++)
						for (int j = 0; j < 4; j++)
							omiga[i][j] = (*_pimg)[r - 1 + i][c - 1 + j];



					for (int k = 0; k < 4; k++)
					{
						for (int l = 0; l < 4; l++)
						{
							beta[k][l] = 0;
							for (int m = 0; m < 4; m++)
							{
								for (int n = 0; n < 4; n++)
								{
									beta[k][l] += controlMatrix[k][m] * controlMatrix[l][n] * omiga[n][m];
								}
							}
						}
					}





					for (int k = 0; k < 4; k++)
					{
						for (int l = 0; l < 4; l++)
						{
							alphas[r][c][k][l] = 0;
							for (int m = 0; m < 4; m++)
							{
								for (int n = 0; n < 4; n++)
								{
									alphas[r][c][k][l] += functionMatrix[k][m] * functionMatrix[l][n] * beta[n][m];
								}
							}
						}
					}




					for (int k = 0; k < 2; k++)
					{
						for (int l = 0; l < 4; l++)
						{
							float temp = alphas[r][c][k][l];
							alphas[r][c][k][l] = alphas[r][c][3 - k][3 - l];
							alphas[r][c][3 - k][3 - l] = temp;
						}
					}



				}
			}

		}
	};

}

#endif