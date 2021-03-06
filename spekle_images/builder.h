#include "image.h"

#include <string>
#include <vector>
#include <algorithm>
#include <utility>



namespace speckle {
	using std::string;
	using std::vector;

	typedef image::image_size image_size;
	typedef image::intensity_t intensity_t;


	struct config {
		image_size size;
		unsigned int number;		/*散斑数量*/
		float radius;				/*散斑半径*/
		float factor_lower;				/*散斑中心灰度/最高灰度*/
		float factor_upper;				/*散斑中心灰度/最高灰度*/
		string filepath;
	};


	inline config default_config() {
		config res;
		res.size = { 512,512 };
		res.radius = 3.1f;
		res.number = 1000;
		res.factor_lower = 0.4f;
		res.factor_upper = 0.8f;

		res.filepath = "C:\\Users\\205\\Desktop\\work\\gen\\";
		return res;
	}



	class image_builder {
		struct pos
		{
			unsigned int r;
			unsigned int c;
		};

	private:
		static vector<pos> _random_positions(int n, const image& img) {
			vector<pos> res;
			res.reserve(n);
			while (res.size() < n) {
				unsigned int r = rand() % img.rows();
				unsigned int c = rand() % img.cols();
				pos p{ r,c };
				res.emplace_back(p);
			}
			return res;
		}

		static vector <intensity_t> _random_intensities(int n, float lower, float upper) {

			vector<intensity_t> res;
			while (res.size() < n) {
				int base = int(lower * image::intensity_max);
				int width = int((upper - lower) * image::intensity_max);
				int val = base + rand() % width;
				res.push_back(val);
			}
			return res;
		}



	public:
		static int build_image(const config& conf, image& res) {
			// 4.1. 2D-DIC for in-plane deformation measurement
			// Zhou Y, Sun C, Chen J. 
			// Adaptive subset offset for systematic error reduction in incremental digital image correlation[J]. 
			// Optics and Lasers in Engineering, 2014, 55: 5-11.

			image temp(conf.size);


			//获得N个位置和灰度强度
			unsigned int N = conf.number;
			vector<pos> positions = _random_positions(N, temp);
			vector<intensity_t> intensities =
				_random_intensities(N, conf.factor_lower, conf.factor_upper);

			float R = conf.radius;

#pragma omp parallel for
			for (int r = 0; r < temp.rows(); ++r) {
				for (int c = 0; c < temp.cols(); ++c) {

					float val = 0.f;

					for (unsigned int j = 0; j < N; ++j) {
						const pos& p = positions[j];
						int center_intensity = static_cast<int>(intensities[j]);
						float coef = (powf((float)p.r - r, 2) + powf((float)p.c - c, 2)) / (2 * R * R);
						coef = exp(-coef);
						val += center_intensity * coef;
					}

					val = image::intensity_max - val;  //取反

					if (val < 0) val = 0;
					if (val > image::intensity_max) val = image::intensity_max;

					temp[r][c] = (intensity_t)val;
				}// col
			}// row


			if (conf.filepath.length() > 0) {

				string path = string(conf.filepath + std::to_string(rand()) + ".bmp");
				temp.setpath(path);
				if (im_write(temp, path)) {
					std::cerr << "[warnning]: fail to save image" << std::endl;
					return 1;
				}
			}

			else return 1;

			res.swap(temp);

			return 0;
		}


	};


}



