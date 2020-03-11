#pragma once
#include <vector>
#include <iostream>

namespace speckle {
	using std::vector;

	class image {
	public:
		struct image_size {
			size_t rows;
			size_t cols;
		};
		typedef int size_t;

		typedef unsigned char intensity_t;
		static const intensity_t intensity_max = 255;

	private:
		vector<vector<intensity_t>> _data;

	public:
		explicit image(image_size size = { 1,1 }) {
			if (size.rows <= 0) size.rows = 1;
			if (size.cols <= 0) size.cols = 1;
			_data.resize(size.rows);
			for (auto& line : _data) line.resize(size.cols);
		}

		size_t rows() const {
			return static_cast<size_t>(_data.size());
		}

		size_t cols() const {
			return static_cast<size_t>(_data[0].size());
		}

		vector<intensity_t>& operator[](size_t row) {
			return _data[row];
		}

		void print(std::ostream& out = std::cout) const {
			int r = rows();
			int c = cols();
			out << std::endl;
			for (int i = 0; i < r; ++i) {
				for (int j = 0; j < c; ++j)	out << static_cast<unsigned int>(_data[i][j]) << " ";
				out << std::endl;
			}
			out << std::endl;
		}

	public:
		void swap(image& another) {
			_data.swap(another._data);
		}
	};


	inline static void test_image_class() {
		image::image_size size{ 3,3 };

		image timg(size);
		timg.print();

		timg[0][0] = unsigned char(1);
		timg[0][1] = unsigned char(2);
		timg[0][2] = unsigned char(-1);
		timg[1][0] = unsigned char(-2);
		timg[1][1] = unsigned char(256);
		timg[1][2] = unsigned char(257);
		timg[2][1] = unsigned char(1);
		timg[2][2] = unsigned char(2);


		image another;

		timg.print();
		another.print();

		another.swap(timg);

		timg.print();
		another.print();
	}

}