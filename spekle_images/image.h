#ifndef _IMAGE_H_SPECKLE_
#define _IMAGE_H_SPECKLE_

#include <vector>
#include <iostream>
#include <opencv2\core.hpp>
#include <opencv2\imgproc.hpp>
#include <opencv2\highgui.hpp>
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
		std::string _path;

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

		image_size size() const {
			image_size res;
			res.cols = cols();
			res.rows = rows();
			return res;
		}

		std::string path() {
			return _path;
		}

		void setpath(std::string path) {
			_path = path;
		}

		void resize(image_size s) {
			_data.resize(s.rows);
			for (auto& line : _data) line.resize(s.cols);
		}

		vector<intensity_t>& operator[](size_t row) {
			return _data[row];
		}

		float_t operator()(size_t x, size_t y) {
			return _data[y][x];
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
			_path.swap(another._path);
		}
	};



	inline int im_write(image& img, std::string path) {

		cv::Mat img_save((int)img.size().rows, (int)img.size().cols, CV_8UC1);

		for (int i = 0; i < img_save.rows; ++i) {
			uchar * p_row = img_save.ptr<uchar>(i);
			for (int j = 0; j < img_save.cols; ++j) {
				p_row[j] = img[i][j];
			}
		}
		if (!cv::imwrite(path, img_save)) {
			std::cerr << "[warnning]: fail to save image" << path << std::endl;
			return 1;
		}
		std::cout << "[IO] save image to: " << path << std::endl;
		return 0;
	}



	inline int im_read(std::string path, image& img) {
		cv::Mat in = cv::imread(path, cv::IMREAD_GRAYSCALE);
		if (!in.data) return 1;

		img.resize({ (size_t)in.rows, (size_t)in.cols });
		img.setpath(path);

		for (int i = 0; i < in.rows; ++i) {
			uchar * p_row = in.ptr<uchar>(i);
			for (int j = 0; j < in.cols; ++j) {
				img[i][j] = p_row[j];
			}
		}
		return 0;
	}




	inline void test_image_class() {
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

#endif