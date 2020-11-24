#pragma once
#ifndef _UTIL_H_SPEKLE_
#define _UTIL_H_SPEKLE_
#include <iostream>
#include <iomanip>

namespace speckle {
	const float_t pi = float(3.141592653f);

	typedef float float_t;

	class point2d {
	public:
		float_t x, y;
	};


	class point3d {
	public:
		float_t x, y, z;
	};


	class matrix3d {
	public:
		float_t data[3][3];
	};


	inline point2d operator-(const point2d& p) {
		return { -p.x,-p.y };
	}

	inline point2d operator-(const point2d& lhs, const point2d& rhs) {
		return { lhs.x - rhs.x, lhs.y - rhs.y };
	}

	inline point2d operator+(const point2d& lhs, const point2d& rhs) {
		return { lhs.x + rhs.x, lhs.y + rhs.y };
	}

	inline bool operator==(const point2d& lhs, const point2d& rhs) {
		return lhs.x == rhs.x && lhs.y == rhs.y;
	}

	inline bool operator!=(const point2d& lhs, const point2d& rhs) {
		return !(lhs == rhs);
	}



	inline point3d operator-(const point3d& p) {
		return { -p.x, -p.y, -p.z };
	}

	inline point3d operator-(const point3d& lhs, const point3d& rhs) {
		return { lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
	}

	inline point3d operator+(const point3d& lhs, const point3d& rhs) {
		return { lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
	}

	inline bool operator==(const point3d& lhs, const point3d& rhs) {
		return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
	}

	inline bool operator!=(const point3d& lhs, const point3d& rhs) {
		return !(lhs == rhs);
	}


	inline std::ostream& operator<<(std::ostream& out, const point2d & p) {
		out << p.x << "," << p.y;
		return out;
	}

	inline std::ostream& operator<<(std::ostream& out, const point3d & p) {
		out << p.x << "," << p.y << "," << p.z;
		return out;
	}


}




#endif // !_UTIL_H_SPEKLE_
