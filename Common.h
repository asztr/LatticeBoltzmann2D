
#ifndef COMMON_H
#define COMMON_H

#include<cmath>

enum Direction {
	East = 1, North, West, South, 
	NorthEast, NorthWest, SouthWest, SouthEast
};

template <typename T> struct vec {
	T x;
	T y;

	vec() { }
	vec(T _x, T _y) {
		x = _x;
		y = _y;
	}

	inline void define(T _x, T _y) { x = _x; y = _y; }
	inline double norm() const { return sqrt(x*x + y*y); }
	inline double norm2() const { return x*x + y*y; }
	inline double angle() const { return atan2(y, x); }
	inline T operator*(const vec<T>& vec2) const { return x*vec2.x + y*vec2.y; }

	inline vec<T> operator*(const double r) const {
		vec<T> ret = {r*x, r*y};
		return ret;
	}
	
	inline vec<T> operator%(const vec<T>& vec2) const {
		vec<T> ret = {(T)x/vec2.x, (T)y/vec2.y};
		return ret;
	}

	inline vec<T> operator+(const vec<T>& vec2) const {
		vec<T> ret = {x+vec2.x, y+vec2.y};
		return ret;
	}
};

#endif