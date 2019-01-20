// https://habr.com/ru/post/436790/
// https://github.com/ssloy/tinyraytracer


#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// Шаблонный класс вектора
template <size_t DIM, typename T> struct vec {
    vec() { 
		for (size_t i=DIM; i--; _data[i] = T()); 
	}
    T& operator[](const size_t i) { 
		assert(i<DIM); 
		return _data[i]; 
	}
    const T& operator[](const size_t i) const { 
		assert(i<DIM); 
		return _data[i]; 
	}
private:
    T _data[DIM];
};

typedef vec<2, float> Vec2f;
typedef vec<3, float> Vec3f;
typedef vec<3, int  > Vec3i;
typedef vec<4, float> Vec4f;

// Шаблонный класс вектора 2
template <typename T> struct vec<2,T> {
    vec(): 
		_x(T()), 
		_y(T()) {
	}
    vec(T X, T Y): 
		_x(X), 
		_y(Y) {		
	}
    template <class U> vec<2,T>(const vec<2,U> &v);
    T& operator[](const size_t i) { 
		assert(i<2); 
		return i<=0 ? _x : _y; 
	}
    const T& operator[](const size_t i) const { 
		assert(i<2); 
		return i<=0 ? _x : _y; 
	}
    T _x,
	T _y;
};

// Шаблонный класс вектора 3
template <typename T> struct vec<3,T> {
    vec(): 
		_x(T()), 
		_y(T()), 
		_z(T()) {
	}
    vec(T X, T Y, T Z): 
		_x(X), 
		_y(Y), 
		_z(Z) {
	}
    T& operator[](const size_t i){ 
		assert(i<3); 
		return i<=0 ? _x : (1==i ? _y : _z); 
	}
    const T& operator[](const size_t i) const { 
		assert(i<3); return i<=0 ? _x : (1==i ? _y : _z); 
	}
    float norm() { 
		return std::sqrt(_x*_x+_y*_y+_z*_z); 
	}
    vec<3,T>& normalize(T l=1) { 
		*this = (*this)*(l/norm()); 
		return *this; 
	}
    T _x;
	T _y;
	T _z;
};

// Шаблонный класс вектора 4
template <typename T> struct vec<4,T> {
    vec(): 
		_x(T()), 
		_y(T()), 
		_z(T()), 
		_w(T()) {
	}
    vec(T X, T Y, T Z, T W): 
		_x(X), 
		_y(Y), 
		_z(Z), 
		_w(W) {
	}
    T& operator[](const size_t i){ 
		assert(i<4); 
		return i<=0 ? _x : (1==i ? _y : (2==i ? _z : _w)); 
	}
    const T& operator[](const size_t i) const { 
		assert(i<4); 
		return i<=0 ? _x : (1==i ? _y : (2==i ? _z : _w)); 
	}
    T _x;
	T _y;
	T _z;
	T _w;
};

// Оператор умножения векторов
template<size_t DIM,typename T> T operator*(const vec<DIM,T>& lhs, const vec<DIM,T>& rhs) {
    T ret = T();
    for (size_t i=DIM; i--; ret+=lhs[i]*rhs[i]);
    return ret;
}

// Оператор сложения векторов
template<size_t DIM,typename T>vec<DIM,T> operator+(vec<DIM,T> lhs, const vec<DIM,T>& rhs) {
    for (size_t i=DIM; i--; lhs[i]+=rhs[i]);
    return lhs;
}

// Оператор вычитания векторов
template<size_t DIM,typename T>vec<DIM,T> operator-(vec<DIM,T> lhs, const vec<DIM,T>& rhs) {
    for (size_t i=DIM; i--; lhs[i]-=rhs[i]);
    return lhs;
}

// Оператор умножения векторов
template<size_t DIM,typename T,typename U> vec<DIM,T> operator*(const vec<DIM,T> &lhs, const U& rhs) {
    vec<DIM,T> ret;
    for (size_t i=DIM; i--; ret[i]=lhs[i]*rhs);
    return ret;
}

// Оператор вычитания векторов
template<size_t DIM,typename T> vec<DIM,T> operator-(const vec<DIM,T> &lhs) {
    return lhs*T(-1);
}

// Векторное произведение векторов
template <typename T> vec<3,T> cross(vec<3,T> v1, vec<3,T> v2) {
    return vec<3,T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

template <size_t DIM, typename T> std::ostream& operator<<(std::ostream& out, const vec<DIM,T>& v) {
    for(unsigned int i=0; i<DIM; i++) {
        out << v[i] << " " ;
    }
    return out ;
}

#endif //__GEOMETRY_H__

