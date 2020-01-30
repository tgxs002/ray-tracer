/* 
	Name: 	elements.h
	Author: Blakey Wu
	Date:	10/5/2019
	Usage:
		Elements is a group of classes that defines points, lines, and polygons, etc.
		This file will only cover the mathmatical properties of them, as well as some transform functions.
		Drawing will be not implemented in this file.
	Purpose:
		It is designed to finish the homework, so clarity is the first proior, not efficiency.
*/

#ifndef ELEMENTS
#define ELEMENTS

#include <vector>
#include <list>
#include <iostream>
#include <assert.h>
#include <cmath>

// base class for all kinds of mats needed
class mat{

protected:
	int m;
	int n;
	int size;
	float* buffer;

public:
	mat(int m, int n);
	virtual ~mat();
	mat(const mat& other);
	mat operator * (const mat& other)const;
	mat operator * (float scale) const;
	mat operator / (float scale) const;
	mat operator + (const mat& other)const;
	mat operator - (const mat& other)const;
	bool operator < (const mat& other)const;
	bool operator > (const mat& other)const;
	bool operator <= (const mat& other)const;
	bool operator >= (const mat& other)const;
	void operator += (const mat& other);
	void operator -= (const mat& other);
	void operator *= (float scale);
	void operator /= (float scale);
	void operator = (const mat& other);
	float get(int index)const;
	float &operator[] (int); 
	float norm2()const;
	float dot(const mat& other)const;
	int getM()const {return m;};
	int getN()const {return n;};
	void print()const;
};

class mat11: public mat{

public:
	mat11() = default;
	mat11(const mat& other);
	float asFloat()const{return this->buffer[0];}
};

class mat21: public mat{

public:
	mat21(float x = 0, float y = 0);
	mat21(const mat& other);
	~mat21() = default;
	float inline getX()const {return this->buffer[0];};
	float inline getY()const {return this->buffer[1];};
	float inline cross(const mat21 & other) const { return this->buffer[0] * other.buffer[1] - this->buffer[1] * other.buffer[0]; }
};

class mat31: public mat{

public:
	mat31(float x = 0, float y = 0, float z = 0);
	mat31(const mat& other);
	~mat31() = default;
	float inline getX()const {return this->buffer[0];};
	float inline getY()const {return this->buffer[1];};
	float inline getZ()const {return this->buffer[2];};
	mat31 inline cross(const mat31 & other) const { return mat31(buffer[1] * other.buffer[2] - buffer[2] * other.buffer[1], buffer[2] * other.buffer[0] - buffer[0] * other.buffer[2], buffer[0] * other.buffer[1] - buffer[1] * other.buffer[0]); }
	mat31 inline change(int i, float val)const{assert(i < 3); mat31 ret(buffer[0], buffer[1], buffer[2]); ret[i] = val; return ret; }
};

class mat12: public mat{

public:
	mat12(float x = 0, float y = 0);
	mat12(const mat& other);
	~mat12() = default;
};

class mat22: public mat{

public:
	mat22(float x11 = 0, float x12 = 0, float x21 = 0, float x22 = 0);
	mat22(const mat& other);
	~mat22() = default;
};

class mat33: public mat{

public:
	mat33(float x11 = 0, float x12 = 0, float x13 = 0, 
		float x21 = 0, float x22 = 0, float x23 = 0, 
		float x31 = 0, float x32 = 0, float x33 = 0);
	mat33(const mat& other);
	mat33(const mat31 & m1, const mat31 & m2, const mat31 & m3);
	mat33 inverse()const;
	~mat33() = default;
};

class eye22: public mat22{
public:
	eye22();
	~eye22() = default;
};

class eye33: public mat33{
public:
	eye33();
	~eye33() = default;
};

class point;
class point_3d;

// base class for all the drawable elements
class element{

public:
	element() = default;
	element(element&) = default;
	virtual ~element() = default;
	virtual point getCentroid() const = 0;
	virtual void transport(const mat21 & v) = 0;
	virtual void linear_transform(const mat22 & r, const point & origin) = 0;
	virtual void print()const = 0;
};

// base class for all the 3d drawable elements
class element_3d{

public:
	element_3d() = default;
	element_3d(element_3d&) = default;
	virtual ~element_3d() = default;
	virtual point_3d getCentroid() const = 0;
	virtual void transport(const mat31 & v) = 0;
	virtual void linear_transform(const mat33 & r, const point_3d & origin) = 0;
	virtual void print()const = 0;
};

class point: public element, public mat21{

public:
	float r, g, b;
	point(float x = 0, float y = 0, float r = 1.0, float g = 1.0, float b = 1.0);
	point(const mat21 & other);
	point(const point& other): mat21(other), r(other.r), g(other.g), b(other.b){};
	~point() = default;
	point getCentroid() const;
	void transport(const mat21 & v);
	void linear_transform(const mat22 & r, const point & origin);
	void print()const;
};

class line: public element{

protected:
	point p1, p2;
public:
	line(const point& p1, const point& p2);
	line(const line& other);
	~line() = default;
	float getLength()const;
	mat21 getDirection()const;
	point getCentroid() const;
	void transport(const mat21 & v);
	void linear_transform(const mat22 & r, const point & origin);
	bool intersect_with(const line & other)const;
	point getP1()const { return this->p1; };
	point getP2()const { return this->p2; };
	void print()const;
};

class polygon: public element{

protected:
	std::vector<point> points;
public:
	polygon(const std::vector<point>& p_set);
	polygon(const polygon & other);
	~polygon() = default;
	point getCentroid() const;
	void transport(const mat21 & v);
	void linear_transform(const mat22 & r, const point & origin);
	void print()const;
	std::vector<point> getPoints()const{ return points; }
};

class point_3d: public element_3d, public mat31{
	public:
	float r, g, b;
	point_3d(float x = 0, float y = 0, float z = 1.0f, float r = 1.0, float g = 1.0, float b = 1.0);
	point_3d(const point_3d & other): mat31(other), r(other.r), g(other.g), b(other.b){}
	point_3d(const mat31 & other);
	~point_3d() = default;
	point_3d getCentroid() const;
	void transport(const mat31 & v);
	void linear_transform(const mat33 & r, const point_3d & origin);
	void print()const;
};

class line_3d: public element_3d{
protected:
	point_3d p1, p2;
public:
	line_3d(const point_3d& p1, const point_3d& p2);
	line_3d(const line_3d& other);
	~line_3d() = default;
	float getLength()const;
	mat31 getDirection()const;
	point_3d getCentroid() const;
	void transport(const mat31 & v);
	void linear_transform(const mat33 & r, const point_3d & origin);
	point_3d getP1() { return this->p1; };
	point_3d getP2() { return this->p2; };
	void print()const;
};

class polyhedron: public element_3d{

public:
	std::vector<point_3d> points;
	std::vector<int> starts;
	std::vector<int> ends;
public:
	polyhedron(const std::vector<point_3d> points,  std::vector<int> starts, std::vector<int> ends);
	polyhedron(const polyhedron & other);
	~polyhedron() = default;
	point_3d getCentroid() const;
	void transport(const mat31 & v);
	void linear_transform(const mat33 & r, const point_3d & origin);
	void print()const;
};

#endif