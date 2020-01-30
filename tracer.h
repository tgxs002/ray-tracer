#ifndef TRACER
#define TRACER

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include "elements.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::sort;

class trace_result;

class face{
public:
	float ita, k_refl, k_refr;
	float sp;
	face(float ita = 0, float refl = 0.3, float refr = 0, float sp = 0.2);
	face(const face & other);
	virtual mat31 normal_(float u, float v)const = 0;
	virtual mat31 interpolate_p(float u, float v)const = 0;
	virtual mat31 interpolate_c(float u, float v)const = 0;
	virtual mat31 getub()const = 0;
	virtual mat31 getlb()const = 0;
	virtual trace_result intersects(const mat31 & from_point, const mat31 & dir)const = 0;
	virtual void print()const = 0;
};

class triangle: public face{
private:
	mat31 p0, p1, p2;
	mat31 c0, c1, c2;
	mat31 normal0, normal1, normal2;
public:
	triangle(const mat31 & p0, const mat31 & p1, const mat31 & p2, 
		const mat31 & c0 = mat31(1,0,0), const mat31 & c1 = mat31(1,0,0), 
		const mat31 & normal0_ = mat31(0,0,0), const mat31 & normal1_ = mat31(0,0,0), const mat31 & normal2_ = mat31(0,0,0), 
		const mat31 & c2 = mat31(1,0,0), 
		float ita = 0, float refl = 0.3, float refr = 0, float sp = 2);
	triangle(const triangle & other);
	mat31 normal_(float u, float v)const;
	mat31 interpolate_p(float u, float v)const;
	mat31 interpolate_c(float u, float v)const;
	mat31 getub()const;
	mat31 getlb()const;
	trace_result intersects(const mat31 & from_point, const mat31 & dir)const;
	void print()const;
};

class ball:public face{
private:
	mat31 center;
	float r;
	mat31 color;
public:
	ball(const mat31 & _center, float _r, const mat31 & _color, 
		float ita = 0, float refl = 0.3, float refr = 0, float sp = 2);
	ball(const ball & other);
	mat31 normal_(float u, float v)const;
	mat31 interpolate_p(float u, float v)const;
	mat31 interpolate_c(float u, float v)const;
	mat31 getub()const;
	mat31 getlb()const;
	trace_result intersects(const mat31 & from_point, const mat31 & dir)const;
	void print()const;
};

struct trace_result{
	const face* f;
	float u;
	float v;
	float t;
	bool front;
	bool valid;
	trace_result(): f(NULL), u(0), v(0), t(0), front(false), valid(false){};
	bool within_triangle()const { return u >= 0 && v >= 0 && (u + v) <= 1; }
	bool within_cube()const { return u >= 0 && v >= 0 && u <= 1 && v <= 1; }
	void print()const;
};

class tracer_node{
private:
	mat31 lb, ub;
	vector<const face*> faces;
	vector<tracer_node> child_nodes;
	vector<face*> cuboid_face;
public:
	tracer_node(const mat31 &  lb, const mat31 & ub);
	void insert(const face * f);
	trace_result trace(const mat31 & from_point, const mat31 & dir, float* best)const;
	void squeeze();
	void print()const;
};

class tracer{
private:
public:
	tracer_node root;
public:
	tracer(mat31 lb, mat31 ub);
	trace_result trace(const mat31 & from_point, const mat31 & dir)const;
	bool test_continuity(const mat31 & start, const mat31 & end)const;
	void load_from_file(string filename);
};

#endif