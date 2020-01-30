#include "tracer.h"

inline float min(float x, float y){
	return x > y? y : x;
}

inline float max(float x, float y){
	return x < y? y : x;	
}

face::face(float ita, float refl, float refr, float sp): ita(ita), k_refl(refl), k_refr(refr), sp(sp){};

triangle::triangle(const mat31 & _p0, const mat31 & _p1, const mat31 & _p2, const mat31 & _c0, const mat31 & _c1, const mat31 & _c2,
	const mat31 & normal0_, const mat31 & normal1_, const mat31 & normal2_,  
	float _ita, float refl, float refr, float _sp): face(_ita, refl, refr, _sp), 
	p0(_p0), p1(_p1), p2(_p2), c0(_c0), c1(_c1), c2(_c2), normal0(normal0_), normal1(normal1_), normal2(normal2_){

	mat31 normal = mat31(_p0 - _p1).cross(mat31(_p1 - _p2));
	normal /= normal.norm2();
	if (normal0_.norm2() == 0){
		normal2 = normal;
		normal1 = normal;
		normal0 = normal;
	}
}

triangle::triangle(const triangle & other): p0(other.p0), p1(other.p1), p2(other.p2), c0(other.c0), c1(other.c1), c2(other.c2), 
	normal0(other.normal0), normal1(other.normal1), normal2(other.normal2), face(other.ita, other.k_refl, other.k_refr, other.sp){

}

mat31 triangle::interpolate_p(float u, float v)const{
	return (p1 - p0) * u + (p2 - p0) * v + p0;
}

mat31 triangle::interpolate_c(float u, float v)const{
	return (c1 - c0) * u + (c2 - c0) * v + c0;
}

mat31 triangle::normal_(float u, float v)const{
	return (normal1 - normal0) * u + (normal2 - normal0) * v + normal0;
}

mat31 triangle::getub()const{
	float x, y, z;
	x = max(p0.getX(), max(p1.getX(), p2.getX()));
	y = max(p0.getY(), max(p1.getY(), p2.getY()));
	z = max(p0.getZ(), max(p1.getZ(), p2.getZ()));
	return mat31(x, y, z);
}

mat31 triangle::getlb()const{
	float x, y, z;
	x = min(p0.getX(), min(p1.getX(), p2.getX()));
	y = min(p0.getY(), min(p1.getY(), p2.getY()));
	z = min(p0.getZ(), min(p1.getZ(), p2.getZ()));
	return mat31(x, y, z);
}

trace_result triangle::intersects(const mat31 & from_point, const mat31 & dir)const{
	trace_result ret;
	ret.f = this;
	mat33 R(mat31()-dir, p1 - p0, p2 - p0);
	R = R.inverse();
	if (R[0] == 0 && R[1] == 0 && R[2] == 0){
		ret.u = 0;
		ret.v = 0;
		ret.front = false;
		ret.valid = false;
	}
	else{
		mat res = mat(R) * mat(from_point - p0);
		ret.u = res[1];
		ret.v = res[2];
		ret.t = res[0];
		ret.valid = true;
		ret.front = dir.dot(normal_(ret.u, ret.v)) < 0;
	}
	return ret;
}

void triangle::print()const{
	cout << "face: ita: " << ita << " k_refr: " << k_refr << " k_refl: " << k_refl << " sp: " << sp << endl;
	cout << "    p0: ";
	p0.print();
	cout << "p1: ";
	p1.print();
	cout << "p2: ";
	p2.print();
	cout << "normal: ";
	normal0.print();
	normal1.print();
	normal2.print();
	cout << endl;
	cout << "    color: ";
	c0.print();
	c1.print();
	c2.print();
	cout << endl;
}

ball::ball(const mat31 & _center, float _r, const mat31 & _color, 
	float ita, float refl, float refr, float sp): face(ita, refl, refr, sp),
	center(_center), r(_r), color(_color){
	}

ball::ball(const ball & other): face(other.ita, other.k_refl, other.k_refr, other.sp), center(other.center), color(other.color), r(other.r){}

using std::sin;
using std::cos;

mat31 ball::normal_(float u, float v)const{
	float alpha = u * 2 * 3.1415926;
	float beta = v * 4 * 3.1415926;
	return mat31(sin(alpha) * cos(beta), sin(alpha) * sin(beta), cos(alpha));
}

mat31 ball::interpolate_p(float u, float v)const{
	return center + normal_(u, v) * r;
}

mat31 ball::interpolate_c(float u, float v)const{
	return color;
}

mat31 ball::getub()const{
	return center + mat31(r,r,r);
}

mat31 ball::getlb()const{
	return center - mat31(r,r,r);
}

trace_result ball::intersects(const mat31 & from_point, const mat31 & dir)const{
	float scale = dir.norm2();
	mat31 n = dir / scale;
	trace_result ret;
	ret.f = this;
	if (n.cross(from_point - center).norm2() >= r){
		ret.u = 0;
		ret.v = 0;
		ret.front = false;
		ret.valid = false;
	}
	else{
		mat31 diff = from_point - center;
		float s2 = n.dot(diff) * 2;
		float a3 = diff.norm2();
		a3 *= a3;
		a3 -= r * r;
		float delta = std::sqrt(s2 * s2 - 4 * a3);
		float x1 = (-s2 + delta) / 2;
		float x2 = (-s2 - delta) / 2;
		if (x1 > 0){
			if (x2 < 0){
				x2 = x1;
				ret.front = false;
			}
			else{
				ret.front = true;
			}
			mat31 inter = diff + n * x2;
			ret.u = std::acos(inter[2] / r) / 3.1415926 / 2;
			ret.v = std::atan2(inter[1], inter[0]);
			if (ret.v < 0)
				ret.v += 2 * 3.1415926;
			ret .v /= (3.1415926 * 4);
			ret.valid = true;
			ret.t = x2 / scale;
		}
	}
	return ret;
}

void ball::print()const{

}

void trace_result::print()const{
	cout << "trace result: " << "u: " << u << " v: " << v << " t: " << t << endl;
	if (f != NULL){
		f->print();
		cout << endl;
	}
}

tracer_node::tracer_node(const mat31 & _lb, const mat31 & _ub): lb(_lb), ub(_ub){
	for (int i = 0; i < 3; i++){
		int next = (i + 1) % 3;
		cuboid_face.push_back(new triangle(_lb, _lb.change(next, _ub.get(next)), _lb.change(i, _ub.get(i))));
		cuboid_face.push_back(new triangle(_ub, _ub.change(next, _lb.get(next)), _ub.change(i, _lb.get(i))));
	}
}

void tracer_node::insert(const face * f){
	mat31 f_lb = f->getlb();
	mat31 f_ub = f->getub();
	for (int i = 0; i < child_nodes.size(); i++){
		if (child_nodes[i].lb < f_lb && child_nodes[i].ub > f_ub){
			child_nodes[i].insert(f);
			return;
		}
	}
	for (int i = 0; i < 3; i++){
		if (!(ub[i] - lb[i] > 0.05)){
			continue;
		}
		if (f_lb[i] > (lb[i] + ub[i]) / 2){
			mat31 bound(lb);
			bound[i] = (lb[i] + ub[i]) / 2;
			tracer_node child(bound, ub);
			child.insert(f);
			child_nodes.push_back(child);
			return;	
		}
		if (f_ub[i] < (lb[i] + ub[i]) / 2){
			mat31 bound(ub);
			bound[i] = (lb[i] + ub[i]) / 2;
			tracer_node child(lb, bound);
			child.insert(f);
			child_nodes.push_back(child);
			return;
		}
	}
	faces.push_back(f);
}

bool cmp_trace_result(const trace_result r1, const trace_result r2)
{
    return (r1.t < r2.t); 
} 

trace_result tracer_node::trace(const mat31 & from_point, const mat31 & dir, float* best)const{
	trace_result ret;
	vector<trace_result> res;

	for (int i = 0; i < child_nodes.size(); i++){
		for (int j = 0; j < 6; j++){
			trace_result temp = child_nodes[i].cuboid_face[j]->intersects(from_point, dir);
			if (temp.valid &&  temp.within_cube() && temp.t <= *best){
				temp.u = i;
				temp.front = true;
				res.push_back(temp);
				break;
			}
		}
	}
	sort(res.begin(), res.end(), cmp_trace_result);
	for (int i = 0; i < res.size(); i++){
		trace_result temp = child_nodes[res[i].u].trace(from_point, dir, best);
		if (temp.valid && temp.within_triangle() && temp.t > 0 && temp.t <= *best){
			ret = temp;
			*best = temp.t;
		}
	}
	for (int i = 0; i < faces.size(); i++){
		trace_result temp = faces[i]->intersects(from_point, dir);
		if (temp.valid && temp.within_triangle() && temp.t > 0 && temp.t <= *best){
			ret = temp;
			*best = temp.t;
		}
	}
	return ret;
}

void tracer_node::squeeze(){
	for (int i = 0; i < child_nodes.size(); i++){
		child_nodes[i].squeeze();
		if (child_nodes[i].child_nodes.size() < 2 && child_nodes[i].faces.size() <= 10){
			auto end = child_nodes[i].faces.end();
			for (auto p = child_nodes[i].faces.begin(); p < end; p++){
				faces.push_back(*p);
			}
			if (child_nodes[i].child_nodes.size() == 1){
				child_nodes.push_back(child_nodes[i].child_nodes[0]);
			}
			child_nodes.erase(child_nodes.begin() + i);
			i--;
		}
	}
}

void tracer_node::print()const{
	cout << "ub: ";
	ub.print();
	cout << endl << "lb: ";
	lb.print();
	cout << endl;
	cout << endl << "childs" << endl;
	for (auto i : child_nodes){
		i.print();
	}
}

tracer::tracer(mat31 lb, mat31 ub): root(lb, ub){}

trace_result tracer::trace(const mat31 & from_point, const mat31 & dir)const{
	float * best = new float;
	*best = 99999999;
	trace_result ret = root.trace(from_point + dir * 0.0001, dir, best);
	delete best;
	return ret;
}

bool tracer::test_continuity(const mat31 & start, const mat31 & end)const{
	float * best = new float;
	*best = 99999999;
	trace_result ret = root.trace(start + (end - start) * 0.0001, end - start, best);
	delete best;
	if (!ret.valid || ret.t >= 1){
		return true;
	}
	else{
		return false;
	}
}

void tracer::load_from_file(string filename){
	std::ifstream input(filename);
	try{
		if (!input){
			std::cout << "can not find input file" << std::endl;
			return;
		}
		else{
			int number;
			input >> number;
			for (int i = 0; i < number; i++){
				int vertex_number;
				float ita;
				float refl;
				float refr;
				input >> ita;
				input >> refl;
				input >> refr;
				input >> vertex_number;
				std::vector<mat31> points;
				for (int j = 0; j < vertex_number; j++){
					float x, y, z;
					input >> x;
					input >> y;
					input >> z;
					points.push_back(mat31(x, y, z));
				}
				vector<mat31> colors;
				std::vector<mat31> normals;
				vector<int> s;
                for (int j = 0; j < vertex_number; j++){
                    float x, y, z;
                    input >> x;
                    x /= 255.0;
                    input >> y;
                    y /= 255.0;
                    input >> z;
                    z /= 255.0;
                    x = 0;
                    y = 0;
                    z = 0;
                    colors.push_back(mat31(x, y, z));
                    normals.push_back(mat31());
                    s.push_back(0);
                }
				int edge_number;
				input >> edge_number;
				std::vector<int> starts;
                std::vector<int> middles;
				std::vector<int> ends;
				for (int j = 0; j < edge_number; j++){
					int start, middle, end;
					input >> start;
                    input >> middle;
					input >> end;
					starts.push_back(start);
                    middles.push_back(middle);
					ends.push_back(end);
				}
                std::vector<float> specularity;
                for (int j = 0; j < edge_number; j++){
                    float ss;
                    input >> ss;
                    specularity.push_back(ss);

                    int start = starts[j] - 1;
					int middle = middles[j] - 1;
					int end = ends[j] - 1;

					mat31 a = points[start];
					mat31 b = points[middle];
					mat31 c = points[end];

					mat31 normal = mat31(a - b).cross(mat31(b - c));
					normal /= normal.norm2();
					normals[start] += normal;
					normals[middle] += normal;
					normals[end] += normal;
					s[start] ++;
					s[middle] ++;
					s[end] ++;
                }

                for (int i = 0; i < vertex_number; i++){
                	normals[i] /= s[i];
                }

				for (int i = 0; i < edge_number; i++){
					int start = starts[i] - 1;
					int middle = middles[i] - 1;
					int end = ends[i] - 1;


					mat31 a = points[start];
					mat31 b = points[middle];
					mat31 c = points[end];

					root.insert(new triangle(a, b, c, colors[start], colors[middle], colors[end], normals[start], normals[middle], normals[end],
						ita, refl, refr, specularity[i]));
				}
			}
			root.squeeze();
			input.close();
		}		
	}
	catch (char* error){
		std::cout << error << std::endl;
	}
}
