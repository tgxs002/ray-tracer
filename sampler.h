#ifndef SAMPLER
#define SAMPLER

#include "tracer.h"
#include <cmath>
#include <random>

class sampler{
private:
	mat31 eye;
	mat31 up_left;
	mat31 dir_x, dir_y;
	int resolution;
	mat31 light_source;
	mat31 color_light_source;
	mat31 ambient;
	float K;
	tracer* scene;
	mat31 sample(const mat31 & from_point, const mat31 & dir, int iter) const;
public:
	sampler(const mat31 & _eye, const mat31 & dir, const mat31 & up, float _angle, int _resolution, const mat31 & _light_source, const mat31 & _color_light_source,
		const mat31 & _ambient, float _K, tracer* scene);
	mat31 sample(int x, int y, int iter) const ;
	void draw(std::string filename);
};

#endif