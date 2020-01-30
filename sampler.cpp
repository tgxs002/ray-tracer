#include "sampler.h"

using std::sin;
using std::cos;
using std::tan;

sampler::sampler(const mat31 & _eye, const mat31 & dir, const mat31 & up, float _angle, int _resolution, const mat31 & _light_source, const mat31 & _color_light_source,
	const mat31 & _ambient, float _K, tracer* _scene):
	eye(_eye), resolution(_resolution), light_source(_light_source), color_light_source(_color_light_source), ambient(_ambient), K(_K), scene(_scene){
	mat31 left = up.cross(dir);
	mat31 real_up = dir.cross(left);
	real_up /= real_up.norm2();
	mat31 real_dir = dir / dir.norm2();
	left /= left.norm2();
	float half = tan(_angle / 180 * 3.1415926);
	up_left = _eye + real_dir + (left + real_up) * half;
	dir_x = (mat31() - left) * (half * 2 / _resolution);
	dir_y = (mat31() - real_up) * (half * 2 / _resolution);
	up_left += ((dir_x + dir_y) / 2);
}

mat31 sampler::sample(int x, int y, int iter)const{
	if (scene == NULL){
		cout << "Tracer lost" << endl;
		exit(0);
	}
	mat31 accu;
	int times = K;
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distrb;
	for (int i = 0; i < times; i++){
		mat31 target = up_left + dir_x * (x + distrb(generator) - 0.5) + dir_y * (y + distrb(generator) - 0.5);	
		accu += sample(eye, target - eye, iter);
	}
	return accu / times;
}	

mat31 sampler::sample(const mat31 & from_point, const mat31 & dir, int iter) const{
	trace_result target = scene->trace(from_point, dir);
	if (!target.valid){
		return mat31(0,0,0);
	}

	mat31 p_t = target.f->interpolate_p(target.u, target.v);
	mat31 c_t = target.f->interpolate_c(target.u, target.v);

	mat31 ls = color_light_source;
	mat31 eye_dis = from_point - p_t;
	mat31 lsp = light_source - p_t;
	float coe0 = (eye_dis.norm2() + lsp.norm2());
	float t_ = lsp.norm2();
	if (t_ == 0){
		std::cout << "Divide by zero(0), try other parameters" << std::endl;
		std::cout << "Light Source:" << std::endl;
		return mat31();
	}
	lsp /= t_;
	mat31 n = target.f->normal_(target.u, target.v);
	t_ = n.norm2();
	if (t_ == 0){
		std::cout << "Divide by zero(1), try other parameters" << std::endl;
		return mat31();
	}
	// n /= t_;
	float coe1 = lsp.dot(n);
	coe1 = coe1 > 0 ? coe1 : 0;
	float cos = n.dot(lsp);
	mat31 r = n * (2 * cos) - lsp;
	t_ = eye_dis.norm2();
	if (t_ == 0){
		std::cout << "Divide by zero(2), try other parameters" << std::endl;
		std::cout << "Eye:" << std::endl;
		return mat31();
	}
	float coe2 = (eye_dis / t_).dot(r);

	coe2 = coe2 > 0 ? std::pow(coe2, target.f->sp) : 0;
	float inte = ls.norm2();

	if (!scene->test_continuity(p_t, light_source)){
		coe1 = 0;
		coe2 = 0;
	}

	mat31 color = ambient + c_t * (1.0 / coe0 * coe1) + ls * (1 * coe2 / coe0);

	if (target.f->k_refl > 0 & iter > 1){
		mat31 _f = dir / (dir.norm2() * (-1));
		mat31 _n = n;
		_n /= (_n.norm2());
		_n *= (2 * _n.dot(_f));
		mat31 _o = _n - _f;
		color += sample(p_t, _o, iter - 1) * target.f->k_refl;
	}

	if (target.f->k_refr > 0 && iter > 1){
		mat31 d_l = n * (dir.dot(n));
		mat31 d_r = dir - d_l;
		if (target.front){
			d_r /= target.f->ita;
		}
		else{
			d_r *= target.f->ita;
		}
		float len = d_r.norm2();
		if (len > 1.0){
			return color;
		}
		float s = std::sqrt(1 - len * len);
		mat31 dir_ = d_l * (s / d_l.norm2()) + d_r;
		color += sample(p_t, dir_, iter - 1) * target.f->k_refr;
	}

	// if (color.norm2() < 0.18){
	// 	cout << "here" << endl;
	// }

	return color;
}

void sampler::draw(std::string filename){

}