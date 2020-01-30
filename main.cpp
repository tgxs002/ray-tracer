#include "sampler.h"
#include <random>
#include <ctime>
#include <chrono>
#include <string>

using std::swap;

int main(){

	tracer tr(mat31(0,0,-0.5), mat31(1,1,1));

	float ux = 0.65, uy = 0.75, uz = 0.8, lx = 0.35, ly = 0.55, lz = -0.4;
	swap(ux, lx);
	swap(uy, ly);
	swap(uz, lz);
	mat31 up_right_back(lx, ly, lz), curb(1,0,0);
	mat31 down_left_front(ux, uy, uz), cdlf(1,1,1);
	mat31 up_right_front(lx, ly, uz), curf(1,1,1);
	mat31 down_left_back(ux, uy, lz), cdlb(0,0,1);
	mat31 down_right_back(lx, uy, lz), cdrb(0,1,0);
	mat31 up_left_front(ux, ly, uz), culf(1,1,1);
	mat31 down_right_front(lx, uy, uz), cdrf(1,1,1);
	mat31 up_left_back(ux, ly, lz), culb(1,0,1);

	/* a ball */
	tr.root.insert(new ball(mat31((lx + ux) / 2 + 0.05, (ly + uy) / 2 - 0.05, (lz + uz) / 2 + 0.05), 0.05, mat31(0,0,0), 1.8, 0.2, 0.7, 2));

	/* load the rabbit */
	tr.load_from_file(string("./bunny_scene"));

	/* left and right */
	tr.root.insert(new triangle(up_left_back, down_left_back, up_left_front, culb, cdlb, culf)); 
	tr.root.insert(new triangle(down_left_front, up_left_front, down_left_back, cdlf, culf, cdlb));  
	tr.root.insert(new triangle(up_right_back, up_right_front, down_right_back, curb, curf, cdrb)); 
	tr.root.insert(new triangle(down_right_front, down_right_back, up_right_front, cdrf, cdrb, curf));  
	/* back and front*/
	tr.root.insert(new triangle(up_left_back, up_right_back, down_left_back, culb, curb, cdlb)); 
	tr.root.insert(new triangle(up_right_back, down_right_back, down_left_back, curb, cdrb, cdlb)); 
	tr.root.insert(new triangle(up_left_front, down_left_front, up_right_front, culf, cdlf, curf));
	tr.root.insert(new triangle(up_right_front, down_left_front, down_right_front, curf, cdlf, cdrf));
	/* up and down */
	tr.root.insert(new triangle(up_right_back, up_left_back, up_left_front, curb, culb, culf));
	tr.root.insert(new triangle(up_right_front, up_right_back, up_left_front, curf, curb, culf));
	tr.root.insert(new triangle(down_right_back, down_left_front, down_left_back, cdrb, cdlf, cdlb));
	tr.root.insert(new triangle(down_right_front, down_left_front, down_right_back, cdrf, cdlf, cdrb));

	std::ifstream config("./config");
	mat31 from_point;
	config >> from_point[0];
	config >> from_point[1];
	config >> from_point[2];
	mat31 at_point;
	config >> at_point[0];
	config >> at_point[1];
	config >> at_point[2];
	mat31 up_vector;
	config >> up_vector[0];
	config >> up_vector[1];
	config >> up_vector[2];
	float view_angle;
	config >> view_angle;
	int reso;
	mat31 light_source;
	config >> light_source[0];
	config >> light_source[1];
	config >> light_source[2];
	config >> reso;
	int refl;
	config >> refl;
	int k;
	config >> k;

	sampler s(from_point, at_point, up_vector, view_angle, reso, light_source, mat31(1,1,1),
		mat31(0.1,0.1,0.1), k, &tr);

	int i = 0;
	vector<mat31> store;
	for (int i = 0; i < reso; i++){
		store.push_back(mat31());
	}
	while (true){
		std::ifstream index_in("./render/index");
		index_in >> i;
		index_in.close();
		if (i >= reso){
			break;
		}
		std::ofstream index_out("./render/index");
		index_out << (i + 1);
		index_out.close();
		auto begin = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for
		for (int j = 0; j < reso; j++){
			store[j] = s.sample(i, j, refl);
		}
		auto end = std::chrono::high_resolution_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	    cout << ms / 1000.0 << " seconds used " << i << " rows complete." << endl;
		/* write to file */
	    std::ofstream output("./render/" + std::to_string(i));
		for (int i = 0; i < reso; i++){
			output << store[i][0] << " " << store[i][1] << " " << store[i][2] << endl;	
		}
		output.close();
	}

	return 0;
}