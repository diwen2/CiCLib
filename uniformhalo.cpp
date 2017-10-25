#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]){
	ofstream out("uniform_0.05.bin",std::ios_base::binary);
	cout << "Creating uniform halo positions ...." << endl;
	vector<vector<float> > coordinates;
	float resolution = 0.05;
	float id = 0;
	float mass = 100;
	for (float x = resolution/2; x < 1;  x = x + resolution) {
		for (float y = resolution/2; y < 1; y = y + resolution) {
			for (float z = resolution/2; z < 1; z = z + resolution) {
				vector<float> coordinate;
				coordinate.push_back(id);
				coordinate.push_back(mass);
				coordinate.push_back(x); // append number count 0 to coordinate
				coordinate.push_back(y);
				coordinate.push_back(z);
				if(out.good())
				{
					cout << "Writing halo: " << id << " " << mass << " " << x << " " << y << " " << z << endl;
					out.write((char *)&coordinate[0],5*sizeof(float));
				}
				coordinates.push_back(coordinate); //append all coordinates
				id = id + 1;
			}
		}
	}
	out.close();
	cout << "Uniform positions saved as uniform_0.05.bin" << endl;
	return 0;
}

vector<vector<double> > initcoordinate(double resolution, double xmin, double xmax) {
	vector<vector<double> > coordinates;
	for (double z = resolution / 2; z < 1; z = z + resolution) {
		for (double y = resolution / 2; y < 1; y = y + resolution) {
			for (double x = xmin; x < xmax + resolution / 2;
					x = x + resolution) {
				vector<double> coordinate;
				coordinate.push_back(x); // append number count 0 to coordinate
				coordinate.push_back(y);
				coordinate.push_back(z);
				coordinate.push_back(0.0);
				coordinates.push_back(coordinate); //append all coordinates
				// cout << coordinate[0] << " ";
				// cout << coordinate[1] << " ";
				// cout << coordinate[2] << " ";
				// cout << coordinate[3] << endl;
			}
		}
	}
	return coordinates;
}



