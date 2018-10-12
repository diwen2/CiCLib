#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

const unsigned int ARRAYSIZE = 512; //number of grids on one side of a cubic box
const double CELL = 25;  //10.052941;  //in Mpc/h
const double unit_l = 0.556657800583054e+27;  //in cm
const double h = 0.720000000000000;
const double Mpc = 3.085677581e+24;  // cm per Mpc
const double r = CELL/(unit_l/Mpc*h);
const int numtasks = 2048;

std::vector<std::vector<double> > initcoordinate(double resolution, int chunksize, int rank) {
	// generate x and y with rank of process. linearly divide xy plane into chunks
	std::vector<std::vector<double> > coordinates;
	//printf("Task %d initializing coordinates....\n", rank);
	for (double z = resolution / 2; z < 1; z = z + resolution) {
		for (int i = 0; i < chunksize; i++) {
			std::vector<double> coordinate{(rank*chunksize+i)%ARRAYSIZE*resolution+resolution/2,
				(rank*chunksize+i)/ARRAYSIZE*resolution+resolution/2, z};
			coordinates.push_back(coordinate);
		}
	}
	return coordinates;
}

int main(int argc, char *argv[]){
	std::ifstream in("counts.bin", std::ios_base::binary);
	in.seekg(0, std::ios::end);
	unsigned int num = in.tellg()/sizeof(unsigned int);
	if (num == 512*512*512)
		std::cout << "Reading " << num << " counts from counts.bin...." << std::endl;
	else
		std::cout << "Warning: number of counts doesn't match 512^3." << std::endl;
	in.seekg(0, std::ios::beg);
	std::vector<unsigned int> counts(num/2);
	std::vector<unsigned int> counts2(num/2);
	if(in.good()){
    	in.read((char *)&counts[0], num/2*sizeof(unsigned int));
		in.read((char *)&counts2[0], num/2*sizeof(unsigned int));
	}
	std::copy(counts2.begin(), counts2.end(), std::back_inserter(counts));
	std::cout << "Found " << counts.size() << " counts." << std::endl;

	double resolution = 1 / (double)ARRAYSIZE;
	int chunksize = ARRAYSIZE*ARRAYSIZE/numtasks;
	std::vector<unsigned int> counts_cut;
	for (unsigned int taskid = 0; taskid < numtasks; taskid++) {
		std::vector<std::vector<double> > co = initcoordinate(resolution, chunksize, taskid);
		for (unsigned int i = 0; i < co.size(); i++) {
			if (co[i][0] > r && co[i][1] > r && co[i][2] > r &&
				co[i][0] < 1.0-r && co[i][1] < 1.0-r && co[i][2] < 1.0-r){
			counts_cut.push_back(counts[taskid*chunksize*ARRAYSIZE+i]);
			}
		}
	}
	std::cout << counts_cut.size() << " counts left after cutting edges." << std::endl;
	std::cout << (double)counts_cut.size()/(double)counts.size()*100 << "% counts left." << std::endl;

	std::ofstream countfile;
	std::string filename = "freq_cell_"+std::to_string(CELL)+"_size_"
			+std::to_string(ARRAYSIZE)+"_cut.txt";
	countfile.open(filename);

	unsigned int sum = 0;
	unsigned int max_count = *std::max_element(counts.begin(), counts.end());
	printf("Global maximum count is %d.\n", max_count);
	for (unsigned int i = 0; i < max_count +1; ++i){
		sum += counts_cut[i];
		countfile << counts_cut[i] << " ";
	}

	countfile.close();

	printf("Total counts is %u.\n", sum);
	std::cout << "Counts saved in file " << filename << "." << std::endl;

	return 0;
}
