#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include "FOFReaderLib/FOFReaderLib.h"
//#include <chrono>
#include <mpi.h>
//#include <cmath>
//#include <math.h>
#include "machar.h"
// Need run parameters
//#define  ARRAYSIZE	20
//#define  CELL		25.0
#define  MASTER		0

const int ARRAYSIZE = 40; //number of grids on one side of a cubic box
const int CELL = 25;  //in Mpc/h
const double unit_l = 0.277801161516035e+28;  //in cm
const double h = 0.720000000000000;
const double Mpc = 3.085677581e+24;  // cm per Mpc
const double radius = CELL/(unit_l/Mpc*h);
//const float radius = 0.05; // >0.0499999 produces wrong counts
const double epsilon = 1.0e-7;

using namespace std;

// initialize a vector of vectors that each hold a halo's parameter from FOF
vector<vector<float> > genfromFOF(int num_files) {

	vector<vector<float> > halo_positions;

	// read 512 FOF halo position files using FOFMasst
	char filename[128];
	for (int j = 0; j < num_files; j++) {
		// need to modify filename prefix for different halo catalog
		sprintf(filename, "%s_%05d",
				"list_503/fof_boxlen648_n2048_lcdmw7_masst", j);
		FOFMasst masst(filename);
		for (int i = 0; i < masst.nHalos(); i++) {
			vector<float> row;
			float id = masst.halos(i)->id();
			float mass = masst.halos(i)->mass();
			float x = masst.halos(i)->x();
			float y = masst.halos(i)->y();
			float z = masst.halos(i)->z();
			row.push_back(id);
			row.push_back(mass);
			row.push_back(x);
			row.push_back(y);
			row.push_back(z);
			halo_positions.push_back(row);
			//cout <<halo_positions[0][0]<<" "<<halo_positions[0][1]<<" "<<halo_positions[0][2]
			//<<" "<<halo_positions[0][3]<<" "<<halo_positions[0][4]<<endl;
		}
	}
	//cout << halo_positions[0][2]<<" "<<halo_positions[6460][2] << endl;
	return halo_positions;
}


// initialize a vector of vectors that each hold a halo's parameter from text file
// should avoid using CSV files because writing to file changes data precision
vector<vector<float> > genfromtxt() {

	vector<vector<float> > halo_positions;

	// read the text file output of FOFReader
	ifstream file("halos_output.csv");

	string line;
	// read the input file a line at a time until file hits end-of-file
	while (getline(file, line)) {
		vector<float> row;
		istringstream iss(line);  // bind iss to the line read
		for (int n = 0; n < 5; ++n) {
			string val;
			iss >> val;
			float num = atof(val.c_str());
			row.push_back(num);
			// cout << num << " ";
		}
		// cout << endl;
		halo_positions.push_back(row);
	}
	// cout << halo_positions[0][2];
	return halo_positions;
}

vector<vector<float> > genfrombin(){
	ifstream in("halo_positions.bin",ios_base::binary);
	//ifstream in("uniform_0.05.bin",ios_base::binary);
	in.seekg(0,ios::end);
	int length = in.tellg();
	in.seekg(0,ios::beg);
	cout << "Input halo positions data size is " << length << " bytes."<< endl;
	vector<vector <float> > read_positions;
	if(in.good())
	{
		vector<float> row(5);
		for (unsigned int i=0; i<length/(5*sizeof(float)); ++i){
			in.read((char *)&row[0],5*sizeof(float));
			read_positions.push_back(row);
		}
	}
	cout << "1st halo is " << read_positions[0][0] <<" "<<read_positions[0][1]<<" "
			<<read_positions[0][2]<<" "<<read_positions[0][3]<<" "<<read_positions[0][4]<<endl;
	cout << "2nd halo is " << read_positions[1][0] <<" "<<read_positions[1][1]<<" "
			<<read_positions[1][2]<<" "<<read_positions[1][3]<<" "<<read_positions[1][4]<<endl;
	cout << "...." << endl;
	cout << length/(5*sizeof(float)) <<"th halo is " << read_positions[length/(5*sizeof(float))-1][0] <<" "
			<<read_positions[length/(5*sizeof(float))-1][1]<<" "<<read_positions[length/(5*sizeof(float))-1][2]<<" "
			<<read_positions[length/(5*sizeof(float))-1][3]<<" "<<read_positions[length/(5*sizeof(float))-1][4]<<endl;
	return read_positions;
}

// function to initialize coordiante grid and count

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

// function to initialize coordiante grid and count
/*
vector<vector<float> > initcoordinate(float resolution, float xmin, float xmax) {
	vector<vector<float> > coordinates;
	for (float z = resolution / 2; z < 1; z = z + resolution) {
		for (float y = resolution / 2; y < 1; y = y + resolution) {
			for (float x = xmin; x < xmax + resolution / 2;
					x = x + resolution) {
				vector<float> coordinate;
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
*/

// ARRAYSIZE/numtasks has to an integer
// i.e. number of grid on an axis must be divisible by number of processes
double xarray[ARRAYSIZE];
//float xarray[ARRAYSIZE];
unsigned int counts[ARRAYSIZE*ARRAYSIZE*ARRAYSIZE];
unsigned int uncertainty[ARRAYSIZE*ARRAYSIZE*ARRAYSIZE];
//vector<vector<float> > halo_positions = genfromtxt();
//vector<vector<float> > halo_positions = genfromFOF(512);
vector<vector<float> > halo_positions = genfrombin();

int main(int argc, char *argv[]) {
	int numtasks, taskid, dest, offset, tag1, tag2, source, chunksize;
	double resolution = 1 / (double)ARRAYSIZE;
	//float resolution = 1/(float)ARRAYSIZE;
	printf("Linear resolution is %.7e.\n",resolution);
	float starttime, endtime;
	struct Machar machine_parameter;
	//machine_parameter.report();
	double eps = machine_parameter.eps;
	double epsneg = machine_parameter.epsneg;

	MPI_Status status;

	/***** Initializations *****/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	printf("MPI task %d has started...\n", taskid);
	chunksize = ARRAYSIZE / numtasks;
	tag2 = 1;
	tag1 = 2;
	starttime = MPI_Wtime();

	/***** Master task only ******/
	if (taskid == MASTER) {

		/* Initialize the array */
		for (unsigned int i = 0; i < ARRAYSIZE; i++) {
			xarray[i] = resolution / 2 + i * resolution;
		}
		printf("Initialized array element = %e %e .... %e\n", xarray[0], xarray[1],
				xarray[ARRAYSIZE - 1]);

		/* Send each task its portion of the array - master keeps 1st part */
		offset = chunksize;
		for (dest = 1; dest < numtasks; dest++) {
			MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
			MPI_Send(&xarray[offset], chunksize, MPI_DOUBLE, dest, tag2,MPI_COMM_WORLD);
			//MPI_Send(&xarray[offset], chunksize, MPI_FLOAT, dest, tag2,MPI_COMM_WORLD);
			printf("Sent %d elements to task %d offset= %d\n", chunksize, dest,
					offset);
			offset = offset + chunksize;
		}

		/* Master does its part of the work */
		offset = 0;
		vector<vector<double> > coord_count = initcoordinate(resolution,
				xarray[0], xarray[chunksize - 1]);
		//vector<vector<float> > coord_count = initcoordinate(resolution,
		//		xarray[0], xarray[chunksize - 1]);
		printf("Task %d coordinates start from %e, %e, %e, %e\n", taskid,
		coord_count[0][0], coord_count[0][1], coord_count[0][2], coord_count[0][3]);

		printf("Task 0 has %d coordinates.\n", (int)coord_count.size());

		for (unsigned long int i = 0; i < coord_count.size(); ++i){
			for (unsigned int j = 0; j < halo_positions.size(); ++j){
				double ratio = ((halo_positions[j][2]-coord_count[i][0])*(halo_positions[j][2]-coord_count[i][0])
						   	   +(halo_positions[j][3]-coord_count[i][1])*(halo_positions[j][3]-coord_count[i][1])
							   +(halo_positions[j][4]-coord_count[i][2])*(halo_positions[j][4]-coord_count[i][2]))
							   /(radius*radius);
				if (ratio <= double(1.0))
				{
					if (double(1.0)-ratio > epsneg)
					{
						coord_count[i][3] = coord_count[i][3] + 1;
					}
					else
					{
						coord_count[i][4] = coord_count[i][4] + 1;
						printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
								halo_positions[j][2], halo_positions[j][3], halo_positions[j][4],
								coord_count[i][0], coord_count[i][1], coord_count[i][2]);
					}
				}
				else
				{
					if (ratio-double(1.0) < eps)
					{
						coord_count[i][4] = coord_count[i][4] + 1;
						printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
								halo_positions[j][2], halo_positions[j][3], halo_positions[j][4],
								coord_count[i][0], coord_count[i][1], coord_count[i][2]);
					}
				}
			}
		}
		/*
		for (unsigned long int row = 0; row < coord_count.size(); row++) {
			for (int col = 0; col < 4; col++) {
				cout << coord_count[row][col] << " ";
			}
			printf("\n");
		}
		*/
		//printf("Task 0 counts are %d %d ....\n", counts[0], counts[1]);
		for (unsigned long int row = 0; row < coord_count.size(); ++row) {
			counts[row] = coord_count[row][3];
			uncertainty[row] = coord_count[row][4];
			//printf("%d ", counts[row]);
		}
		printf("\n");
		printf("Task 0 counts are %d %d ....\n", counts[0], counts[1]);

		/* Wait to receive results from each task */
		for (int i = 1; i < numtasks; i++) {
			source = i;
			MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
			MPI_Recv(&counts[offset*ARRAYSIZE*ARRAYSIZE], coord_count.size(), MPI_UNSIGNED, source, tag2, MPI_COMM_WORLD, &status);
			MPI_Recv(&uncertainty[offset*ARRAYSIZE*ARRAYSIZE], coord_count.size(), MPI_UNSIGNED, source, tag2, MPI_COMM_WORLD, &status);
		}

		printf("\nPrinting all counts:\n");
		for (unsigned long int i = 0; i < ARRAYSIZE*ARRAYSIZE*ARRAYSIZE; ++i) {
			printf("%d ", counts[i]);
		}
		printf("\n");

		unsigned int max_count = *max_element(counts, counts+ARRAYSIZE*ARRAYSIZE*ARRAYSIZE);
		printf("\nMaximum count is %d.\n", max_count);

		ofstream countfile;
		string filename = "freq_cell_"+to_string(CELL)+"_size_"+to_string(ARRAYSIZE)+"_p_"+to_string(numtasks)+".txt";
		countfile.open(filename);

		unsigned long int sum = 0;
		vector<int> Nfreq(max_count+1);
		printf("\nPrinting counts for every N ....\n");
		for (unsigned int i = 0; i < max_count +1; ++i){
			// loop to count the occurrences in vector
			unsigned int icount = count(counts, counts+ARRAYSIZE*ARRAYSIZE*ARRAYSIZE, i);
			sum = sum + icount;
			cout << icount << " ";
			countfile << icount << " ";
			Nfreq[i] = icount;
		}

		/*
		printf("\nPrinting N frequencies for every N:\n");
		for (unsigned int i = 0; i < max_count+1; ++i){
			cout << Nfreq[i] << " ";
		}
		*/
		countfile.close();

		printf("\nTotal counts is %lu.\n", sum);
		cout << "Counts saved in file " << filename << "." << endl;

		unsigned int num_uncertain = 0;
		for (unsigned long int i = 0; i < ARRAYSIZE*ARRAYSIZE*ARRAYSIZE; ++i) {
			num_uncertain += uncertainty[i];
		}
		printf("\nTotal number of uncertain cases is %d.\n", num_uncertain);

		// treat all uncertain cases as inside the cells and add uncertainty to counts
		for (unsigned long int i = 0; i < ARRAYSIZE*ARRAYSIZE*ARRAYSIZE; ++i) {
			// original counts overwritten by counts including uncertain cases
			counts[i] = counts[i] + uncertainty[i];
			printf("%d ", counts[i]);
		}
		printf("\n");

		unsigned int max_uncertain_count = *max_element(counts,
				counts + ARRAYSIZE * ARRAYSIZE * ARRAYSIZE);
		printf("\nMaximum uncertain count is %d.\n", max_uncertain_count);

		ofstream uncertaincountfile;
		string fileuncertain = "uncertain_cell_" + to_string(CELL) + "_size_"
				+ to_string(ARRAYSIZE) + "_p_" + to_string(numtasks) + ".txt";
		uncertaincountfile.open(fileuncertain);

		unsigned long int uncertain_sum = 0;
		vector <int> Nfreq_uncertain(max_uncertain_count + 1);
		printf("\nPrinting counts for every N ....\n");
		for (unsigned int i = 0; i < max_uncertain_count + 1; ++i) {
			// loop to count the occurrences in vector
			unsigned int icount = count(counts, counts + ARRAYSIZE * ARRAYSIZE * ARRAYSIZE, i);
			uncertain_sum = uncertain_sum + icount;
			cout << icount << " ";
			uncertaincountfile << icount << " ";
			Nfreq_uncertain[i] = icount;
		}
		/*
		printf("\nPrinting N uncertain frequencies for every N ....\n");
		for (unsigned int i = 0; i < max_uncertain_count+1; ++i){
			cout << Nfreq_uncertain[i] << " ";
		}
		*/
		uncertaincountfile.close();

		printf("\nTotal counts after including uncertain cases is %lu.\n", uncertain_sum);
		cout << "Counts with uncertain cases saved in file " << fileuncertain << "." << endl;

		// store plus/minus error in 2*N long array
		ofstream plusminusfile;
		string plusminusfilename = "plusminus_cell_" + to_string(CELL) + "_size_"
				+ to_string(ARRAYSIZE) + "_p_" + to_string(numtasks) + ".txt";
		plusminusfile.open(plusminusfilename);
		int plusminus[2*max_uncertain_count+2];
		int totaluncertain = 0;
		printf("\nPrinting error bars of every N ....\n");
		if (max_count < max_uncertain_count){
			// append zeros to Nfreq if Nfreq is shorter than Nfreq_uncertain
			vector <int> zeros(max_uncertain_count-max_count, 0);
			Nfreq.insert(Nfreq.end(), zeros.begin(), zeros.end());
			if (Nfreq.size() == Nfreq_uncertain.size()){
				printf("\nAdded zeros after maximum count of N.\n");
				plusminus[0] = 0;
				plusminus[1] = Nfreq[0]-Nfreq_uncertain[0];
				plusminusfile << plusminus[0] << " " << plusminus[1] << " ";
				for (unsigned int i = 1; i < max_uncertain_count +1; ++i){
					plusminus[2*i] = plusminus[2*i-1];
					plusminus[2*i+1] = Nfreq[i] + plusminus[2*i] - Nfreq_uncertain[i];
					plusminusfile << plusminus[2*i] << " " << plusminus[2*i+1] << " ";
					totaluncertain += plusminus[2*i];
				}
				printf("\nTotal number of uncertain cases is %d.\n", totaluncertain);
			}
		}
		else if (max_count == max_uncertain_count){
			plusminus[0] = 0;
			plusminus[1] = Nfreq[0]-Nfreq_uncertain[0];
			plusminusfile << plusminus[0] << " " << plusminus[1] << " ";
			for (unsigned int i = 1; i < max_uncertain_count +1; ++i){
				plusminus[2*i] = plusminus[2*i-1];
				plusminus[2*i+1] = Nfreq[i] + plusminus[2*i] - Nfreq_uncertain[i];
				plusminusfile << plusminus[2*i] << " " << plusminus[2*i+1] << " ";
				totaluncertain += plusminus[2*i];
			}
			printf("\nTotal number of uncertain cases is %d.\n", totaluncertain);
		}
		else {
			printf("\nCounting error occured!\n");
		}
		plusminusfile.close();
		cout << "\nError bars saved in file " << plusminusfile << "." << endl;
		endtime = MPI_Wtime();
		printf("\nTime is %e seconds.\n", endtime-starttime);
	} /* end of master section */

	/***** Non-master tasks only *****/

	if (taskid > MASTER) {

		/* Receive my portion of array from the master task */
		source = MASTER;
		MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
		MPI_Recv(&xarray[offset], chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
		//MPI_Recv(&xarray[offset], chunksize, MPI_FLOAT, source, tag2, MPI_COMM_WORLD, &status);

		vector<vector<double> > coord_count = initcoordinate(resolution,
				xarray[offset], xarray[offset + chunksize - 1]);
		//vector<vector<float> > coord_count = initcoordinate(resolution,
		//		xarray[offset], xarray[offset + chunksize - 1]);
		printf("Task %d coordinates start from %e, %e, %e, %e\n", taskid,
				coord_count[0][0], coord_count[0][1], coord_count[0][2], coord_count[0][3]);
		printf("Task %d has %lu coordinates.\n", taskid, coord_count.size());

		for (unsigned long int i = 0; i < coord_count.size(); ++i){
			for (unsigned int j = 0; j < halo_positions.size(); ++j){
				double ratio = ((halo_positions[j][2]-coord_count[i][0])*(halo_positions[j][2]-coord_count[i][0])
						   	   +(halo_positions[j][3]-coord_count[i][1])*(halo_positions[j][3]-coord_count[i][1])
							   +(halo_positions[j][4]-coord_count[i][2])*(halo_positions[j][4]-coord_count[i][2]))
							   /(radius*radius);
				if (ratio <= double(1.0))
				{
					if (double(1.0)-ratio > epsneg)
					{
						coord_count[i][3] = coord_count[i][3] + 1;
					}
					else
					{
						coord_count[i][4] = coord_count[i][4] + 1;
						printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
								halo_positions[j][2], halo_positions[j][3], halo_positions[j][4],
								coord_count[i][0], coord_count[i][1], coord_count[i][2]);
					}
				}
				else
				{
					if (ratio-double(1.0) < eps)
					{
						coord_count[i][4] = coord_count[i][4] + 1;
						printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
								halo_positions[j][2], halo_positions[j][3], halo_positions[j][4],
								coord_count[i][0], coord_count[i][1], coord_count[i][2]);
					}
				}
			}
		}
		/*
		for (unsigned long int row = 0; row < coord_count.size(); ++row) {
			for (int col = 0; col < 4; col++) {
				cout << coord_count[row][col] << " ";
			}
			printf("\n");
		}
		*/

		unsigned long int counts_offset = offset*ARRAYSIZE*ARRAYSIZE;
		//printf("Task %d counts are %d %d....\n", taskid, counts[counts_offset], counts[counts_offset+1]);
		for (unsigned long int row = 0; row < coord_count.size(); ++row) {
			counts[counts_offset+row] = coord_count[row][3];
			uncertainty[counts_offset+row] = coord_count[row][4];
			//printf("%d ", counts[counts_offset+row]);
		}
		printf("\n");
		printf("Task %d counts are %d %d ....\n", taskid, counts[counts_offset], counts[counts_offset+1]);

		/* Send my results back to the master task */
		dest = MASTER;
		MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
		MPI_Send(&counts[counts_offset], coord_count.size(), MPI_UNSIGNED, MASTER, tag2, MPI_COMM_WORLD);
		MPI_Send(&uncertainty[counts_offset], coord_count.size(), MPI_UNSIGNED, MASTER, tag2, MPI_COMM_WORLD);

	} /* end of non-master */

	MPI_Finalize();
	return 0;

} /* end of main */
