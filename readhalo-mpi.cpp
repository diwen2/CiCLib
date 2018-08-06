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
//#include <omp.h>
//#include <cmath>
//#include <math.h>
#include "machar.h"
// Need run parameters
//#define  ARRAYSIZE	20
//#define  CELL		25.0
#define  MASTER		0

const unsigned int ARRAYSIZE = 512; //number of grids on one side of a cubic box
const double CELL = 10;  //10.052941;  //in Mpc/h
const double unit_l = 0.277801161516035e+28;  //in cm
const double h = 0.720000000000000;
const double Mpc = 3.085677581e+24;  // cm per Mpc
const double radiussq = CELL*CELL/(unit_l/Mpc*h)/(unit_l/Mpc*h);
const unsigned int chknum  = 10;    // number of checkpoint saves before job completion
const unsigned int cacheblk = 4096;    // optimal cache blocking size
unsigned int chkpt = 0;    // Restarted job should start from the last saved checkpoint.

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

vector<vector<double> > genfrombin(){
	ifstream in("halo_positions.bin",ios_base::binary);
	//ifstream in("uniform_0.05.bin",ios_base::binary);
	in.seekg(0,ios::end);
	unsigned int length = in.tellg();
	in.seekg(0,ios::beg);
	//cout << "Input halo positions data size is " << length << " bytes."<< endl;
	vector<vector <double> > read_positions;
	if(in.good())
	{
		vector<double> row(3);
		for (unsigned int i=0; i<length/(3*sizeof(double)); ++i){
			in.read((char *)&row[0],3*sizeof(double));
			read_positions.push_back(row);
		}
	}
	/*
	cout << "1st halo is " << read_positions[0][0] <<" "<<read_positions[0][1]<<" "
			<<read_positions[0][2]<<" "<<read_positions[0][3]<<" "<<read_positions[0][4]<<endl;
	cout << "2nd halo is " << read_positions[1][0] <<" "<<read_positions[1][1]<<" "
			<<read_positions[1][2]<<" "<<read_positions[1][3]<<" "<<read_positions[1][4]<<endl;
	cout << "...." << endl;
	cout << length/(5*sizeof(float)) <<"th halo is " << read_positions[length/(5*sizeof(float))-1][0] <<" "
			<<read_positions[length/(5*sizeof(float))-1][1]<<" "<<read_positions[length/(5*sizeof(float))-1][2]<<" "
			<<read_positions[length/(5*sizeof(float))-1][3]<<" "<<read_positions[length/(5*sizeof(float))-1][4]<<endl;
	*/
	return read_positions;
}

// function to initialize coordiante grid and count

vector<vector<double> > initcoordinate(double resolution, int chunksize, int rank) {
	// generate x and y with rank of process. linearly divide xy plane into chunks
	vector<vector<double> > coordinates;
	printf("Task %d initializing coordinates....\n", rank);
	for (double z = resolution / 2; z < 1; z = z + resolution) {
		for (int i = 0; i < chunksize; i++) {
			vector<double> coordinate{(rank*chunksize+i)%ARRAYSIZE*resolution+resolution/2,
				(rank*chunksize+i)/ARRAYSIZE*resolution+resolution/2, z};
			//coordinate.push_back(0.0);
			//coordinate.push_back(0.0);
			coordinates.push_back(coordinate);
			//cout << coordinate[0] << " ";
			//cout << coordinate[1] << " ";
			//cout << coordinate[2] << " ";
			//cout << coordinate[3] << " ";
			//cout << coordinate[4] << endl;
		}
	}
	return coordinates;
}
/*
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

int main(int argc, char *argv[]) {
	int numtasks, taskid, chunksize, bufsize;
	unsigned int max_count, max_uncertain_count, task_max, task_uncertain_max, num_uncertain;
	double resolution = 1 / (double)ARRAYSIZE;
	//float resolution = 1/(float)ARRAYSIZE;
	float starttime, endtime;
	struct Machar machine_parameter;
	//machine_parameter.report();
	double eps = machine_parameter.eps;
	double epsneg = machine_parameter.epsneg;

	/***** Initializations *****/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	starttime = MPI_Wtime();

	chunksize = ARRAYSIZE*ARRAYSIZE/numtasks; // ARRAYSIZE^2 is divisible by numtasks
	unsigned int counts[chunksize*ARRAYSIZE] = {0};
	unsigned int uncertainty[chunksize*ARRAYSIZE] = {0};
	unsigned int uncertain_counts[chunksize*ARRAYSIZE] = {0};

	printf("MPI task %d has started...\n", taskid);
	/* Load halo catalog */
	printf("Task %d is loading halo catalog...\n", taskid);
	//vector<vector<float> > halo_positions = genfromtxt();
	//vector<vector<float> > halo_positions = genfromFOF(512);
	vector<vector<double> > halo_positions = genfrombin();
	printf("Task %d finished loading %ld halos.\n", taskid, halo_positions.size());

	vector<vector<double> > coord_count = initcoordinate(resolution, chunksize, taskid);
	//vector<vector<double> > coord_count = initcoordinate(resolution,
	//		xarray[offset], xarray[offset + chunksize - 1]);
	//vector<vector<float> > coord_count = initcoordinate(resolution,
	//		xarray[offset], xarray[offset + chunksize - 1]);
	printf("Task %d coordinates start with %e, %e, %e\n", taskid,
			coord_count[0][0], coord_count[0][1], coord_count[0][2]);
	printf("Task %d has %lu coordinates.\n", taskid, coord_count.size());

	unsigned int coordnum = coord_count.size();
	unsigned int halonum = halo_positions.size();
	unsigned int temp = halonum / (chknum * cacheblk);
        unsigned int chkblk = temp * cacheblk;
        bufsize = chunksize * ARRAYSIZE * sizeof(unsigned int);
	
        if (chkpt > 0)    // read existing checkpoint data if chkpt > 0
        {
	        MPI_File fh_counts;
        	MPI_File fh_uncertain;
        	MPI_Status status;
		string cstr = "count_" + to_string(chkpt) + ".bin";
		string ustr = "uncertainty_" + to_string(chkpt) + ".bin";
		MPI_File_open(MPI_COMM_WORLD, cstr.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_counts);
        	MPI_File_open(MPI_COMM_WORLD, ustr.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_uncertain);
                if (taskid == MASTER)
			printf("Checkpoint %d data files found.\nLoading checkpoint %d data ....\n", chkpt, chkpt);
		MPI_File_read_at(fh_counts, taskid*bufsize, counts, coordnum, MPI_UNSIGNED, & status);
                MPI_File_read_at(fh_uncertain, taskid*bufsize, uncertainty, coordnum, MPI_UNSIGNED, & status);
        	MPI_File_close(&fh_counts);
		MPI_File_close(&fh_uncertain);
		if (taskid == MASTER)
			printf("Successfully read checkpoint %d data.\nContinue computation after checkpoint %d ....\n", chkpt, chkpt);
	}
	//omp_set_nested(1);
	//#pragma omp parallel for
	for (unsigned int j = chkpt*chkblk; j < halonum; j += cacheblk){
	    for (unsigned int i = 0; i < coordnum; ++i){
		unsigned int jj = 0;
		while (jj < cacheblk && j + jj < halonum){
			double ratio = ((halo_positions[j+jj][0]-coord_count[i][0])*(halo_positions[j+jj][0]-coord_count[i][0])
				       +(halo_positions[j+jj][1]-coord_count[i][1])*(halo_positions[j+jj][1]-coord_count[i][1])
				       +(halo_positions[j+jj][2]-coord_count[i][2])*(halo_positions[j+jj][2]-coord_count[i][2]))
				       /radiussq;
			
			if (ratio <= double(1.0))
			{
				if (double(1.0)-ratio > epsneg)
				{
					counts[i] = counts[i] + 1;
				}
				else
				{
					uncertainty[i] = uncertainty[i] + 1;
					printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
							halo_positions[j+jj][0], halo_positions[j+jj][1], halo_positions[j+jj][2],
							coord_count[i][0], coord_count[i][1], coord_count[i][2]);
				}
			}
			else
			{
				if (ratio-double(1.0) < eps)
				{
					uncertainty[i] = uncertainty[i] + 1;
					printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
							halo_positions[j+jj][0], halo_positions[j+jj][1], halo_positions[j+jj][2],
							coord_count[i][0], coord_count[i][1], coord_count[i][2]);
				}
			}
			++jj;
			/*
			if (double(1.0)-ratio > epsneg)
			{
				coord_count[i][3] = coord_count[i][3] + 1;
			}
			else if (double(1.0)-ratio <= epsneg && ratio <= double(1.0))
			{
				coord_count[i][4] = coord_count[i][4] +1;
				printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
						halo_positions[j][2], halo_positions[j][3], halo_positions[j][4],
						coord_count[i][0], coord_count[i][1], coord_count[i][2]);
			}
			else if (ratio-double(1.0) < eps && ratio > double(1.0))
			{
				coord_count[i][4] = coord_count[i][4] + 1;
				printf("Can't decide if halo (%e, %e, %e) is in cell (%e, %e, %e).\n",
						halo_positions[j][2], halo_positions[j][3], halo_positions[j][4],
						coord_count[i][0], coord_count[i][1], coord_count[i][2]);
			}
			*/
		}
	    }
		if ((j + cacheblk) % chkblk == 0)
		{
			MPI_File fh_counts;
                	MPI_File fh_uncertain;
                	MPI_Status status;
			++chkpt;
                	string cstr = "count_" + to_string(chkpt) + ".bin";
                	string ustr = "uncertainty_" + to_string(chkpt) + ".bin";
			MPI_File_open(MPI_COMM_WORLD, cstr.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_counts);
			MPI_File_write_at(fh_counts, taskid*bufsize, counts, coordnum, MPI_UNSIGNED, & status);
			MPI_File_close(&fh_counts);
			MPI_File_open(MPI_COMM_WORLD, ustr.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_uncertain);
			MPI_File_write_at(fh_uncertain, taskid*bufsize, uncertainty, coordnum, MPI_UNSIGNED, & status);
			MPI_File_close(&fh_uncertain);
			MPI_Barrier(MPI_COMM_WORLD);    // Barrier makes sure printf executes only after all tasks finish.
			if (taskid == MASTER)
				printf("Checkpoint %d data saved in binary files.\n", chkpt);
		}
	}
	
        // save counts from the last checkpoint to the end
        MPI_File fh_counts;
        MPI_File fh_uncertain;
        MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, "counts.bin", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_counts);
        MPI_File_write_at(fh_counts, taskid*bufsize, counts, coordnum, MPI_UNSIGNED, & status);
        MPI_File_close(&fh_counts);
        MPI_File_open(MPI_COMM_WORLD, "uncertainties.bin", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_uncertain);
        MPI_File_write_at(fh_uncertain, taskid*bufsize, uncertainty, coordnum, MPI_UNSIGNED, & status);
        MPI_File_close(&fh_uncertain);
	MPI_Barrier(MPI_COMM_WORLD);
	if (taskid == MASTER)
                printf("Last checkpoint reached.\ncounts.bin and uncertainties.bin saved.\n");

	unsigned int task_uncertainty = 0;
	for (unsigned int row = 0; row < coordnum; ++row) {
		//counts[row] = coord_count[row][3];
		//uncertainty[row] = coord_count[row][4];
		// treat all uncertain cases as inside the cells and add uncertainty to counts
		uncertain_counts[row] = counts[row] + uncertainty[row];
		task_uncertainty += uncertainty[row];
		//printf("%d ", counts[row]);
	}

	printf("Task %d counts are %d %d ....\n", taskid, counts[0], counts[1]);
	printf("Task %d uncertain counts are %d %d ....\n",
			taskid, uncertain_counts[0], uncertain_counts[1]);

	task_max = *max_element(counts, counts+chunksize*ARRAYSIZE);
	printf("Task %d maximum count is %d.\n", taskid, task_max);
	task_uncertain_max = *max_element(uncertain_counts,
			uncertain_counts+chunksize*ARRAYSIZE);
	printf("Task %d maximum uncertain count is %d.\n", taskid, task_uncertain_max);

	MPI_Reduce(&task_uncertainty, &num_uncertain, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&task_max, &max_count, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&max_count, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Reduce(&task_uncertain_max, &max_uncertain_count, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&max_uncertain_count, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	unsigned int task_sum = 0;
	unsigned int task_Nfreq[max_uncertain_count+1];
	unsigned int task_uncertain_Nfreq[max_uncertain_count+1];
	for (unsigned int i = 0; i < max_count +1; ++i){
		// loop to count the occurrences in vector
		task_Nfreq[i] = count(counts, counts+chunksize*ARRAYSIZE, i);
		task_uncertain_Nfreq[i] =
				count(uncertain_counts, uncertain_counts+chunksize*ARRAYSIZE, i);
		task_sum = task_sum + task_Nfreq[i];
		//cout << task_Nfreq[i] << " ";
	}

	unsigned int Nfreq[max_uncertain_count+1];
	MPI_Reduce(&task_Nfreq, &Nfreq, max_uncertain_count+1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
	unsigned int Nfreq_uncertain[max_uncertain_count + 1];
	MPI_Reduce(&task_uncertain_Nfreq, &Nfreq_uncertain, max_uncertain_count+1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);


	/***** Master task only ******/
	if (taskid == MASTER) {

		/* Print machine's smallest float larger or smaller than 1 */
		printf("Machine eps equals to %.5e.\n", eps);
		printf("Machine epsneg equals to %.5e.\n", epsneg);
		printf("Linear resolution is %.7e.\n",resolution);

		printf("\nTotal number of uncertain cases is %d.\n", num_uncertain);
		printf("Global maximum count is %d.\n", max_count);
		printf("Global maximum uncertain count is %d.\n", max_uncertain_count);

		ofstream countfile;
		string filename = "freq_cell_"+to_string(CELL)+"_size_"+to_string(ARRAYSIZE)+"_p_"+to_string(numtasks)+".txt";
		countfile.open(filename);

		unsigned int sum = 0;
		printf("\nPrinting counts for every N ....\n");
		for (unsigned int i = 0; i < max_uncertain_count +1; ++i){
			printf("%d ", Nfreq[i]);
			sum += Nfreq[i];
			countfile << Nfreq[i] << " ";
		}

		countfile.close();

		printf("\nTotal counts is %u.\n", sum);
		cout << "Counts saved in file " << filename << "." << endl;

		ofstream uncertaincountfile;
		string fileuncertain = "uncertain_cell_" + to_string(CELL) + "_size_"
				+ to_string(ARRAYSIZE) + "_p_" + to_string(numtasks) + ".txt";
		uncertaincountfile.open(fileuncertain);

		unsigned int uncertain_sum = 0;
		printf("\nPrinting uncertain counts for every N ....\n");
		for (unsigned int i = 0; i < max_uncertain_count + 1; ++i) {
			printf("%d ", Nfreq_uncertain[i]);
			uncertain_sum += Nfreq_uncertain[i];
			uncertaincountfile << Nfreq_uncertain[i] << " ";
		}

		uncertaincountfile.close();

		printf("\nTotal counts after including uncertain cases is %d.\n", uncertain_sum);
		cout << "Counts with uncertain cases saved in file " << fileuncertain << "." << endl;

		// store plus/minus error in 2*N long array
		ofstream minusplusfile;
		string minusplusfilename = "minusplus_cell_" + to_string(CELL) + "_size_"
				+ to_string(ARRAYSIZE) + "_p_" + to_string(numtasks) + ".txt";
		minusplusfile.open(minusplusfilename);
		int minusplus[2*max_uncertain_count+2];
		int totaluncertain = 0;
		//printf("\nPrinting error bars of every N ....\n");
		if (max_count == max_uncertain_count){
			//for every N, write lower error, then upper error
			minusplus[0] = Nfreq[0]-Nfreq_uncertain[0];
			minusplus[1] = 0;
			minusplusfile << minusplus[0] << " " << minusplus[1] << " ";
			for (unsigned int i = 1; i < max_uncertain_count +1; ++i){
				//i'th lower error = i'th count + (i-1)th lower error - i'th uncertain count
				minusplus[2*i] = Nfreq[i] + minusplus[2*i-2] - Nfreq_uncertain[i];
				//i'th upper error = (i-1)th lower error
				minusplus[2*i+1] = minusplus[2*i-2];
				minusplusfile << minusplus[2*i] << " " << minusplus[2*i+1] << " ";
				totaluncertain += minusplus[2*i+1];
			}
			printf("\nTotal number of uncertain cases is %d.\n", totaluncertain);
		}
		else {
			printf("\nCounting error occured!\n");
		}
		minusplusfile.close();
		cout << "\nError bars saved in file " << minusplusfilename << "." << endl;
		endtime = MPI_Wtime();
		printf("\nTime taken is %e seconds.\n", endtime-starttime);
	}
	/* end of master section */

	MPI_Finalize();
	return 0;

} /* end of main */
