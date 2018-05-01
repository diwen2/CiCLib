/*
 * ranhalo.cpp
 *
 *  Created on: Apr 24, 2018
 *      Author: Di Wen
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

int main (void)
{

	ofstream out("halo_positions_poisson.bin",std::ios_base::binary);
	cout << "Creating Poisson halo numbers ...." << endl;
	vector<vector<float> > coordinates;
	/*s = no. of division per side = 1/((4/3*pi)^1/3*radius) = 40.21 */
	unsigned int s = 40;
	unsigned int nhalo = 0;  /*sum of total number of halos*/
	float resolution = 1.0/s;
	float id = 0;
	float mass = 100;
	vector<unsigned int> halo(s*s*s, 0);

	/* Generate a Poisson random number as
	 * the no. of halos in a cubic cell */

	const gsl_rng_type * T;
	gsl_rng * r;

	//unsigned int icell = 0;
	unsigned int ncell = 0;
	double mu = 3045306/(40*40*40);

	/* create a generator chosen by the
     environment variable GSL_RNG_TYPE */

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	/* print n random variates chosen from
     the poisson distribution with mean
     parameter mu */

	for (unsigned int icell = 0; icell < s*s*s; ++icell)
	{
		halo[icell] = gsl_ran_poisson (r, mu);
		nhalo = nhalo + halo[icell];
		printf (" %u", halo[icell]);
	}

	printf ("\nTotal halo number is %d.\n", nhalo);
	gsl_rng_free (r);

	for (unsigned int i = 0; i < s; ++i) {
		for (unsigned int j = 0; j < s; ++j) {
			for (unsigned int k = 0; k < s; ++k) {

				float xa = resolution*i;
				float xb = resolution*(i+1);
				float ya = resolution*j;
				float yb = resolution*(j+1);
				float za = resolution*k;
				float zb = resolution*(k+1);

				for (unsigned int n = 0; n < halo[ncell]; ++n) {
					const gsl_rng_type * T;
					gsl_rng * rx, * ry, * rz;
					gsl_rng_env_setup();
					T = gsl_rng_default;
					rx = gsl_rng_alloc (T);
					ry = gsl_rng_alloc (T);
					rz = gsl_rng_alloc (T);
					//printf ("Time: %ld\n", chrono::high_resolution_clock::now().time_since_epoch().count());
					gsl_rng_set(rx, chrono::high_resolution_clock::now().time_since_epoch().count());
					gsl_rng_set(ry, chrono::high_resolution_clock::now().time_since_epoch().count());
					gsl_rng_set(rz, chrono::high_resolution_clock::now().time_since_epoch().count());

					float x, y, z;
					x = gsl_ran_flat(rx, xa, xb);
					y = gsl_ran_flat(ry, ya, yb);
					z = gsl_ran_flat(rz, za, zb);

					gsl_rng_free(rx);
					gsl_rng_free(ry);
					gsl_rng_free(rz);

					vector<float> coordinate;
					coordinate.push_back(id);
					coordinate.push_back(mass);
					coordinate.push_back(x); // append number count 0 to coordinate
					coordinate.push_back(y);
					coordinate.push_back(z);
					if(out.good())
					{
						//cout << "Writing halo: " << id << " " << mass << " " << x << " " << y << " " << z << endl;
						out.write((char *)&coordinate[0],5*sizeof(float));
					}
					coordinates.push_back(coordinate); //append all coordinates
					id = id + 1;
				}
				ncell = ncell +1;
			}
		}
	}
	out.close();
	cout << "Uniform positions saved as halo_positions_poisson.bin" << endl;


	return 0;
}

