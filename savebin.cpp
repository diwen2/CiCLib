#include <fstream>
#include <iostream>
#include <string>
#include "FOFReaderLib/FOFReaderLib.h"

using namespace std;

int main(int argc, char *argv[])
{
	string databin = "halo_positions_439.bin";
	ofstream out(databin,std::ios_base::binary);
	vector<vector<float> > halo_positions;
	char filename[128];
	int totalhalo = 0;
	for (int j = 0; j < 512; j++) {
		// need to modify filename prefix for different halo catalog
		sprintf(filename, "%s_%05d",
				"list_439/fof_boxlen648_n2048_wcdmw7_masst", j);
		FOFMasst masst(filename);
		cout << masst.nHalos() << " halos in " << filename << endl;
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
			if(out.good())
			{
				//if(mass > 119.6) // >119.6 for rpcdm, >106.9 for lcdm
				{
					//cout << "Writing halo: " << id << " " << mass << " " << x << " " << y << " " << z << endl;
					out.write((char *)&row[0],5*sizeof(float));
					halo_positions.push_back(row);
				}
			}
		}
		totalhalo += masst.nHalos();
	}
	out.close();
	cout << "Total number of halos is " << totalhalo << endl;
	cout << "Halo positions saved as halo_positions.bin" <<endl;
	ifstream in(databin,ios_base::binary);
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
	int nparmin = read_positions[0][1];
	for (unsigned int i = 0; i < read_positions.size(); i++)
	{
		if (read_positions[i][1] < nparmin)
		{
			nparmin = read_positions[i][1];
		}
	}
	cout << "nparmin is " << nparmin <<endl;
	int nparmax = read_positions[0][1];
	for (unsigned int i = 0; i < read_positions.size(); i++)
	{
		if (read_positions[i][1] > nparmax)
		{
			nparmax = read_positions[i][1];
		}
	}
	cout << "nparmax is " << nparmax <<endl;
	cout << "1st halo is " << read_positions[0][0] <<" "<<read_positions[0][1]<<" "
			<<read_positions[0][2]<<" "<<read_positions[0][3]<<" "<<read_positions[0][4]<<endl;
	cout << "2nd halo is " << read_positions[1][0] <<" "<<read_positions[1][1]<<" "
			<<read_positions[1][2]<<" "<<read_positions[1][3]<<" "<<read_positions[1][4]<<endl;
	cout << "...." << endl;
	cout << length/(5*sizeof(float)) <<"th halo is " << read_positions[length/(5*sizeof(float))-1][0] <<" "
			<<read_positions[length/(5*sizeof(float))-1][1]<<" "<<read_positions[length/(5*sizeof(float))-1][2]<<" "
			<<read_positions[length/(5*sizeof(float))-1][3]<<" "<<read_positions[length/(5*sizeof(float))-1][4]<<endl;
}
