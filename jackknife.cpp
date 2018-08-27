#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>

const unsigned int numsub = 64;

int main(int argc, char *argv[]) {
	for (unsigned int i = 0; i < numsub; ++i){
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
		std::cout << "Erasing " << counts[num/numsub*i] << " " << counts[num/numsub*i+1] << " " << counts[num/numsub*i+2]
			 	<< " .... " << counts[num/numsub*(i+1)-3] << " " << counts[num/numsub*(i+1)-2] 
				<< " " << counts[num/numsub*(i+1)-1] << std::endl;
		counts.erase(counts.begin()+num/numsub*i, counts.begin()+num/numsub*(i+1));
		std::cout << "Jackknife subsample " << i << " has " << counts.size() <<" counts." << std::endl;

		unsigned int max_count = *std::max_element(counts.begin(), counts.end());
		std::cout << "Jackknife subsample " << i << " maximum count is " << max_count << std::endl;

		std::ofstream countfile;
		std::string filename = "freq_jackknife" + std::to_string(numsub) +"_"+std::to_string(i)+".txt";
		countfile.open(filename);
		for (unsigned int n = 0; n < max_count +1; ++n){
			countfile << std::count(counts.begin(), counts.end(), n) << " ";
		}
		countfile.close();
		std::cout << "Jackknife subsample " << i << " count frequency saved in file " << filename << "." << std::endl;
	}
	return 0;
}
