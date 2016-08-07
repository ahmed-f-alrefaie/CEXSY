#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Util.h"
#pragma once
using namespace std;

enum ProfileType{
	GAUSSIAN,
	DOPPLER,
	VOIGT
	
};

enum FileType{
	HITRAN_TYPE,
	EXOMOL_TYPE
};

struct TemperaturePressure{
	double temperature;
	double pressure;
	
};

struct BroadenerInput{
	std::string filename;
	double mixture;
};

class Input{
	private:
		double half_width;
		double mean_mass;
		double lorentz;
		double temperature;
		double partition;
		double pressure;
		double nu_start;
		double nu_end;
		double cutoff;
		vector<TemperaturePressure> temp_press_grid;
		vector<double> temperatures;
		vector<double> partitions;
		int Npoints;
		int num_files;
		double gamma_air;
		double n_air;
		double hitran_mixture_air;
		size_t memory;
		int num_threads;
		int max_points;
		FileType which_file;
		vector<string> trans_files;
		vector<BroadenerInput> broadener_files;
		string state_file;
		string hitran_file;
		ProfileType StringToProfile(string & pString);
		//istream file_stream;
		
		
	public:
		Input();
		~Input();
		ProfileType profile;
		void ReadStream();
		void ReadInput();
		void ReadInputII();
		double GetHalfWidth(){return half_width;};
		double GetMeanMass(){return mean_mass;};
		double GetTemperature(){return temperature;};
		double GetNuStart(){return nu_start;};
		double GetNuEnd(){return nu_end;};
		double GetPartition(){return partition;};
		double GetPressure(){return pressure;}
		double GetHitranMixture(){return hitran_mixture_air;}
		double GetMaxPoints(){return max_points;}
		int GetNumThreads(){return num_threads;}
		size_t GetMemory(){return memory;}
		double GetGamma(){return gamma_air;};
		double GetGammaN(){return n_air;};
		
		double GetLorentzCutoff(){return cutoff;};
		
		FileType GetFileType(){return which_file;}

		int GetNpoints(){return Npoints;};
		string GetStateFile(){return state_file;};
		vector<TemperaturePressure> GetTemperaturePressureGrid(){return temp_press_grid;}
		vector<double> GetTemperatures(){return temperatures;}
		vector<string> GetTransFiles(){return trans_files;};
		vector<BroadenerInput> GetBroadeners(){return broadener_files;}
		








};
