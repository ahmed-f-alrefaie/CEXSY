#include "../StateReader.h"
#include "exomol_objects.h"
#include "exomol_functions.h"
#include "../../common/Util.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
//#include <bzlib.h>
#pragma once

#define BUFFER_SIZE 100241

class ExomolStateReader : public StateReader{
	protected:
		//FILE* trans_file;
		std::ifstream file_stream;
		std::istream * stream;
		exomol_states exomol;
		//bz2_stream zip_stream;
		//bool use_broadeners;
		//bool use_compression;
		//double m_pressure;
		//double default_gam;
		//double default_n;
		//bool GetGammaN(state* state_i,state* state_f,double & gamma, double & n);
		//void InitializeBroadener(std::string filename);
	public:
		ExomolStateReader(std::string pFilename,std::vector<TemperaturePressure> grid);
		~ExomolStateReader();
		virtual bool OpenFile(std::string pFilename);
		virtual bool CloseFile();
		//virtual void AddBroadeners(Broadeners* broad)=0;
		bool ReadNextState(double & nu,int & gns,double & e_i, double & aif, int & Ji);
};
