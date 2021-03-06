#include <string>
#include <vector>

#include "../common/Input.h"
#pragma once


class StateReader{
	protected:
		std::vector<double> partitions;	
		//std::vector<Broadeners*> m_broadeners;
		//std::vector<std::vector<double> > m_gammaL;
		std::vector<TemperaturePressure> m_grid;
		int maxJ;
		bool open_file;	
	public: 
		StateReader(std::vector<TemperaturePressure> grid): m_grid(grid){};
		~StateReader(){};
		virtual bool OpenFile(std::string pFilename)=0;
		virtual bool CloseFile()=0;
		//virtual void AddBroadeners(Broadeners* broad)=0;
		virtual bool ReadNextState(double & nu,int & gns,double & e_i, double & aif, int & Ji)=0;
		std::vector<double> GetPartitions(){return partitions;};
		//std::vector< std::vector<double> > GetGammaL(){return gammaL;};
		int GetMaxJ(){return maxJ;};
};
		
		
