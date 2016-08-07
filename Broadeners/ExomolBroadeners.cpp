#pragma once
#include "ExomolBroadeners.h"
#include <fstream>
#include <iostream>
#include <cmath>

	
ExomolBroadeners::ExomolBroadeners(double default_gam,double default_n,double mixture, int maxJ) : Broadeners(default_gam,default_n,296.0,1.0,mixture,maxJ){

	//Initialize the vector
	for(int i = 0; i <= maxJ; i++){
		gamma.push_back(broadinfo());
		gamma.back().gamma = default_gam;
		gamma.back().n = default_n;	
	}	
	
	
	
	






}	

void ExomolBroadeners::InitializeBroadener(std::string filename){
	std::ifstream i_broad(filename.c_str());
	if(!i_broad){
		std::cout<<"Broadener file with name "<<filename<<" not found!"<<std::endl;
		exit(0);
	}
	
	double gam;
	double n;
	std::string line;
	int Ji,Jf;
	while(getline(i_broad,line)){
		std::vector<std::string> split_line = split(line);
		Ji=atoi(trim(split_line[3]).c_str());
		gam=atof(trim(split_line[1]).c_str());
		n=atof(trim(split_line[2]).c_str());
		if(trim(split_line[0])!="a0")
			continue;
			
			
		gamma[Ji].gamma = gam;
		gamma[Ji].n = n;

	}		
	
	

}


double ExomolBroadeners::GetGamma(int Ji,double temperature,double pressure){
	
	return gamma[Ji].gamma*m_mixture*pow(m_ref_temp/temperature,gamma[Ji].n)*(pressure/m_ref_press)/BAR_TO_ATM;
	

}



