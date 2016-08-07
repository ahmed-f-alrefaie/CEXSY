#pragma once
#include "Broadeners.h"
#include <string>
#include <map>
class ExomolBroadeners : public Broadeners{

private:
	
	
	
protected:

	
	
public:
	//TODO: switch to more generic quanta
	double GetGamma(int Ji,double temperature,double pressure);
	void InitializeBroadener(std::string filename);
	ExomolBroadeners(double default_gam,double default_n,double mixture,int maxJ);
};


