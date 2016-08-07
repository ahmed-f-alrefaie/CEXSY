#pragma once
#include <vector>
#include "../common/defines.h"

struct broadinfo{
	double gamma;
	double n;	
};

class Broadeners{

private:
	
	
	
protected:

	double m_def_gam;
	double m_def_n;
	
	double m_ref_temp;
	double m_ref_press;
	
	double m_mixture;
	
	std::vector<broadinfo> gamma;
	
	int m_maxJ;
	
	
	Broadeners(double default_gam,double default_n,double ref_temp,double ref_press,double mixture,int maxJ) : m_def_gam(default_gam), m_def_n(default_n), 
													  m_ref_temp(ref_temp), m_ref_press(ref_press), m_mixture(mixture), m_maxJ(maxJ) {};
	 

public:
	//TODO: switch to more generic quanta
	virtual double GetGamma(int Ji,double temperature,double pressure)=0;
	virtual void InitializeBroadener(std::string filename)=0;

};


