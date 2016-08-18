#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <map>
#include "profiles.h"


class VoigtKampff{

public:



	VoigtKampff(std::vector<double> pGammaL,double pGammaD,double pRes,double pLorenzCutoff);
	

	//void Initialize()

	double ComputeVoigt(double dfreq,double gammaL);
	double ComputeDoppler(double dfreq,double gammaD, double nu);
	double ComputeVoigt(double dfreq,int gammaL,double nu);
	double HumlicekTest(double dfreq,double gammaL,double nu);
	void ComputeVoigtVectorized(const double* __restrict freq,double* __restrict intens,const double abscoef,const int ib,const int ie,const int gammaL,const double start_nu,const double nu);
	void DoDopplerVectorized(const double* __restrict freq, double* __restrict intens, const double abscoef, const int ib, const int ie, const double gammaD,const double nu);

private:
	int middle_point;
	void DoVectorized(double* __restrict intens,const std::vector<double> & gammaVoigt,const double abscoef,const int start,const int end);
	int m_Npoints;
	double m_res;
	double m_lorentz_cutoff;
	double m_gammaD;
	double m_mag;
	std::vector< std::vector<double> > m_voigt_grid;
	std::map<double,int> m_gamma_mapping;
	std::vector<double> m_gammaL;
};




