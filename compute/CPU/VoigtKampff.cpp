#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include "VoigtKampff.h"
#include <cstdlib>
#include <cstdio>
const double PI = 3.14159265359;
const double SQRTLN2 = 0.832554611;
const double ISQRTPI = 0.564189584; 

const double DISTANCE_MAGIC_NUMBER = 1.9;
const double REFERNECE_NU = 1.0;
double CubicInterpolate(
   double y0,double y1,
   double y2,double y3,
   double mu)
{
   double a0,a1,a2,a3,mu2;

   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;

   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
};

VoigtKampff::VoigtKampff(std::vector<double> pGammaL,double pGammaD,double pRes,double pLorenzCutoff){

	m_Npoints = 2*pLorenzCutoff/pRes + 1;
	middle_point = m_Npoints/2;
	m_res=pRes;
	m_gammaD = pGammaD;
	m_lorentz_cutoff = pLorenzCutoff;
	
	std::cout<<"Initializing Voigt....."<<std::endl;
	std::cout<<"Npoints: "<<m_Npoints<<" middle: "<<middle_point<<"Humlecik points: "<<DISTANCE_MAGIC_NUMBER/m_res<<" range : "<< 0.0 << " - "<< double(m_Npoints)*m_res<<std::endl;
	m_gammaL = pGammaL;
	double x,y;
	double nu = 0.0;
	for(int i = 0; i < pGammaL.size(); i++){
		nu = -m_lorentz_cutoff;
		double gammaL = pGammaL[i];
		m_voigt_grid.push_back(std::vector<double>());	
		m_gamma_mapping[gammaL]=i;
	
		for(int j = 0; j < m_Npoints; j++){	
			
			x = SQRTLN2*fabs(nu)/(pGammaD*REFERNECE_NU);
			y = SQRTLN2*gammaL/(pGammaD*REFERNECE_NU);
		
			m_voigt_grid.back().push_back(humlic(x,y)*SQRTLN2*ISQRTPI/(pGammaD*REFERNECE_NU));
		
			nu+=m_res;
		}
		//std::cout<<"size:"<<m_voigt_grid.back().size()<<std::endl;
	}

	std::cout<<"Done!!!"<<std::endl;



}

double VoigtKampff::ComputeVoigt(double dfreq,int gammaL,double nu){
	
	double dfreq_=fabs(dfreq);
	if(dfreq_ > m_lorentz_cutoff)
		return 0.0;
	int point = dfreq/m_res + middle_point;
	//printf("point:%d\n",point);
	//point = std::max(point,0c);
	//point = std::min(point,m_Npoints);
	
		
		//TODO: Use intepolation 
		/*int i1 = std::max(point-2,0);
		int i2 = std::max(point-1,0);
		int i3 =  std::min(point+2,m_Npoints);
		int i4 = std::min(point+2,m_Npoints);
		double p1=m_voigt_grid[gammaL][i1];
		double p2=m_voigt_grid[gammaL][i2];
		double p3=m_voigt_grid[gammaL][i3];
		double p4=m_voigt_grid[gammaL][i4];
		return CubicInterpolate(p1,p2,p3,p4,dfreq);///ComputeDoppler(dfreq,m_gammaD,nu);
		*/
	//double dist = 2.0*m_gammaD*nu/m_gammaL[gammaL];
	//if(dfreq_ < dist) return HumlicekTest(dfreq,m_gammaL[gammaL],nu);
	return m_voigt_grid.at(gammaL).at(point);
		//;
	//return 0.0;
	
}	

void VoigtKampff::DoVectorized(double* __restrict intens,const std::vector<double> & gammaVoigt,const double abscoef,const int start,const int end){
	//int count = end-start;
	
	for(int i = start; i < end; i++){
		intens[i]+=gammaVoigt[i]*abscoef;
	}


}

void VoigtKampff::ComputeVoigtVectorized(const double* __restrict freq,double* __restrict intens,const double abscoef,const int ib,const int ie,const int gammaL,const double start_nu,const double nu){


	
	int center_point = (nu-start_nu)/m_res;
	int middle_shift = center_point - middle_point;
	double gammaL_ = m_gammaL[gammaL];
	
	//int dist = DISTANCE_MAGIC_NUMBER/m_res;//(m_gammaD*nu/m_gammaL[gammaL])/m_res;
	int dist;
	if (gammaL_ == 0.0)
		dist = m_Npoints;
	else
		dist = (m_gammaD*nu / m_gammaL[gammaL]) / m_res;
	int ib_rel = ib - middle_shift;
	int ie_rel = ie - middle_shift;
	
	double gammaD_ = (m_gammaD*nu);

	

	//int start_dist 
	//printf("Middle calc %d %d\n",start_dist,end_dist);
	
	
	
	double dfreq;
	//DO PURE DOPPLER
	if (gammaL_ == 0.0) {
		//USe vectorized doppler
		DoDopplerVectorized(freq, intens, abscoef, ib, ie, gammaD_, nu);



	}//OTHERWISE WE VOIGT BABEEEEEEEEEEEEEEE
	else {


		//int ib_rel = (freq[ib] - nu)/m_res  + middle_point;
		//int ie_rel = (freq[ie] - nu)/m_res  + middle_point;

		//printf("ib: %d ie: %d start:%d end:%d dist: %d center:%d nu: %12.6f\n",ib,ie,ib_rel,ie_rel,dist,center_point,nu);
		int left_start = ib_rel;
		int left_end = std::min(middle_point - dist, ie_rel);
		int right_start = std::max(middle_point + dist, ib_rel);
		int right_end = ie_rel;
		//If we have a calculation on the left
		if (ib < center_point) {


			DoVectorized(intens + ib - ib_rel, m_voigt_grid[gammaL], abscoef, left_start, left_end);
			//printf("Left calc %d %d\n",left_start,left_end);


		}
		if (ie >= center_point) {

			DoVectorized(intens + ib - ib_rel, m_voigt_grid[gammaL], abscoef, right_start, right_end);
			//printf("Right calc %d %d\n",right_start,right_end);

		}

		//Do the middle
		int start_dist = std::max(center_point - dist, ib);
		int end_dist = std::min(center_point + dist, ie);


		for (int i = start_dist; i < end_dist; i++) {
			//int index = i + ib - ib_rel;
			//if(index < ib) continue;
			//if(index >= ie) continue;
			dfreq = freq[i] - nu;
			intens[i] += HumlicekTest(dfreq, gammaL_, nu)*abscoef;
		}
	}
	


	//int dist_hum = dist;
	
	//for(int i = dist+1; i < dist
	
	
	
	//exit(0);
	
	
	
	



}

void VoigtKampff::DoDopplerVectorized(const double * __restrict freq, double * __restrict intens, const double abscoef, const int ib, const int ie, const double gammaD, const double nu)
{
	double fact = SQRTLN2 / gammaD;
	double x0 = fact*m_res*0.5;
	
	for (int i = ib; i < ie; i++) {
	//	double dfreq = freq[i] - nu;
		//double xp = fact*(dfreq)+x0;
		//double xm = fact*(dfreq)-x0;
	//	double de = erf(xp) - erf(xm);
		intens[i] += abscoef*0.5/m_res*(erf(fact*(freq[i] - nu) + x0)
			- erf(fact*(freq[i] - nu) - x0));
			//intens[i] += exp(i);
			/*intens[i] +=
				//DE
				( erf(   fact*(freq[i] - nu) + x0   )
				- erf(  fact*(freq[i] - nu) - x0  ) )

				*abscoef*0.5 /
				(freq[i] - nu);
		}*/
	}
}


double VoigtKampff::ComputeDoppler(double dfreq,double gammaD,double nu){
	
    	 //double dop = 0.03912;//(3155.87/nu)*4.523*(gammaD*nu);
         /*  
         double xp = SQRTLN2/dop*(freq+ m_res/2.0 - nu);
         double xm = SQRTLN2/dop*(freq- m_res/2.0 - nu);
             
         double de = erf(xp)-erf(xm);
        
         return 0.5/m_res*de;
	*/
	
	//double alpha = -1.0/(dop);
	//double de = exp(alpha*dfreq_*dfreq_);
	//double de = (dfreq*dfreq + (dop)*(dop))/(dop*dop*nu) ;
	//return de;
	return (-3.295e-05*dfreq*dfreq*dfreq*dfreq + 2.916e-08*dfreq*dfreq*dfreq 
		+ 0.2924*dfreq*dfreq -6.335e-06*dfreq + 0.02171);


}


	//void Initialize()
double VoigtKampff::ComputeVoigt(double dfreq,double gammaL){
	
	double dfreq_=fabs(dfreq);
	if(dfreq_ > m_lorentz_cutoff)
		return 0.0;
	int point = fabs(dfreq_)/m_res;
	//point = std::max(point,0);
	
	
	if(m_gamma_mapping.count(gammaL) > 0){
		int gamma_idx = m_gamma_mapping[gammaL];
		
		//TODO: Use intepolation 
		/*int i1 = std::max(point-1,0);
		int i2 = point;
		int i3 =  std::min(point+2,m_Npoints);
		int i4 = std::min(point+2,m_Npoints);
		double p1=m_voigt_grid[gamma_idx][i1];
		double p2=m_voigt_grid[gamma_idx][i2];
		double p3=m_voigt_grid[gamma_idx][i3];
		double p4=m_voigt_grid[gamma_idx][i4];
		return CubicInterpolate(p1,p2,p3,p4,dfreq);
		*/
		
		return m_voigt_grid[gamma_idx][point];
		//;
	}
	return 0.0;
	
}



double VoigtKampff::HumlicekTest(double dfreq,double gammaL,double nu){
	double x,y;
	double dfreq_=fabs(dfreq);
	double gammaD = m_gammaD*nu;
	if(dfreq_ > m_lorentz_cutoff)
		return 0.0;
	
	x = SQRTLN2*dfreq_/(gammaD);
	y = SQRTLN2*gammaL/gammaD;
	return humlic(x,y)*SQRTLN2*ISQRTPI/gammaD;
}



