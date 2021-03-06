#include <cmath>
#include <cstdlib>
#include "../common/Input.h"
#include "../common/BaseManager.h"
#pragma once

class BaseWorker : public BaseManager{

protected:
		size_t N_trans;
		double m_lorentz_cutoff;
		double m_dfreq;
		int m_Npoints;
		int m_grid_size;
		std::vector<std::vector<double> > m_gammaL;
		std::vector<double> m_gammaD;
		std::vector<double> m_partitions;
		std::vector<TemperaturePressure> m_grid;
		BaseWorker(std::vector<TemperaturePressure> grid,std::vector<double> pPartitions);
		
public:
		virtual void InitializeLorentzian(std::vector< std::vector<double> > pgammaL)=0;
		virtual void InitializeDoppler(double meanmass)=0;
		virtual void InitializeVectors(int Npoints)=0;
		//virtual void InitializeVectors(int Npoints,int N_intens)=0;
		virtual void InitializeConstants(double dfreq,double lorentz_cutoff)=0;
		
		virtual void PrepareForComputation()=0;
		
		virtual void TransferFreq(double* h_freq,double* h_intens,int N)=0;
		//virtual void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma,double * h_n)=0;
		virtual void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,int* gamma)=0;
		

		virtual void ExecuteCrossSection(int N_ener)=0;
		virtual void TransferResults(double* h_freq,double* h_intens,int N,int TP)=0;
		virtual void Cleanup()=0;
		virtual bool ReadyForWork()=0;
		size_t GetNtrans();
		

};

