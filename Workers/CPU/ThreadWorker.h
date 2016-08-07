#include <cmath>
#include <thread>
#include "../BaseWorker.h"
#include <vector>
#include "../../common/defines.h"
#include "../../compute/CPU/VoigtKampff.h"
#include "../../common/Input.h"

#pragma once

class ThreadWorker : public BaseWorker{

private:
		bool* threads_done;
		bool master_thread_done;
		int working_threads;
		int total_threads;
		double m_dfreq;
		std::vector<std::thread> workers;
		std::vector<double**> t_intens;
		double* g_freq;
		double* g_energies;
		double* g_nu;
		double* g_aif;
		int* g_gns;
		int * g_gamma;
		std::vector<double*> g_tp_intens;
		std::vector<double> m_beta;
		std::vector<VoigtKampff*> voigt_calc;
		std::mutex g_intens_mutex;	
		double start_nu;	
		void DistributeWork();
		void JoinAllThreads();
		void ComputeVoigt(int id,int Nener);
		//bool ReadyForThreadWork();		

public:
		ThreadWorker(std::vector<TemperaturePressure> grid,std::vector<double> pPartitions,int num_threads);
		void InitializeLorentzian(std::vector< std::vector<double> > pgammaL);
		void InitializeDoppler(double meanmass);
		void InitializeVectors(int Npoints);
		//virtual void InitializeVectors(int Npoints,int N_intens)=0;
		void InitializeConstants(double dfreq,double lorentz_cutoff);
		
		void PrepareForComputation();
		
		void TransferFreq(double* h_freq,double* h_intens,int N);
		//virtual void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma,double * h_n)=0;
		void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,int* gamma);
		

		void ExecuteCrossSection(int N_ener);
		void TransferResults(double* h_freq,double* h_intens,int N,int TP);
		void Cleanup();
		bool ReadyForWork();

		

};
