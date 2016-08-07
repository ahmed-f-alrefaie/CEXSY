//##include "MultiGPUManager.h"
#include "../common/BaseManager.h"
#include "../Readers/StateReader.h"
#include "../Readers/Exomol/ExomolStateReader.h"
#include "../common/Input.h"
#include "../common/Timer.h"
#include "../Workers/CPU/ThreadWorker.h"
#include "../Broadeners/ExomolBroadeners.h"
#pragma once


class BaseProfile{

	private:
		BaseWorker* worker;
		StateReader* state_reader;
		Input* input;
		double* freq;
		double* intens;
		double dfreq;
		int num_intens;
		double lorentz_cutoff;
		double* h_energies; 
		double* h_nu;
		double * h_aif;
		int* h_gns;
		int* h_Ji;
		int Npoints;
		
		size_t Ntrans;
		size_t num_trans_fit;
		double start_nu,end_nu,min_nu,max_nu,max_hw;
		
		
		
		std::vector< std::vector<double> > gamma_grid;
		
		vector<string> transition_files;
		vector<BroadenerInput> broadener_files;
		vector<TemperaturePressure> TP_grid;
		
		double m_beta;
		double total_intens;
		
		
		std::vector<Broadeners* > broadeners;
		
		void ComputeTotalIntensity();
		void SetupBroadening();
		void PerformCrossSection();
	protected:

	public:
		BaseProfile(Input* pInput);
		~BaseProfile();
		void Initialize();
		void ExecuteCrossSection();
		void OutputProfile();
};
