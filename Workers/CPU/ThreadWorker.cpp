#include "ThreadWorker.h"
#include "../../common/Timer.h"
#include "fmath.hpp"
#include <algorithm>
/*		bool* threads_done;
		bool master_thread_done;
		int total_threads;
		std::vector<std::thread> workers;
		std::thread* master_thread;
		double** t_intens;
		double* g_freq;
		double* g_intens;
		double* g_energies;
		double* g_nu;
		double* g_aif;
		int * g_Ji;
		double * g_n;
		int* g_gns;
		std::vector<double*> g_tp_intens;
		std::vector<VoigtKampff*> voigt_calc;
		std::mutex g_intens_mutex;	
		double start_nu;	
		void DistributeWork();
		void JoinAllThreads();
		//bool ReadyForThreadWork();		
*/
//public:


ThreadWorker::ThreadWorker(std::vector<TemperaturePressure> grid,std::vector<double> pPartitions, int num_threads) : BaseWorker(grid,pPartitions), total_threads(num_threads) {
	printf("Utilizing ThreadWorker with %i cores\n",total_threads);
	//Initialize memory

	//Initialize our threads
	//workers = new std::thread[total_threads];
	//t_intens = new double*[total_threads];
	//threads_done = new bool[total_threads];
	//master_thread_done=true;
	for(int i = 0; i < total_threads; i++)
		threads_done.push_back(true);
	//master_thread = new std::thread(&ThreadWorker::DistributeWork,this);



}

void ThreadWorker::InitializeLorentzian(std::vector< std::vector<double> > pgammaL){

	//Simple
	m_gammaL = pgammaL;


}
void ThreadWorker::InitializeDoppler(double meanmass){
	for(int i = 0; i < m_grid.size(); i++){
		double temperature = m_grid[i].temperature;
		double gammaD = sqrt(2.0*BOLTZ*LN2*temperature*AVOGNO/((meanmass)))/VELLGT;
		m_gammaD.push_back(gammaD);
		//
		voigt_calc.push_back(new VoigtKampff(m_gammaL[i],m_gammaD.back(),m_dfreq,m_lorentz_cutoff));
		m_beta.push_back(SEC_RAD/temperature);
	}

}
void ThreadWorker::InitializeVectors(int Npoints){

	m_Npoints = Npoints;
	size_t memory_for_work = BaseManager::GetAvailableGlobalMemory()/2l;

	InitializeMemory(memory_for_work); //Half the memory we specify
	//Remove it from global memory
	 BaseManager::TrackGlobalMemory(memory_for_work);
	
	
	g_freq = new double[Npoints];
	TrackMemory(size_t(Npoints)*sizeof(double));
	for(int i = 0; i < m_grid.size(); i++){
		g_tp_intens.push_back(NULL);
		g_tp_intens.back() = new double[Npoints];
		TrackMemory(sizeof(double)*size_t(Npoints));
	}
	
	
	for(int i = 0; i < m_grid_size; i++){
		t_intens.push_back(new double*[total_threads]);
	//We create seperate intensities for each thread
		for(int j = 0; j < total_threads; j++){
			t_intens.back()[j] = new double[Npoints];
			TrackMemory(size_t(Npoints)*sizeof(double));
			for(int k = 0; k < Npoints; k++)
				t_intens.back()[j][k]=0.0;
		}
	}	
	N_trans = GetAvailableMemory()/(3l*sizeof(double)+2l*sizeof(int));
	printf("Thread Worker = %zu\n",N_trans);
	g_energies = new double[N_trans];
	TrackMemory(sizeof(double)*size_t(N_trans));
	g_nu = new double[N_trans];
	TrackMemory(sizeof(double)*size_t(N_trans));
	g_aif = new double[N_trans];
	TrackMemory(sizeof(double)*size_t(N_trans));
	g_gns = new int[N_trans];	
	TrackMemory(sizeof(int)*size_t(N_trans));	
	g_gamma = new int[N_trans];	
	TrackMemory(sizeof(int)*size_t(N_trans));	

}
		//virtual void InitializeVectors(int Npoints,int N_intens)=0;
void ThreadWorker::InitializeConstants(double dfreq,double lorentz_cutoff){
	m_lorentz_cutoff = lorentz_cutoff;
	m_dfreq = dfreq;

}
		
void ThreadWorker::PrepareForComputation(){



}
		
void ThreadWorker::TransferFreq(double* h_freq,double* h_intens,int N){
	memcpy(g_freq,h_freq,size_t(N)*sizeof(double));
	start_nu = g_freq[0];
}
		//virtual void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma,double * h_n)=0;
void ThreadWorker::TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,int* h_gamma){
	JoinAllThreads();	
	memcpy(g_energies,h_energies,size_t(Nener)*sizeof(double));
	memcpy(g_nu,h_nu,size_t(Nener)*sizeof(double));
	memcpy(g_aif,h_aif,size_t(Nener)*sizeof(double));
	memcpy(g_gns,h_gns,size_t(Nener)*sizeof(int));
	memcpy(g_gamma,h_gamma,size_t(Nener)*sizeof(int));

}
		

void ThreadWorker::ExecuteCrossSection(int N_ener){
		int threads_to_work = total_threads;
		if(N_ener==N_trans) threads_to_work--;
		working_threads = threads_to_work;
		for(int i = 0; i < working_threads; i++){
			// printf("PUSH");
			//fflush(0);
			while(!threads_done[i]){}; //Wait for our threads to finish
			 threads_done[i]=false;
			// printf("HARD");	
			//fflush(0);		
			 workers.push_back(std::thread(&ThreadWorker::ComputeVoigt,this,i,N_ener));
		}

}
void ThreadWorker::TransferResults(double* h_freq,double* h_intens,int N,int TP){
	printf("Wait for the jobs to finish\n");
	Timer::getInstance().StartTimer("Joining Threads");
	JoinAllThreads(); // If we aren't done wait
	Timer::getInstance().EndTimer("Joining Threads");
	printf("Transferring from all threads into host!\n");
	Timer::getInstance().StartTimer("Applying Results");
	for(int i = 0; i < total_threads; i++){
		for(int j = 0; j < N; j++){
			h_intens[j]+=t_intens[TP][i][j];
		}
	}
	Timer::getInstance().EndTimer("Applying Results");
}



void ThreadWorker::Cleanup(){};
bool ThreadWorker::ReadyForWork(){
	bool status=true;
	for(int i = 0; i < total_threads; i++){
			status &= threads_done[i];
	}
	return status;

}

void ThreadWorker::ComputeVoigt(int id,int Nener){
	
	int thread_id = id;
	double* intens_ptr;
	VoigtKampff* v_ptr;
	fflush(0);
	for(int i = thread_id; i < Nener; i+= working_threads){
		//int i = 
		double nu_if = g_nu[i];
		int ib =std:: max(round( ( nu_if-m_lorentz_cutoff-start_nu)/m_dfreq ),0.0);
		int ie =  std::min(round( ( nu_if+m_lorentz_cutoff-start_nu)/m_dfreq ),(double)m_Npoints);
		
		if((ib + ie)==0.0) continue;
		
		//if(thread_id==18 && i %1000 == 0)printf ("nu_if:%12.6f lorenz:%12.6f dfreq:%12.6f ib: %d, ie:%d \n",nu_if,m_lorentz_cutoff,m_dfreq,ib,ie);
		if(nu_if==0) continue;
		
		double aif = g_aif[i];
		double gns = (double)g_gns[i];
		double energy = g_energies[i];
		int gammaIdx = g_gamma[i];
		for(int temp = 0; temp < m_grid_size; temp++){
			//Which intensity are we doing?
			double beta = m_beta[temp];
			 intens_ptr = t_intens[temp][thread_id];
			 v_ptr = voigt_calc[temp];
			double abscoef= CMCOEF*aif*gns
				//*fmath::expd(-beta*energy)*(1.0-fmath::expd(-beta*nu_if))/
				*exp(-beta*energy)*(1.0-exp(-beta*nu_if))/
				(nu_if*nu_if*m_partitions[temp]);	
		
		
		 

			if(abscoef==0.0) continue;
			//TODO: Create vectorised version
			//for(int j = ib; j < ie; j++){
				//printf("freq: %12.6f j:%d ib:%d ie:%d\n",g_freq[j],j,ib,ie);
				
				//double dfreq_ = nu_if - g_freq[j];
				
				//intens_ptr[j]+=abscoef*v_ptr->ComputeVoigt(dfreq_,gammaIdx,nu_if);
				v_ptr->ComputeVoigtVectorized(g_freq,intens_ptr,abscoef,ib,ie,gammaIdx,start_nu,nu_if);
			//}
		}

		
	}
	//printf("Thread %i DONE!!!!\n",thread_id);
	//fflush(0);
	/*g_intens_mutex.lock();
	for(int i = start_idx; i < Npoints + start_idx; i++){
		g_intens[i]+= intens_ptr[i];
		intens_ptr[i] = 0.0;
	}
	g_intens_mutex.unlock();
	*/
	threads_done[thread_id]=true;	
	


}

void ThreadWorker::JoinAllThreads(){
	printf("Joining threads...");
	for(int i = 0; i < workers.size();i++){
		printf("..%d..",i);
		fflush(0);
		//while (!threads_done[i]);;
		workers[i].join();
	}
	workers.clear(); // Destroy the threads
	printf("done!\n");
}
