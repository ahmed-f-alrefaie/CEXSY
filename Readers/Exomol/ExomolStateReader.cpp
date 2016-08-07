#include "ExomolStateReader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
using namespace std;


ExomolStateReader::ExomolStateReader(std::string pFilename,std::vector<TemperaturePressure> grid) : StateReader(grid){
	printf("Using Exomol File format\n");
	initialize_states(pFilename.c_str(), &exomol);
	int N_states = exomol.N;
	
	maxJ = 0.0;
	for(int i = 0; i < N_states; i++){
		maxJ = std::max(exomol.states[i].J,maxJ);
	
	}
	std::cout<<"Max J: "<<maxJ<<std::endl;
	//Initializwe the broadeners
	//InitializeBroadener(pbroadener);
	cout<<"Temperature (K)\tPartition\n";
	for(int i = 0; i < grid.size(); i++){
		double temp= grid[i].temperature;
		double Q = compute_partition( &exomol,temp);
		partitions.push_back(Q);
		
		cout<<temp<<"\t"<<partitions.back()<<endl;
	
	
	
	
	}
	


}

ExomolStateReader::~ExomolStateReader(){
	CloseFile();
}

bool ExomolStateReader::OpenFile(std::string pFilename){
	stream = &file_stream;
	file_stream.open(pFilename.c_str());
	return file_stream.is_open();
	
	
}
bool ExomolStateReader::CloseFile(){
	//fclose(trans_file);
	//if(!use_compression){
		
	if(file_stream.is_open())
		file_stream.close();
	return true;
	//}else{
	//	close_transition_file_bz2_();
	//	use_compression = false;
	//}
}
bool ExomolStateReader::ReadNextState(double & nu,int & gns,double & e_i, double & aif, int & Ji){
	int id_f,id_i;	


	

	if ( ((*stream)>>id_f>>id_i>>aif)==false)
		return false;

	

	id_f-=1;
	id_i-=1;
	nu = exomol.states[id_f].energy-exomol.states[id_i].energy;
	e_i = exomol.states[id_i].energy;
	gns = (int)exomol.states[id_f].gns;
	Ji = exomol.states[id_i].J;
	return true;	

}



