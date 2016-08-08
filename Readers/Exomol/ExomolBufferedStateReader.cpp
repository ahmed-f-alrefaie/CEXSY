#include "ExomolBufferedStateReader.h"
#include "../../common/Util.h"
#include "../../common/BaseManager.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <algorithm>
using namespace std;


ExomolBufferedStateReader::ExomolBufferedStateReader(std::string pFilename,std::vector<TemperaturePressure> grid,std::vector<std::string> filename_list) : ExomolStateReader(pFilename,grid){
	printf("Using buffered read\n");
	size_t max_memory = 0;	
	for(int i = 0; i < filename_list.size(); i++){
	
		max_memory = std::max(max_memory,GetFilenameSize(filename_list[i]));
		
	}
	
	printf("Maximum memory required is %zu bytes\n",max_memory);
	
	
	BaseManager::TrackGlobalMemory(max_memory);


}

ExomolBufferedStateReader::~ExomolBufferedStateReader(){
	CloseFile();
}

bool ExomolBufferedStateReader::OpenFile(std::string pFilename){
	printf("Buffering read...\n");
	file_stream.open(pFilename.c_str());
	if(!file_stream.is_open())
		return false;
	//stream.open(pFilename.c_str());
	buffer << file_stream.rdbuf();
	stream = &buffer;
	printf("...done!\n");
	file_stream.close();
	return true;
	
	
}
bool ExomolBufferedStateReader::CloseFile(){
	//fclose(trans_file);
	//if(!use_compression){
	buffer.clear();
	//if(stream.is_open())
	//	stream.close();
	return true;
	//}else{
	//	close_transition_file_bz2_();
	//	use_compression = false;
	//}
}


