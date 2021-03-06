#include "Util.h"
#include <string>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <assert.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <functional>
#include <cctype>
#include <locale>
#include <cmath>
#include <ctime>

#ifdef __linux
#include <sys/time.h>
#include <unistd.h>
#elif _WIN32
#include <windows.h>
#endif

using namespace std;



bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile.is_open();
}

void destroy_arr_valid(void** ptr)
{
	if(*ptr!= NULL)
	{
		delete[] *ptr;
		*ptr = 0;
	}
}




double GenRandNumber(double LO, double HI){
return LO + (double) (rand()) /( (double) (RAND_MAX/(HI-LO)));
}

size_t GetFilenameSize(std::string name)
{
	FILE* pFile;
	pFile = fopen(name.c_str(),"r");
	size_t or_pos = ftell(pFile);
	fseek(pFile,0,SEEK_END);
	size_t file_size = ftell(pFile);
	//Return to original position
	fseek(pFile,or_pos,SEEK_SET);
	fclose(pFile);
	return file_size;
};

// trim from start
std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

std::vector<std::string> split(std::string const &input) { 
    std::istringstream buffer(input);
    std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)), 
                                 std::istream_iterator<std::string>());
    return ret;
}


void ReadFortranRecord(FILE* IO, void* data_ptr) {
        unsigned int file_record_size = 0;
        //Read the first 4 bytes
        fread(&file_record_size, 4, 1, IO);
        //Read the output to the data pointer
        fread(data_ptr, file_record_size, 1, IO);
        //Read last 4 bytes of record
        fread(&file_record_size, 4, 1, IO);
}

void wrapvalue(unsigned long int & var, unsigned long int min,unsigned long int max){if(var >= max) var=min;}


void assertdouble(double & d1, double & d2, double tol)
{
	assert(abs(d1-d2) < tol);

}

void CreateRandomVector(double** vector,size_t count){
	(*vector) = new double[count];
	for(size_t i = 0; i < count; i++){
		(*vector)[i] = GenRandNumber(1e-18, 1.0);
	}


}

void CreateZeroVector(double** vector,size_t count){
	(*vector) = new double[count];
	for(size_t i = 0; i < count; i++){
		(*vector)[i] = 0.0;
	}
}


double readd(){
	char* line_ptr;
	line_ptr=strtok(NULL," ,"); //Get the temperature
	if(line_ptr!=NULL)return strtod(line_ptr,NULL);
	else{
		printf("Read error[ readd]\n");
                exit(0);
		return -1.0;
	}
};	

int readi(){
	char* line_ptr;
	line_ptr=strtok(NULL," ,"); //Get the temperature
	if(line_ptr!=NULL)return strtol(line_ptr,NULL,0);
	else{
		printf("Read error[ readi]\n");
                exit(0);
		return -1.0;
	}
};

char* readc(){
	char* line_ptr;
	line_ptr=strtok(NULL," ,"); //Get the temperature
	return line_ptr;
};

int64 GetTimeMs64()
{

#ifdef WIN32
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	* to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	long int ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

	return ret;
#else
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	uint ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
	ret /= 1000;

	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);

	return ret;
#endif

};


unsigned long long getTotalSystemMemory()
{
	printf("Get memory\n");
#ifdef __linux 
    printf("With linux!!\n");
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#elif _WIN32
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
#endif
}

void DistributeWork(int & start, int & end,int N,int id, int workers){
	int chunk_size = N/workers;
	int remainder = N % workers;
	if(id < (workers-remainder)){
		start = id*chunk_size;

	}else{
		start = (workers-remainder)*chunk_size;
		chunk_size++;
		start += (id-(workers-remainder))*chunk_size;
	}
	end = min(start+ chunk_size,N);
}	
	//int start = thread_id*chunk_size;

