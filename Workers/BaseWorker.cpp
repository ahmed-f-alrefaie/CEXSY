
#include "BaseWorker.h"
#pragma once

BaseWorker::BaseWorker(std::vector<TemperaturePressure> grid,std::vector<double> pPartitions) : BaseManager(), m_grid(grid), m_partitions(pPartitions), m_grid_size(grid.size())
{ 

}; 
//BaseWorker::~BaseWorker(){};
size_t BaseWorker::GetNtrans(){return N_trans;};
