#ifndef F_MATCH_FPGA_H
#define F_MATCH_FPGA_H

#include "../lib/xcl2.cpp"
#include "../lib/xcl2.hpp"
#include <vector>
#include <map>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <unordered_map>
#include <sys/time.h>
#include <string>
#include "krnl_intersection.cpp"

#ifndef END_FLAG
#define END_FLAG 4294967290
#endif

#ifndef OFFSET_LEN
#define OFFSET_LEN 25600000
#endif

#ifndef ADJ_LEN
#define ADJ_LEN 2560000
#endif


#ifndef ADJ_NUM_UP_BOUND
#define ADJ_NUM_UP_BOUND 25600
#endif

#ifndef OUT_LEN
#define OUT_LEN 64
#endif

#ifndef FLATTEN_NUM
#define FLATTEN_NUM 4
#endif

using id_t = unsigned;

class Host {
public:
    vector<cl::Context> context_list;
    vector<cl::Device> device_list;
    vector<cl::Program> program_list;
    vector<cl::CommandQueue> q_list;
    vector<cl::Event*> sm_events;
    int outputSize;
    vector<vector<id_t>> *deviceResults;

    Host(){};
    ~Host();
    void init(string binary_file);

};

class FPGAControl{
public:
	LPCSR* graph_data;
	vector<vector<id_t,aligned_allocator<id_t>>> final_output;
	int pe_num = 4;

	FPGAControl(LPCSR* graph,, int pe);
	uint64_t get_duration_ns(const cl::Event *event);
	void runOnFPGA(Host* host, 
                            unsigned int** query_plan, 
                            vector<unsigned>& join_order,
                            Graph& query,
                            vector<unsigned>* candidate,
                            vector<char>* candidate_bitmap);
};

#endif
