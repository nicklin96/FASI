
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
#include "FPGA.h"
#include <string>
#include "krnl_intersection.cpp"

Host::~Host(){
    for(auto e : sm_events){
        delete e;
    }
}

void Host::init(string binary_file){
    cl_int err;
    int device_id = 0;
    cout << "initial " << binary_file << std::endl;

    for(device_id=0; device_id<1;device_id++){
        auto devices = xcl::get_xil_devices();
        auto tmp_device = devices[device_id];

        cl::Context context(tmp_device, NULL, NULL, NULL, &err);
        
        auto fileBuf = xcl::read_binary_file(binary_file);
        cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

        std::swap(devices[0],devices[device_id]);
        devices.resize(1);
        cl::Program program(context, devices, bins, NULL, &err);

        cl::CommandQueue q(context, tmp_device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);

        context_list.push_back(context);
        device_list.push_back(tmp_device);
        program_list.push_back(program);
        q_list.push_back(q);
        sm_events.emplace_back();
    }
     
    std::cout << "size of context " << context_list.size() << std::endl;
    std::cout << "size of device " << device_list.size() << std::endl;
    std::cout << "size of program " << program_list.size() << std::endl;
    std::cout << "size of q " << q_list.size() << std::endl;
    std::cout << "size of events " << ll_events.size() << " " << hl_events.size() << " " << hh_events.size() << std::endl;
}

FPGAControl::FPGAControl(LPCSR* graph, int pe):graph_data(graph),pe_num(pe){
    for(int j = 0; j < pe; j++){
        final_output.emplace_back();
        for(int k = 0; k < OUT_LEN; k++){
            final_output[j].push_back(0);
        }
    }


}

FPGAControl::~FPGAControl(){

}

uint64_t FPGAControl::get_duration_ns(const cl::Event *event) {
    cl_int err;
    uint64_t nstimestart, nstimeend;
    OCL_CHECK(err,
                err = event->getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_START, &nstimestart));
    OCL_CHECK(err,
                err = event->getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_END, &nstimeend));
    return (nstimeend - nstimestart) ;
}

void FPGAControl::runOnFPGA(Host* host, 
                            unsigned int** query_plan, 
                            vector<unsigned>& join_order,
                            Graph& query,
                            vector<unsigned>* candidate,
                            vector<char>* candidate_bitmap
                            ){

    std::vector<cl::Kernel> sm_kernels(pe_num);
    cl_int sm_err;

    for(int i = 0 ; i < pe_num ; i++){
        OCL_CHECK(sm_err,sm_kernels[i] = cl::Kernel(host->program_list[0], "multiIntersect", &hh_err));
    }

    cl_mem_flags read = CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY;
    cl_mem_flags write = CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE;

    std::cout << "*** init data buffers ***" << std::endl;

    cl::Buffer *neighborhood_buffer = new cl::Buffer(host->context_list[0], read, sizeof(unsigned long long) * this->graph_data->neighborhood_structure.size(),
    this->graph_data->neighborhood_structure.data(),
    &sm_err);

    cl::Buffer *index_buffer = new cl::Buffer(host->context_list[0], read, sizeof(id_t) * this->graph_data->index_list.size(),
    this->graph_data->index_list.data(),
    &sm_err);

    //load LPCSR
    std::vector<cl::Buffer*> offset_buffers(6,nullptr);
    std::vector<cl::Buffer*> adjacency_buffers(6,nullptr);
    std::vector<cl::Buffer*> candidate_buffers(6,nullptr);
    std::vector<cl::Buffer*> candidate_bitmap_buffers(6,nullptr);
    
    int valid_vertex_count = 0;
    int query_plan_count = 0;
    for(int valid_vertex_count = 0; valid_vertex_count < 6;valid_vertex_count++){
        if(query.vertex_list[join_order[valid_vertex_count]]->degree==1){
            break;
        }
        query_plan_count += query_plan[i][0];
        offset_buffers[valid_vertex_count] = new cl::Buffer(host->context_list[0], read, 
                                             sizeof(id_t) * this->graph_data->offset_list[join_order[valid_vertex_count]].size(),
                                             this->graph_data->offset_list[join_order[valid_vertex_count]].data(),
                                             &sm_err);

        adjacency_buffers[valid_vertex_count] = new cl::Buffer(host->context_list[0], read, 
                                                sizeof(id_t) * this->graph_data->adjacency_list[join_order[valid_vertex_count]].size(),
                                                this->graph_data->adjacency_list[join_order[valid_vertex_count]].data(),
                                                &sm_err);

        candidate_buffers[valid_vertex_count] = new cl::Buffer(host->context_list[0], read, 
                                                sizeof(id_t) * candidate[join_order[valid_vertex_count]].size(),
                                                candidate[join_order[valid_vertex_count]].data(),
                                                &sm_err);

        candidate_bitmap_buffers[valid_vertex_count] = new cl::Buffer(host->context_list[0], read, 
                                                sizeof(char) * candidate_bitmap[join_order[valid_vertex_count]].size(),
                                                candidate_bitmap[join_order[valid_vertex_count]].data(),
                                                &sm_err);

		
    }

	for(int k = 0; k < this->pe_num; k++){
		final_output.emplace_back();
		output_buffers.append(new cl::Buffer(host->context_list[0], write, 
                                                sizeof(id_t) * final_output[k].size(),
                                                final_output[k].data(),
                                                &sm_err));  
	}

    cl::Buffer *plan_buffer = new cl::Buffer(host->context_list[0], read, sizeof(id_t) * query_plan_count,
    query_plan,
    &sm_err);

    cl::Buffer *order_buffer = new cl::Buffer(host->context_list[0], read, sizeof(id_t) * valid_vertex_count,
    join_order.data(),
    &sm_err);

    int arg_count = 0;

    for(int i = 0; i < this->pe_num; i++){
        arg_count = 0;
        OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *neighborhood_buffer));
        OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *index_buffer));
		for(int j = 0; j < valid_vertex_count; j++){
			OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *offset_buffers[j]));
			OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *adjacency_buffers[j]));
			OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *candidate_buffers[j]));
			OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *candidate_bitmap_buffers[j]));
		}
		OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *plan_buffer));
		OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, *order_buffer));
		OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, valid_vertex_count));
		OCL_CHECK(sm_err, sm_err = sm_kernels[i].setArg(arg_count++, query_plan_count));
		
    }

    OCL_CHECK(sm_err,sm_err = host->q_list[0].enqueueMigrateMemObjects(
                        {*neighborhood_buffer,
                        *index_buffer,
                        *plan_buffer,
                        *order_buffer},
                        0));

    for(int i = 0; i < valid_vertex_count; i++){
        //std::cout << "writing hh task no." << i << std::endl;
        OCL_CHECK(sm_err,sm_err = host->q_list[0].enqueueMigrateMemObjects(
                                    {
                                    *offset_buffers[i],
                                    *adjacency_buffers[i]},
                                    0));
    }
	for(int i = 0; i < this->pe_num; i++){
		OCL_CHECK(sm_err,sm_err = host->q_list[0].enqueueMigrateMemObjects(
                                    {
                                    *output_buffers[i],},
                                    0));

	}
    
    host->q_list[0].finish();

    struct timeval start{}, end{};
    vector<cl::Event*> sm_events(hh_kernel_num);

    for(int i = 0; i < this->pe_num; i++){
        sm_events[i] = new cl::Event();
    }

    gettimeofday(&start,nullptr);
    for(int i = 0; i < this->pe_num; i++){
        OCL_CHECK(sm_err, sm_err = host->q_list[0].enqueueTask(sm_kernels[i], NULL, sm_events[i]));
    }
    host->q_list[0].finish();

    gettimeofday(&end,nullptr);
    
    double time_cost =(double) (end.tv_sec - start.tv_sec) * 1000 +
                        (end.tv_usec - start.tv_usec) / 1000.0;
    //std::cout << "Time cost = " << time_cost << std::endl;
    uint64_t kernel_time = 0;

    std::cout << "*** getting results ***" << endl;
    for(int i = 0; i < this->pe_num; i++){
        OCL_CHECK(sm_err,sm_err = host->q_list[0].enqueueMigrateMemObjects({*output_buffers[i]}, CL_MIGRATE_MEM_OBJECT_HOST));
    }

    host->q_list[0].finish();

    gettimeofday(&end,nullptr);
    time_cost =(double) (end.tv_sec - start.tv_sec) * 1000 +
                                (end.tv_usec - start.tv_usec) / 1000.0;

    //std::cout << "Kernel time = " << kernel_time << std::endl;
    std::cout << "Time cost = " << time_cost << std::endl;

}

delete neighborhood_buffer;
delete index_buffer;
delete plan_buffer;
delete order_buffer;
for(int i = 0; i < valid_vertex_count; i++){
	delete offset_buffers;
	delete adjacency_buffers;
	delete candidate_buffers;
	delete candidate_bitmap_buffers;
}
}


#endif
