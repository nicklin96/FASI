#include <iostream>
#include <vector>
#include <unordered_map>
#include "Graph.h"

#define UINT_MAX 0xffffffff

using namespace std;

class LPCSR
{
public:
	unsigned label_num = 0;
	unsigned vertex_num = 0;
	unsigned port_num = 0;
	unordered_map<unsigned, unsigned> vertex_label_mapping;
	vector<unsigned> label_index_list;
	vector<unsigned long long> neighborhood_structure; // 32bits(initial address) + 32bits(neighbor labels)
	vector<unsigned> index_list;
	vector<unsigned>** candidate_filter_list;
	vector<unsigned>* offset_list;
	vector<unsigned>* adjacency_list;
	LPCSR();
	~LPCSR();
	void init(Graph G, unsigned pnum, unsigned D); //on FPGA
	void init(Graph G); //on CPU
	void export_LPCSR(const string& filename); //export LPCSR
	void import_LPCSR(const string& filename); //import LPCSR
};