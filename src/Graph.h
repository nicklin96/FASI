#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>
#include <unordered_map>

using namespace std;

class vertex
{
public:
	unsigned vid = 0;
	unsigned label = 0;
	unsigned degree = 0;
	unsigned partition_max_degree = 0;
	vector<vertex*> neighbor_list;
	vector<vertex*>* partition_neighbor_list;
	unsigned* partition_neighbor_degree_sum;
	vertex() {}
	~vertex() {}
};

class Graph
{
public:
	unsigned label_num = 0;
	unsigned vertex_num = 0;
	unsigned edge_num = 0;
	unsigned max_degree = 0;
	unsigned average_degree = 0;
	unordered_map<unsigned, unsigned> vertex_label_mapping;
	vector<vertex*> vertex_list;
	vector<vertex*>* label_list;
	vector<unsigned> label_index_list;
	Graph();
	~Graph();
	void init(const string& filename);
	void init(const string& filename, unordered_map<unsigned, unsigned>& mapping);
	void partition(bool flag = false);
	unsigned count(unsigned D);
	void export_data_graph(unsigned D, const string& filename);
	void export_data_graph_without_v(unsigned v, const string& filename);
	void export_vertex_label_mapping(const string& filename);
	bool has_edge(unsigned id1, unsigned id2);
};
