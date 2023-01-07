#include <iostream>
#include <set>
#include <map>
#include <unordered_map>
#include <queue>
#include <string>
#include <sys/time.h>
#include <omp.h>
#include <ctime>
#include "LPCSR.h"
#include <algorithm>
using namespace std;

clock_t begin_time, end_time;
struct timeval begin_tv, end_tv;
unsigned* eight_bit_sum = new unsigned[256];
unsigned long long limit_res_len = 100000000;
unsigned bitmap_size = 1024 * 64;
unsigned limit_query_node_number = 6;

bool vid_cmp(vertex* v1, vertex* v2) {
	return v1->vid < v2->vid;
}

vector<unsigned> determine_join_order(Graph& query_graph, unsigned start_node, vector<unsigned>* candidate_list, double* score)
{
	vector<unsigned> join_order;
	unsigned var_num = query_graph.vertex_num;
	set<unsigned> connective_node;
	bool* visited = new bool[var_num];
	for (unsigned i = 0; i < var_num; ++i)
		visited[i] = false;
	unsigned next_node = start_node;
	double next_node_score = score[start_node];
	for (unsigned i = 0; i < var_num; ++i)
	{
		join_order.push_back(next_node);
		visited[next_node] = true;
		for (unsigned i = 0; i < query_graph.vertex_list[next_node]->degree; ++i)
		{
			if (visited[query_graph.vertex_list[next_node]->neighbor_list[i]->vid])
				continue;
			else
				connective_node.insert(query_graph.vertex_list[next_node]->neighbor_list[i]->vid);
		}
		if (connective_node.empty())
			continue;
		next_node = *connective_node.begin();
		next_node_score = score[next_node];
		for (set<unsigned>::iterator it = connective_node.begin(); it != connective_node.end(); it++)
		{
			if (score[*it] < next_node_score)
			{
				next_node = *it;
				next_node_score = score[*it];
			}
		}
		connective_node.erase(next_node);
	}
	delete[] visited;
	return join_order;
}

unsigned** generate_query_plan(Graph& query_graph, vector<unsigned>& join_order) //query plan
{
	cout << "query_plan: " << endl;
	unsigned** query_plan = new unsigned*[join_order.size() + 1];
	unsigned width = 1 + 3 * join_order.size();
	query_plan[join_order.size()] = new unsigned[width];
	query_plan[join_order.size()][0] = width;
	width = 1;
	for (unsigned i = 0; i < join_order.size(); ++i)
	{
		if (i == 0)
		{
			cout << "the 0th iteration:";
			width += 3;
			query_plan[i] = new unsigned[width];
			query_plan[i][0] = width;
			query_plan[i][1] = join_order[i];
			query_plan[join_order.size()][1 + i * 3] = join_order[i];
			query_plan[i][2] = query_graph.vertex_list[join_order[i]]->label;
			query_plan[join_order.size()][2 + i * 3] = query_graph.vertex_list[join_order[i]]->label;
			if (query_graph.vertex_list[join_order[i]]->degree == 1)
				query_plan[i][3] = 0;
			else
			{
				query_plan[i][3] = 1;
				query_plan[join_order.size()][3 + i * 3] = 1000000000;
			}
			for (unsigned j = 1; j < query_plan[i][0]; ++j)
				cout << " " << query_plan[i][j];
			cout << endl;
			continue;
		}
		cout << "the " << i << "th iteration:";
		width += 3;
		query_plan[i] = new unsigned[width];
		query_plan[i][0] = width;
		unsigned deg = 0;
		for (unsigned j = 0; j < (width - 4) / 3; ++j)
		{
			query_plan[i][j * 3 + 1] = join_order[j];
			query_plan[i][j * 3 + 2] = query_graph.vertex_list[join_order[j]]->label;
			if (query_graph.has_edge(join_order[j], join_order[i]))
			{
				query_plan[i][j * 3 + 3] = 1;
				if (query_graph.vertex_list[join_order[i]]->degree == 1)
					query_plan[join_order.size()][3 + i * 3] = j;
				deg++;
			}
			else
				query_plan[i][j * 3 + 3] = 0;
		}
		query_plan[i][width - 3] = join_order[i];
		query_plan[join_order.size()][1 + i * 3] = join_order[i];
		query_plan[i][width - 2] = query_graph.vertex_list[join_order[i]]->label;
		query_plan[join_order.size()][2 + i * 3] = query_graph.vertex_list[join_order[i]]->label;
		if (query_graph.vertex_list[join_order[i]]->degree == 1)
		{
			query_plan[i][width - 1] = 0;
			width -= 3;
		}
		else
		{
			query_plan[i][width - 1] = deg;
			query_plan[join_order.size()][3 + i * 3] = 1000000000;
		}
		for (unsigned j = 1; j < query_plan[i][0]; ++j)
			cout << " " << query_plan[i][j];
		cout << endl;
	}
	cout << endl;
	return query_plan;
}

unsigned set_intersection(unsigned* begin_loc, unsigned* end_loc, unsigned* join_list, vector<unsigned>& adj, unsigned k, vector<unsigned>& same_label_vertices, unsigned lower_boundary, vector<char>& candidate_bitmap, vector<unsigned>& inter_result, unsigned res_begin_loc, vector<unsigned>& res_out)
{
	unsigned ret = 0;
	unsigned min_size = 1000000000;
	unsigned start_set_id;
	for (unsigned i = 0; i < k; ++i)
	{
		if (join_list[i] == 1)
		{
			if ((end_loc[i] - begin_loc[i]) < min_size)
			{
				min_size = end_loc[i] - begin_loc[i];
				start_set_id = i;
			}
		}
	}
	for (unsigned i = begin_loc[start_set_id]; i < end_loc[start_set_id]; ++i)
	{
		unsigned common = adj[i];
		bool find = true;
		for (unsigned j = 0; j < k; ++j)
		{
			bool flag = false;
			if (j == start_set_id)
				continue;
			if (join_list[j] == 0)
				continue;
			unsigned low = begin_loc[j], high = end_loc[j], middle = 0;
			while (low < high) {
				middle = (low + high) / 2;
				if (common == adj[middle]) {
					flag = true;
					break;
				}
				else if (common < adj[middle]) {
					high = middle;
				}
				else if (common > adj[middle]) {
					low = middle + 1;
				}
			}
			if (flag == false)
			{
				find = false;
				break;
			}
		}
		for (unsigned j = 0; j < same_label_vertices.size(); ++j)
		{
			if (common == same_label_vertices[j])
			{
				find = false;
				break;
			}
		}
		if (!find)
			continue;
		unsigned idx = (common - lower_boundary) / 8;
		char off = 1 << (7 - ((common - lower_boundary) % 8));
		if ((candidate_bitmap[idx] & off) != 0)
		{
			for (unsigned j = res_begin_loc; j < res_begin_loc + k; ++j)
				res_out.push_back(inter_result[j]);
			res_out.push_back(common);
			ret++;
		}
	}
	return ret;
}

vector<unsigned> sets_intersection(vector<unsigned>* adj, unsigned size)
{
	if (size == 1)
		return adj[0];
	unsigned i;
	vector<unsigned> common_adj;
	for (i = 0; i < size; ++i)
	{
		if (adj[i].size() == 0)
			continue;
		else
		{
			common_adj = adj[i];
			i++;
			break;
		}
	}
	for (; i < size; ++i)
	{
		if (adj[i].size() == 0)
			continue;
		else
		{
			vector<unsigned> temp_adj = common_adj;
			//common_adj = set_intersection(temp_adj, 0, 0, adj[i]);
		}
	}
	return common_adj;
}

void multi_join(LPCSR& csr, unsigned* query_plan, vector<unsigned>& inter_result, vector<char>& candidate_bitmap, vector<unsigned>& res_out, unsigned** isolated_node_begin_loc, unsigned** isolated_node_end_loc, bool** isolated_node_visited)
{
	struct timeval btv1, etv1;
	struct timeval btv2, etv2;
	struct timeval btv3, etv3;
	struct timeval btv4, etv4;
	double tv1 = 0;
	double tv2 = 0;
	double tv3 = 0;
	double tv4 = 0;
	unsigned len = query_plan[0];
	unsigned k = (len - 4) / 3; //the width of a result in the inter_result
	unsigned label = query_plan[len - 2];
	unsigned lower_boundary;
	unsigned batch_size = 8;
	vector<unsigned> same_label_idx;
	for (unsigned i = 0; i < len - 3; i = i + 3)
	{
		if (query_plan[i + 2] == label)
			same_label_idx.push_back(i / 3);
	}
	unsigned ret;
	if (query_plan[len - 1] == 0)
	{
		unsigned idx;
		for (unsigned i = 3; i < len - 3; i = i + 3)
		{
			if (query_plan[i] == 1)
			{
				idx = i / 3 - 1;
				lower_boundary = csr.label_index_list[query_plan[i - 1]];
				break;
			}
		}
		for (unsigned i = 0; i < inter_result.size(); i = i + k)
		{
			if (isolated_node_visited[query_plan[len - 3]][inter_result[i + idx] - lower_boundary])
				continue;
			unsigned long long ns = csr.neighborhood_structure[inter_result[i + idx]];
			if (((ns >> label) & 1) == 0)
				continue;
			unsigned start_address = ns >> 32;
			unsigned bit = label;
			while (bit >= 8)
			{
				start_address += eight_bit_sum[ns & 255];
				bit -= 8;
				ns = ns >> 8;
			}
			start_address += eight_bit_sum[ns & ((unsigned long long)pow(2, bit) - 1)];
			unsigned index = csr.index_list[start_address];
			isolated_node_begin_loc[query_plan[len - 3]][inter_result[i + idx] - lower_boundary] = csr.offset_list[label][index];
			isolated_node_end_loc[query_plan[len - 3]][inter_result[i + idx] - lower_boundary] = csr.offset_list[label][index + 1];
			isolated_node_visited[query_plan[len - 3]][inter_result[i + idx] - lower_boundary] = true;
		}
		res_out = inter_result;
	}
	else if (query_plan[len - 1] == 1)
	{
		unsigned* begin_loc = new unsigned[k];
		unsigned* end_loc = new unsigned[k];
		lower_boundary = csr.label_index_list[label];
		unsigned idx;
		for (unsigned i = 3; i < len - 3; i = i + 3)
		{
			if (query_plan[i] == 1)
			{
				idx = i / 3 - 1;
				break;
			}
		}
		for (unsigned i = 0; i < inter_result.size(); i = i + k)
		{
			if (i > 0)
			{
				if (inter_result[i + idx] == inter_result[i + idx - k])
				{
					unsigned common;
					for (unsigned j = begin_loc[idx]; j < end_loc[idx]; ++j)
					{
						common = csr.adjacency_list[label][j];
						//inter_common++;
						bool find = true;
						for (unsigned t = 0; t < same_label_idx.size(); ++t)
						{
							if (common == inter_result[i + same_label_idx[t]])
							{
								find = false;
								break;
							}
						}
						if (!find)
							continue;
						unsigned loc = (common - lower_boundary) / 8;
						char off = 1 << (7 - ((common - lower_boundary) % 8));
						if ((candidate_bitmap[loc] & off) != 0)
						{
							for (unsigned t = i; t < i + k; ++t)
								res_out.push_back(inter_result[t]);
							res_out.push_back(common);
						}
					}
					continue;
				}
			}
			unsigned long long ns = csr.neighborhood_structure[inter_result[i + idx]];
			unsigned start_address = ns >> 32;
			unsigned bit = label;
			while (bit >= 8)
			{
				start_address += eight_bit_sum[ns & 255];
				bit -= 8;
				ns = ns >> 8;
			}
			start_address += eight_bit_sum[ns & ((unsigned long long)pow(2, bit) - 1)];
			unsigned index = csr.index_list[start_address];
			begin_loc[idx] = csr.offset_list[label][index];
			end_loc[idx] = csr.offset_list[label][index + 1];
			unsigned common;
			for (unsigned j = begin_loc[idx]; j < end_loc[idx]; ++j)
			{
				common = csr.adjacency_list[label][j];
				//inter_common++;
				bool find = true;
				for (unsigned t = 0; t < same_label_idx.size(); ++t)
				{
					if (common == inter_result[i + same_label_idx[t]])
					{
						find = false;
						break;
					}
				}
				if (!find)
					continue;
				unsigned loc = (common - lower_boundary) / 8;
				char off = 1 << (7 - ((common - lower_boundary) % 8));
				if ((candidate_bitmap[loc] & off) != 0)
				{
					for (unsigned t = i; t < i + k; ++t)
						res_out.push_back(inter_result[t]);
					res_out.push_back(common);
				}
			}
		}
	}
	else
	{
		unsigned* begin_loc = new unsigned[batch_size * k];
		unsigned* end_loc = new unsigned[batch_size * k];
		unsigned* left = new unsigned[k];
		unsigned* right = new unsigned[k];
		lower_boundary = csr.label_index_list[label];
		vector<unsigned> idx;
		for (unsigned i = 3; i < len - 3; i = i + 3)
		{
			if (query_plan[i] == 1)
				idx.push_back(i / 3 - 1);
		}
		
		unsigned batch_num = inter_result.size() / (batch_size * k);
		unsigned rem = inter_result.size() % (batch_size * k);
		for (unsigned bt = 0; bt < batch_num; ++bt)
		{
			for (unsigned j = 0; j < idx.size(); ++j)
			{
				for (unsigned i = 0; i < batch_size; ++i)
				{
					if (bt + i > 0)
					{
						if (inter_result[(bt * batch_size + i) * k + idx[j]] == inter_result[(bt * batch_size + i) * k + idx[j] - k])
						{
							begin_loc[i * k + idx[j]] = begin_loc[((i - 1 + batch_size) % batch_size) * k + idx[j]];
							end_loc[i * k + idx[j]] = end_loc[((i - 1 + batch_size) % batch_size) * k + idx[j]];
							continue;
						}
					}
					unsigned long long ns = csr.neighborhood_structure[inter_result[(bt * batch_size + i) * k + idx[j]]];
					unsigned start_address = ns >> 32;
					unsigned bit = label;
					while (bit >= 8)
					{
						start_address += eight_bit_sum[ns & 255];
						bit -= 8;
						ns = ns >> 8;
					}
					start_address += eight_bit_sum[ns & ((unsigned long long)pow(2, bit) - 1)];
					unsigned index = csr.index_list[start_address];
					begin_loc[i * k + idx[j]] = csr.offset_list[label][index];
					end_loc[i * k + idx[j]] = csr.offset_list[label][index + 1];
				}
			}

			for (unsigned i = 0; i < batch_size; ++i)
			{
				unsigned start_set_idx1 = idx[0];
				unsigned start_set_idx2 = idx[1];
				if (k > 2)
				{
					unsigned min_size = 1000000000;
					unsigned second_min_size = 1000000000;
					for (unsigned j = 0; j < idx.size(); ++j)
					{
						if ((end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]]) < min_size)
						{
							min_size = end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]];
							start_set_idx1 = idx[j];
							second_min_size = min_size;
							start_set_idx2 = start_set_idx1;
						}
						else
						{
							if ((end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]]) < second_min_size)
							{
								second_min_size = end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]];
								start_set_idx2 = idx[j];
							}
						}
					}
				}

				for (unsigned j = 0; j < idx.size(); ++j)
				{
					left[idx[j]] = begin_loc[i * k + idx[j]];
					right[idx[j]] = end_loc[i * k + idx[j]];
				}

				unsigned left1 = left[start_set_idx1];
				unsigned left2 = left[start_set_idx2];
				unsigned right1 = right[start_set_idx1];
				unsigned right2 = right[start_set_idx2];
				while (true)
				{
					if (left1 == right1)
						break;
					if (left2 == right2)
						break;
					if (csr.adjacency_list[label][left1] > csr.adjacency_list[label][left2])
					{
						left2++;
						continue;
					}
					else if (csr.adjacency_list[label][left1] < csr.adjacency_list[label][left2])
					{
						left1++;
						continue;
					}
					else
					{
						unsigned common = csr.adjacency_list[label][left1];
						bool find = true;
						if (k > 2)
						{
							for (unsigned t = 0; t < idx.size(); ++t)
							{
								bool flag = false;
								if ((idx[t] == start_set_idx1) || (idx[t] == start_set_idx2))
									continue;
								unsigned low = left[idx[t]], high = right[idx[t]], middle = 0;
								while (low < high) {
									middle = (low + high) / 2;
									if (common == csr.adjacency_list[label][middle]) {
										flag = true;
										left[idx[t]] = middle + 1;
										break;
									}
									else if (common < csr.adjacency_list[label][middle]) {
										high = middle;
									}
									else if (common > csr.adjacency_list[label][middle]) {
										low = middle + 1;
									}
								}
								if (flag == false)
								{
									find = false;
									break;
								}
							}
						}
						for (unsigned t = 0; t < same_label_idx.size(); ++t)
						{
							if (common == inter_result[(bt * batch_size + i) * k + same_label_idx[t]])
							{
								find = false;
								break;
							}
						}
						if (!find)
						{
							left1++;
							left2++;
							continue;
						}
						//inter_common++;
						unsigned loc = (common - lower_boundary) / 8;
						char off = 1 << (7 - ((common - lower_boundary) % 8));
						if ((candidate_bitmap[loc] & off) != 0)
						{
							for (unsigned t = 0; t < k; ++t)
								res_out.push_back(inter_result[(bt * batch_size + i) * k + t]);
							res_out.push_back(common);
						}
						left1++;
						left2++;
						continue;
					}
				}
			}
		}
		for (unsigned j = 0; j < idx.size(); ++j)
		{
			for (unsigned i = 0; i < rem / k; ++i)
			{
				if (batch_num + i > 0)
				{
					if (inter_result[(batch_num * batch_size + i) * k + idx[j]] == inter_result[(batch_num * batch_size + i) * k + idx[j] - k])
					{
						begin_loc[i * k + idx[j]] = begin_loc[((i - 1 + batch_size) % batch_size) * k + idx[j]];
						end_loc[i * k + idx[j]] = end_loc[((i - 1 + batch_size) % batch_size) * k + idx[j]];
						continue;
					}
				}
				unsigned long long ns = csr.neighborhood_structure[inter_result[(batch_num * batch_size + i) * k + idx[j]]];
				unsigned start_address = ns >> 32;
				unsigned bit = label;
				while (bit >= 8)
				{
					start_address += eight_bit_sum[ns & 255];
					bit -= 8;
					ns = ns >> 8;
				}
				start_address += eight_bit_sum[ns & ((unsigned long long)pow(2, bit) - 1)];
				unsigned index = csr.index_list[start_address];
				begin_loc[i * k + idx[j]] = csr.offset_list[label][index];
				end_loc[i * k + idx[j]] = csr.offset_list[label][index + 1];
			}
		}

		for (unsigned i = 0; i < rem / k; ++i)
		{
			unsigned start_set_idx1 = idx[0];
			unsigned start_set_idx2 = idx[1];
			if (k > 2)
			{
				unsigned min_size = 1000000000;
				unsigned second_min_size = 1000000000;
				for (unsigned j = 0; j < idx.size(); ++j)
				{
					if ((end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]]) < min_size)
					{
						min_size = end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]];
						start_set_idx1 = idx[j];
						second_min_size = min_size;
						start_set_idx2 = start_set_idx1;
					}
					else
					{
						if ((end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]]) < second_min_size)
						{
							second_min_size = end_loc[i * k + idx[j]] - begin_loc[i * k + idx[j]];
							start_set_idx2 = idx[j];
						}
					}
				}
			}

			for (unsigned j = 0; j < idx.size(); ++j)
			{
				left[idx[j]] = begin_loc[i * k + idx[j]];
				right[idx[j]] = end_loc[i * k + idx[j]];
			}

			unsigned left1 = left[start_set_idx1];
			unsigned left2 = left[start_set_idx2];
			unsigned right1 = right[start_set_idx1];
			unsigned right2 = right[start_set_idx2];
			while (true)
			{
				if (left1 == right1)
					break;
				if (left2 == right2)
					break;
				if (csr.adjacency_list[label][left1] > csr.adjacency_list[label][left2])
				{
					left2++;
					continue;
				}
				else if (csr.adjacency_list[label][left1] < csr.adjacency_list[label][left2])
				{
					left1++;
					continue;
				}
				else
				{
					unsigned common = csr.adjacency_list[label][left1];
					bool find = true;
					if (k > 2)
					{
						for (unsigned t = 0; t < idx.size(); ++t)
						{
							bool flag = false;
							if ((idx[t] == start_set_idx1) || (idx[t] == start_set_idx2))
								continue;
							unsigned low = left[idx[t]], high = right[idx[t]], middle = 0;
							while (low < high) {
								middle = (low + high) / 2;
								if (common == csr.adjacency_list[label][middle]) {
									flag = true;
									left[idx[t]] = middle + 1;
									break;
								}
								else if (common < csr.adjacency_list[label][middle]) {
									high = middle;
								}
								else if (common > csr.adjacency_list[label][middle]) {
									low = middle + 1;
								}
							}
							if (flag == false)
							{
								find = false;
								break;
							}
						}
					}
					for (unsigned t = 0; t < same_label_idx.size(); ++t)
					{
						if (common == inter_result[(batch_num * batch_size + i) * k + same_label_idx[t]])
						{
							find = false;
							break;
						}
					}
					if (!find)
					{
						left1++;
						left2++;
						continue;
					}
					//inter_common++;
					unsigned loc = (common - lower_boundary) / 8;
					char off = 1 << (7 - ((common - lower_boundary) % 8));
					if ((candidate_bitmap[loc] & off) != 0)
					{
						for (unsigned t = 0; t < k; ++t)
							res_out.push_back(inter_result[(batch_num * batch_size + i) * k + t]);
						res_out.push_back(common);
					}
					left1++;
					left2++;
					continue;
				}
			}
		}


/*
		for (unsigned i = 0; i < inter_result.size(); i = i + k)
		{
			for (unsigned j = 0; j < idx.size(); ++j)
			{
				if (i > 0)
				{
					if (inter_result[i + idx[j]] == inter_result[i + idx[j] - k])
						continue;
				}
				unsigned long long ns = csr.neighborhood_structure[inter_result[i + idx[j]]];
				unsigned start_address = ns >> 32;
				unsigned bit = label;
				while (bit >= 8)
				{
					start_address += eight_bit_sum[ns & 255];
					bit -= 8;
					ns = ns >> 8;
				}
				start_address += eight_bit_sum[ns & ((unsigned long long)pow(2, bit) - 1)];
				unsigned index = csr.index_list[start_address];
				begin_loc[idx[j]] = csr.offset_list[label][index];
				end_loc[idx[j]] = csr.offset_list[label][index + 1];
			}
			unsigned min_size = 1000000000;
			unsigned start_set_idx;
			for (unsigned j = 0; j < idx.size(); ++j)
			{
				if ((end_loc[idx[j]] - begin_loc[idx[j]]) < min_size)
				{
					min_size = end_loc[idx[j]] - begin_loc[idx[j]];
					start_set_idx = idx[j];
				}
			}
			for (unsigned j = begin_loc[start_set_idx]; j < end_loc[start_set_idx]; ++j)
			{
				unsigned common = csr.adjacency_list[label][j];
				bool find = true;
				for (unsigned t = 0; t < idx.size(); ++t)
				{
					bool flag = false;
					if (idx[t] == start_set_idx)
						continue;
					unsigned low = begin_loc[idx[t]], high = end_loc[idx[t]], middle = 0;
					while (low < high) {
						middle = (low + high) / 2;
						if (common == csr.adjacency_list[label][middle]) {
							flag = true;
							break;
						}
						else if (common < csr.adjacency_list[label][middle]) {
							high = middle;
						}
						else if (common > csr.adjacency_list[label][middle]) {
							low = middle + 1;
						}
					}
					if (flag == false)
					{
						find = false;
						break;
					}
				}
				for (unsigned t = 0; t < same_label_idx.size(); ++t)
				{
					if (common == inter_result[i + same_label_idx[t]])
					{
						find = false;
						break;
					}
				}
				if (!find)
					continue;
				unsigned loc = (common - lower_boundary) / 8;
				char off = 1 << (7 - ((common - lower_boundary) % 8));
				if ((candidate_bitmap[loc] & off) != 0)
				{
					for (unsigned t = i; t < i + k; ++t)
						res_out.push_back(inter_result[t]);
					res_out.push_back(common);
				}
			}
		}*/
		//cout << "search time: " << tv1 << "ms" << endl;
		//cout << "locate time: " << tv2 << "ms" << endl;
	}
	return;
}

unsigned CPU_match(LPCSR& csr, Graph& query_graph, vector<unsigned>& join_order, unsigned** query_plan, vector<unsigned>* high_degree_candidate_list, vector<char>* candidate_bitmap, vector<unsigned>* inter_result, unsigned** isolated_node_begin_loc, unsigned** isolated_node_end_loc, bool** isolated_node_visited, vector<unsigned>& res)
{
	unsigned var_num = query_graph.vertex_num;
	/*unsigned start_node = query_plan[it + 1];
	unsigned start_node_label = query_plan[it + 2];
	vector<char>* temp_candidate_bitmap = new vector<char>[var_num];
	for (unsigned i = 0; i < var_num; ++i)
	{
		if (i == start_node)
			continue;
		temp_candidate_bitmap[i] = candidate_bitmap[i];
	}
	for (unsigned i = start_node + 1; i < var_num; ++i)
	{
		unsigned boundary = csr.label_index_list[query_graph.vertex_list[i]->label];
		for (unsigned j = 0; j < high_degree_candidate_list[i].size(); ++j)
		{
			unsigned idx = (high_degree_candidate_list[i][j] - boundary) / 8;
			unsigned off = 1 << (7 - ((high_degree_candidate_list[i][j] - boundary) % 8));
			temp_candidate_bitmap[i][idx] += off;
		}
	}*/

	struct timeval btv, etv;
	double CPU_time;
	for (unsigned i = 0; i < var_num; ++i)
	{
		gettimeofday(&btv, NULL);
		unsigned len = query_plan[i][0];
		unsigned next_node = query_plan[i][len - 3];
		if (i == 0)
		{
			for (unsigned j = 0; j < high_degree_candidate_list[next_node].size(); ++j)
				inter_result[i].push_back(high_degree_candidate_list[next_node][j]);	
			gettimeofday(&etv, NULL);
			CPU_time = (double)(etv.tv_sec - btv.tv_sec) * 1000.0 + (double)(etv.tv_usec - btv.tv_usec) / 1000.0;
			cout << "the " << i << "th iteration time: " << CPU_time << "ms" << endl;
			continue;
		}
		//unsigned inter_common = 0;
		multi_join(csr, query_plan[i], inter_result[i - 1], candidate_bitmap[next_node], inter_result[i], isolated_node_begin_loc, isolated_node_end_loc, isolated_node_visited);
		gettimeofday(&etv, NULL);
		CPU_time = (double)(etv.tv_sec - btv.tv_sec) * 1000.0 + (double)(etv.tv_usec - btv.tv_usec) / 1000.0;
		cout << "the " << i << "th iteration time: " << CPU_time << "ms" << endl;
		//cout << "the " << i << "th iteration intersection res num: " << inter_common << endl;
		cout << "the " << i << "th iteration output num: " << inter_result[i].size()/(i + 1) << endl;
	}
	unsigned k = 0;
	unsigned len = query_plan[var_num][0];
	for (unsigned i = 1; i < len; i = i + 3)
	{
		if (query_plan[var_num][i + 2] == 1000000000)
			k++;
		else
			break;
	}
	if (k == var_num)
	{
		res = inter_result[var_num - 1];
		return inter_result[var_num - 1].size() / var_num;
	}
	unsigned long long res_num = 0;
	for (unsigned long long i = 0; i < inter_result[var_num - 1].size(); i = i + k)
	{
		unsigned long long temp_num = 1;
		for (unsigned long long j = k * 3 + 1; j < len; j = j + 3)
		{
			unsigned isolated_node = query_plan[var_num][j];
			unsigned idx = query_plan[var_num][j + 2];
			unsigned lower_boundary = csr.label_index_list[query_plan[var_num][idx * 3 + 2]];
			temp_num *= isolated_node_end_loc[isolated_node][inter_result[var_num - 1][i + idx] - lower_boundary] - isolated_node_begin_loc[isolated_node][inter_result[var_num - 1][i + idx] - lower_boundary];
		}
		for (unsigned j = 0; j < k; ++j)
		{
			for (unsigned long long t = res_num; t < res_num + temp_num; ++t)
				res[t * var_num + j] = inter_result[var_num - 1][i + j];
		}
		unsigned cycle_num = 1;
		unsigned repeat_num = temp_num;
		for (unsigned j = k; j < var_num; ++j)
		{
			unsigned isolated_node = query_plan[var_num][j * 3 + 1];
			unsigned label = query_plan[var_num][j * 3 + 2];
			unsigned idx = query_plan[var_num][j * 3 + 3];
			unsigned lower_boundary = csr.label_index_list[query_plan[var_num][idx * 3 + 2]];
			unsigned neighbor_num = isolated_node_end_loc[isolated_node][inter_result[var_num - 1][i + idx] - lower_boundary] - isolated_node_begin_loc[isolated_node][inter_result[var_num - 1][i + idx] - lower_boundary];
			repeat_num /= neighbor_num;
			unsigned long long x = res_num;
			for (unsigned t = isolated_node_begin_loc[isolated_node][inter_result[var_num - 1][i + idx] - lower_boundary]; t < isolated_node_end_loc[isolated_node][inter_result[var_num - 1][i + idx] - lower_boundary]; ++t)
			{
				for (unsigned z = 0; z < cycle_num; ++z)
				{
					for (unsigned y = z * temp_num / cycle_num; y < z * temp_num / cycle_num + repeat_num; ++y)
					{
						res[(x + y) * var_num + j] = csr.adjacency_list[label][t];
					}
				}
				x += repeat_num;
			}
			cycle_num *= neighbor_num;
		}
		res_num += temp_num;
	}
	return res_num;
}

unsigned high_degree_subgraph_match(LPCSR& csr, vector<unsigned> query_plan, unsigned start_candidate, vector<char>* candidate_bitmap)
{
	unsigned it = 0;
	unsigned var_num = query_plan[it++];
	unsigned step_size = 2;
	vector<vector<unsigned>>* inter_result = new vector<vector<unsigned>>[2];
	for (unsigned i = 0; i < var_num; ++i)
	{
		if (i == 0)
		{
			vector<unsigned> tuple;
			tuple.push_back(start_candidate);
			inter_result[i % 2].push_back(tuple);
			it += step_size;
			step_size++;
		}
		else
		{
			inter_result[i % 2].clear();
			unsigned boundary = csr.label_index_list[query_plan[it + i + 1]];
			for (unsigned j = 0; j < inter_result[(i + 1) % 2].size(); ++j)
			{
				vector<unsigned>* adj = new vector<unsigned>[i];
				for (unsigned k = 0; k < i; ++k)
				{
					if (query_plan[it + k] == 0)
						continue;
					unsigned long long ns;
					if (k == 0)
						ns = csr.neighborhood_structure[start_candidate];
					else
						ns = csr.neighborhood_structure[inter_result[(i + 1) % 2][j][k - 1]];
					unsigned start_address = ns >> 32;
					unsigned index_off = 0;
					for (unsigned t = 0; t < query_plan[it + i + 1]; ++t)
					{
						if (ns & 1 == 1)
							index_off++;
						ns >> 1;
					}
					unsigned index = csr.index_list[start_address + index_off];
					unsigned begin_offset = csr.offset_list[query_plan[it + i + 1]][index];
					unsigned end_offset = csr.offset_list[query_plan[it + i + 1]][index + 1];
					for (unsigned t = begin_offset; t < end_offset; ++t)
						adj[k].push_back(csr.adjacency_list[query_plan[it + i + 1]][t]);
				}
				vector<unsigned> common_adj = sets_intersection(adj, i);
				for (unsigned k = 0; k < common_adj.size(); ++k)
				{
					unsigned common_id = common_adj[k];
					unsigned idx = (common_id - boundary) / 8;
					unsigned off = 1 << (7 - ((common_id - boundary) % 8));
					if ((candidate_bitmap[query_plan[it + i]][idx] & (char)off) == 0)
						continue;
					vector<unsigned> tuple;
					for (unsigned t = 0; t < i - 1; ++t)
						tuple.push_back(inter_result[(i + 1) % 2][j][t]);
					tuple.push_back(common_id);
					inter_result[i % 2].push_back(tuple);
				}
				delete[] adj;
			}
			it += step_size;
			step_size++;
		}
	}
	return inter_result[(var_num - 1) % 2].size();
}

unsigned CPU_handler(LPCSR& csr, Graph& query_graph, vector<unsigned> join_order, vector<unsigned>* high_degree_candidate_list, vector<char>* candidate_bitmap)
{/*
	vector<unsigned> query_plan = generate_query_plan(query_graph, join_order);
	unsigned it = 0;
	unsigned var_num = query_plan[it];
	unsigned start_node = query_plan[it + 1];
	unsigned start_node_label = query_plan[it + 2];
	vector<char>* temp_candidate_bitmap = new vector<char>[var_num];
	for (unsigned i = 0; i < var_num; ++i)
	{
		if (i == start_node)
			continue;
		temp_candidate_bitmap[i] = candidate_bitmap[i];
	}
	for (unsigned i = start_node + 1; i < var_num; ++i)
	{
		unsigned boundary = csr.label_index_list[query_graph.vertex_list[i]->label];
		for (unsigned j = 0; j < high_degree_candidate_list[i].size(); ++j)
		{
			unsigned idx = (high_degree_candidate_list[i][j] - boundary) / 8;
			unsigned off = 1 << (7 - ((high_degree_candidate_list[i][j] - boundary) % 8));
			temp_candidate_bitmap[i][idx] += off;
		}
	}
	unsigned result_num = 0;
	cout << "high_degree_vertex_num:" << high_degree_candidate_list[start_node].size() << endl;
	for (unsigned i = 0; i < high_degree_candidate_list[start_node].size(); ++i)
		result_num += high_degree_subgraph_match(csr, query_plan, high_degree_candidate_list[start_node][i], temp_candidate_bitmap);
	return result_num;*/
}

void print(Graph& G)
{
	cout << endl;
	cout << "Each vertex's neighbor_list:" << endl;
	for (unsigned i = 0; i < G.vertex_num; ++i)
	{
		cout << G.vertex_list[i]->vid << ": ";
		for (unsigned j = 0; j < G.vertex_list[i]->degree; ++j)
			cout << G.vertex_list[i]->neighbor_list[j]->vid << " ";
		cout << endl;
	}
}

int main(int argc, char* argv[])
{
	string op = argv[1];
	if (op == "-p")
	{
		//read data graph
		string filename = argv[2];
		Graph data_graph;
		begin_time = clock();
		data_graph.init(filename);
		end_time = clock();
		double read_time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
		cout << "data_graph = " << filename << endl;
		cout << "label_num = " << data_graph.label_num << endl;
		cout << "vertex_num = " << data_graph.vertex_num << endl;
		cout << "edge_num = " << data_graph.edge_num << endl;
		cout << "max_degree = " << data_graph.max_degree << endl;
		cout << "average_degree = " << data_graph.average_degree << endl;
		cout << "read_time = " << read_time * 1000 << "ms" << endl;
		//print(data_graph);

		//graph partition
		unsigned cut_degree = 10000;
		begin_time = clock();
		data_graph.partition();
		end_time = clock();
		double partition_time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
		cout << endl;
		cout << "partition_max_degree = " << data_graph.max_degree << endl;
		cout << "high_degree_vertex_num = " << data_graph.count(cut_degree) << endl;
		cout << "partition_time = " << partition_time * 1000 << "ms" << endl;

		//export the data graph without high-degree vertices
		//string export_filename = "graph_without_high_degree_vertices10000";
		//string export_filename = "graph_without_vertex296205";
		//string vertex_label_mapping = "vertex_label_mapping";
		//begin_time = clock();
		//data_graph.export_data_graph(cut_degree, export_filename);
		//data_graph.export_data_graph_without_v(296205, export_filename);
		//data_graph.export_vertex_label_mapping(vertex_label_mapping);
		//end_time = clock();
		//double export_time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
		//cout << endl;
		//cout << "export_time = " << export_time * 1000 << "ms" << endl;

		/*
		//graph reordering
		unsigned decay_factor = 0.8;
		begin_time = clock();
		GRO(data_graph, decay_factor, cut_degree);
		end_time = clock();
		double reorder_time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
		cout << endl;
		cout << "reorder_time = " << reorder_time * 1000 << "ms" << endl;*/

		//construct LPCSR
		unsigned port_num = 2;
		LPCSR csr_on_CPU, csr_on_FPGA;
		begin_time = clock();
		csr_on_CPU.init(data_graph);
		//csr_on_FPGA.init(data_graph, port_num, cut_degree);
		end_time = clock();
		double construct_time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
		cout << endl;
		cout << "construct_time = " << construct_time * 1000 << "ms" << endl;

		//export LPCSR
		string LPCSR_filename = filename + ".LPCSR";
		begin_time = clock();
		csr_on_CPU.export_LPCSR(LPCSR_filename);
		end_time = clock();
		double export_LPCSR_time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
		cout << endl;
		cout << "export_LPCSR_time = " << export_LPCSR_time * 1000 << "ms" << endl;
		return 0;
	}

	//import LPCSR
	string LPCSR_filename = argv[2];
	LPCSR csr_on_CPU, csr_on_FPGA;
	begin_time = clock();
	csr_on_CPU.import_LPCSR(LPCSR_filename);
	end_time = clock();
	double import_LPCSR_time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
	cout << endl;
	cout << "import_LPCSR_time = " << import_LPCSR_time * 1000 << "ms" << endl;

	//compute 8-bits sum
	for (unsigned i = 0; i < 256; ++i)
	{
		unsigned sum = 0;
		unsigned num = i;
		for (unsigned j = 0; j < 8; ++j)
		{
			if ((num & 1) == 1)
				sum++;
			num = num >> 1;
		}
		eight_bit_sum[i] = sum;
	}

	unsigned long long limit_res_num = atoi(argv[4]);
	limit_res_len = limit_res_num * (unsigned long long)limit_query_node_number;
	vector<unsigned> result(limit_res_len);
	vector<char>* candidate_bitmap = new vector<char>[limit_query_node_number];
	while (true)
	{
		if (bitmap_size < (csr_on_CPU.vertex_num * 2) / csr_on_CPU.label_num)
			bitmap_size *= 2;
		else
			break;
	}
	for (unsigned i = 0; i < limit_query_node_number; ++i)
		candidate_bitmap[i].resize(bitmap_size, 0);


	cout << "preprocessing finish!" << endl;
	

		//read query graph
		string query_filename = argv[3];
		Graph query_graph;
		query_graph.init(query_filename, csr_on_CPU.vertex_label_mapping);
		cout << endl;
		cout << "query_graph = " << query_filename << endl;
		cout << "query_node_num = " << query_graph.vertex_num << endl;
		cout << "edge_num = " << query_graph.edge_num << endl;
		cout << "max_degree = " << query_graph.max_degree << endl;

		//generate candidate lists
		gettimeofday(&begin_tv, NULL);
		unsigned var_num = query_graph.vertex_num;
		vector<unsigned>* candidate_list = new vector<unsigned>[var_num];
		unsigned* min_size = new unsigned[var_num]; //the minimum candidate size
		unsigned* min_label = new unsigned[var_num]; //neighbor label with the minimum candidate size
		unsigned long long* query_node_neighborhood_structure = new unsigned long long[var_num];
		for (unsigned i = 0; i < var_num; ++i)
		{
			vertex* query_node = query_graph.vertex_list[i];
			unsigned degree = query_node->degree;
			query_node_neighborhood_structure[i] = 0;
			min_size[i] = 1000000000;	//the minimum candidate size
			for (unsigned j = 0; j < degree; ++j)
			{
				unsigned nlabel = query_node->neighbor_list[j]->label;
				if (csr_on_CPU.candidate_filter_list[nlabel][query_node->label].size() < min_size[i])
				{
					min_size[i] = csr_on_CPU.candidate_filter_list[nlabel][query_node->label].size();
					min_label[i] = nlabel;
				}
				query_node_neighborhood_structure[i] |= (unsigned long long)(1 << nlabel);
			}
			if (min_size[i] == 0)
			{
				gettimeofday(&end_tv, NULL);
				double generate_candidate_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;
				cout << "embeddings = 0" << endl;
				cout << "CPU_time = " << generate_candidate_time << "ms" << endl;
				return 0;
			}
		}
		unsigned** isolated_node_begin_loc = new unsigned*[var_num];
		unsigned** isolated_node_end_loc = new unsigned*[var_num];
		bool** isolated_node_visited = new bool*[var_num];
		for (unsigned i = 0; i < var_num; ++i)
		{
			vertex* query_node = query_graph.vertex_list[i];
			unsigned degree = query_node->degree;
			if (degree == 1)
			{
				//unsigned len = csr_on_FPGA.label_index_list[min_label[i] + 1] - csr_on_FPGA.label_index_list[min_label[i]];
				unsigned len = csr_on_CPU.label_index_list[min_label[i] + 1] - csr_on_CPU.label_index_list[min_label[i]];
				isolated_node_begin_loc[i] = new unsigned[len];
				isolated_node_end_loc[i] = new unsigned[len];
				isolated_node_visited[i] = new bool[len];
				for (unsigned j = 0; j < len; ++j)
					isolated_node_visited[i][j] = false;
				continue;
			}
			for (unsigned j = 0; j < min_size[i]; ++j)
			{
				unsigned vid = csr_on_CPU.candidate_filter_list[min_label[i]][query_node->label][j];
				if ((csr_on_CPU.neighborhood_structure[vid] & query_node_neighborhood_structure[i]) == query_node_neighborhood_structure[i])
				{
					candidate_list[i].push_back(vid);
				}
			}
		}
		gettimeofday(&end_tv, NULL);
		double generate_candidate_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;
		cout << endl;
		for (unsigned i = 0; i < var_num; ++i)
			cout << "node: " << query_graph.vertex_list[i]->vid << ", degree: " << query_graph.vertex_list[i]->degree << ", label: " << query_graph.vertex_list[i]->label << ", candidates: " << candidate_list[i].size() << endl;
		cout << "generate_candidate_time = " << generate_candidate_time << "ms" << endl;

		//high degree candidates
		//vector<unsigned>* high_degree_candidate_list = new vector<unsigned>[var_num];

		/*	//generate candidate bitmaps without high degree candidates
			vector<char>* candidate_bitmap = new vector<char>[var_num];
			unsigned bitmap_size = 1024 * 1024 * 8;
			for (unsigned i = 0; i < var_num; ++i)
			candidate_bitmap[i].resize(bitmap_size);
			for (unsigned i = 0; i < var_num; ++i)
			{
			unsigned boundary = csr_on_FPGA.label_index_list[query_graph.vertex_list[i]->label];
			for (unsigned j = 0; j < candidate_list[i].size(); ++j)
			{
			if ((csr_on_FPGA.neighborhood_structure[candidate_list[i][j]] << 32) == 0)
			{
			high_degree_candidate_list[i].push_back(candidate_list[i][j]);
			continue;
			}
			unsigned idx = (candidate_list[i][j] - boundary) / 8;
			unsigned off = 1 << (7 - ((candidate_list[i][j] - boundary) % 8));
			candidate_bitmap[i][idx] += off;
			}
			}*/

		//generate candidate bitmaps with high degree candidates

		gettimeofday(&begin_tv, NULL);
		for (unsigned i = 0; i < var_num; ++i)
		{
			if (query_graph.vertex_list[i]->degree == 1)
				continue;
			//unsigned boundary = csr_on_FPGA.label_index_list[query_graph.vertex_list[i]->label];
			unsigned boundary = csr_on_CPU.label_index_list[query_graph.vertex_list[i]->label];
			for (unsigned j = 0; j < candidate_list[i].size(); ++j)
			{
				unsigned idx = (candidate_list[i][j] - boundary) / 8;
				unsigned off = 1 << (7 - ((candidate_list[i][j] - boundary) % 8));
				candidate_bitmap[i][idx] += off;
			}
		}
		gettimeofday(&end_tv, NULL);
		double generate_bitmap_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;
		cout << endl;
		cout << "generate_bitmap_time = " << generate_bitmap_time << "ms" << endl;

		//determine FPGA join order
		gettimeofday(&begin_tv, NULL);
		double* score = new double[var_num];
		for (unsigned i = 0; i < var_num; ++i)
		{
			if (query_graph.vertex_list[i]->degree == 1)
				score[i] = 1000000000.0;
			else
				score[i] = double(candidate_list[i].size()) / double(query_graph.vertex_list[i]->degree);
		}
		unsigned next_node = 0;
		double next_node_score = score[0];
		for (unsigned i = 1; i < var_num; ++i)
		{
			if (score[i] < next_node_score)
			{
				next_node = i;
				next_node_score = score[i];
			}
		}
		vector<unsigned> FPGA_join_order = determine_join_order(query_graph, next_node, candidate_list, score);
		gettimeofday(&end_tv, NULL);
		double determine_FPGA_join_order_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;

		/*cout << endl;
		cout << "computed_join_order: [";
		for (unsigned i = 0; i < var_num; ++i)
			cout << " " << FPGA_join_order[i];
		cout << " ]" << endl;

		for (unsigned i = 0; i < var_num; ++i)
			FPGA_join_order[i] = atoi(argv[i + 5]);*/

		cout << endl;
		cout << "join_order: [";
		for (unsigned i = 0; i < var_num; ++i)
			cout << " " << FPGA_join_order[i];
		cout << " ]" << endl;
		cout << "determine_FPGA_join_order_time = " << determine_FPGA_join_order_time << "ms" << endl;

		/*//generate FPGA input
		vector<unsigned> FPGA_query_plan = generate_query_plan(query_graph, FPGA_join_order);
		vector<unsigned> initial_node_candidates;
		vector<char>* other_node_candidates = new vector<char>[var_num - 1];
		for (unsigned i = 0; i < candidate_list[FPGA_join_order[0]].size(); ++i)
		{
		if ((csr_on_FPGA.neighborhood_structure[candidate_list[FPGA_join_order[0]][i]] << 32) == 0)
		continue;
		initial_node_candidates.push_back(candidate_list[FPGA_join_order[0]][i]);
		}
		for (unsigned i = 1; i < var_num; ++i)
		other_node_candidates[i - 1] = candidate_bitmap[FPGA_join_order[i]];*/

		//subgraph match on CPU
		gettimeofday(&begin_tv, NULL);
		vector<unsigned>* inter_result = new vector<unsigned>[var_num];
		unsigned** query_plan = generate_query_plan(query_graph, FPGA_join_order);
		double CPU_time = 0;
		unsigned result_num = CPU_match(csr_on_CPU, query_graph, FPGA_join_order, query_plan, candidate_list, candidate_bitmap, inter_result, isolated_node_begin_loc, isolated_node_end_loc, isolated_node_visited, result);
		/*	unsigned result_num = 0;
			for (unsigned i = 0; i < var_num; ++i)
			{
			cout << i << endl;
			begin_time = clock();
			if (high_degree_candidate_list[i].size() == 0)
			continue;
			vector<unsigned> join_order = determine_join_order(query_graph, i, candidate_list, score);
			cout << "join_order: [";
			for (unsigned j = 0; j < var_num; ++j)
			cout << " " << join_order[j];
			cout << " ]" << endl;
			//unsigned res = CPU_handler(csr_on_CPU, query_graph, join_order, high_degree_candidate_list, candidate_bitmap);
			unsigned res = CPU_match(csr_on_CPU, query_graph, high_degree_candidate_list, candidate_bitmap, join_order);
			cout << "i: " << res << endl;
			result_num += res;
			end_time = clock();
			double time = (double)(end_time - begin_time) / CLOCKS_PER_SEC;
			cout << "start_from_" << i << "_CPU_time = " << time * 1000 << "ms" << endl;
			CPU_time += time;
			}*/
		gettimeofday(&end_tv, NULL);
		CPU_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;
		double total_query_time = generate_candidate_time + generate_bitmap_time + determine_FPGA_join_order_time + CPU_time;
		cout << endl;
		cout << "embeddings = " << result_num << endl;
		cout << "CPU_time = " << CPU_time << "ms" << endl;
		cout << "total_query_time = " << total_query_time << "ms" << endl;

		delete[] eight_bit_sum;
		delete[] candidate_list;
		delete[] min_size;
		delete[] min_label;
		delete[] query_node_neighborhood_structure;
		for (unsigned i = 0; i < var_num; ++i)
		{
			vertex* query_node = query_graph.vertex_list[i];
			unsigned degree = query_node->degree;
			if (degree == 1)
			{
				delete[] isolated_node_begin_loc[i];
				delete[] isolated_node_end_loc[i];
				delete[] isolated_node_visited[i];
			}
		}
		delete[] isolated_node_begin_loc;
		delete[] isolated_node_end_loc;
		delete[] isolated_node_visited;
		delete[] candidate_bitmap;
		delete[] score;
		delete[] inter_result;
		for (unsigned i = 0; i <= var_num; ++i)
			delete[] query_plan[i];
		delete[] query_plan;
	return 0;
}