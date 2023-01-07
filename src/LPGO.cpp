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

struct timeval begin_tv, end_tv;

void InsertSort(vector<pair<int, double>>& vec, int left, int right)
{
	for (int i = left + 1; i <= right; ++i)
	{
		int j = i;
		while ((j - 1 >= left) && (vec[j].second > vec[j - 1].second))
		{
			swap(vec[j - 1], vec[j]);
			--j;
		}
	}
}

int FindMidMid(vector<pair<int, double>>& vec, int left, int right)
{
	if (right - left + 1 <= 5)
	{
		InsertSort(vec, left, right);
		return (left + right) >> 1;
	}
	int j = left - 1;
	for (int i = left; i <= right; i += 5)
	{
		if (right - i + 1 <= 5)
		{
			InsertSort(vec, i, right);
			swap(vec[++j], vec[(i + right) >> 1]);
			continue;
		}
		InsertSort(vec, i, i + 4);
		swap(vec[++j], vec[i + 2]);
	}
	return FindMidMid(vec, left, j);
}

int Partiton(vector<pair<int, double>>& vec, int left, int right, int pivot_index)
{
	swap(vec[pivot_index], vec[right]);
	int j = left;
	for (int i = left; i < right; ++i)
	{
		if (vec[i].second >= vec[right].second)
			swap(vec[j++], vec[i]);
	}
	swap(vec[j], vec[right]);
	return j;
}

void TopK(vector<pair<int, double>>& vec, int left, int right, int k)
{
	if (left == right)
		return;
	int pivot_index = FindMidMid(vec, left, right);
	int index = Partiton(vec, left, right, pivot_index);
	int num = index - left + 1;
	if (num == k)
		return;
	else if (num > k)
		TopK(vec, left, index - 1, k);
	else
		TopK(vec, index + 1, right, k - num);
}

bool compare(pair<unsigned, double> p1, pair<unsigned, double> p2)
{
	return (p1.second < p2.second); 
}

void LPGO(Graph& G, int w, string filename)
{
	cout << endl;
	vector<int> vertex_order(G.vertex_num);
	vector<int> inverted_vertex_order(G.vertex_num);
	vector<bool> vertex_visited(G.vertex_num);
	vector<double> vertex_score(G.vertex_num);
	vector<vector<vector<double>>> common_neighbors(w - 1);
	vector<vector<vector<int>>> timestamp(w - 1);
	for (int i = 0; i < w - 1; ++i)
	{
		common_neighbors[i].resize(G.vertex_num);
		timestamp[i].resize(G.vertex_num);
		for (int j = 0; j < G.vertex_num; ++j)
		{
			common_neighbors[i][j].resize(G.label_num);
			timestamp[i][j].resize(G.label_num);
		}
	}
	int iter = 0;
	int visited_iter = 0;
	for (int i = 0; i < G.label_num; ++i)
	{
		cout << "label: " << i << endl;
		vector<int> high_degree_vertex;
		int j;

		//find the initial node
		for (j = 0; j < G.label_list[i].size(); ++j)
		{
			if (G.label_list[i][j]->partition_max_degree > 10000)
			{
				high_degree_vertex.push_back(G.label_list[i][j]->vid);
				vertex_visited[G.label_list[i][j]->vid] = true;
				continue;
			}
			else
			{
				vertex_order[iter] = G.label_list[i][j]->vid;
				inverted_vertex_order[G.label_list[i][j]->vid] = iter;
				iter++;
				vertex_visited[G.label_list[i][j]->vid] = true;
				++j;
				break;
			}
		}

		//reorder
		vector<int> priority_queue(w - 1);
		vertex* node;
		vertex* neighbor;
		vertex* brother_node;
		for (; j < G.label_list[i].size(); ++j)
		{
			if (j == G.label_list[i].size() - 1)
			{
				while (true)
				{
					if (vertex_visited[visited_iter])
						visited_iter++;
					else
					{
						vertex_order[iter] = visited_iter;
						inverted_vertex_order[visited_iter] = iter;
						iter++;
						vertex_visited[visited_iter] = true;
						break;
					}
				}
				visited_iter = G.label_list[i].size();
			}
			if (visited_iter == G.label_list[i].size())
				break;
			node = G.vertex_list[vertex_order[iter - 1]];
			int vmax_id;
			set<int> update;

			for (int k = 0; k < node->degree; ++k)
			{
				neighbor = node->neighbor_list[k];
				if (neighbor->partition_max_degree > 10000)
					continue;
				for (int t = 0; t < neighbor->partition_neighbor_list[i].size(); ++t)
				{
					brother_node = neighbor->partition_neighbor_list[i][t];
					if (vertex_visited[brother_node->vid])
						continue;
					else
					{
						if (brother_node->partition_max_degree > 10000)
						{
							high_degree_vertex.push_back(brother_node->vid);
							vertex_visited[brother_node->vid] = true;
							continue;
						}
						auto succ = update.insert(brother_node->vid);
						if (succ.second)
							common_neighbors[(iter - G.label_index_list[i] - 1) % (w - 1)][brother_node->vid][neighbor->label] = 0;
						common_neighbors[(iter - G.label_index_list[i] - 1) % (w - 1)][brother_node->vid][neighbor->label] += 1;
						timestamp[(iter - G.label_index_list[i] - 1) % (w - 1)][brother_node->vid][neighbor->label] = iter - G.label_index_list[i];
					}
				}
			}

			vector<pair<int, double>> vec;
			for (set<int>::iterator k = update.begin(); k != update.end(); k++)
			{
				vertex_score[*k] = 0;
				vector<double> read_length(G.label_num);
				for (int l = 0; l < G.label_num; ++l)
					read_length[l] = G.vertex_list[*k]->partition_neighbor_list[l].size();
				for (int t = iter - G.label_index_list[i] - 1; t >= max(0, iter - (int)G.label_index_list[i] - w + 1); t--)
				{
					for (int l = 0; l < G.label_num; ++l)
					{
						if (timestamp[t % (w - 1)][*k][l] <= iter - G.label_index_list[i] - w + 1)
							continue;
						read_length[l] += G.vertex_list[vertex_order[G.label_index_list[i] + t]]->partition_neighbor_list[l].size();
						vertex_score[*k] += common_neighbors[t % (w - 1)][*k][l] / (read_length[l] * (double)G.label_num);
					}
				}
				vec.push_back(pair<int, double>(*k, vertex_score[*k]));
			}
			for (int k = 0; k < w - 1; ++k)
			{
				if (priority_queue[k] == 0)
					continue;
				if (update.find(priority_queue[k]) != update.end())
					continue;
				vertex_score[priority_queue[k]] = 0;
				vector<double> read_length(G.label_num);
				for (int l = 0; l < G.label_num; ++l)
					read_length[l] = G.vertex_list[priority_queue[k]]->partition_neighbor_list[l].size();
				for (int t = iter - G.label_index_list[i] - 1; t >= max(0, iter - (int)G.label_index_list[i] - w + 1); t--)
				{
					for (int l = 0; l < G.label_num; ++l)
					{
						if (timestamp[t % (w - 1)][priority_queue[k]][l] <= iter - G.label_index_list[i] - w + 1)
							continue;
						read_length[l] += G.vertex_list[vertex_order[G.label_index_list[i] + t]]->partition_neighbor_list[l].size();
						vertex_score[priority_queue[k]] += common_neighbors[t % (w - 1)][priority_queue[k]][l] / (read_length[l] * (double)G.label_num);
					}
				}
				vec.push_back(pair<int, double>(priority_queue[k], vertex_score[priority_queue[k]]));
			}
			if (vec.empty())
			{
				while (true)
				{
					if (vertex_visited[visited_iter])
						visited_iter++;
					else
					{
						vmax_id = visited_iter;
						break;
					}
				}
			}
			else
			{
				if (vec.size() < 20)
					InsertSort(vec, 0, vec.size() - 1);
				else
				{
					TopK(vec, 0, vec.size() - 1, w);
					InsertSort(vec, 0, w - 1);
				}
				vmax_id = vec[0].first;
			}

			if (G.vertex_list[vmax_id]->partition_max_degree > 10000)
			{
				high_degree_vertex.push_back(vmax_id);
				vertex_visited[vmax_id] = true;
			}
			else
			{
				vertex_order[iter] = vmax_id;
				inverted_vertex_order[vmax_id] = iter;
				iter++;
				vertex_visited[vmax_id] = true;
			}

			if (vec.empty())
			{
				for (int k = 0; k < w - 1; ++k)
					priority_queue[k] = 0;
			}
			else
			{
				if (vec.size() >= w)
				{
					for (int k = 0; k < w - 1; ++k)
						priority_queue[k] = vec[k + 1].first;
				}
				else
				{
					int k;
					for (k = 0; k < vec.size() - 1; ++k)
						priority_queue[k] = vec[k + 1].first;
					for (k = vec.size() - 1; k < w - 1; ++k)
						priority_queue[k] = 0;
				}
			}
		}
		for (j = 0; j < high_degree_vertex.size(); ++j)
		{
			vertex_order[iter] = high_degree_vertex[j];
			inverted_vertex_order[high_degree_vertex[j]] = iter;
			iter++;
		}
	}
	cout << vertex_order.size() << endl;
	vector<unsigned> vertex_label(G.label_num);
	for (unordered_map<unsigned, unsigned>::iterator it = G.vertex_label_mapping.begin(); it != G.vertex_label_mapping.end(); it++)
		vertex_label[it->second] = it->first;

	ofstream fout(filename.c_str());
	fout << "t 1 " << G.vertex_num << endl;
	for (unsigned i = 0; i < G.vertex_num; ++i)
		fout << "v " << i << " " << vertex_label[G.vertex_list[vertex_order[i]]->label] << endl;
	vertex* node;
	for (unsigned i = 0; i < G.vertex_num; ++i)
	{
		node = G.vertex_list[vertex_order[i]];
		vector<unsigned> neighbors;
		for (unsigned j = 0; j < node->degree; ++j)
		{
			unsigned new_vid = inverted_vertex_order[node->neighbor_list[j]->vid];
			if (i < new_vid)
				neighbors.push_back(new_vid);
		}
		sort(neighbors.begin(), neighbors.end());
		for (unsigned it = 0; it < neighbors.size(); it++)
			fout << "e " << i << " " << neighbors[it] << " 0" << endl;
	}
	fout.close();
}

/*void LPGO(Graph& G, int w, string filename)
{
	cout << endl;
	vector<unsigned> vertex_order(G.vertex_num);
	vector<unsigned> inverted_vertex_order(G.vertex_num);
	unsigned iter = 0;
	for (unsigned i = 0; i < G.label_num; ++i)
	{
		cout << "label: " << i << endl;

		vector<pair<unsigned, double>> vec;
		for (unsigned j = 0; j < G.label_list[i].size(); ++j)
			vec.push_back(pair<unsigned, double>(G.label_list[i][j]->vid, G.label_list[i][j]->degree));
		sort(vec.begin(), vec.end(), compare);

		for (unsigned j = 0; j < G.label_list[i].size(); ++j)
		{
			vertex_order[iter] = vec[j].first;
			inverted_vertex_order[vec[j].first] = iter;
			iter++;
		}
	}

	vector<unsigned> vertex_label(G.label_num);
	for (unordered_map<unsigned, unsigned>::iterator it = G.vertex_label_mapping.begin(); it != G.vertex_label_mapping.end(); it++)
		vertex_label[it->second] = it->first;

	ofstream fout(filename.c_str());
	fout << "t 1 " << G.vertex_num << endl;
	for (unsigned i = 0; i < G.vertex_num; ++i)
		fout << "v " << i << " " << vertex_label[G.vertex_list[vertex_order[i]]->label] << endl;
	cout<<"vertex output finish"<<endl;
	vertex* node;
	for (unsigned i = 0; i < G.vertex_num; ++i)
	{
		if(i%1000000==0)
			cout<<"vertex "<<i<<" finish"<<endl;
		node = G.vertex_list[vertex_order[i]];
		vector<unsigned> neighbors;
		for (unsigned j = 0; j < node->degree; ++j)
		{
			unsigned new_vid = inverted_vertex_order[node->neighbor_list[j]->vid];
			if (i < new_vid)
				neighbors.push_back(new_vid);
		}
		sort(neighbors.begin(),neighbors.end());
		for (unsigned it = 0; it < neighbors.size(); it++)
			fout << "e " << i << " " << neighbors[it] << " 0" << endl;
	}
	fout.close();
}*/

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
	//read data graph
	string filename = argv[1];
	Graph data_graph;
	gettimeofday(&begin_tv, NULL);
	data_graph.init(filename);
	gettimeofday(&end_tv, NULL);
	double read_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;
	cout << "data_graph = " << filename << endl;
	cout << "label_num = " << data_graph.label_num << endl;
	cout << "vertex_num = " << data_graph.vertex_num << endl;
	cout << "edge_num = " << data_graph.edge_num << endl;
	cout << "max_degree = " << data_graph.max_degree << endl;
	cout << "average_degree = " << data_graph.average_degree << endl;
	cout << "read_time = " << read_time << "ms" << endl;
	//print(data_graph);
	
	//graph partition
	unsigned cut_degree = 10000;
	gettimeofday(&begin_tv, NULL);
	data_graph.partition(true);
	gettimeofday(&end_tv, NULL);
	double partition_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;
	cout << endl;
	cout << "partition_max_degree = " << data_graph.max_degree << endl;
	cout << "high_degree_vertex_num = " << data_graph.count(cut_degree) << endl;
	cout << "partition_time = " << partition_time << "ms" << endl;
	
	//graph reordering
	string output_filename = argv[2];
	int w = 8;
	gettimeofday(&begin_tv, NULL);
	LPGO(data_graph, w, output_filename);
	gettimeofday(&end_tv, NULL);
	double reorder_time = (double)(end_tv.tv_sec - begin_tv.tv_sec) * 1000.0 + (double)(end_tv.tv_usec - begin_tv.tv_usec) / 1000.0;
	cout << endl;
	cout << "reorder_time = " << reorder_time << "ms" << endl;

	return 0;
}
