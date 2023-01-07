#include "Graph.h"

Graph::Graph()
{

}

Graph::~Graph()
{

}

void Graph::init(const string& filename) // init the data graph
{
	char type;
	unsigned vid;
	unsigned src_vid;
	unsigned dst_vid;
	unsigned label;
	unsigned max_label = 0;
	ifstream fin(filename.c_str());
	while (fin>>type)
	{
		if (type == 't')
		{
			char temp;
			fin >> temp;
			fin >> temp;
		}
		else if (type == 'v')
		{
			fin >> vid >> label;
			vertex* node = new vertex();
			node->vid = vid;
			unordered_map<unsigned, unsigned>::iterator it = this->vertex_label_mapping.find(label);
			if (it == this->vertex_label_mapping.end())
			{
				this->vertex_label_mapping.insert(pair<unsigned, unsigned>(label, max_label));
				node->label = max_label;
				max_label++;
			}
			else
				node->label = it->second;
			this->vertex_list.push_back(node);
		}
		else if (type == 'e')
		{
			fin >> src_vid >> dst_vid >> label;
			this->vertex_list[src_vid]->neighbor_list.push_back(this->vertex_list[dst_vid]); 
			this->vertex_list[src_vid]->degree++;
			this->vertex_list[dst_vid]->neighbor_list.push_back(this->vertex_list[src_vid]);
			this->vertex_list[dst_vid]->degree++;
			this->edge_num++;
		}
	}
	this->label_num = max_label;
	this->vertex_num = this->vertex_list.size();
	fin.close();
	for (unsigned i = 0; i < this->vertex_num; ++i)
	{
		unsigned degree = this->vertex_list[i]->degree;
		this->average_degree += degree;
		if (degree > this->max_degree)
			this->max_degree = degree;
	}
	this->average_degree = this->average_degree / this->vertex_num;
}

void Graph::init(const string& filename, unordered_map<unsigned, unsigned>& mapping) //init the query graph
{
	char type;
	unsigned vid;
	unsigned src_vid;
	unsigned dst_vid;
	unsigned label;
	ifstream fin(filename.c_str());
	while (fin >> type)
	{
		if (type == 't')
		{
			char temp;
			fin >> temp;
			fin >> temp;
		}
		else if (type == 'v')
		{
			fin >> vid >> label;
			vertex* node = new vertex();
			node->vid = vid;
			node->label = mapping.find(label)->second;
			this->vertex_list.push_back(node);
		}
		else if (type == 'e')
		{
			fin >> src_vid >> dst_vid >> label;
			this->vertex_list[src_vid]->neighbor_list.push_back(this->vertex_list[dst_vid]);
			this->vertex_list[src_vid]->degree++;
			this->vertex_list[dst_vid]->neighbor_list.push_back(this->vertex_list[src_vid]);
			this->vertex_list[dst_vid]->degree++;
			this->edge_num++;
		}
	}
	this->vertex_num = this->vertex_list.size();
	fin.close();
	for (unsigned i = 0; i < this->vertex_num; ++i)
	{
		unsigned degree = this->vertex_list[i]->degree;
		if (degree > this->max_degree)
			this->max_degree = degree;
	}
}

//partition G based on vertex label
void Graph::partition(bool flag)
{
	this->label_list = new vector<vertex*>[this->label_num];
	for (unsigned i = 0; i < this->vertex_num; ++i)
		this->label_list[this->vertex_list[i]->label].push_back(this->vertex_list[i]);
	cout << endl;
	unsigned label_index = 0;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		this->label_index_list.push_back(label_index);
		cout << "label_" << i << "_num = " << this->label_list[i].size() << endl;
		label_index += this->label_list[i].size();
	}
	this->label_index_list.push_back(label_index);
	unsigned new_vid = 0;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			this->label_list[i][j]->vid = new_vid;
			if (flag)
				this->vertex_list[new_vid] = this->label_list[i][j];
			new_vid++;
		}
	}
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			this->label_list[i][j]->partition_neighbor_list = new vector<vertex*>[this->label_num];
			for (unsigned k = 0; k < this->label_list[i][j]->degree; ++k)
			{
				vertex* neighbor = this->label_list[i][j]->neighbor_list[k];
				this->label_list[i][j]->partition_neighbor_list[neighbor->label].push_back(neighbor);
				if (this->label_list[i][j]->partition_neighbor_list[neighbor->label].size() > this->label_list[i][j]->partition_max_degree)
					this->label_list[i][j]->partition_max_degree = this->label_list[i][j]->partition_neighbor_list[neighbor->label].size();
			}
		}
	}
	this->max_degree = 0;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			if (this->label_list[i][j]->partition_max_degree > this->max_degree)
				this->max_degree = this->label_list[i][j]->partition_max_degree;
			this->label_list[i][j]->partition_neighbor_degree_sum = new unsigned[this->label_num];
			for (unsigned k = 0; k < this->label_num; ++k)
				this->label_list[i][j]->partition_neighbor_degree_sum[k] = 0;
			for (unsigned k = 0; k < this->label_list[i][j]->degree; ++k)
			{
				vertex* neighbor = this->label_list[i][j]->neighbor_list[k];
				for (unsigned t = 0; t < this->label_num; ++t)
					this->label_list[i][j]->partition_neighbor_degree_sum[t] += neighbor->partition_neighbor_list[t].size();
			}
		}
	}
	ofstream f1("long.txt");
	for(unsigned i=0;i<this->vertex_num;++i)
		f1<<this->vertex_list[i]->vid<<endl;
	ofstream f2("short.txt");
	for(unsigned i=0;i<this->label_num;++i)
		for(unsigned j=0;j<this->label_list[i].size();++j)
		f2<<this->label_list[i][j]->vid<<endl;
	f1.close();
	f2.close();
}

//count the number of vertices with partition_max_degree greater than D
unsigned Graph::count(unsigned D)
{
	unsigned num = 0;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			if (this->label_list[i][j]->partition_max_degree >= D)
				num++;
		}
	}
	return num;
}

void Graph::export_data_graph(unsigned D, const string& filename) //export the data graph without high-degree vertices
{
	vector<unsigned> vertex_label(this->label_num);
	for (unordered_map<unsigned, unsigned>::iterator it = this->vertex_label_mapping.begin(); it != this->vertex_label_mapping.end(); it++)
		vertex_label[it->second] = it->first;

	ofstream fout(filename.c_str());
	ofstream f("high_degree_vertices");
	fout << "t 1 " << this->vertex_num  << endl;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			fout << "v " << this->label_list[i][j]->vid << " " << vertex_label[i] << endl;
		}
	}
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			if (this->label_list[i][j]->partition_max_degree < D)
			{
				for (unsigned k = 0; k < this->label_num; ++k)
				{
					for (unsigned t = 0; t < this->label_list[i][j]->partition_neighbor_list[k].size(); ++t)
					{
						vertex* neighbor = this->label_list[i][j]->partition_neighbor_list[k][t];
						if (neighbor->partition_max_degree >= D)
							continue;
						if (neighbor->vid < this->label_list[i][j]->vid)
							continue;
						fout << "e " << this->label_list[i][j]->vid << " " << neighbor->vid << " 0" << endl;
					}
				}
			}
			else
				f << "v " << this->label_list[i][j]->vid << " " << vertex_label[i] << endl;
		}
	}
	fout.close();
	f.close();
}

void Graph::export_data_graph_without_v(unsigned v, const string& filename) //export the data graph without vertex v
{
	vector<unsigned> vertex_label(this->label_num);
	for (unordered_map<unsigned, unsigned>::iterator it = this->vertex_label_mapping.begin(); it != this->vertex_label_mapping.end(); it++)
		vertex_label[it->second] = it->first;

	ofstream fout(filename.c_str());
	fout << "t 1 " << this->vertex_num << endl;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			fout << "v " << this->label_list[i][j]->vid << " " << vertex_label[i] << endl;
		}
	}
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_list[i].size(); ++j)
		{
			if (this->label_list[i][j]->vid != v)
			{
				for (unsigned k = 0; k < this->label_num; ++k)
				{
					for (unsigned t = 0; t < this->label_list[i][j]->partition_neighbor_list[k].size(); ++t)
					{
						vertex* neighbor = this->label_list[i][j]->partition_neighbor_list[k][t];
						if (neighbor->vid == v)
							continue;
						if (neighbor->vid < this->label_list[i][j]->vid)
							continue;
						fout << "e " << this->label_list[i][j]->vid << " " << neighbor->vid << " 0" << endl;
					}
				}
			}
		}
	}
	fout.close();
}

void Graph::export_vertex_label_mapping(const string& filename)
{
	vector<unsigned> vertex_label(this->label_num);
	for (unordered_map<unsigned, unsigned>::iterator it = this->vertex_label_mapping.begin(); it != this->vertex_label_mapping.end(); it++)
		vertex_label[it->second] = it->first;
	ofstream fout(filename.c_str());
	for (unsigned i = 0; i < this->label_num; ++i)
		fout << i << " " << vertex_label[i] << endl;
	fout.close();
}

bool Graph::has_edge(unsigned id1, unsigned id2) //check if there exists an edge between two vertices
{
	int len1 = this->vertex_list[id1]->degree;
	int len2 = this->vertex_list[id2]->degree;
	unsigned vid1 = this->vertex_list[id1]->vid;
	unsigned vid2 = this->vertex_list[id2]->vid;
	if (len1 < len2)
	{
		int low = 0;
		int high = len1 - 1;
		while (low <= high)
		{
			int mid = (low + high) / 2;
			if (this->vertex_list[id1]->neighbor_list[mid]->vid < vid2)
				low = mid + 1;
			else if (this->vertex_list[id1]->neighbor_list[mid]->vid > vid2)
				high = mid - 1;
			else
				return true;
		}
	}
	else
	{
		int low = 0;
		int high = len2 - 1;
		while (low <= high)
		{
			int mid = (low + high) / 2;
			if (this->vertex_list[id2]->neighbor_list[mid]->vid < vid1)
				low = mid + 1;
			else if (this->vertex_list[id2]->neighbor_list[mid]->vid > vid1)
				high = mid - 1;
			else
				return true;
		}
	}
	return false;
}
