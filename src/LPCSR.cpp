#include "LPCSR.h"


LPCSR::LPCSR()
{
}


LPCSR::~LPCSR()
{
}

void LPCSR::init(Graph G, unsigned pnum, unsigned D)
{
	this->label_num = G.label_num;
	this->vertex_num = G.vertex_num;
	this->port_num = pnum;
	this->vertex_label_mapping = G.vertex_label_mapping;
	this->label_index_list = G.label_index_list;
	this->offset_list = new vector<unsigned>[this->label_num];
	this->adjacency_list = new vector<unsigned>[this->label_num];
	unsigned* index = new unsigned[this->label_num];
	unsigned* offset = new unsigned[this->label_num];
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		index[i] = 0;
		offset[i] = 0;
	}
	unsigned long long initial_address = 0;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < G.label_list[i].size(); ++j)
		{
			unsigned long long vertex_neighborhood_structure = initial_address << 32;
			if (G.label_list[i][j]->partition_max_degree >= D)
			{
				neighborhood_structure.push_back(vertex_neighborhood_structure);
				continue;
			}
			for (unsigned k = 0; k < this->label_num; ++k)
			{
				if (G.label_list[i][j]->partition_neighbor_list[k].size() > 0)
				{
					unsigned len = G.label_list[i][j]->partition_neighbor_list[k].size();
					for (unsigned t = 0; t < G.label_list[i][j]->partition_neighbor_list[k].size(); ++t)
					{
						vertex* neighbor = G.label_list[i][j]->partition_neighbor_list[k][t];
						if (neighbor->partition_max_degree >= D)
						{
							len--;
							continue;
						}
						this->adjacency_list[k].push_back(neighbor->vid);
					}
					if (len == 0)
						continue;
					this->offset_list[k].push_back(offset[k]);
					offset[k] += len;
					this->index_list.push_back(index[k]);
					index[k]++;
					vertex_neighborhood_structure |= (unsigned long long)(1 << k);
					initial_address++;
				}
			}
			neighborhood_structure.push_back(vertex_neighborhood_structure);
		}
	}
	for (unsigned i = 0; i < this->label_num; ++i)
		this->offset_list[i].push_back(offset[i]);
}

void LPCSR::init(Graph G)
{
	this->label_num = G.label_num;
	this->vertex_num = G.vertex_num;
	this->vertex_label_mapping = G.vertex_label_mapping;
	this->label_index_list = G.label_index_list;
	this->candidate_filter_list = new vector<unsigned>*[this->label_num];
	this->offset_list = new vector<unsigned>[this->label_num];
	this->adjacency_list = new vector<unsigned>[this->label_num];
	unsigned* index = new unsigned[this->label_num];
	unsigned* offset = new unsigned[this->label_num];
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		index[i] = 0;
		offset[i] = 0;
		this->candidate_filter_list[i] = new vector<unsigned>[this->label_num];
	}
	unsigned long long initial_address = 0;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < G.label_list[i].size(); ++j)
		{
			unsigned long long vertex_neighborhood_structure = initial_address << 32;
			for (unsigned k = 0; k < this->label_num; ++k)
			{
				if (G.label_list[i][j]->partition_neighbor_list[k].size() > 0)
				{
					unsigned len = G.label_list[i][j]->partition_neighbor_list[k].size();
					for (unsigned t = 0; t < G.label_list[i][j]->partition_neighbor_list[k].size(); ++t)
					{
						vertex* neighbor = G.label_list[i][j]->partition_neighbor_list[k][t];
						this->adjacency_list[k].push_back(neighbor->vid);
					}
					this->offset_list[k].push_back(offset[k]);
					this->candidate_filter_list[k][i].push_back(G.label_list[i][j]->vid);
					offset[k] += len;
					this->index_list.push_back(index[k]);
					index[k]++;
					vertex_neighborhood_structure |= (unsigned long long)(1 << k);
					initial_address++;
				}
			}
			neighborhood_structure.push_back(vertex_neighborhood_structure);
		}
	}
	for (unsigned i = 0; i < this->label_num; ++i)
		this->offset_list[i].push_back(offset[i]);
}

void LPCSR::export_LPCSR(const string& filename)
{
	ofstream fout(filename.c_str());
	fout << this->vertex_num << endl;

	for (unsigned i = 0; i < this->vertex_num; ++i)
	{
		if (i == 0)
			fout << this->neighborhood_structure[i];
		else
			fout << " " << this->neighborhood_structure[i];
	}
	fout << endl;

	fout << this->index_list.size() << endl;

	for (unsigned i = 0; i < this->index_list.size(); ++i)
	{
		if (i == 0)
			fout << this->index_list[i];
		else
			fout << " " << this->index_list[i];
	}
	fout << endl;

	fout << this->label_num << endl;

	for (unordered_map<unsigned, unsigned>::iterator i = this->vertex_label_mapping.begin(); i != this->vertex_label_mapping.end(); i++)
		fout << i->first << " " << i->second << " ";
	fout << endl;

	for (unsigned i = 0; i <= this->label_num; ++i)
	{
		if (i == 0)
			fout << this->label_index_list[i];
		else
			fout << " " << this->label_index_list[i];
	}
	fout << endl;

	for (unsigned i = 0; i < this->label_num; ++i)
	{
		for (unsigned j = 0; j < this->label_num; ++j)
		{
			fout << this->candidate_filter_list[i][j].size() << endl;
			for (unsigned k = 0; k < this->candidate_filter_list[i][j].size(); ++k)
			{
				if (k == 0)
					fout << this->candidate_filter_list[i][j][k];
				else
					fout << " " << this->candidate_filter_list[i][j][k];
			}
			fout << endl;
		}
	}

	for (unsigned i = 0; i < this->label_num; ++i)
	{
		fout << this->offset_list[i].size() << endl;
		for (unsigned j = 0; j < this->offset_list[i].size(); ++j)
		{
			if (j == 0)
				fout << this->offset_list[i][j];
			else
				fout << " " << this->offset_list[i][j];
		}
		fout << endl;
	}

	for (unsigned i = 0; i < this->label_num; ++i)
	{
		fout << this->adjacency_list[i].size() << endl;
		for (unsigned j = 0; j < this->adjacency_list[i].size(); ++j)
		{
			if (j == 0)
				fout << this->adjacency_list[i][j];
			else
				fout << " " << this->adjacency_list[i][j];
		}
		fout << endl;
	}
	fout.close();
}

void LPCSR::import_LPCSR(const string& filename)
{
	ifstream fin(filename.c_str());
	fin >> this->vertex_num;

	this->neighborhood_structure.resize(this->vertex_num);
	for (unsigned i = 0; i < this->vertex_num; ++i)
		fin >> this->neighborhood_structure[i];

	unsigned size;
	fin >> size;
	this->index_list.resize(size);
	for (unsigned i = 0; i < size; ++i)
		fin >> this->index_list[i];

	fin >> this->label_num;

	unsigned a, b;
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		fin >> a >> b;
		this->vertex_label_mapping.insert(pair<unsigned, unsigned>(a, b));
	}

	this->label_index_list.resize(label_num + 1);
	for (unsigned i = 0; i <= this->label_num; ++i)
		fin >> this->label_index_list[i];

	this->candidate_filter_list = new vector<unsigned>*[this->label_num];
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		this->candidate_filter_list[i] = new vector<unsigned>[this->label_num];
		for (unsigned j = 0; j < this->label_num; ++j)
		{
			fin >> size;
			this->candidate_filter_list[i][j].resize(size);
			for (unsigned k = 0; k < size; ++k)
				fin >> this->candidate_filter_list[i][j][k];
		}
	}

	this->offset_list = new vector<unsigned>[this->label_num];
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		fin >> size;
		this->offset_list[i].resize(size);
		for (unsigned j = 0; j < size; ++j)
			fin >> this->offset_list[i][j];
	}

	this->adjacency_list = new vector<unsigned>[this->label_num];
	for (unsigned i = 0; i < this->label_num; ++i)
	{
		fin >> size;
		this->adjacency_list[i].resize(size);
		for (unsigned j = 0; j < size; ++j)
			fin >> this->adjacency_list[i][j];
	}
	fin.close();
}