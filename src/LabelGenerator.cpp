#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <dirent.h>
#include <stdio.h>
using namespace std;

int main(int argc, char* argv[])
{
	string graph_name = argv[1];
	int label_num = atoi(argv[2]);
	DIR *dir = opendir("../query/query/");
	if (opendir(("../query/" + graph_name + "_query/").c_str()) == NULL)
	{
		string cmd = "mkdir ../query/" + graph_name + "_query/";
		system(cmd.c_str());
	}
	if (opendir(("../query/" + graph_name + "_DAF_query/").c_str()) == NULL)
	{
		string cmd = "mkdir ../query/" + graph_name + "_DAF_query/";
		system(cmd.c_str());
	}
	if (opendir(("../query/" + graph_name + "_CECI_query/").c_str()) == NULL)
	{
		string cmd = "mkdir ../query/" + graph_name + "_CECI_query/";
		system(cmd.c_str());
	}	
	struct dirent *ptr;
	while ((ptr = readdir(dir)) != NULL)
	{
		string filename = ptr->d_name;
		if ((filename == ".") || (filename == ".."))
			continue;

		ifstream fin1("../query/query/" + filename);
		ifstream fin2("../query/DAF_query/" + filename);
		ifstream fin3("../query/CECI_query/" + filename);
		ofstream fout1("../query/" + graph_name + "_query/" + filename);
		ofstream fout2("../query/" + graph_name + "_DAF_query/" + filename);
		ofstream fout3("../query/" + graph_name + "_CECI_query/" + filename);
		set<int> labels;
		char type;
		while (fin1 >> type)
		{
			if (type == 't')
			{
				fin2 >> type;
				fin3 >> type;
				int node_num, edge_num;
				fin1 >> node_num >> edge_num;
				int temp;
				fin2 >> temp;
				fin2 >> temp;
				fin2 >> temp;
				fin3 >> temp;
				fin3 >> temp;
				fout1 << "t " << node_num << " " << edge_num << endl;
				fout2 << "t 1 " << node_num << " " << edge_num * 2 << endl;
				fout3 << "t " << node_num << " " << edge_num << endl;
			}
			else if (type == 'v')
			{
				fin3 >> type;
				int vid, label;
				fin1 >> vid >> label;
				fin2 >> vid >> label;
				fin3 >> vid >> label;
				int neighbor_num;
				fin2 >> neighbor_num;
				fin3 >> neighbor_num;
				while(true)
				{
					label = rand() % label_num;
					auto succ = labels.insert(label);
					if(succ.second)
						break;
				}
				fout1 << "v " << vid << " " << label << endl;
				fout2 << vid << " " << label << " " << neighbor_num;
				fout3 << "v " << vid << " " << label << " " << neighbor_num;
				int temp;
				for (int i = 0; i < neighbor_num; ++i)
				{
					fin2 >> temp;
					fout2 << " " << temp;
				}
				fout2 << endl;
			}
			else if (type == 'e')
			{
				fin3 >> type;
				int src_vid, dst_vid, label;
				fin1 >> src_vid >> dst_vid >> label;
				fin3 >> src_vid >> dst_vid;
				fout1 << "e " << src_vid << " " << dst_vid << " " << label << endl;
				fout3 << "e " << src_vid << " " << dst_vid << endl;
			}
		}
		fout3 << "t # -1" << endl;
		fin1.close();
		fin2.close();
		fin3.close();
		fout1.close();
		fout2.close();
		fout3.close();
	}
	closedir(dir);
	return 0;
}