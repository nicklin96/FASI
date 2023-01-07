#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

int main(int argc, char* argv[])
{
	int threshold = 60;
	int file_num = atoi(argv[1]);
	int node_num = atoi(argv[2]);
	for (int i = 1; i <= file_num; ++i)
	{
		int edge_num = 0;
		vector<int>* neighbors = new vector<int>[node_num];
		for (int j = 1; j < node_num; ++j)
		{
			if (j == 1)
			{
				neighbors[0].push_back(1);
				neighbors[1].push_back(0);
				edge_num++;
				continue;
			}
			bool flag = true;
			while (flag)
			{
				for (int k = 0; k < j; ++k)
				{
					int prob = rand() % 100;
					if (prob < threshold)
					{
						neighbors[k].push_back(j);
						neighbors[j].push_back(k);
						edge_num++;
						flag = false;
					}
				}
			}
		}
		string filename1 = "../query/query/query_v" + to_string(node_num) + "_" + to_string(i) + ".graph";
		string filename2 = "../query/DAF_query/query_v" + to_string(node_num) + "_" + to_string(i) + ".graph";
		string filename3 = "../query/CECI_query/query_v" + to_string(node_num) + "_" + to_string(i) + ".graph";
		ofstream fout1(filename1);
		ofstream fout2(filename2);
		ofstream fout3(filename3);
		fout1 << "t " << node_num << " " << edge_num << endl;
		for (int j = 0; j < node_num; ++j)
			fout1 << "v " << j << " 0" << endl;
		for (int j = 0; j < node_num - 1; ++j)
		{
			for (int k = 0; k < neighbors[j].size(); ++k)
			{
				if (neighbors[j][k] > j)
					fout1 << "e " << j << " " << neighbors[j][k] << " 0" << endl;
			}
		}
		fout1.close();
		fout2 << "t 1 " << node_num << " " << edge_num*2 << endl;
		for (int j = 0; j < node_num; ++j)
		{
			fout2 << j << " 0 " << neighbors[j].size();
			for (int k = 0; k < neighbors[j].size(); ++k)
			{
				fout2 << " " << neighbors[j][k];
			}
			fout2 << endl;
		}
		fout2.close();
		fout3 << "t " << node_num << " " << edge_num << endl;
		for (int j = 0; j < node_num; ++j)
			fout3 << "v " << j << " 0 " << neighbors[j].size() << endl;
		for (int j = 0; j < node_num - 1; ++j)
		{
			for (int k = 0; k < neighbors[j].size(); ++k)
			{
				if (neighbors[j][k] > j)
					fout3 << "e " << j << " " << neighbors[j][k] << endl;
			}
		}
		fout3 << "t # -1" << endl;
		fout3.close();		
	}
	return 0;
}