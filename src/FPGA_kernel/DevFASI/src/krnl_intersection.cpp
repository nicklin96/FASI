#include "common.h"
#include <ap_int.h>
#include <hls_stream.h>
//#include <iostream>

static unsigned int l1[SLOT_NUM][DATA_SIZE];
static unsigned int l2[SLOT_NUM][DATA_SIZE];
static unsigned int l3[SLOT_NUM][DATA_SIZE];
static unsigned int l4[SLOT_NUM][DATA_SIZE];

static unsigned int offset_buffer1[OFFSET_BUFFER_SIZE];
static unsigned int adj_buffer1[ADJ_BUFFER_SIZE];

static unsigned int offset_buffer2[OFFSET_BUFFER_SIZE];
static unsigned int adj_buffer2[ADJ_BUFFER_SIZE];

static unsigned int offset_buffer3[OFFSET_BUFFER_SIZE];
static unsigned int adj_buffer3[ADJ_BUFFER_SIZE];

static unsigned int offset_buffer4[OFFSET_BUFFER_SIZE];
static unsigned int adj_buffer4[ADJ_BUFFER_SIZE];

static unsigned int off_buf_begin[4], off_buf_end[4];

static unsigned short slot_lock[SLOT_NUM];

static unsigned int chunk[8][CHUNK_SIZE];

static char CBIndex[CBI_SIZE];
static unsigned int selected[V_LIMIT];
static hls::stream<unsigned int> readOut("output_read");
static hls::stream<unsigned int> mergeOut("output_merge");
static hls::stream<unsigned int> bsearchOut_1("output_bsearch_1");
static hls::stream<unsigned int> bsearchOut_2("output_bsearch_2");
static int s1;
static int s2;
static int s3;
static int s4;

static void mergeTwo(){
	int size_left,size_right,chunk_size;
	unsigned int source;
	source = readOut.read();
	while(source!=END_FLAG){
		size_left = l1[source][0] + 1;
		size_right = l2[source][0] + 1;
		merge:for(int l = 1, r = 1; l < size_left && r < size_right ; ){
			//std::cout << "merging " << l << " " << r << std::endl;
			if(l1[source][l] == l2[source][r]){
				mergeOut << l1[source][l];
				mergeOut << source;
				l++;
				r++;
			}else if(l1[source][l] > l2[source][r]){
				r++;
			}else{
				l++;
			}
		}
		source = readOut.read();
	}
	mergeOut << END_FLAG;
	mergeOut << END_FLAG;

}

static void binarySearch(unsigned int list[4][DATA_SIZE], hls::stream<unsigned int>& inStream, hls::stream<unsigned int>& outStream){
	int mid;
	unsigned int target,source;
	target = inStream.read();
	source = inStream.read();
	while(target!=END_FLAG){
		bsearch:for(int l=1,r=list[source][0]; l <= r; ){
			mid = l + ((r - l) >> 1);
			//std::cout << "searching for " << target << std::endl;
			if(list[source][mid] == target){
				l = r+1;
				outStream << target;
				outStream << source;
			}else if(list[source][mid] > target){
				r = mid - 1;
			}else{
				l = mid + 1;
			}
		}
		target = inStream.read();
		source = inStream.read();
	}
	outStream << END_FLAG;
	outStream << END_FLAG;

}

static void write_result(unsigned int* res_list, hls::stream<unsigned int>& outStream, int size){
	unsigned int res,source;
	int i;
	res = outStream.read();
	source = outStream.read();
	//std::cout << res << std::endl;
	for(i = 1; i <= size && res!=END_FLAG; i++){
		res_list[i] = res;
		res = outStream.read();
		source = outStream.read();
	}
	res_list[0] = i;
	slot_lock[source] = 0;
}

static void init_buffer(char* cbi, unsigned int* join_list,
						unsigned int cbi_size){
	unsigned int local_cbi_size = cbi_size;
	for(unsigned int i = 0; i < local_cbi_size; i++){
		CBIndex[i] = cbi[i];
	}
	for(unsigned int j = 0; j < 4; j++){
		selected[j] = join_list[j];
	}
}

static void read_chunk(unsigned int* cand1, unsigned int k, unsigned int begin){
	unsigned int limit = CHUNK_SIZE * k;
	for(int i = 0; i < limit; i++){
		chunk[i] = cand1[i + begin];
	}
}

static void read_fpcsr(unsigned int* offset,
					   unsigned int* adj,
					   unsigned int k){
	unsigned int begin,end,source;
	unsigned int target[4];
	source = 0;

	off_buf_begin[0] = END_FLAG;off_buf_end[0] = END_FLAG;
	off_buf_begin[1] = END_FLAG;off_buf_end[1] = END_FLAG;
	off_buf_begin[2] = END_FLAG;off_buf_end[2] = END_FLAG;
	off_buf_begin[3] = END_FLAG;off_buf_end[3] = END_FLAG;

	for(unsigned int t = 0; t < CHUNK_SIZE; t++){
		//find a slot
		for(;slot_lock[source] == 1;source++){
			source = source % SLOT_NUM;
		}
		slot_lock[source] = 1;
		target[0] = chunk[t * k + selected[0]];
		target[1] = chunk[t * k + selected[1]];
		target[2] = chunk[t * k + selected[2]];
		target[3] = chunk[t * k + selected[3]];
		//update buffer


	}


}

static void read_csr(unsigned int* list1, unsigned int* list2,
        			 unsigned int* list3, unsigned int* list4,
					 unsigned int size1, unsigned int size2,
					 unsigned int size3, unsigned int size4){

	unsigned int source;
	source = 0;
	s1 = l1[source][0] = size1;
	s2 = l2[source][0] = size2;
	s3 = l3[source][0] = size3;
	s4 = l4[source][0] = size4;
	for(int i1 = 1; i1 <= s1; i1++){
		l1[source][i1] = list1[i1 - 1];
	}
	for(int i2 = 1; i2 <= s2; i2++){
		l2[source][i2] = list2[i2 - 1];
	}
	for(int i3 = 1; i3 <= s3; i3++){
		l3[source][i3] = list3[i3 - 1];
	}
	for(int i4 = 1; i4 <= s4; i4++){
		l4[source][i4] = list4[i4 - 1];
	}
	readOut << source;
	readOut << END_FLAG;
}

extern "C"{
void kerneltwo(
		IN_TYPE offset1[OFFSET_LEN],IN_TYPE adj1[ADJ_LEN],
		IN_TYPE offset2[OFFSET_LEN],IN_TYPE adj2[ADJ_LEN],
		IN_TYPE offset3[OFFSET_LEN],IN_TYPE adj3[ADJ_LEN],
		IN_TYPE offset4[OFFSET_LEN],IN_TYPE adj4[ADJ_LEN],
		IN_TYPE cand1[CAN_LEN],CAN_TYPE cand2[CAN2_LEN],CAN_TYPE findex[INDEX_LEN],
		IN_TYPE size1,IN_TYPE size2,IN_TYPE index_size,
		IN_TYPE pid1,IN_TYPE pid2,
		IN_TYPE out_buf[OUT_LEN]
		)
{
#pragma HLS INTERFACE m_axi port=cand1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=cand2 offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=out_buf offset=slave bundle=gmem2

#pragma HLS INTERFACE m_axi port=offset1 offset=slave bundle=gmem3
#pragma HLS INTERFACE m_axi port=adj1 offset=slave bundle=gmem4
#pragma HLS INTERFACE m_axi port=offset2 offset=slave bundle=gmem5
#pragma HLS INTERFACE m_axi port=adj2 offset=slave bundle=gmem6

#pragma HLS INTERFACE m_axi port=offset3 offset=slave bundle=gmem7
#pragma HLS INTERFACE m_axi port=adj3 offset=slave bundle=gmem8
#pragma HLS INTERFACE m_axi port=offset4 offset=slave bundle=gmem9
#pragma HLS INTERFACE m_axi port=adj4 offset=slave bundle=gmem10
#pragma HLS INTERFACE m_axi port=findex offset=slave bundle=gmem11

#pragma HLS INTERFACE s_axilite port=offset1 bundle=control
#pragma HLS INTERFACE s_axilite port=adj1 bundle=control
#pragma HLS INTERFACE s_axilite port=offset2 bundle=control
#pragma HLS INTERFACE s_axilite port=adj2 bundle=control
#pragma HLS INTERFACE s_axilite port=offset3 bundle=control
#pragma HLS INTERFACE s_axilite port=adj3 bundle=control
#pragma HLS INTERFACE s_axilite port=offset4 bundle=control
#pragma HLS INTERFACE s_axilite port=adj4 bundle=control

#pragma HLS INTERFACE s_axilite port=cand1 bundle=control
#pragma HLS INTERFACE s_axilite port=cand2 bundle=control
#pragma HLS INTERFACE s_axilite port=out_buf bundle=control
#pragma HLS INTERFACE s_axilite port=findex bundle=control
#pragma HLS INTERFACE s_axilite port=size1 bundle=control
#pragma HLS INTERFACE s_axilite port=size2 bundle=control
#pragma HLS INTERFACE s_axilite port=index_size bundle=control
#pragma HLS INTERFACE s_axilite port=pid1 bundle=control
#pragma HLS INTERFACE s_axilite port=pid2 bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
	IN_TYPE size_c1 = 0;
	IN_TYPE size_c2 = 0;
	IN_TYPE size_i = 0;

	IN_TYPE pre_adj_count = 0;
	IN_TYPE pre_head = END_FLAG,pre_value = END_FLAG;
	IN_TYPE pre_end = 0,pre_begin = 0,pre_i=0,pre_size=0;

	IN_TYPE pre_adj_count2 = 0;
	IN_TYPE pre_head2 = END_FLAG,pre_value2 = END_FLAG;
	IN_TYPE pre_end2 = 0,pre_begin2 = 0,pre_i2=0,pre_size2=0;

	static IN_TYPE pre_off_buffer[32];
	static IN_TYPE pre_adj_buffer[40960];

	static IN_TYPE pre_off_buffer2[32];
	static IN_TYPE pre_adj_buffer2[81920];

	//static CAN_TYPE block_buffer[512];

	static IN_TYPE preid1;
	static IN_TYPE preid2;
	static IN_TYPE cand1_buffer[1024][16];
	static IN_TYPE cand1_count = 0;
	static IN_TYPE cand1_bf_count[1024];
#pragma HLS BIND_STORAGE variable=cand1_bf_count type=ram_2p
	static IN_TYPE a1_count=0;
	static IN_TYPE a2_count=0;
	static IN_TYPE a1_cache[4096];
	static IN_TYPE a2_cache[4096];
	static CAN_TYPE index_buffer[212256];
#pragma HLS BIND_STORAGE variable=index_buffer type=ram_2p impl=uram
	static IN_TYPE cand2_bf_count = 0;
	static IN_TYPE out_count=0;
	static IN_TYPE cand2_count = 0;
	static IN_TYPE out_iter = 0;
	IN_TYPE out_index,out_cache_count,oi;
	IN_TYPE out_array[10];
#pragma HLS ARRAY_PARTITION variable=out_array dim=1 complete
	IN_TYPE out_count_array[10];
#pragma HLS ARRAY_PARTITION variable=out_count_array dim=1 complete
	IN_TYPE out_cache[102400];
#pragma HLS BIND_STORAGE variable=out_cache type=ram_t2p impl=uram

	IN_TYPE i,j,tmp,tmp2,begin,end,begin2,end2,adj_size,adj_size2,iter,jj,kk,index,a1c,a2c,former1,former2,efcount,efoff,bi,blc_begin,ci,csize;
	int bid,off,bbid,boff,pre_bbid,cbid,coff;
	CAN_TYPE tmpblock,tmpblock2,tmpblock3,tmptype,tmpblock4,tmpblock5,tmpblock6,tmpblock7,tmpblock8,tmpblock9;
	size_c1=0;
	size_c2=0;
	pre_bbid = 0;
	preid1=13;
	preid2=5;
	a1_count=0;
	a2_count=0;
	a1c=0;
	a2c=0;
	out_count=0;
	cand1_count=0;
	cand2_bf_count=0;
	out_iter=0;
	i=0;
	j=0;
	begin=0;
	end=0;
	adj_size=0;
	iter=0;
	jj=0;
	kk=0;
	index=0;
	former1=END_FLAG;
	former2=END_FLAG;
	size_c1 = size1;
	size_c2 = size2;
	size_i = index_size;
	preid1 = pid1;
	preid2 = pid2;
	blc_begin = 0;
	out_index=0;
	out_cache_count=0;
	out_array[0] = 0;
	out_array[1] = 0;
	out_array[2] = 0;
	out_array[3] = 0;
	out_array[4] = 0;
	out_array[5] = 0;
	out_array[6] = 0;
	out_array[7] = 0;
	out_array[8] = 0;
	out_array[9] = 0;
	out_count_array[0] = 0;
	out_count_array[1] = 1;
	out_count_array[2] = 1;
	out_count_array[3] = 1;
	out_count_array[4] = 1;
	out_count_array[5] = 1;
	out_count_array[6] = 1;
	out_count_array[7] = 1;
	out_count_array[8] = 1;
	out_count_array[9] = 1;

	INDEX_LOAD:for(i=0;i<size_i;i++){
		index_buffer[i] = findex[i];
	}

	for(;cand1_count<size1;){

		INITIAL2:for(i=0;i<1024;i++){
			cand1_bf_count[i]=0;
		}

		CAND1_READ:for(ci=0;ci<1024 && cand1_count<size_c1;ci++){
			cand1_buffer[ci][0] = cand1[cand1_count];
			cand1_buffer[ci][1] = cand1[cand1_count + 1];
			cand1_count += 2;
		}

		CAND1_LOAD_OUTTER:for(i=0;i<ci;i++){
			tmp = cand1_buffer[i][0];
			cand1_bf_count[i]++;
			tmp2 = cand1_buffer[i][1];
			cand1_bf_count[i]++;
			if(preid1==13){
				if(pre_head==END_FLAG || tmp - pre_head > 30){
					pre_head = tmp;
					for(pre_i=0;pre_i<32;pre_i++){
						pre_off_buffer[pre_i] = offset2[tmp+pre_i];
					}
					pre_begin = pre_off_buffer[0];
					pre_end = pre_off_buffer[31];
					pre_size = pre_end - pre_begin;
					for(pre_i=0;pre_i<pre_size;pre_i++){
						pre_adj_buffer[pre_i] = adj2[pre_begin+pre_i];
					}
				}
				a1_count=0;
				begin = pre_off_buffer[tmp - pre_head] - pre_begin;
				end = pre_off_buffer[tmp - pre_head + 1] - pre_begin;
				adj_size = end - begin;
				READ_ADJ1:for(iter=0;iter<adj_size;iter++){
					a1_cache[a1_count++] = pre_adj_buffer[iter+begin];
				}

				if(pre_head2 == END_FLAG || tmp2 - pre_head2 > 0){
					pre_head2 = tmp2;
					for(pre_i2=0;pre_i2<2;pre_i2++){
						pre_off_buffer2[pre_i2] = offset1[tmp2+pre_i2];
					}
					pre_begin2 = pre_off_buffer2[0];
					pre_end2 = pre_off_buffer2[1];
					pre_size2 = pre_end2 - pre_begin2;
					for(pre_i2=0;pre_i2<pre_size2;pre_i2++){
						pre_adj_buffer2[pre_i2] = adj1[pre_begin2+pre_i2];
					}
				}
				a2_count=0;
				begin2 = pre_off_buffer2[tmp2 - pre_head2] - pre_begin2;
				end2 = pre_off_buffer2[tmp2 - pre_head2 + 1] - pre_begin2;
				adj_size2 = end2 - begin2;
				READ_ADJ2:for(iter=0;iter<adj_size2;iter++){
					a2_cache[a2_count++] = pre_adj_buffer2[iter+begin2];
				}
			}else{
				if(pre_head==END_FLAG || tmp - pre_head > 10){
					pre_head = tmp;
					for(pre_i=0;pre_i<12;pre_i++){
						pre_off_buffer[pre_i] = offset4[tmp+pre_i];
					}
					pre_begin = pre_off_buffer[0];
					pre_end = pre_off_buffer[11];
					pre_size = pre_end - pre_begin;
					for(pre_i=0;pre_i<pre_size;pre_i++){
						pre_adj_buffer[pre_i] = adj4[pre_begin+pre_i];
					}
				}
				a1_count=0;
				begin = pre_off_buffer[tmp - pre_head] - pre_begin;
				end = pre_off_buffer[tmp - pre_head + 1] - pre_begin;
				adj_size = end - begin;
				READ_ADJ3:for(iter=0;iter<adj_size;iter++){
					a1_cache[a1_count++] = pre_adj_buffer[iter+begin];
				}

				if(pre_head2 == END_FLAG || tmp2 - pre_head2 > 30){
					pre_head2 = tmp2;
					for(pre_i2=0;pre_i2<32;pre_i2++){
						pre_off_buffer2[pre_i2] = offset3[tmp2+pre_i2];
					}
					pre_begin2 = pre_off_buffer2[0];
					pre_end2 = pre_off_buffer2[31];
					pre_size2 = pre_end2 - pre_begin2;
					for(pre_i2=0;pre_i2<pre_size2;pre_i2++){
						pre_adj_buffer2[pre_i2] = adj3[pre_begin2+pre_i2];
					}
				}
				a2_count=0;
				begin2 = pre_off_buffer2[tmp2 - pre_head2] - pre_begin2;
				end2 = pre_off_buffer2[tmp2 - pre_head2 + 1] - pre_begin2;
				adj_size2 = end2 - begin2;
				READ_ADJ4:for(iter=0;iter<adj_size2;iter++){
					a2_cache[a2_count++] = pre_adj_buffer2[iter+begin2];
				}
			}

			for(a1c = 0, a2c = 0; a1c < a1_count && a2c < a2_count;){
				if(a1_cache[a1c]>a2_cache[a2c]){
					a2c++;
				}else if(a1_cache[a1c]<a2_cache[a2c]){
					a1c++;
				}else{
					cand1_buffer[i][cand1_bf_count[i]] = a1_cache[a1c];
					a1c++;
					a2c++;
					cand1_bf_count[i]++;
				}
			}
		}
		//out_cache_count = out_cache_count * out_count_array[out_index];


		for(i=0;i<1024;i++){
			for(jj=2;jj<cand1_bf_count[i];jj++){

				bbid = cand1_buffer[i][jj] / 4096;
				//cbid = cand1_buffer[i][jj] / 1024;
				boff = (cand1_buffer[i][jj] / 512) % 8;
				//coff = cand1_buffer[i][jj] % 8;
				tmpblock4 = index_buffer[bbid];
				//tmpblock7 = cindex_buffer[cbid];
				tmpblock5 = 1<<(7-boff);
				//tmpblock8 = 1<<(7-coff);
				tmpblock6 = tmpblock4&tmpblock5;
				//tmpblock9 = tmpblock7&tmpblock8;
				if(tmpblock5==tmpblock6){
					bid = cand1_buffer[i][jj] / 8;
					off = cand1_buffer[i][jj] % 8;
					tmpblock = cand2[bid];
					tmpblock2 = 1<<(7-off);
					tmpblock3 = tmpblock&tmpblock2;
					if(tmpblock3 == tmpblock2){
						//out_buf[out_count++] = cand1_buffer[i][0];
						//out_buf[out_count++] = cand1_buffer[i][1];
						//out_buf[out_count++] = cand1_buffer[i][jj];
						out_cache[out_cache_count++] = cand1_buffer[i][0];
						out_cache[out_cache_count++] = cand1_buffer[i][1];
						out_cache[out_cache_count++] = cand1_buffer[i][jj];
					}
				}

			}
		}
		out_array[9] = out_cache_count;
		for(oi=0;oi<out_array[out_index];oi++){
			out_buf[out_count] = out_cache[oi];
			out_count++;
		}
		out_index = (out_index + 1) % 10;
		out_cache_count = out_cache_count * out_count_array[out_index];
	}
	for(oi=0;oi<out_cache_count;oi++){
		out_buf[out_count] = out_cache[oi];
		out_count++;
	}
	efcount = out_count % 64;
	efoff = 64 - efcount + 1;
	for(i=0;i<efoff;i++){
		out_buf[out_count] = END_FLAG;
		out_count++;
	}
}
}

extern "C"{
void kernelone(
		IN_TYPE offset1[OFFSET_LEN],IN_TYPE adj1[ADJ_LEN],
		IN_TYPE offset2[OFFSET_LEN],IN_TYPE adj2[ADJ_LEN],
		IN_TYPE cand1[CAN_LEN],CAN_TYPE cand2[CAN2_LEN],CAN_TYPE findex[INDEX_LEN],
		IN_TYPE size1,IN_TYPE size2,IN_TYPE index_size,
		IN_TYPE pid1,CAN_TYPE type1,
		IN_TYPE out_buf[OUT_LEN]
		)
{
#pragma HLS INTERFACE m_axi port=cand1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=cand2 offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=out_buf offset=slave bundle=gmem2

#pragma HLS INTERFACE m_axi port=offset1 offset=slave bundle=gmem3
#pragma HLS INTERFACE m_axi port=adj1 offset=slave bundle=gmem4
#pragma HLS INTERFACE m_axi port=offset2 offset=slave bundle=gmem5
#pragma HLS INTERFACE m_axi port=adj2 offset=slave bundle=gmem6

#pragma HLS INTERFACE m_axi port=findex offset=slave bundle=gmem11

#pragma HLS INTERFACE s_axilite port=offset1 bundle=control
#pragma HLS INTERFACE s_axilite port=adj1 bundle=control
#pragma HLS INTERFACE s_axilite port=offset2 bundle=control
#pragma HLS INTERFACE s_axilite port=adj2 bundle=control

#pragma HLS INTERFACE s_axilite port=cand1 bundle=control
#pragma HLS INTERFACE s_axilite port=cand2 bundle=control
#pragma HLS INTERFACE s_axilite port=out_buf bundle=control
#pragma HLS INTERFACE s_axilite port=findex bundle=control
#pragma HLS INTERFACE s_axilite port=size1 bundle=control
#pragma HLS INTERFACE s_axilite port=size2 bundle=control
#pragma HLS INTERFACE s_axilite port=pid1 bundle=control
#pragma HLS INTERFACE s_axilite port=type1 bundle=control
#pragma HLS INTERFACE s_axilite port=index_size bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

	IN_TYPE size_c1 = 0;
	IN_TYPE size_c2 = 0;
	IN_TYPE size_i = 0;
	IN_TYPE pre_adj_count = 0;
	IN_TYPE pre_head = END_FLAG,pre_value = END_FLAG;
	IN_TYPE pre_end = 0,pre_begin = 0,pre_i=0,pre_size=0;

	static IN_TYPE preid1;

	static IN_TYPE cand1_buffer[1024][16];
	static IN_TYPE cand1_count = 0;
	static IN_TYPE cand1_bf_count[1024];
#pragma HLS BIND_STORAGE variable=cand1_bf_count type=ram_2p

	static IN_TYPE pre_off_buffer[64];
	static IN_TYPE pre_adj_buffer[512];

	//static CAN_TYPE block_buffer[512];

	static CAN_TYPE index_buffer[212256];
#pragma HLS BIND_STORAGE variable=index_buffer type=ram_2p impl=uram
	static IN_TYPE cand2_bf_count = 0;
	static IN_TYPE out_count=0;
	static IN_TYPE cand2_count = 0;
	static IN_TYPE out_iter = 0;
	IN_TYPE out_index,out_cache_count,oi;
	IN_TYPE out_array[10];
#pragma HLS ARRAY_PARTITION variable=out_array dim=1 complete
	IN_TYPE out_count_array[10];
#pragma HLS ARRAY_PARTITION variable=out_count_array dim=1 complete
	IN_TYPE out_cache[102400];
#pragma HLS BIND_STORAGE variable=out_cache type=ram_t2p impl=uram

	IN_TYPE i,j,tmp,begin,end,adj_size,iter,jj,kk,index,efcount,efoff,bi,blc_begin,ci,csize;
	int bid,off,bbid,boff,pre_bbid,cbid,coff;
	CAN_TYPE tmpblock,tmpblock2,tmpblock3,tmptype,tmpblock4,tmpblock5,tmpblock6,tmpblock7,tmpblock8,tmpblock9;

	size_c1 = 0;
	size_c2 = 0;
	pre_bbid = 0;
	preid1 = 13;
	cand1_count = 0;
	cand2_bf_count = 0;
	out_count=0;
	cand2_count=0;
	out_iter = 0;
	i=0;
	j=0;
	begin=0;
	end=0;
	adj_size=0;
	iter=0;
	jj=0;
	kk=0;
	index = 0;
	size_c1 = size1;
	size_c2 = size2;
	size_i = index_size;
	preid1 = pid1;
	tmptype = type1;
	blc_begin = 0;
	out_index=0;
	out_cache_count=0;
	out_array[0] = 0;
	out_array[1] = 0;
	out_array[2] = 0;
	out_array[3] = 0;
	out_array[4] = 0;
	out_array[5] = 0;
	out_array[6] = 0;
	out_array[7] = 0;
	out_array[8] = 0;
	out_array[9] = 0;
	out_count_array[0] = 0;
	out_count_array[1] = 1;
	out_count_array[2] = 1;
	out_count_array[3] = 1;
	out_count_array[4] = 1;
	out_count_array[5] = 1;
	out_count_array[6] = 1;
	out_count_array[7] = 1;
	out_count_array[8] = 1;
	out_count_array[9] = 1;


	INDEX_LOAD:for(i=0;i<size_i;i++){
		index_buffer[i] = findex[i];
	}

	for(;cand1_count<size_c1;){
		//cout<<"a new iteration"<<endl;
		INITIAL2:for(i=0;i<1024;i++){
			cand1_bf_count[i]=0;
		}

		CAND1_READ:for(ci=0;ci<1024 && cand1_count<size_c1;ci++){
			cand1_buffer[ci][0] = cand1[cand1_count];
			cand1_count++;
		}

		CAND1_LOAD_OUTTER:for(i=0;i<ci;i++){
			tmp = cand1_buffer[i][0];
			cand1_bf_count[i]++;
			if(preid1==3){
				if(pre_head==END_FLAG || tmp - pre_head > 62){
					pre_head = tmp;
					for(pre_i=0;pre_i<64;pre_i++){
						pre_off_buffer[pre_i] = offset1[tmp+pre_i];
					}
					pre_begin = pre_off_buffer[0];
					pre_end = pre_off_buffer[63];
					pre_size = pre_end - pre_begin;
					for(pre_i=0;pre_i<pre_size;pre_i++){
						pre_adj_buffer[pre_i] = adj1[pre_begin+pre_i];
					}
				}
				begin = pre_off_buffer[tmp - pre_head] - pre_begin;
				end = pre_off_buffer[tmp - pre_head + 1] - pre_begin;
				adj_size = end - begin;
				READ_ADJ1:for(iter=0;iter<adj_size;iter++){
					cand1_buffer[i][cand1_bf_count[i]] = pre_adj_buffer[iter+begin];
					cand1_bf_count[i]++;
				}
			}else{
				if(pre_head==END_FLAG || tmp - pre_head > 62){
					pre_head = tmp;
					for(pre_i = 0;pre_i < 64;pre_i++){
						pre_off_buffer[pre_i] = offset2[tmp+pre_i];
					}
					pre_begin = pre_off_buffer[0];
					pre_end = pre_off_buffer[63];
					pre_size = pre_end - pre_begin;
					for(pre_i = 0;pre_i < pre_size;pre_i++){
						pre_adj_buffer[pre_i] = adj2[pre_begin+pre_i];
					}
				}
				begin = pre_off_buffer[tmp - pre_head] - pre_begin;
				end = pre_off_buffer[tmp - pre_head + 1] - pre_begin;
				adj_size = end - begin;
				READ_ADJ2:for(iter=0;iter<adj_size;iter++){
					cand1_buffer[i][cand1_bf_count[i]] = pre_adj_buffer[iter+begin];
					cand1_bf_count[i]++;
				}
			}
		}

		//out_cache_count = out_cache_count * out_count_array[out_index];

		SEARCH_OUTTER:for(i=0;i<1024;i++){
			SEARCH_INNER:for( jj = 1;jj < cand1_bf_count[i] ; jj++ ){
				bbid = cand1_buffer[i][jj] / 4096;
				//cbid = cand1_buffer[i][jj] / 1024;
				boff = (cand1_buffer[i][jj] / 512) % 8;
				//coff = cand1_buffer[i][jj] % 8;
				tmpblock4 = index_buffer[bbid];
				//tmpblock7 = cindex_buffer[cbid];
				tmpblock5 = 1<<(7-boff);
				//tmpblock8 = 1<<(7-coff);
				tmpblock6 = tmpblock4&tmpblock5;
				//tmpblock9 = tmpblock7&tmpblock8;
				if(tmpblock6==tmpblock5){
					bid = cand1_buffer[i][jj] / 8;
					off = cand1_buffer[i][jj] % 8;
					tmpblock = cand2[bid];
					tmpblock2 = 1<<(7-off);
					tmpblock3 = tmpblock&tmpblock2;
					if(tmpblock3 == tmpblock2){
						//out_buf[out_count++] = cand1_buffer[i][0];
						//out_buf[out_count++] = cand1_buffer[i][jj];
						out_cache[out_cache_count++] = cand1_buffer[i][0];
						out_cache[out_cache_count++] = cand1_buffer[i][jj];
					}
				}
			}
		}

		out_array[9] = out_cache_count;
		for(oi=0;oi<out_array[out_index];oi++){
			out_buf[out_count] = out_cache[oi];
			out_count++;
		}
		//memcpy(out_buf_1 + out_1_count,(const IN_TYPE*)out_cache,out_array[out_index]*sizeof(IN_TYPE));
		//out_1_count += out_array[out_index];
		out_index = (out_index + 1) % 10;
		out_cache_count = out_cache_count * out_count_array[out_index];
	}
	for(oi=0;oi<out_cache_count;oi++){
		out_buf[out_count] = out_cache[oi];
		out_count++;
	}
	efcount = out_count % 64;
	efoff = 64 - efcount + 1;
	for(i=0;i<efoff;i++){
		out_buf[out_count] = END_FLAG;
		out_count++;
	}
}
}


extern "C"{
void multiIntersect(unsigned int* list1, unsigned int* list2,
		            unsigned int* list3, unsigned int* list4,
					int size1, int size2, int size3, int size4,
					unsigned int* out){
#pragma HLS INTERFACE m_axi port = list1 offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi port = list2 offset = slave bundle = gmem2
#pragma HLS INTERFACE m_axi port = list3 offset = slave bundle = gmem3
#pragma HLS INTERFACE m_axi port = list4 offset = slave bundle = gmem4
#pragma HLS INTERFACE m_axi port = out offset = slave bundle = gmem5
#pragma HLS INTERFACE s_axilite port = list1 bundle = control
#pragma HLS INTERFACE s_axilite port = list2 bundle = control
#pragma HLS INTERFACE s_axilite port = list3 bundle = control
#pragma HLS INTERFACE s_axilite port = list4 bundle = control
#pragma HLS INTERFACE s_axilite port = out bundle = control
#pragma HLS INTERFACE s_axilite port = size1 bundle = control
#pragma HLS INTERFACE s_axilite port = size2 bundle = control
#pragma HLS INTERFACE s_axilite port = size3 bundle = control
#pragma HLS INTERFACE s_axilite port = size4 bundle = control

#pragma HLS STREAM variable = mergeOut depth = 64
#pragma HLS STREAM variable = readOut depth = 32
#pragma HLS STREAM variable = bsearchOut_1 depth = 32
#pragma HLS STREAM variable = bsearchOut_2 depth = 32

#pragma HLS dataflow
	read_csr(list1,list2,list3,list4,size1,size2,size3,size4);
	mergeTwo();
	binarySearch(l3, mergeOut, bsearchOut_1);
	binarySearch(l4, bsearchOut_1, bsearchOut_2);
	write_result(out, bsearchOut_2, 1000000);

}
}
