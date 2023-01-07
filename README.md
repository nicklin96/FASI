# ICDE23_SUBMISSION_905
FASI: FPGA-friendly Subgraph Isomorphism on Massive Graphs

## Input File Format
The graph file format for data graphs and query graphs is a text format to store an vertex-labeled undirected graph.
* The first line of the file should be "t #vertices #edges".
- Following lines of "v vertex-ID vertex-label" indicate the vertices in the graph.

- The vertices should be written in the file in ascending order of their IDs, and a vertex ID should be in [0, #vertices - 1].
* Following lines of "e src-ID dst-ID edge-ID" after the vertices indicate the undirected edges in the graph.

For example:

```
t 6 7
v 0 0
v 1 1
v 2 2
v 3 3
v 4 4
v 5 5
e 0 1 0
e 1 2 0
e 1 3 0
e 1 4 0
e 1 5 0
e 2 3 0
e 2 4 0
```

## make
```
cd src
make
```

## Data Preprocessing

### Graph Reordering
````./LPGO [original_graph] [reordered_graph]```` will reorder a ````[original_graph]```` with above format into a ````[reordered_graph]````.

For example:
```
cd src
./LPGO ../data/patents.nt ../data/patents.LPGO.nt
```

### Generate LPCSR
Call FASI with option ````-p```` to generate the LPCSR structure.

For example:
```
cd src
./FASI -p ../data/patents.LPGO.nt
```
After preprocessing, we get the LPCSR structure ````patents.LPCSR.nt````.

## Generate query
````./QueryGenerator [query_num] [vertex_num]```` is responsible for generating ````[query_num]```` queries with ````[vertex_num]```` query vertices and ````./LabelGenerator [graph_name] [label_num]```` takes charge of assigning labels within ````[label_num]```` for the graph ````[graph_name]````.

For example:
```
cd src
./QueryGenerator 5 6 
./LabelGenerator patents 20
```

## Compiling the FPGA processing units
Under ./src/FPGA_kernel is the C++-based project export from Vitis.

To compile the binary container containing FPGA processing units, first you should import the project into Vitis.

Afterwards, project configuration should be loaded and you can immediately launch the HLS and Vivado compiler.

The following Xilinx runtime and platforms are required (Avaliable at https://www.xilinx.com/products/boards-and-kits/alveo/u200.html#gettingStarted):

Xilinx runtime : xrt_202210.2.13.466_7.8.2003

Deployment target platform : xilinx-u200-gen3x16-xdma_2022.1_2022_0415_2123

Development target Platform : xilinx-u200-gen3x16-xdma-2-202110-1-dev-1-3514848

## Run
Call FASI with option ````-q```` to run a query using the LPCSR structure.

For example:
```
cd src
./FASI -q ../data/patents.LPCSR.nt ../query/[query_name] [output_limit]
```
where ````[query_name]```` is a query file and ````[output_limit]```` is the upper limit of the output results.

## Compared algorithms
Five of our compared algorithms are open-sourced on github, you can find them with the following links. Thanks for the authors sharing!:

DAF: https://github.com/SNUCSE-CTA/DAF

CECI: https://github.com/iHeartGraph/ceci-release

RAPIDMATCH: https://github.com/RapidsAtHKUST/RapidMatch

GSI: https://github.com/pkumod/GSI

GPSM: https://github.com/bookug/GpSM

The FPGA-based compared algorithms, FAST, does not have public code. We emailed the authors and managed to get the code. You can contact the authors via the email address in this paper: 

Jin, X., Yang, Z., Lin, X., Yang, S., Qin, L., & Peng, Y. (2021, April). Fast: Fpga-based subgraph matching on massive graphs. In 2021 IEEE 37th International Conference on Data Engineering (ICDE) (pp. 1452-1463). IEEE.
