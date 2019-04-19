/*
 * Copyright 2011, 2012 Medha Atre
 * 
 * This file is part of BitMat.
 * 
 * BitMat is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * BitMat is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with BitMat.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Report any bugs or feature requests to <firstname>.<lastname>@gmail.com
 * ("firstname" and "lastname" of the first copyright holder above).
 */

#ifndef _BITMAT_H_
#define _BITMAT_H_

#include <sys/types.h>
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <unordered_set>
#include <stack>
#include <algorithm>
#include <fstream>
#include <list>
#include <set>
#include <cstring>
#include <string>
#include <cmath>
#include <fcntl.h>
#include <stdio.h>
//#include <math.h>
#include <sys/mman.h>
#include <assert.h>
//#include <bits/stdc++.h>
#include <thread>

#include <sys/types.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <time.h>

using namespace std;

#define REVERSEMAP			1
#define REVERSEMAP_THR			0
#define REVERSEMAP_POST			!REVERSEMAP_THR

#define SORT_MEM_GB			2

#define SUB_EDGE			1001
#define PRED_EDGE			1002
#define OBJ_EDGE			1003
#define SO_EDGE				1004
#define NORM_BRACE			1005
#define OPT_BRACE			1006

/*
 * ROW_VEC - BitMat's row vec only uptodate
 * ROW_UNFOLD - column vec out of date (BitMat
 * 				and row vec uptodate)
 * Similar terminology for COLUMN_VEC/UNFOLD
 * UPTODATE - everything uptodate
 */
#define ROW_VEC				2001 //only row vec uptodate
#define ROW_UNFOLD			2002 //row vec uptodate and unfolded
#define COLUMN_VEC			2003 //only column vec uptodate
#define COLUMN_UNFOLD		2004 //column vec uptodate and unfolded
#define BMUPTODATE			2005 //both row, col vec uptodate and unfolded, except num-triples
#define ALLUPTODATE			2006 //everything uptodate

#define BM_ROW_SIZE			4 //sizeof(unsigned char)
#define ROW_SIZE_BYTES 		4
#define GAP_SIZE_BYTES		4
#define TOTAL_ROWSIZE(c)	(((c-1)*GAP_SIZE_BYTES)+1)
#define ROW_SIZE_HEADER(c)		(((c-1)/GAP_SIZE_BYTES)+1)
#define IS_ROW_SIZE_SHORT(c)			((unsigned int)(((c-1)/GAP_SIZE_BYTES)+1) < 0xffff)
//#define ROW_SIZE(c) ((c[0] & 0x80) ? 4:2)
//#define GAP_SIZE_BYTES		sizeof(unsigned int)
//#define ROWSIZE(c)			(c-1)*GAP_SIZE_BYTES+1


#define SUB_DIMENSION 		3001
#define PRED_DIMENSION 		3002
#define OBJ_DIMENSION 		3003
#define SPO_BITMAT			3004
#define OPS_BITMAT			3005
#define PSO_BITMAT			3006
#define POS_BITMAT			3007
#define SOP_BITMAT			3008
#define OSP_BITMAT			3009
#define SO_DIMENSION		3010
#define UNDECIDED			3011

#define JVAR_NODE			4001
#define TP_NODE				4002
#define	NON_JVAR_NODE		4003

#define TRIPLE_STR_SPACE    50
#define ROW					10012
#define COLUMN				10013

#define MMAPFILES			1

#define SELECTIVITY_THRESHOLD	128

#define USE_MORE_BYTES		0

#define MAX_TPS_IN_QUERY		32														  
#define MAX_JVARS_IN_QUERY		32
#define	MAX_JVARS			-MAX_JVARS_IN_QUERY
#define MAX_NODES_IN_GRAPH		MAX_JVARS_IN_QUERY + MAX_TPS_IN_QUERY

#define	ASSERT_ON			0
#define	REALLOC_ON			0
#define DEBUG				0

#define WHITE				10
#define GREY				11
#define BLACK				12

#define NUMTHREADS			8

#define RESBITVEC			0

#define MICROSEC 		1000000.0

extern map <std::string, std::string> config;
extern unsigned int gnum_subs, gnum_preds, gnum_objs, gnum_comm_so,
gsubject_bytes, gobject_bytes, gpredicate_bytes, gcommon_so_bytes, grow_size;
extern unsigned int row_size_bytes, gap_size_bytes;
extern unsigned int comp_folded_arr;
extern unsigned int table_col_bytes;
extern unsigned char *grow;

extern unsigned int graph_tp_nodes, graph_tp_nodes_new;
extern unsigned int graph_jvar_nodes, graph_var_nodes;
extern unsigned char **subjmapping, **predmapping, **objmapping;
extern unsigned char *distresult;
extern unsigned long distresult_size;
extern map<int, struct node *> selvars;
extern map<int, struct node *> seljvars;
extern map<string, char *> mmap_table;
//extern set<struct node *> unvisited;
//extern bool resbitvec, isCycle;
extern set<string> null_perm;
////////////////////////////////////////////////////////////
// Data structures for multi-join graph

typedef unordered_set<std::string> ResultSet;
extern ResultSet resset;

struct row {
	unsigned int rowid; //nodeid from the orig graph
	unsigned char *data;
	bool operator<(const struct row &a) const
	{
		return (this->rowid < a.rowid);
	}
};

typedef struct BitMat {
	list<struct row> bm;
	vector<struct row> vbm;
//	bool single_row;
//	unsigned int num_subs, num_preds, num_objs, num_comm_so;
	unsigned int num_rows, num_totalBMs, num_columns, num_comm_so;
//	unsigned int row_size_bytes, gap_size_bytes;
//	unsigned int subject_bytes, predicate_bytes, object_bytes, common_so_bytes, dimension;//, last_unfold;
	unsigned int row_bytes, totalBMs_bytes, column_bytes, common_so_bytes, dimension;//, last_unfold;
	unsigned long num_triples;
	unsigned char *rowfold; //row_bytearr
	unsigned char *colfold; //column_bytearr
	unsigned int last_op;
//	unsigned char *subfold_prev;
//	unsigned char *objfold_prev;

	void freebm(bool list = true, bool vector = true)
	{
		if (list) {
			for (std::list<struct row>::iterator it = bm.begin(); it != bm.end(); ){
				free((*it).data);
				it = bm.erase(it);
			}
		}
		if (vector) {
			for (std::vector<struct row>::iterator it = vbm.begin(); it != vbm.end(); ){
				free((*it).data);
				it = vbm.erase(it);
			}
		}

		if (rowfold != NULL) {
			free(rowfold);
			rowfold = NULL;
		}

		if (colfold != NULL) {
			free(colfold);
			colfold = NULL;
		}
		num_triples = 0;
	}

//	void freerows(void)
//	{
//		for (std::list<struct row>::iterator it = bm.begin(); it != bm.end(); ){
//			free((*it).data);
//			it = bm.erase(it);
//		}
//	}

	void reset(void)
	{
		for (std::list<struct row>::iterator it = bm.begin(); it != bm.end(); ){
			free((*it).data);
			it = bm.erase(it);
		}
		if (rowfold != NULL) {
			memset(rowfold, 0, row_bytes);
		}
		if (colfold != NULL) {
			memset(colfold, 0, column_bytes);
		}
		num_triples = 0;
	}

	void clone_to(BitMat *second)
	{
		second->num_rows = num_rows;
		second->num_totalBMs = num_totalBMs;
		second->num_columns = num_columns;
		second->num_comm_so = num_comm_so;
		second->dimension = dimension;
		second->row_bytes = row_bytes;
		second->totalBMs_bytes = totalBMs_bytes;
		second->column_bytes = column_bytes;
		second->common_so_bytes = common_so_bytes;

		if (second->rowfold != NULL) free(second->rowfold);
		if (second->colfold != NULL) free(second->colfold);

		second->rowfold =  (unsigned char *) malloc(row_bytes * sizeof(unsigned char));
		memcpy(second->rowfold, rowfold, row_bytes);
		second->colfold = (unsigned char *) malloc(column_bytes * sizeof(unsigned char));
		memcpy(second->colfold, colfold, column_bytes);
		second->num_triples = num_triples;

		for (std::list<struct row>::iterator it = bm.begin(); it != bm.end(); it++){
			unsigned char *old_data = (*it).data;

			unsigned int rw_size = 0;
			memcpy(&rw_size, old_data, BM_ROW_SIZE);

			unsigned char *data = (unsigned char *) malloc(rw_size);
			memcpy(data, old_data, rw_size);
			struct row r = {(*it).rowid, data};

			second->bm.push_back(r);
		}
		for (std::vector<struct row>::iterator it = vbm.begin(); it != vbm.end(); it++){
			unsigned char *old_data = (*it).data;

			unsigned int rw_size = 0;
			memcpy(&rw_size, old_data, BM_ROW_SIZE);

			unsigned char *data = (unsigned char *) malloc(rw_size);
			memcpy(data, old_data, rw_size);
			struct row r = {(*it).rowid, data};

			second->vbm.push_back(r);
		}
	}

	BitMat()
	{
		num_rows = num_totalBMs = num_columns = num_comm_so = num_triples = 0;
		rowfold = colfold = NULL;
	}

	~BitMat()
	{
		for (std::list<struct row>::iterator it = bm.begin(); it != bm.end(); ){
			free((*it).data);
			it = bm.erase(it);
		}
		if (rowfold != NULL)
			free(rowfold);
		if (colfold != NULL)
			free(colfold);
//		free(subfold_prev);
//		free(objfold_prev);
	}

} BitMat;

extern BitMat bmorig;

struct triple {
	unsigned int sub;
	unsigned int pred;
	unsigned int obj;
};

struct twople {
	unsigned int sub;
	unsigned int obj;
};

struct level {
	int level;
	char ch;
};

typedef struct blevel {
	unsigned char *joinres;
	unsigned int joinres_size;
	int joinres_dim;
	bool empty;
//	int m_level, s_level;
//	stack<struct strlevel> level;
	string strlevel;
} BitVecJoinres;

typedef struct triplepattern {
	int nodenum; // unique in the graph
	int sub;
	int pred;
	int obj;
	string substr;
	string predstr;
	string objstr;
//	deque<struct level> level;
//	deque<string> pseudomaster;
//	deque<string> real_masters;
	string strlevel;
//	string old_level;
	BitMat bitmat;
//	BitMat bitmat2;
	int unfolded;
	int last_unfold_dim;
//	bool rollback;
	int numpass;
	unsigned int triplecnt;
	unsigned short numjvars;
	unsigned short numvars;
	struct node *gtp;

	bool containsJvar(int jvar) {
		if (sub > MAX_JVARS && sub == jvar)
			return true;
		if (obj > MAX_JVARS &&  obj == jvar)
			return true;
		if (pred > MAX_JVARS && pred == jvar)
			return true;
		return false;
	}

//	string deque_to_string()
//	{
//		strlevel.clear();
//		for (deque<struct level>::iterator it = this->level.begin(); it != this->level.end(); it++) {
//			char str[20];
//			sprintf(str, "%c%d", (*it).ch, (*it).level);
//			strlevel.append(str);
//		}
//		return strlevel;
//	}
	bool isSlave(string that) {
		if (this->strlevel.find(that) == 0) {
			return true;
		}
		return false;
	}
	bool isSlaveOfAnyone(set<string> thatvec)
	{
		for (set<string>::iterator it = thatvec.begin(); it != thatvec.end(); it++) {
			if (this->isSlave(*it))
				return true;
		}
		return false;
	}
	std::string toString()
	{
		char str[1024]; sprintf(str, "[%d %d %d] [%s %s %s] %s bitmat.dim %d", this->sub, this->pred,
				this->obj, this->substr.c_str(), this->predstr.c_str(), this->objstr.c_str(),
					this->strlevel.c_str(), this->bitmat.dimension);
		return string(str);
	}

	std::string plainString()
	{
		char str[1024]; sprintf(str, "[%d %d %d]", this->sub, this->pred, this->obj);
		return string(str);
	}

} TP;

//extern map<struct node*, struct triple> q_to_gr_node_map;  // Maps query node to the graph triple.

typedef struct var {
	int nodenum;
//	unsigned int val;
	struct node *gjvar;
//	bool sojoin;
	//this is for s-o join
	//in case of s-o join if S dimension gets 
	//evaluated first, then joinres points to
	//subject array o/w to compressed obj array
	//joinres_so_dim indicates the dimension of joinres_so
//	unsigned char *joinres;
//	int joinres_dim;
//	unsigned char *joinres_so;
	vector<BitVecJoinres*> joinres;
//	int joinres_so_dim;
//	unsigned int joinres_size;
//	unsigned int joinres_so_size;
	unsigned int upperlim, lowerlim, range;

	std::string toString()
	{
		char str[1024];
		sprintf(str, "%d %d", nodenum);

		return string(str);
	}

} JVAR;

typedef struct linkedlist {
	struct node *gnode;
	int edgetype;
	struct linkedlist *next;
} LIST;

struct node {
	int type;
	int level;
	//following 2 vars are for cycle detection algo
	int color;
	bool isNeeded;
	void *data;
	LIST *nextadjvar;
	unsigned int numadjvars;
	unsigned int numadjjvars;
	LIST *nextTP;
	unsigned int numTPs;

	bool isNeighbor(struct node *other)
	{
		LIST *neighbor = NULL;

		if (other->type == TP_NODE) {
			neighbor = nextTP;
		} else if (other->type == JVAR_NODE || other->type == NON_JVAR_NODE) {
			neighbor = nextadjvar;
		} else
			assert(0);
		for (; neighbor; neighbor = neighbor->next) {
			if (neighbor->gnode == other) {
				return true;
			}
		}
		return false;
	}

	bool isNeighborSameLevel (struct node *other)
	{
		if (type == TP_NODE && other->type == TP_NODE) {
			if (((TP*)data)->strlevel == ((TP*)other->data)->strlevel) {
				return isNeighbor(other);
			}
		}
		return false;
	}
};

extern struct node graph[MAX_NODES_IN_GRAPH];
extern struct node *jvarsitr2[MAX_JVARS_IN_QUERY];
extern map<struct node*, vector <struct node*> > tree_edges_map;  // for tree edges in the query graph.

inline unsigned int set_size(unsigned int i, char *s)
{
	assert(i <= 0x80000000);
	if (i > 0x8000) {
		//large size
		s[0] = (i >> 24) & 0xFF;
		s[1] = (i >> 16) & 0xFF;
		s[2] = (i >> 8) & 0xFF;
		s[3] = i & 0xFF;
		s[0] |= 0x80;
//		printf("2 -- %x %x %x %x\n", s[0], s[1], s[2], s[3]);
		return 4;
	} else {
		//small size
		s[0] = (i >> 8) & 0xFF;
		s[1] = i & 0xFF;
		return 2;
//		printf("2.2 -- %x %x\n", s[0], s[1]);
	}
}

inline unsigned int get_size(char *s)
{
	unsigned int j = 0;
	if (s[0] & 0x80) {
		//large size
		s[0] &= 0x7f;
		j = (unsigned int) s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
	} else {
		//small size
		j = (unsigned short) s[0] << 8 | s[1];
	}
	return j;
}

void parse_config_file(char *fname);

void init_bitmat(BitMat *bitmat, unsigned int snum, unsigned int pnum, unsigned int onum, unsigned int commsonum, int dimension);

void shallow_init_bitmat(BitMat *bitmat, unsigned int snum, unsigned int pnum, unsigned int onum, unsigned int commsonum, int dimension);

unsigned int dgap_compress(unsigned char *in, unsigned int size, unsigned char *out);

unsigned int dgap_compress_new(unsigned char *in, unsigned int size, unsigned char *out);

unsigned long count_bits_in_row(unsigned char *in, unsigned int size);

void dgap_uncompress(unsigned char *in, unsigned int insize, unsigned char *out, unsigned int outsize);

unsigned int concatenate_comp_arr(unsigned char *in1, unsigned int size1,
								unsigned char *in2, unsigned int size2);

void map_to_row_wo_dgap(unsigned char **in, unsigned int pbit, unsigned int obit,
				unsigned int spos, bool cflag, bool start);

unsigned long count_triples_in_row(unsigned char *in, unsigned int size);

unsigned long count_triples_in_bitmat(BitMat *bitmat);

unsigned int fold(unsigned char **ts, int ret_dimension, unsigned char *foldarr);

void simple_fold(BitMat *bitmat, int ret_dimension, unsigned char *foldarr, unsigned int size);

void simple_unfold(BitMat *bitmat, unsigned char *maskarr, unsigned int maskarr_size, int ret_dimension);

unsigned int count_size_of_bitmat(unsigned char **bitmat);
unsigned long get_size_of_bitmat(int dimension, unsigned int node);

void list_triples_in_bitmat(unsigned char **bitmat, unsigned int dimension, unsigned long num_triples, char *triplefile);

void go_over_triples_in_row(unsigned char *in, unsigned int rownum, struct node *n,
							int curr_node, struct node *parent, int in_dimension, ofstream& outfile);

unsigned long list_enctriples_in_bitmat(unsigned char **bitmat, unsigned int dimension, unsigned int num_triples, char *triplefile);

void init_TP_nodes(TP **tplist, unsigned int listsize);

bool cycle_detect(struct node *node);

void build_graph(TP **tplist, unsigned int tplist_size,
					JVAR **jvarlist, unsigned int jvarlist_size);

void build_graph_new();

void build_graph_and_init(char *fname);

bool parse_query_and_build_graph(char *fname);

unsigned int load_from_dump_file(char *fname_dump, unsigned int bmnum, BitMat *bitmat,
		bool readtcnt, bool readarray, unsigned char *maskarr, unsigned int mask_size, int maskarr_dim,
		char *fpos, int fd, bool fold_objdim);

unsigned int load_from_dump_file2(char *fname_dump, unsigned int bmnum, BitMat *bitmat, bool readtcnt,
			bool readarray, unsigned char *maskarr, unsigned int maskarr_size, int maskarr_dim,
			char *fpos, int fd, bool loadvbm);

void spanning_tree(struct node *n, int curr_node, int nodetype);

vector<struct node *>& spanning_tree_tpnodes(void);

void spanning_tree_bfs(struct node *n, int nodetype);

struct node *optimize_and_get_start_node();

struct node *get_start_node_subgraph_matching();

void print_mem_usage();

void clear_rows(BitMat *bitmat, bool clearbmrows, bool clearfoldarr, bool optimize);

void load_data_vertically(char *file, vector<struct twople> &triplelist, BitMat *bitmat, char *fname_dump, bool ondisk, bool invert, bool loadlist, char *tmpfile);

bool add_row(BitMat *bitmat, char *fname, unsigned int bmnum, unsigned int readdim, unsigned int rownum, bool load_objfold);

bool filter_and_load_bitmat(BitMat *bitmat, int fd, char *fpos, unsigned char *and_array, unsigned int and_array_size);

bool init_tp_nodes_new(bool bushy);

unsigned int get_and_array(BitMat *bitmat, unsigned char *and_array, unsigned int bit);

unsigned long get_offset(char *fname, unsigned int bmnum);
	
void print_node_info(struct node *n);

void print_spanning_tree(int nodetype);

void build_jvar_tree(bool best_match_reqd, bool iscylic, bool allslaves_with_one_jvar);

bool prune_triples_new(bool bushy, bool opt, vector<vector<struct node *>> &cycles,
						bool best_match_reqd, bool allslaves_one_jvar);

bool populate_all_tp_bitmats();

//void match_query_graph_new(FILE **outfile, FILE **intfile, FILE **intfile_preds, struct node **bfsarr, int sidx, int eidx, char *null_pad_str, bool bestm,
//		map<struct node *, struct triple> &q_to_gr_node_map, map<struct node *, struct triple> &copy_of_resmap);

void list_enctrips_bitmat_new(BitMat *bitmat, unsigned int bmnum, vector<twople> &twoplist, FILE *outfile);

void list_enctrips_bitmat2(BitMat *bitmat, unsigned int bmnum, vector<twople> &twoplist, FILE *outfile);

void load_mappings(char *subjfile, char *predfile, char *objfile);

void cleanup_graph(void);

void fix_all_tp_bitmats(void);

void fix_all_tp_bitmats2(vector<struct node *> &tplist);

void init_bitmat_rows(BitMat *bitmat, bool initbm, bool initfoldarr);

unsigned long count_size_of_bitmat(BitMat *bitmat);

void list_all_data(int dimension, unsigned int, char *outfile);

void print_graph(void);

bool update_graph();

void print_results(FILE *);

bool mmap_all_files();
bool munmap_all_files();

bool update_graph_for_final(void);

unsigned int wrapper_load_from_dump_file(char *fname_dump, unsigned int bmnum, struct node *, bool readtcnt, bool readarray);

unsigned int wrapper_load_from_dump_file2(BitMat *bitmat, unsigned int bmnum);

void get_all_triples_in_query(char *file);

//void identify_true_jvars(void);

void print_stats(char *fname_dump, unsigned int numbms, bool compfold, unsigned int numsubs, unsigned int numobjs);

void get_all_triples_in_query2(char *file);

void test_new_simple_fold();

unsigned int print_set_bits_in_row(unsigned char *in, unsigned int size);

void enumerate_cycles(struct node *jvar, struct node *parent,
			vector<struct node *> &path, set<struct node *> &visited, vector<vector<struct node*>> &cycles);

bool get_longest_cycles(vector<vector<struct node*>> &cycles);

void mark_needed_tps(vector<vector<struct node*>> &cycles);

void test_cycle_detection();

void best_match_subsumption(char *file);

void best_match_postproc(char *file, char *outfile);

bool is_allslave_bottomup_pass(void);

bool slaves_with_more_than_one_jvar();

bool init_tp_triplecnt();

void count_all_triples_bm();

void set_all_bm_lastop(unsigned int lastop);

unsigned int set_null_padding(vector<struct node *> &tplist);

void start_multiway_join(char *outfile, vector<struct node *> &tplist, bool bestm);

bool build_semijoins_levelwise(map<string, vector<struct node *>> &);

bool prune_levelwise(map<string, vector<struct node *>> levelwise_semij);

void dict_initialize(int forward = 1, int reverse = 1, char *path = NULL, bool nodes = true);

void dict_close(int forward = 1, int reverse = 1, bool nodes = true);

void dict_forward_lookup(vector<string> &input, vector<unsigned int> &output, bool nodes);

void dict_reverse_lookup(vector<unsigned int> &input, unsigned int minint, unsigned int start_idx,
		unsigned int end_idx, unsigned char **output, bool nodes, bool file);

void conditional_fold(BitMat *bitmat, int ret_dimension, unsigned char *foldarr, unsigned int foldarr_size);

bool add_row2(BitMat *bitmat, unsigned int dimension, unsigned int bmnum, unsigned int rownum);

void read_bitmat_triplecnt(BitMat *bitmat, unsigned int dimension, unsigned int bmnum);

bool parse_query_new(char *fname);

unsigned char * get_maskbitarr(BitMat *bitmat, int dim, unsigned int *size);

#endif

//void cumulative_dgap(unsigned char *in, unsigned int size, unsigned char *out);
//void de_cumulative_dgap(unsigned char *in, unsigned int size, unsigned char *out);
//void bitmat_cumulative_dgap(unsigned char **bitmat);
//void bitmat_de_cumulative_dgap(unsigned char **bitmat);
//unsigned int dgap_AND(unsigned char *in1, unsigned int size1,
//					unsigned char *in2, unsigned int size2,
//					unsigned char *res);
//unsigned int fixed_p_fixed_o (unsigned char *in, unsigned int size, unsigned char *out,
//								unsigned int ppos, unsigned int opos);
//unsigned int fixed_p_var_o(unsigned char *in, unsigned int size, unsigned char *out,
//							unsigned int ppos);
//unsigned int var_p_fixed_o(unsigned char *in, unsigned int size, unsigned char *out,
//							unsigned int opos, unsigned char *and_array, unsigned int and_array_size);
//unsigned char **filter(unsigned char **ts, unsigned int sub, unsigned int pred, unsigned int obj);
//void unfold(BitMat *bitmat, unsigned char *maskarr, unsigned int maskarr_size, int ret_dimension);
//void map_to_row_wo_dgap_ops(unsigned char **in, unsigned int pbit, unsigned int obit,
//				unsigned int spos, bool cflag, bool start);
//void match_query(struct node* n, int curr_node, struct node *parent, ofstream& outfile);
//unsigned char **load_ops_bitmat(char *fname, BitMat *ops_bm, BitMat *spo_bm);
//bool do_fold_unfold_ops(struct node *gnode, unsigned int gnode_idx, unsigned int loopcnt);
//void multi_join(char *fname);
//void load_data(char *file);
//vector<struct node *> final_res_gen_tp_seq();
//unsigned int count_intersect(unsigned int bm1, unsigned int bm2);
//unsigned int count_intersect_oo(unsigned int bm1, unsigned int bm2);
//unsigned int count_intersect_ss(unsigned int bm1, unsigned int bm2);
//void match_query_graph(FILE *outfile, struct node **tparr, int, int, set<string> &null_tmp, bool opt);
