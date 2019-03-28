/*
 * Copyright 2010-2014 Medha Atre
 * 
 * This file is part of BitMat.
 * 
 * BitMat is a free software: you can redistribute it and/or modify
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

#include "bitmat.hpp"

unsigned int gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, gsubject_bytes, gobject_bytes, gpredicate_bytes, gcommon_so_bytes;//, grow_size;
unsigned int row_size_bytes, gap_size_bytes, bm_row_size;
unsigned int comp_folded_arr;
map <std::string, std::string> config;
bool resbitvec = false;

////////////////////////////////////////////////////////////
int main(int args, char **argv)
{
	struct timeval start_time, stop_time, start_prune, start_final, start_init;
//	clock_t  t_start, t_end;
	double curr_time;
	double st, en, prune, init;
	int c = 0;
	unsigned int bmnum = 0;
	char qfile[25][124], outfile[25][124];
	bool loaddata = false, querydata = true, pruning = false;
	char *dump_file = NULL;

	int q_count = 0, op_count = 0;

	if (args < 11) {
		printf("Copyright 2011, 2012 Medha Atre\n\n");
//		printf("Usage: bitmat -l [y/n] -Q [y/n] -f config-file -q query-file -o res-output-file -t resultsettype(b: bitvec, h:hashmap)\n");
		printf("Usage: bitmat -l [y/n] -Q [y/n] -f config-file -q query-file -o res-output-file -b bitmatnum\n\n");
		exit (-1);
	}

	printf("Copyright 2011, 2012 Medha Atre\n\n");

	while((c = getopt(args, argv, "t:l:Q:f:q:o:b:")) != -1) {
		switch (c) {
			case 'f':
				parse_config_file(optarg);
				break;
			case 'l':
				if (!strcmp(optarg, "y")) {
					loaddata = true;
				}
				break;
			case 'Q':
				if (!strcmp(optarg, "n")) {
					cout << "Query option" << endl;
					querydata = false;
				}
				break;
			case 'q':
				strcpy(qfile[q_count], optarg);
				q_count++;
				break;
			case 'o':
				strcpy(outfile[op_count], optarg);
				op_count++;
				break;
			case 't':
				if (!strcmp(optarg, "b")) {
					resbitvec = true;
				} else if (!strcmp(optarg, "h")) {
					resbitvec = false;
				} else
					assert(0);
				assert(resbitvec == RESBITVEC);
				break;
//			case 'd':
//				dump_file = (char *)malloc (strlen(optarg));
//				strcpy(dump_file, optarg);
//				break;
			case 'b':
				bmnum = strtoul(optarg, NULL, 10);
				break;
			default:
				printf("Usage: bitmat -f config-file -q query-file -o res-output-file\n");
				exit (-1);

		}
	}

	printf("Process id = %d\n", (unsigned)getpid());


	if (row_size_bytes != ROW_SIZE_BYTES || gap_size_bytes != GAP_SIZE_BYTES || bm_row_size != BM_ROW_SIZE) {
		cerr << "**** ERROR: Descrepancy in the row/gap_size_bytes values" << endl;
		cerr << "row_size_bytes " << row_size_bytes << " ROW_SIZE_BYTES " << ROW_SIZE_BYTES <<
			" gap_size_bytes " << gap_size_bytes << " GAP_SIZE_BYTES " << GAP_SIZE_BYTES <<
			" bm_row_size " << bm_row_size << " BM_ROW_SIZE " << BM_ROW_SIZE << endl;
		exit(-1);
	}
	
	if (loaddata) {

		gettimeofday(&start_time, (struct timezone *)0);

		vector<struct twople> triplelist;

		cout << "*********** BUILDING BITMATS ***********" << endl;
		BitMat bmorig_spo;
		init_bitmat(&bmorig_spo, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
		cout << "Building SPO bitmat" << endl;
		load_data_vertically((char *)config[string("RAWDATAFILE_SPO")].c_str(), triplelist, &bmorig_spo, (char *)config[string("BITMATDUMPFILE_SPO")].c_str(), true, false, true, (char *)config[string("TMP_STORAGE")].c_str());
		bmorig_spo.freebm();
//		clear_rows(&bmorig_spo, true, true, false);

//		BitMat bmorig_ops;
//		init_bitmat(&bmorig_ops, gnum_objs, gnum_preds, gnum_subs, gnum_comm_so, OPS_BITMAT);
//		cout << "Loading vertically for OPS bitmat" << endl;
//		load_data_vertically((char *)config[string("RAWDATAFILE_OPS")].c_str(), triplelist, &bmorig_ops, (char *)config[string("BITMATDUMPFILE_OPS")].c_str(), true, false, true, (char *)config[string("TMP_STORAGE")].c_str());
//		bmorig_ops.freebm();
////		clear_rows(&bmorig_ops, true, true, false);
//
//		BitMat bmorig_pso;
//		init_bitmat(&bmorig_pso, gnum_preds, gnum_subs, gnum_objs, gnum_comm_so, PSO_BITMAT);
//		cout << "Loading vertically for PSO bitmat" << endl;
//		load_data_vertically((char *)config[string("RAWDATAFILE_PSO")].c_str(), triplelist, &bmorig_pso, (char *)config[string("BITMATDUMPFILE_PSO")].c_str(), true, false, true, NULL);
//		bmorig_pso.freebm();
////		clear_rows(&bmorig_pso, true, true, false);
//
//		BitMat bmorig_pos;
//		init_bitmat(&bmorig_pos, gnum_preds, gnum_objs, gnum_subs, gnum_comm_so, POS_BITMAT);
//		cout << "Loading vertically for POS bitmat" << endl;
//		load_data_vertically((char *)config[string("RAWDATAFILE_POS")].c_str(), triplelist, &bmorig_pos, (char *)config[string("BITMATDUMPFILE_POS")].c_str(), true, false, true, NULL);
//		bmorig_pos.freebm();
////		clear_rows(&bmorig_pos, true, true, false);


		gettimeofday(&stop_time, (struct timezone *)0);
		en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
		st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
		curr_time = en-st;

		printf("Time for building BitMats: %f\n", curr_time);
	}

	/////////////////////////////////////////////////////////////

///////////////////
// BitMat ops
///////////////////
	
	if (querydata) {
		BitMat bmorig_spo;
		init_bitmat(&bmorig_spo, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
	#if MMAPFILES
		mmap_all_files();
	#endif

		/*
		 * For Preetam: bmnum is the edge-label, which in turn decides
		 * which adjancency matrix (BitMat) to use. You can change this number,
		 * and appropriate BitMat will be picked up.
		 */
		unsigned int bmnum = 121;
		bmorig_spo.reset();

		gettimeofday(&start_time, (struct timezone *)0);

		unsigned int ret = wrapper_load_from_dump_file2(&bmorig_spo, bmnum);

		gettimeofday(&stop_time, (struct timezone *)0);
		en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
		st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
		curr_time = en-st;
		printf("Time for loading bitmat %d: %f\n", bmnum, curr_time);

		unsigned int foldarr_size = bmorig_spo.row_bytes;
		unsigned char *foldarr = (unsigned char *) malloc (foldarr_size);

		gettimeofday(&start_time, (struct timezone *)0);

		simple_fold(&bmorig_spo, ROW, foldarr, foldarr_size);

		gettimeofday(&stop_time, (struct timezone *)0);
		en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
		st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
		curr_time = en-st;
		printf("Time for folding (ROW): %f\n", curr_time);

		unsigned int foldarr2_size = bmorig_spo.column_bytes;
		unsigned char *foldarr2 = (unsigned char *) malloc (foldarr2_size);

		gettimeofday(&start_time, (struct timezone *)0);

		simple_fold(&bmorig_spo, COLUMN, foldarr2, foldarr2_size);

		gettimeofday(&stop_time, (struct timezone *)0);
		en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
		st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
		curr_time = en-st;
		printf("Time for folding (COLUMN): %f\n", curr_time);

		/*
		 * For Preetam: write some code to generate a maskarray.
		 * Right now you can do it randomly.
		 */

//		unsigned char *maskarr; //TODO: populate this.
//		unsigned int maskarr_size = 0; //TODO: populate this.
//
//		gettimeofday(&start_time, (struct timezone *)0);
//
//		simple_unfold(&bmorig_spo, maskarr, maskarr_size, ROW);
//
//		gettimeofday(&stop_time, (struct timezone *)0);
//		en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
//		st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
//		curr_time = en-st;
//		printf("Time for unfolding (ROW): %f\n", curr_time);

		/*
		 * For Preetam: You may want to reset the BitMat before doing unfold
		 * on another dimension.
		 */

//		bmorig_spo.reset();
//		wrapper_load_from_dump_file2(&bmorig_spo, bmnum);
//
//		gettimeofday(&start_time, (struct timezone *)0);
//
//		simple_unfold(&bmorig_spo, maskarr, maskarr_size, COLUMN);
//
//		gettimeofday(&stop_time, (struct timezone *)0);
//		en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
//		st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
//		curr_time = en-st;
//		printf("Time for unfolding (COLUMN): %f\n", curr_time);

	#if MMAPFILES
		munmap_all_files();
	#endif
	}

	return 0;

}
///////////////////////////////////////
void parse_config_file(char *fname)
{
	FILE *fp;
	char str[300];
	char key[50];
	char value[249];
	bool rawdata = false;

	fp = fopen(fname, "r");

	if (fp == NULL) {
		fprintf(stderr, "Could not open config file %s\n", fname);
		exit (-1);
	}

	while(!feof(fp)) {
		memset(str, 0, sizeof(str));
		memset(key, 0, sizeof(key));
		memset(value, 0, sizeof(value));
		if (fgets(str, 300, fp) == NULL)
			break;
		
		//Commented line
		if (index(str, '#') == str) {
			continue;
		}
		
		char *delim = index(str, '=');
		strncpy(key, str, delim - str);
		char *newline = index(str, '\n');
		*newline = '\0';
		strcpy(value, delim + 1);
		config[string(key)]= string(value);
	}

	gnum_subs = strtoul(config[string("NUM_SUBJECTS")].c_str(), NULL, 10);
	gnum_preds = strtoul(config[string("NUM_PREDICATES")].c_str(), NULL, 10);
	gnum_objs = strtoul(config[string("NUM_OBJECTS")].c_str(), NULL, 10);
	gnum_comm_so = strtoul(config[string("NUM_COMMON_SO")].c_str(), NULL, 10);
	table_col_bytes = strtoul(config[string("TABLE_COL_BYTES")].c_str(), NULL, 10);
	row_size_bytes = strtoul(config[string("ROW_SIZE_BYTES")].c_str(), NULL, 10);
	bm_row_size = strtoul(config[string("BM_ROW_SIZE")].c_str(), NULL, 10);
	gap_size_bytes = strtoul(config[string("GAP_SIZE_BYTES")].c_str(), NULL, 10);

	comp_folded_arr = strtoul(config[string("COMPRESS_FOLDED_ARRAY")].c_str(), NULL, 10);

	gsubject_bytes = (gnum_subs%8>0 ? gnum_subs/8+1 : gnum_subs/8);
	gpredicate_bytes = (gnum_preds%8>0 ? gnum_preds/8+1 : gnum_preds/8);
	gobject_bytes = (gnum_objs%8>0 ? gnum_objs/8+1 : gnum_objs/8);
	gcommon_so_bytes = (gnum_comm_so%8>0 ? gnum_comm_so/8+1 : gnum_comm_so/8);

//	grow_size = GAP_SIZE_BYTES * gnum_objs;
//	grow = (unsigned char *) malloc (grow_size);

}
