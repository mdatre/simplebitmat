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

#include "bitmat.hpp"

unsigned char *grow;
unsigned int grow_size;
unsigned char *comp_rowfold;
unsigned char *comp_colfold;
unsigned int comp_rowfold_size;
//map<struct node*, struct triple> q_to_gr_node_map;  // Maps query node to the graph triple.
//map<struct node*, struct triple> copy_of_resmap;
unsigned int gmaxrowsize;
unsigned long gtotal_size;
unsigned int gtotal_triples;
unsigned int prev_bitset, prev_rowbit;
unsigned int table_col_bytes;
unsigned char *distresult;
unsigned long distresult_size;
ResultSet resset;
map<string, char *> mmap_table;
vector< pair<int, off_t> > vectfd(8);
//char buf_tmp[0xf000000];
//char buf_dump[0xf000000];
//char buf_table[0xf000000];

char *buf_tmp;
char *buf_dump;
char *buf_table;

//For file buffering
//char buffer[0x8000000];

//BitMat bmorig;
//set of unvisited nodes in the query graph,
//it is useful when you are at a leaf node
//and want to jump to a node which is not a neighbor
//of your parent node in case of following type of Q graph traversal in the
//order of nodes 1, 2, 3, 4, 7,... then you want to jump
//from 7 to 6 or 5.
//
//   1   2   3   4
//   o---o---o---o
//  /   /   /
// o   o   o
// 5   6   7 

//deque<struct node *> unvisited;

/**
 * Comparator (less than) for the map.
 */
//struct lt_triple {
//	bool operator()(const TP& a, const TP& b) const {
//		if((a.sub < b.sub) || (a.sub == b.sub && a.pred < b.pred) || (a.sub == b.sub && a.pred == b.pred && a.obj < b.obj))
//			return true;
//
//		return false;
//	}
//};

void match_query_graph_new(FILE **outfile, FILE **intfile, FILE **intfile_preds, struct node **bfsarr, int sidx, int eidx, bool bestm,
		map<struct node *, struct triple> &q_to_gr_node_map, map<struct node *, struct triple> &copy_of_resmap);

template <class ForwardIterator, class T>
ForwardIterator mybinary_search(ForwardIterator first, ForwardIterator last, const T& value)
{
	first = lower_bound(first,last,value);
	if (first != last && !(value < *first))
		return (first);
	else
		return (last);
}

//////////////////////////////////////////////////////////
void init_bitmat(BitMat *bitmat, unsigned int snum, unsigned int pnum, unsigned int onum, unsigned int commsonum, int dimension)
{
//	bitmat->bm = (unsigned char **) malloc (snum * sizeof (unsigned char *));
//	memset(bitmat->bm, 0, snum * sizeof (unsigned char *));
	bitmat->bm.clear();
	bitmat->num_rows = snum;
	bitmat->num_totalBMs = pnum;
	bitmat->num_columns = onum;
	bitmat->num_comm_so = commsonum;

//	row_size_bytes = bitmat->row_size_bytes;
//	gap_size_bytes = bitmat->gap_size_bytes;
	bitmat->dimension = dimension;

	bitmat->row_bytes = (snum%8>0 ? snum/8+1 : snum/8);
	bitmat->totalBMs_bytes = (pnum%8>0 ? pnum/8+1 : pnum/8);
	bitmat->column_bytes = (onum%8>0 ? onum/8+1 : onum/8);
	bitmat->common_so_bytes = (commsonum%8>0 ? commsonum/8+1 : commsonum/8);

	bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof(unsigned char));
	memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof(unsigned char));
	bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof(unsigned char));
	memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof(unsigned char));
//	bitmat->subfold_prev = NULL;
//	bitmat->objfold_prev = NULL;
//	bitmat->single_row = false;
	bitmat->num_triples = 0;

}

void shallow_init_bitmat(BitMat *bitmat, unsigned int snum, unsigned int pnum,
						unsigned int onum, unsigned int commsonum, int dimension)
{
	bitmat->bm.clear();
	bitmat->num_rows = snum;
	bitmat->num_totalBMs = pnum;
	bitmat->num_columns = onum;
	bitmat->num_comm_so = commsonum;

//	row_size_bytes = bitmat->row_size_bytes;
//	gap_size_bytes = bitmat->gap_size_bytes;
	bitmat->dimension = dimension;

	bitmat->row_bytes = (snum%8>0 ? snum/8+1 : snum/8);
	bitmat->totalBMs_bytes = (pnum%8>0 ? pnum/8+1 : pnum/8);
	bitmat->column_bytes = (onum%8>0 ? onum/8+1 : onum/8);
	bitmat->common_so_bytes = (commsonum%8>0 ? commsonum/8+1 : commsonum/8);

	bitmat->rowfold = NULL;
	bitmat->colfold = NULL;
//	bitmat->subfold_prev = NULL;
//	bitmat->objfold_prev = NULL;
//	bitmat->single_row = false;
	bitmat->num_triples = 0;
}

/////////////////////////////////////////////////
void init_bitmat_rows(BitMat *bitmat, bool initbm, bool initfoldarr)
{
	if (initfoldarr) {
		bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof(unsigned char));
		memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof(unsigned char));
		bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof(unsigned char));
		memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof(unsigned char));
	}
}
unsigned int dgap_compress_new(unsigned char *in, unsigned int size, unsigned char *out)
{
	unsigned int i = 0;
	unsigned int count = 0;
	unsigned int idx = 0;

	bool flag = (in[0] & 0x80 ? 0x01 : 0x00);
	out[0] = (in[0] & 0x80 ? 0x01 : 0x00);
	idx += 1;

	for (i = 0; i < size; i++) {
		if (in[i] == 0x00) {
			//All 0 bits
			if (!flag) {
				count += 8;
			} else {
	//			printf("%u ", count);
				memcpy(&out[idx], &count, GAP_SIZE_BYTES);
				flag = 0;
				idx += GAP_SIZE_BYTES;
				count = 8;
			}
		} else if (in[i] == 0xff) {
			if (flag) {
				count += 8;
			} else {
				memcpy(&out[idx], &count, GAP_SIZE_BYTES);
				flag = 1;
				idx += GAP_SIZE_BYTES;
				count = 8;
			}
		} else {
			//mix of 0s and 1s byte
			for (unsigned short j = 0; j < 8; j++) {
				if (!(in[i] & (0x80 >> j))) {
					//0 bit
					if (!flag) {
						count++;
					} else {
						memcpy(&out[idx], &count, GAP_SIZE_BYTES);
						flag = 0;
						idx += GAP_SIZE_BYTES;
						count = 1;
					}
				} else {
					//1 bit
					if (flag) {
						count++;
					} else {
						memcpy(&out[idx], &count, GAP_SIZE_BYTES);
						flag = 1;
						idx += GAP_SIZE_BYTES;
						count = 1;
					}

				}
				
			}
		}
	}

	//TODO: can omit pushing last 0 bits count
//	if (flag) {
		memcpy(&out[idx], &count, GAP_SIZE_BYTES);
		idx += GAP_SIZE_BYTES;
//	}
//	if (idx >= pow(2, 8*ROW_SIZE_BYTES)) {
//		printf("**** dgap_compress: ERROR size is greater than 2^%u %u\n", 8*ROW_SIZE_BYTES, idx);
//		fflush(stdout);
//		exit(1);
//	}
	return idx;
}

////////////////////////////////////////////////////////////

unsigned long count_bits_in_row(unsigned char *in, unsigned int size)
{
	if (size == 0)
		return 0;

#if USE_MORE_BYTES
	unsigned long count;
#else	
	unsigned int count;
#endif	
	unsigned int i;
//	unsigned char tmp[size];

//	memcpy(tmp, in, size);

	count = 0;
/*	
	for (i = 0; i < size*8; i++) {
		if (tmp[i/8] & 0x01u) {
			count++;
		}
		tmp[i/8] >>= 1;
	}
*/
	for (i = 0; i < size; i++) {
		if (in[i] == 0xff) {
			count += 8;
		} else if (in[i] > 0x00) {
			for (unsigned short j = 0; j < 8; j++) {
				if (in[i] & (0x80 >> j)) {
					count++;
				}
			}
		}
	}

	return count;
}
unsigned int print_set_bits_in_row(unsigned char *in, unsigned int size)
{
	unsigned int count;
	unsigned int i;

	if (in == NULL || size == 0)
		return 0;

	count = 0;
/*	
	for (i = 0; i < size*8; i++) {
		if (tmp[i/8] & 0x01u) {
			count++;
		}
		tmp[i/8] >>= 1;
	}
*/
	for (i = 0; i < size; i++) {
		if (in[i] == 0xff) {
			count += 8;
			for (unsigned short j=1; j<= 8; j++) {
				cout << i*8+j << endl;
			}
		} else if (in[i] > 0x00) {
			for (unsigned short j = 0; j < 8; j++) {
				if (in[i] & (0x80 >> j)) {
					cout << i*8+(j+1) << endl;
					count++;
				}
			}
		}
	}

	return count;
}


////////////////////////////////////////////////////////////
/*
 * The caller function is responsible for appropriate
 * initialization of the out array. It's not memset everytime as
 * in the fold function this property is exploited
 */
void dgap_uncompress(unsigned char *in, unsigned int insize, unsigned char *out, unsigned int outsize)
{
#if USE_MORE_BYTES
	unsigned long tmpcnt = 0, bitcnt = 0;
#else
	unsigned int tmpcnt = 0, bitcnt = 0;
#endif
	unsigned int cnt = 0, total_cnt = 0, bitpos = 0;
	bool flag;
//	unsigned char b;

	total_cnt = (insize-1)/GAP_SIZE_BYTES;
	flag = in[0] & 0x01;
	unsigned char format = in[0] & 0x02;

	while (cnt < total_cnt) {
		memcpy(&tmpcnt, &in[cnt*GAP_SIZE_BYTES+1], GAP_SIZE_BYTES);
//		printf("Inside while loop %u\n", cnt);
		if (format == 0x02) {
			out[(tmpcnt-1)/8] |= (0x80 >> ((tmpcnt-1) % 8));
		} else {
			if (flag) {
				for (bitpos = bitcnt; bitpos < bitcnt+tmpcnt; bitpos++) {
					out[bitpos/8] |= (0x80 >> (bitpos % 8));
				}
			}
		}
		cnt++;
		flag = !flag;
		bitcnt += tmpcnt;
	}
	if (format == 0x02) {
		assert((tmpcnt-1)/8 < outsize);
	} else {
		assert((bitcnt-1)/8 < outsize);
	}
}

////////////////////////////////////////////////////////////

void cumulative_dgap(unsigned char *in, unsigned int size, unsigned char *out)
{
	unsigned int cnt = 0;
	unsigned int total_cnt = (size-1)/GAP_SIZE_BYTES;
#if USE_MORE_BYTES
	unsigned long bitcnt = 0;
	unsigned long tmpcnt = 0;
#else
	unsigned int bitcnt = 0;
	unsigned int tmpcnt = 0;
#endif	

	out[0] = in[0];

	while (cnt < total_cnt) {
		memcpy(&tmpcnt, &in[(cnt)*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
		bitcnt += tmpcnt;

		memcpy(&out[(cnt)*GAP_SIZE_BYTES + 1], &bitcnt, GAP_SIZE_BYTES);

		cnt++;
	}

}
////////////////////////////////////////////////////////////

void de_cumulative_dgap(unsigned char *in, unsigned int size, unsigned char *out)
{
	unsigned int cnt = 1;
	unsigned int total_cnt = (size-1)/GAP_SIZE_BYTES;
#if USE_MORE_BYTES
	unsigned long bitcnt = 0, tmpcnt = 0;
#else	
	unsigned int bitcnt = 0, tmpcnt = 0;
#endif	

	out[0] = in[0];

	memcpy(&tmpcnt, &in[TOTAL_ROWSIZE(cnt)], GAP_SIZE_BYTES);
	cnt = 2;
	while (cnt <= total_cnt) {
		memcpy(&bitcnt, &in[TOTAL_ROWSIZE(cnt)], GAP_SIZE_BYTES);
		tmpcnt = bitcnt - tmpcnt;

		memcpy(&out[TOTAL_ROWSIZE(cnt)], &tmpcnt, GAP_SIZE_BYTES);
		tmpcnt = bitcnt;

		cnt++;
	}

}

////////////////////////////////////////////////////////////
//void bitmat_cumulative_dgap(unsigned char **bitmat)
//{
//	unsigned char *rowptr;
//	unsigned int rowsize = 0;
//	unsigned int i;
//
//	for (i = 0; i < gnum_subs; i++) {
//		if (bitmat[i] == NULL)
//			continue;
//
//		rowptr = bitmat[i] + ROW_SIZE_BYTES;
//		memcpy(&rowsize, bitmat[i], ROW_SIZE_BYTES);
//
//		cumulative_dgap(rowptr, rowsize, rowptr);
//	}
//}

////////////////////////////////////////////////////////////
//void bitmat_de_cumulative_dgap(unsigned char **bitmat)
//{
//	unsigned char *rowptr;
//	unsigned int rowsize = 0;
//	unsigned int i;
//
//	for (i = 0; i < gnum_subs; i++) {
//		if (bitmat[i] == NULL)
//			continue;
//
//		rowptr = bitmat[i] + ROW_SIZE_BYTES;
//		memcpy(&rowsize, bitmat[i], ROW_SIZE_BYTES);
//
//		de_cumulative_dgap(rowptr, rowsize, rowptr);
//	}
//}

////////////////////////////////////////////////////////////
unsigned int dgap_AND_new(unsigned char *in1, unsigned int size1,
					unsigned char *in2, unsigned int size2, 
					unsigned char *res)
{
	//get the initial flags whether the seq starts
	//with a 0 or 1
	bool start1 = in1[0] & 0x01;
	bool start2 = in2[0] & 0x01;
	unsigned int digit_cnt1 = (size1 - 1)/GAP_SIZE_BYTES;
	unsigned int digit_cnt2 = (size2 - 1)/GAP_SIZE_BYTES;

	if (start1 && start2)
		res[0] = 0x01;
	else 
		res[0] = 0x00;

	unsigned int cd1_pos = 0, cd1_pos_prev;
	unsigned int cd2_pos = 0, cd2_pos_prev;
	unsigned int r_pos = 0;

#if USE_MORE_BYTES
	unsigned long digit1 = 0, digit2 = 0;
#else
	unsigned int digit1 = 0, digit2 = 0;
	unsigned int tdigit1 = 0, tdigit2 = 0;
#endif

	bool flag1, flag2, flag_r, start_r;
	bool start = true;
	
	cd1_pos_prev = cd1_pos;
	cd2_pos_prev = cd2_pos;
	memcpy(&digit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
	memcpy(&digit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
	flag1 = start1;
	flag2 = start2;
	tdigit1 = digit1;
	tdigit2 = digit2;
	start_r = res[0];

	while (cd1_pos < digit_cnt1 && cd2_pos < digit_cnt2) {
		// no need to memcpy if cd1/2_pos is same
		if (cd1_pos_prev != cd1_pos) {
			memcpy(&tdigit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);	
			digit1 += tdigit1;
			cd1_pos_prev = cd1_pos;
		}
		if (cd2_pos_prev != cd2_pos) {
			memcpy(&tdigit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);	
			digit2 += tdigit2;
			cd2_pos_prev = cd2_pos;
		}

		flag_r = ((r_pos % 2) == 0) ? start_r : !start_r;

		// case '00'
		if (!flag1 && !flag2) {
//			cout << "case: 00 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (digit1 >= digit2) {

				if (start || !flag_r) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				}
				cd1_pos++;
				flag1 = !flag1;

				if (digit1 == digit2) {
					cd2_pos++;
					flag2 = !flag2;					
				} else {
					cd2_pos++;
					flag2 = !flag2;					
					while (cd2_pos < digit_cnt2) {
						memcpy(&tdigit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
						digit2 += tdigit2;
						if (digit2 > digit1) {
							cd2_pos_prev = cd2_pos;
							break;
						}
						cd2_pos++;
						flag2 = !flag2;					
					}
				}

			} else {
				if (start || !flag_r) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				}
				cd2_pos++;
				flag2 = !flag2;

				cd1_pos++;
				flag1 = !flag1;
				while (cd1_pos < digit_cnt1) {
					memcpy(&tdigit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
					digit1 += tdigit1;
					if (digit1 > digit2) {
						cd1_pos_prev = cd1_pos;
						break;
					}
					cd1_pos++;
					flag1 = !flag1;
				}

			}
		}
		// case '11'
		else if (flag1 && flag2) {
//			cout << "case: 11 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (digit1 <= digit2) {
				if (start) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				}
				cd1_pos++;
				flag1 = !flag1;

				if (digit1 == digit2) {
					cd2_pos++;
					flag2 = !flag2;
				}

			} else {
				if (start) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				}

				cd2_pos++;
				flag2 = !flag2;

			}
		}
		// case '01'
		else if (!flag1 && flag2) {
//			cout << "case: 01 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (start || !flag_r) {
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
			} else {
				r_pos++;
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
			}

			cd1_pos++;
			flag1 = !flag1;

			if (digit1 == digit2) {
				cd2_pos++;
				flag2 = !flag2;
			} else if (digit2 < digit1) {
				cd2_pos++;
				flag2 = !flag2;
				while (cd2_pos < digit_cnt2) {
					memcpy(&tdigit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
					digit2 += tdigit2;
					if (digit2 > digit1) {
						cd2_pos_prev = cd2_pos;
						break;
					}
					cd2_pos++;
					flag2 = !flag2;
				}
			}

		}
		// case '10'
		else if (flag1 && !flag2) {
//			cout << "case: 10 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (start || !flag_r) {
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
			} else {
				r_pos++;
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
			}

			cd2_pos++;
			flag2 = !flag2;

			if (digit1 == digit2) {
				cd1_pos++;
				flag1 = !flag1;
			} else if (digit1 < digit2) {
				cd1_pos++;
				flag1 = !flag1;
				while (cd1_pos < digit_cnt1) {
					memcpy(&tdigit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
					digit1 += tdigit1;
					if (digit1 > digit2) {
						cd1_pos_prev = cd1_pos;
						break;
					}
					cd1_pos++;
					flag1 = !flag1;
				}

			}

		}

		if (start)
			start = false;

	}
	unsigned int ressize = (r_pos+1) * GAP_SIZE_BYTES + 1;
	de_cumulative_dgap(res, ressize, res);
//	return ((r_pos+1) * GAP_SIZE_BYTES + 1);	
	return (ressize);	

}

unsigned int dgap_AND(unsigned char *in1, unsigned int size1,
					unsigned char *in2, unsigned int size2, 
					unsigned char *res)
{
	//get the initial flags whether the seq starts
	//with a 0 or 1
	bool start1 = in1[0] & 0x01;
	bool start2 = in2[0] & 0x01;
	unsigned int digit_cnt1 = (size1 - 1)/GAP_SIZE_BYTES;
	unsigned int digit_cnt2 = (size2 - 1)/GAP_SIZE_BYTES;

	if (start1 && start2)
		res[0] = 0x01;
	else 
		res[0] = 0x00;

	unsigned int cd1_pos = 0, cd1_pos_prev;
	unsigned int cd2_pos = 0, cd2_pos_prev;
	unsigned int r_pos = 0;

#if USE_MORE_BYTES
	unsigned long digit1 = 0, digit2 = 0;
#else
	unsigned int digit1 = 0, digit2 = 0;
#endif

	bool flag1, flag2, flag_r, start_r;
	bool start = true;
	
	cd1_pos_prev = cd1_pos;
	cd2_pos_prev = cd2_pos;
	memcpy(&digit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
	memcpy(&digit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
	flag1 = start1;
	flag2 = start2;
	start_r = res[0];

	while (cd1_pos < digit_cnt1 && cd2_pos < digit_cnt2) {
		// no need to memcpy if cd1/2_pos is same
		if (cd1_pos_prev != cd1_pos) {
			memcpy(&digit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);	
			cd1_pos_prev = cd1_pos;
		}
		if (cd2_pos_prev != cd2_pos) {
			memcpy(&digit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);	
			cd2_pos_prev = cd2_pos;
		}

		flag_r = ((r_pos % 2) == 0) ? start_r : !start_r;

		// case '00'
		if (!flag1 && !flag2) {
//			cout << "case: 00 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (digit1 >= digit2) {

				if (start || !flag_r) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				}
				cd1_pos++;
				flag1 = !flag1;

				if (digit1 == digit2) {
					cd2_pos++;
					flag2 = !flag2;					
				} else {
					cd2_pos++;
					flag2 = !flag2;					
					while (cd2_pos < digit_cnt2) {
						memcpy(&digit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
						if (digit2 > digit1) {
							cd2_pos_prev = cd2_pos;
							break;
						}
						cd2_pos++;
						flag2 = !flag2;					
					}
				}

			} else {
				if (start || !flag_r) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				}
				cd2_pos++;
				flag2 = !flag2;

				cd1_pos++;
				flag1 = !flag1;
				while (cd1_pos < digit_cnt1) {
					memcpy(&digit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
					if (digit1 > digit2) {
						cd1_pos_prev = cd1_pos;
						break;
					}
					cd1_pos++;
					flag1 = !flag1;
				}

			}
		}
		// case '11'
		else if (flag1 && flag2) {
//			cout << "case: 11 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (digit1 <= digit2) {
				if (start) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
				}
				cd1_pos++;
				flag1 = !flag1;

				if (digit1 == digit2) {
					cd2_pos++;
					flag2 = !flag2;
				}

			} else {
				if (start) {
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				} else {
					r_pos++;
					memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
				}

				cd2_pos++;
				flag2 = !flag2;

			}
		}
		// case '01'
		else if (!flag1 && flag2) {
//			cout << "case: 01 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (start || !flag_r) {
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
			} else {
				r_pos++;
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit1, GAP_SIZE_BYTES);
			}

			cd1_pos++;
			flag1 = !flag1;

			if (digit1 == digit2) {
				cd2_pos++;
				flag2 = !flag2;
			} else if (digit2 < digit1) {
				cd2_pos++;
				flag2 = !flag2;
				while (cd2_pos < digit_cnt2) {
					memcpy(&digit2, &in2[cd2_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
					if (digit2 > digit1) {
						cd2_pos_prev = cd2_pos;
						break;
					}
					cd2_pos++;
					flag2 = !flag2;
				}
			}

		}
		// case '10'
		else if (flag1 && !flag2) {
//			cout << "case: 10 digit1 " << digit1 << " digit2 " << digit2 << endl;
			if (start || !flag_r) {
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
			} else {
				r_pos++;
				memcpy(&res[r_pos*GAP_SIZE_BYTES+1], &digit2, GAP_SIZE_BYTES);
			}

			cd2_pos++;
			flag2 = !flag2;

			if (digit1 == digit2) {
				cd1_pos++;
				flag1 = !flag1;
			} else if (digit1 < digit2) {
				cd1_pos++;
				flag1 = !flag1;
				while (cd1_pos < digit_cnt1) {
					memcpy(&digit1, &in1[cd1_pos*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
					if (digit1 > digit2) {
						cd1_pos_prev = cd1_pos;
						break;
					}
					cd1_pos++;
					flag1 = !flag1;
				}

			}

		}

		if (start)
			start = false;

	}
	return ((r_pos+1) * GAP_SIZE_BYTES + 1);	

}

////////////////////////////////////////////////////////////
/*
 * It always concatenates to in1, so make sure to have in1
 * large enough.. and returns new in1 size
 */  
unsigned int concatenate_comp_arr(unsigned char *in1, unsigned int size1,
								unsigned char *in2, unsigned int size2)
{
	bool flag1 = in1[0] & 0x01;
	bool flag2 = in2[0] & 0x01;

	unsigned int total_cnt1 = (size1-1)/GAP_SIZE_BYTES;
	unsigned int total_cnt2 = (size2-1)/GAP_SIZE_BYTES;

	bool last_digit_0 = ((!flag1) && (total_cnt1%2)) || ((flag1) && !(total_cnt1%2));
	bool last_digit_1 = (flag1 && (total_cnt1%2)) || ((!flag1) && !(total_cnt1%2));

	bool first_digit_0 = !flag2;
	bool first_digit_1 = flag2;
								
#if USE_MORE_BYTES
	unsigned long tmpcnt1 = 0;
	unsigned long tmpcnt2 = 0;
#else
	unsigned int tmpcnt1 = 0;
	unsigned int tmpcnt2 = 0;
#endif

	//concatenate in case of these conditions
	if ((last_digit_0 && first_digit_0) || (last_digit_1 && first_digit_1)) {
		memcpy(&tmpcnt1, &in1[TOTAL_ROWSIZE(total_cnt1)], GAP_SIZE_BYTES);
		memcpy(&tmpcnt2, &in2[1], GAP_SIZE_BYTES);

		tmpcnt1 += tmpcnt2;
		memcpy(&in1[TOTAL_ROWSIZE(total_cnt1)], &tmpcnt1, GAP_SIZE_BYTES);
		if (size2 - (GAP_SIZE_BYTES+1) > 0)
			memcpy(&in1[size1], &in2[GAP_SIZE_BYTES+1], size2-(GAP_SIZE_BYTES+1));

		return (size1 + size2 - (GAP_SIZE_BYTES+1));
		
	} else {
		memcpy(&in1[size1], &in2[1], size2-1);
		return (size1 + size2 - 1);
	}

}


////////////////////////////////////////////////////////////
/*
 * cflag is completion flag  
 */  

//TODO: make provision for loading OPS bitmat
/*
void map_to_row_wo_dgap(unsigned char **in, unsigned int pbit, unsigned int obit,
				unsigned int spos, bool cflag, bool start)
{
#ifdef USE_MORE_BYTES	
	unsigned long bitset = 0, tmpcnt = 0;
#else	
	unsigned int bitset = 0, tmpcnt = 0;
#endif	
	unsigned int size = 0;
	unsigned char *rowptr;

//	printf("sbit %u pbit %u obit %u\n", spos, pbit, obit);
	bitset = (pbit - 1) * (gcolumn_bytes << 3) + obit;

//	printf("prev_bitset %u bitset %u\n", prev_bitset, bitset);

	//Complete earlier row and start a new one
	if (cflag) {
		tmpcnt = gnum_preds * (gcolumn_bytes << 3) - prev_bitset;
		memcpy(&size, grow, ROW_SIZE_BYTES);

		//for the last 0s
		if (tmpcnt > 0) {
//			printf("last 0s appended\n");
			unsigned char tmp[GAP_SIZE_BYTES+1];
			tmp[0] = 0x00;
			memcpy(&tmp[1], &tmpcnt, GAP_SIZE_BYTES);
			size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
			if (size >= pow(2, 8*ROW_SIZE_BYTES)) {
				printf("**** map_to_row_wo_dgap: ERROR2 size greater than 2^%u %u -- s:%u p:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, pbit, obit);
				exit (-1);
			}

		}
		
		//TODO: you need different *in
		if (*in != NULL) {
			printf("*** ERROR: something wrong *in not null\n");
			fflush(stdout);
			exit(-1);
		}
		*in = (unsigned char *) malloc (size + ROW_SIZE_BYTES);
		memcpy(*in, &size, ROW_SIZE_BYTES);
		memcpy(*in+ROW_SIZE_BYTES, &grow[ROW_SIZE_BYTES], size);


		//////////////////////////
		//DEBUG
		//////////////////////////
//		cnt = 0;
//		total_cnt = (size-1)/4;
//
//		rowptr = *in + 2;
//		bitcnt = 0;
//		bool flag = rowptr[0];
//		printf("[%u] ", rowptr[0]);
//		for (cnt = 0; cnt < total_cnt; cnt++, flag = !flag) {
//			memcpy(&tmpcnt, &rowptr[cnt*4+1], 4);
////			printf("%u ", tmpcnt);
//			if (flag)
//				bitcnt += tmpcnt;
//		}
//		printf("\n#triples %u\n", bitcnt);

		///////////////////////////

		if (size > gmaxrowsize)
			gmaxrowsize = size;
		
//		printf("rowsize %u\n", size);
		gtotal_size += size + ROW_SIZE_BYTES;
	}
	if (start) {
//		*in = (unsigned char *) malloc (1024);
		prev_bitset = bitset;
		if (bitset - 1 > 0) {
			grow[ROW_SIZE_BYTES] = 0x00;
			tmpcnt = bitset - 1;
			memcpy(&grow[ROW_SIZE_BYTES+1], &tmpcnt, GAP_SIZE_BYTES);
			tmpcnt = 1;
			memcpy(&grow[ROW_SIZE_BYTES + GAP_SIZE_BYTES + 1], &tmpcnt, GAP_SIZE_BYTES);
			size = 2*GAP_SIZE_BYTES + 1;
			memcpy(grow, &size, ROW_SIZE_BYTES);
		} else {
			grow[ROW_SIZE_BYTES] = 0x01;
			tmpcnt = 1;
			memcpy(&grow[ROW_SIZE_BYTES+1], &tmpcnt, GAP_SIZE_BYTES);
			size = GAP_SIZE_BYTES+1;
			memcpy(grow, &size, ROW_SIZE_BYTES);
		}
	} else {
		unsigned char tmp[GAP_SIZE_BYTES+1];
		memcpy(&size, grow, ROW_SIZE_BYTES);
		//append 0s in between 2 set bits
		if (bitset - prev_bitset > 1) {
//			printf("0s in between 2 bitsets\n");
			tmp[0] = 0x00;
			tmpcnt = bitset - prev_bitset - 1;
			memcpy(&tmp[1], &tmpcnt, GAP_SIZE_BYTES);
//			printf("before concate1\n");
			size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
			if (size >= pow(2, 8*ROW_SIZE_BYTES)) {
				printf("**** map_to_row_wo_dgap: ERROR size greater than 2^%u %u -- s:%u p:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, pbit, obit);
				exit (-1);
			}
		} else if (bitset - prev_bitset == 0) {
			//no op
			return;
		}
		//now append the set bit
		tmp[0] = 0x01;
		tmpcnt = 1;
		memcpy(&tmp[1], &tmpcnt, GAP_SIZE_BYTES);
//		memcpy(&size, grow, 2);
//		printf("before concate2\n");
		size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
		if (size >= pow(2, 8*ROW_SIZE_BYTES)) {
			printf("**** map_to_row_wo_dgap: ERROR2 size greater than 2^%u %u -- s:%u p:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, pbit, obit);
			exit (-1);
		}
//		printf("after concate2\n");
		memcpy(grow, &size, ROW_SIZE_BYTES);
		prev_bitset = bitset;
	}
	
//	printf("Exiting map_to_row_wo_dgap\n");
}
*/
////////////////////////////////////////////////////////////
unsigned long count_dgaps_in_row(unsigned char *in, unsigned int size)
{
	if (in == NULL)
		return 0;

	unsigned char format = in[0] & 0x02;
	if (format == 0x02)
		return 0;

	assert((size-1) % GAP_SIZE_BYTES == 0);

	return ((size-1)/GAP_SIZE_BYTES);
}
////////////////////////////////////////////////////////////

FILE* map_to_row_wo_dgap_vertical(BitMat *bitmat, unsigned int spos, unsigned int obit,
				unsigned int sprev, bool cflag, bool start, bool ondisk, bool listload, FILE *tmpdump,
				bool changebm)
{
//	assert(BM_ROW_SIZE <= sizeof (unsigned short));
#if USE_MORE_BYTES
	unsigned long bitset = 0, later_0 = 0, ini_0 = 0, mid_0 = 0;
#else	
	unsigned int bitset = 0, later_0 = 0, ini_0 = 0, mid_0 = 0;
#endif	
	unsigned int size = 0;
	unsigned char *rowptr;

	assert(GAP_SIZE_BYTES == sizeof(ini_0));

	bitset = obit;

//	if (ondisk)
//		assert(tmpdump != NULL);

	if (!changebm && bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT && ondisk) {
		bitmat->colfold[(obit-1)/8] |= (0x80 >> ((obit-1) % 8));
	}

	//Complete earlier row and start a new one
	if (cflag) {
//		later_0 = (bitmat->column_bytes << 3) - prev_bitset;
		size = 0;
		memcpy(&size, grow, ROW_SIZE_BYTES);

		//for the last 0s
		//COMMENT: Optimization: you don't need to store last 0s
//		if (later_0 > 0) {
//			unsigned char tmp[GAP_SIZE_BYTES+1];
//			tmp[0] = 0x00;
//			memcpy(&tmp[1], &later_0, GAP_SIZE_BYTES);
//			size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
//			if (size >= pow((double)2, (double)8*ROW_SIZE_BYTES)) {
//				printf("**** map_to_row_wo_dgap_vertical: ERROR2 size greater than 2^%u %u -- s:%u p:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, pbit, obit);
//				exit (-1);
//			}
//		}
		
		unsigned int total_setbits = 0, total_gaps = 0;
		total_setbits = count_triples_in_row(&grow[ROW_SIZE_BYTES], size);
		total_gaps = count_dgaps_in_row(&grow[ROW_SIZE_BYTES], size);

		if (total_gaps >= 0xffff || total_setbits >= 0xffff) {
			cout << "total_gaps " << total_gaps << " total_setbits " << total_setbits << endl;
			assert(0);
		}
//		assert(total_gaps < 0xffff && total_setbits < 0xffff);

		unsigned char *data = NULL;

		if (total_gaps > total_setbits) {
			//Just store the positions of setbits.
//			vector<unsigned int> setbits;
			data = grow + ROW_SIZE_BYTES;
			unsigned char *tmpdata = (unsigned char *) malloc (1 + total_setbits*GAP_SIZE_BYTES +
																BM_ROW_SIZE);
			bool flag = data[0] & 0x01; unsigned char format = data[0] & 0x02;
			if (format == 0x02) assert(0);

			unsigned int cnt = 0, total_bits = 0, gapcnt = 0;
			tmpdata[BM_ROW_SIZE] = 0x02; //new format
			size = 1;
			while (cnt < total_gaps) {
				memcpy(&gapcnt, &data[cnt*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
				if (flag) {
					for (unsigned int i = total_bits+1; i <= total_bits+gapcnt; i++) {
						memcpy(&tmpdata[size + BM_ROW_SIZE], &i, GAP_SIZE_BYTES);
						size += GAP_SIZE_BYTES;
//						setbits.push_back(i);
					}
				}
				flag = !flag;
				cnt++;
				total_bits += gapcnt;
			}
			assert(size == (1 + total_setbits * GAP_SIZE_BYTES));
			unsigned int sz = total_setbits + 1; //((size-1)/GAP_SIZE_BYTES) + 1;
			memcpy(tmpdata, &sz, BM_ROW_SIZE);
			data = tmpdata;
		} else {
			data = (unsigned char *) malloc (size + BM_ROW_SIZE);
			unsigned int sz = total_gaps + 1;// ((size-1)/GAP_SIZE_BYTES) + 1;
			memcpy(data, &sz, BM_ROW_SIZE);
			memcpy(&data[BM_ROW_SIZE], &grow[ROW_SIZE_BYTES], size);
		}
//		memcpy(data, &size, ROW_SIZE_BYTES);
//		memcpy(&data[ROW_SIZE_BYTES], &grow[ROW_SIZE_BYTES], size);

		if (tmpdump) {
			unsigned int rowsize = 0;
			memcpy(&rowsize, data, BM_ROW_SIZE);
			fwrite(data, TOTAL_ROWSIZE(rowsize) + BM_ROW_SIZE, 1, tmpdump);
			free(data);
		} else {
			struct row r = {sprev, data};
			if (listload) {
				bitmat->bm.push_back(r);
			} else {
				bitmat->vbm.push_back(r);
			}
		}

		if (ondisk) {
			if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
				unsigned int t = sprev - prev_rowbit;
	//			cout << "sprev " << sprev << " prev_rowbit " << prev_rowbit << " comp_rowfold_size " << comp_rowfold_size << endl;
				if (comp_rowfold_size == 0) {
					//Add some zeros before
					if (t > 0) {
						comp_rowfold[ROW_SIZE_BYTES] = 0x00;
						memcpy(&comp_rowfold[1+ROW_SIZE_BYTES], &t, GAP_SIZE_BYTES);
						comp_rowfold_size += 1 + GAP_SIZE_BYTES;
					} else {
						comp_rowfold[ROW_SIZE_BYTES] = 0x01;
						comp_rowfold_size += 1;
					}

					t = 1;
					memcpy(&comp_rowfold[comp_rowfold_size + ROW_SIZE_BYTES], &t, GAP_SIZE_BYTES);
					comp_rowfold_size += GAP_SIZE_BYTES;

				} else {
					//Concatenate middle 0s
					unsigned char tmp[1+GAP_SIZE_BYTES];
					//Not consecutive rows, hence pad 0s
					if (t > 1) {
						tmp[0] = 0x00;
						t -= 1;
						memcpy(&tmp[1], &t, GAP_SIZE_BYTES);

						comp_rowfold_size = concatenate_comp_arr(&comp_rowfold[ROW_SIZE_BYTES],
								comp_rowfold_size, tmp, GAP_SIZE_BYTES+1);
//						assert(comp_rowfold_size < pow(2, 8*ROW_SIZE_BYTES));
					}
					//Concatenate a 1 to set this row-bit
					tmp[0] = 0x01;
					t = 1;
					memcpy(&tmp[1], &t, GAP_SIZE_BYTES);
					comp_rowfold_size = concatenate_comp_arr(&comp_rowfold[ROW_SIZE_BYTES], comp_rowfold_size, tmp, GAP_SIZE_BYTES+1);
//					assert(comp_rowfold_size < pow(2, 8*ROW_SIZE_BYTES));
				}
				prev_rowbit = sprev;
			} else {
				bitmat->rowfold[(sprev-1)/8] |= (0x80 >> ((sprev-1) % 8));
			}
		}

	}
	if (!changebm) {
		if (start) {
	//		cout << "starting here" << endl;
			prev_bitset = bitset;
			if (bitset - 1 > 0) {
				grow[ROW_SIZE_BYTES] = 0x00;
				ini_0 = bitset - 1;
	//			cout << "before memcpy to grow" << endl;
				memcpy(&grow[ROW_SIZE_BYTES+1], &ini_0, GAP_SIZE_BYTES);
	//			cout << "after memcpy to grow" << endl;
				unsigned int tmpcnt = 1;
				memcpy(&grow[ROW_SIZE_BYTES + GAP_SIZE_BYTES + 1], &tmpcnt, GAP_SIZE_BYTES);
				size = 2*GAP_SIZE_BYTES + 1;
				memcpy(grow, &size, ROW_SIZE_BYTES);
			} else {
				grow[ROW_SIZE_BYTES] = 0x01;
				unsigned int tmpcnt = 1;
				memcpy(&grow[ROW_SIZE_BYTES+1], &tmpcnt, GAP_SIZE_BYTES);
				size = GAP_SIZE_BYTES+1;
				memcpy(grow, &size, ROW_SIZE_BYTES);
			}
		} else {
	//		cout << "not starting here" << endl;
			unsigned char tmp[GAP_SIZE_BYTES+1];
			size = 0;
			memcpy(&size, grow, ROW_SIZE_BYTES);
			//append 0s in between 2 set bits
			if (bitset - prev_bitset > 1) {
				tmp[0] = 0x00;
				mid_0 = bitset - prev_bitset - 1;
				memcpy(&tmp[1], &mid_0, GAP_SIZE_BYTES);
				size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
				if (size >= grow_size)
					assert(0);
				if (size >= pow((double)2, (double)8*ROW_SIZE_BYTES)) {
					printf("**** map_to_row_wo_dgap_vertical: ERROR size greater than 2^%u %u -- s:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, obit);
					exit (-1);
				}
			} else if (bitset - prev_bitset == 0) {
				//no op
				return tmpdump;
			}
			//now append the set bit
			tmp[0] = 0x01;
			unsigned int tmpcnt = 1;
			memcpy(&tmp[1], &tmpcnt, GAP_SIZE_BYTES);
			size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
			if (size >= grow_size)
				assert(0);
			if (size >= pow((double)2, (double)8*ROW_SIZE_BYTES)) {
				printf("**** map_to_row_wo_dgap_vertical: ERROR2 size greater than 2^%u %u -- s:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, obit);
				exit (-1);
			}
			memcpy(grow, &size, ROW_SIZE_BYTES);
			prev_bitset = bitset;
		}
	}

	return tmpdump;
	
//	cout << "Exiting map_to_row_wo_dgap_vertical" << endl;
}

////////////////////////////////////////////////////////////
//static unsigned int dbg_counter;

/*
void map_to_row_wo_dgap_ops(unsigned char **in, unsigned int pbit, unsigned int obit,
				unsigned int spos, bool cflag, bool start)
{
#ifdef USE_MORE_BYTES	
	unsigned long bitset = 0, tmpcnt = 0;
#else
	unsigned int bitset = 0, tmpcnt = 0;
#endif
	unsigned int size = 0;
	unsigned char *rowptr;

//	printf("sbit %u pbit %u obit %u\n", spos, pbit, obit);
	bitset = (pbit - 1) * (gsubject_bytes << 3) + obit;

	dbg_counter++;

//	printf("prev_bitset %u bitset %u\n", prev_bitset, bitset);

	//Complete earlier row and start a new one
	if (cflag) {
		tmpcnt = gnum_preds * (gsubject_bytes << 3) - prev_bitset;
		memcpy(&size, grow, ROW_SIZE_BYTES);
//		printf("Completion flag set tmpcnt %u prev_bitset %u size %u\n", tmpcnt, prev_bitset, size);

		//for the last 0s
		if (tmpcnt > 0) {
//			printf("last 0s appended\n");
			unsigned char tmp[GAP_SIZE_BYTES+1];
			tmp[0] = 0x00;
			memcpy(&tmp[1], &tmpcnt, GAP_SIZE_BYTES);
//			printf("Before concatenate_comp_arr\n");
			size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
//			printf("After concatenate_comp_arr\n");
		}
		
		//TODO: you need different *in
		if (*in != NULL) {
			printf("*** ERROR: something wrong *in not null for dbg_counter %d\n", dbg_counter);
			fflush(stdout);
			exit(-1);
		}
		*in = (unsigned char *) malloc (size + ROW_SIZE_BYTES);
		memcpy(*in, &size, ROW_SIZE_BYTES);
		memcpy(*in+ROW_SIZE_BYTES, &grow[ROW_SIZE_BYTES], size);

		if (size > gmaxrowsize)
			gmaxrowsize = size;
		
//		printf("rowsize %u\n", size);
		gtotal_size += size + ROW_SIZE_BYTES;
	}
	if (start) {
//		*in = (unsigned char *) malloc (1024);
		prev_bitset = bitset;
		if (bitset - 1 > 0) {
			grow[ROW_SIZE_BYTES] = 0x00;
			tmpcnt = bitset - 1;
			memcpy(&grow[ROW_SIZE_BYTES+1], &tmpcnt, GAP_SIZE_BYTES);
			tmpcnt = 1;
			memcpy(&grow[ROW_SIZE_BYTES + GAP_SIZE_BYTES + 1], &tmpcnt, GAP_SIZE_BYTES);
			size = 2*GAP_SIZE_BYTES + 1;
			memcpy(grow, &size, ROW_SIZE_BYTES);
		} else {
			grow[ROW_SIZE_BYTES] = 0x01;
			tmpcnt = 1;
			memcpy(&grow[ROW_SIZE_BYTES+1], &tmpcnt, GAP_SIZE_BYTES);
			size = GAP_SIZE_BYTES+1;
			memcpy(grow, &size, ROW_SIZE_BYTES);
		}
	} else {
		unsigned char tmp[GAP_SIZE_BYTES+1];
		memcpy(&size, grow, ROW_SIZE_BYTES);
		//append 0s in between 2 set bits
		if (bitset - prev_bitset > 1) {
//			printf("0s in between 2 bitsets\n");
			tmp[0] = 0x00;
			tmpcnt = bitset - prev_bitset - 1;
			memcpy(&tmp[1], &tmpcnt, GAP_SIZE_BYTES);
//			printf("before concate1\n");
			size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
			if (size >= pow(2, 8*ROW_SIZE_BYTES)) {
				printf("**** map_to_row_wo_dgap_ops: ERROR size greater than 2^%u %u -- s:%u p:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, pbit, obit);
				exit (-1);
			}
		} else if (bitset - prev_bitset == 0) {
			//no op
			return;
		}
		//now append the set bit
		tmp[0] = 0x01;
		tmpcnt = 1;
		memcpy(&tmp[1], &tmpcnt, GAP_SIZE_BYTES);
//		memcpy(&size, grow, 2);
//		printf("before concate2\n");
		size = concatenate_comp_arr(&grow[ROW_SIZE_BYTES], size, tmp, GAP_SIZE_BYTES+1);
		if (size >= pow(2, 8*ROW_SIZE_BYTES)) {
			printf("**** map_to_row_wo_dgap_ops: ERROR2 size greater than 2^%u %u -- s:%u p:%u o:%u\n", 8*ROW_SIZE_BYTES, size, spos, pbit, obit);
			exit (-1);
		}
//		printf("after concate2\n");
		memcpy(grow, &size, ROW_SIZE_BYTES);
		prev_bitset = bitset;
	}
	
}
*/

///////////////////////////////////////////////////////////

unsigned long count_triples_in_row(unsigned char *in, unsigned int size)
{
	if (in == NULL)
		return 0;

	bool flag = in[0] & 0x01; unsigned char format = in[0] & 0x02;
//	bool flag2;
	unsigned int cnt = 0;
	unsigned int total_cnt = (size-1)/GAP_SIZE_BYTES;
#if USE_MORE_BYTES
	unsigned long triplecnt = 0;
	unsigned long tmpcnt = 0;
#else
	unsigned int triplecnt = 0;
	unsigned int tmpcnt = 0;
#endif	

//	printf("[%u] ", flag);
	while (cnt < total_cnt) {
		memcpy(&tmpcnt, &in[cnt*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
		if (format == 0x02) {
			triplecnt++;
		} else if (format == 0x00) {
			if (flag) {
				triplecnt += tmpcnt;
			}
		}

		cnt++;
		flag = !flag;
	}
//	printf("\ntriplecnt %u\n", triplecnt);
	return triplecnt;
}

////////////////////////////////////////////////////////////
/*
unsigned int fixed_p_fixed_o (unsigned char *in, unsigned int size, unsigned char *out,
								unsigned int ppos, unsigned int opos)
{
//	printf("Inside fixed_p_fixed_o\n");
#ifdef USE_MORE_BYTES	
	unsigned long obit = (ppos - 1) * (gobject_bytes << 3) + opos;
	unsigned long ini_0 = obit - 1;
	unsigned long later_0 = gnum_preds * (gobject_bytes << 3) - obit;
	unsigned long tmpcnt, bitcnt = 0;
#else	
	unsigned int obit = (ppos - 1) * (gobject_bytes << 3) + opos;
	unsigned int ini_0 = obit - 1;
	unsigned int later_0 = gnum_preds * (gobject_bytes << 3) - obit;
	unsigned int tmpcnt, bitcnt = 0;
#endif	

	bool flag = in[0];
	unsigned int cnt = 0;
	unsigned int total_cnt = (size-1)/GAP_SIZE_BYTES;
	unsigned long t = 1;

	if (obit > 1) {
		out[0] = 0x00;
	} else {
		//obit is at first position
		out[0] = in[0];
		if (out[0]) {
			//obit is set to 1
			memcpy(&out[1], &t, GAP_SIZE_BYTES);
			memcpy(&out[GAP_SIZE_BYTES+1], &later_0, GAP_SIZE_BYTES);
			return (2*GAP_SIZE_BYTES + 1);
		} else {
			//obit is set to 0
//			free(out);
//			out = NULL;
//			later_0 += 1;
//			memcpy(&out[1], &later_0, 4);
			return 0;
		}
	}
	
	while (cnt < total_cnt) {
		memcpy(&tmpcnt, &in[(cnt)*GAP_SIZE_BYTES+1], GAP_SIZE_BYTES);
		bitcnt += tmpcnt;
		if (bitcnt >= obit && flag) {
			memcpy(&out[1], &ini_0, GAP_SIZE_BYTES);
			memcpy(&out[GAP_SIZE_BYTES+1], &t, GAP_SIZE_BYTES);
			if (later_0) {
				memcpy(&out[2*GAP_SIZE_BYTES + 1], &later_0, GAP_SIZE_BYTES);
				return (3*GAP_SIZE_BYTES + 1);
			} else {
				return (2*GAP_SIZE_BYTES + 1);
			}
		} else if (bitcnt >= obit) {
//			free(out);
//			out = NULL;
//			ini_0 += (later_0 + 1);
//			memcpy(&out[1], &ini_0, 4);
			return 0;
		}
		flag = !flag;
		cnt++;
	}
}
*/
////////////////////////////////////////////////////////////
// has got tested

//unsigned int fixed_p_var_o(unsigned char *in, unsigned int size, unsigned char *out,
//							unsigned int ppos)
//{
////	printf("Inside fixed_p_var_o\n");
//#ifdef USE_MORE_BYTES
//	unsigned long ini_0 = (ppos - 1) * (gobject_bytes << 3);
//	unsigned long later_0 = (gnum_preds - ppos) * (gobject_bytes << 3);
//	unsigned long ppos_start_bit = ini_0 + 1, gap;
//	unsigned long ppos_end_bit = ppos_start_bit + (gobject_bytes << 3) - 1;
//	unsigned long tmpcnt, bitcnt = 0;
//	unsigned long t = gobject_bytes *8;
//#else
//	unsigned int ini_0 = (ppos - 1) * (gobject_bytes << 3);
//	unsigned int later_0 = (gnum_preds - ppos) * (gobject_bytes << 3);
//	unsigned int ppos_start_bit = ini_0 + 1, gap;
//	unsigned int ppos_end_bit = ppos_start_bit + (gobject_bytes << 3) - 1;
//	unsigned int tmpcnt, bitcnt = 0;
//	unsigned int t = gobject_bytes *8;
//#endif
//
//	unsigned int idx = 0;
//	bool flag = in[0] & 0x01;
//	unsigned int cnt = 1;
//	unsigned int total_cnt = (size-1)/GAP_SIZE_BYTES;
//	bool inside = false;
//
//	//TODO: change everywhere
//	//ppos_start_bit = ppos_end_bit + 1;
//
//	if (ppos > 1) {
//		out[0] = 0;
//	} else {
//		out[0] = in[0];
//	}
//
//	for (cnt = 1; cnt <= total_cnt; cnt++, flag = !flag) {
//		memcpy(&tmpcnt, &in[(cnt-1)*GAP_SIZE_BYTES+1], GAP_SIZE_BYTES);
//		bitcnt += tmpcnt;
////		printf("bitcnt %u ppost_start %u ppos_end %u\n", bitcnt, ppos_start_bit, ppos_end_bit);
//		if (bitcnt >= ppos_start_bit && bitcnt <= ppos_end_bit) {
//			inside = true;
//			gap = (bitcnt - ppos_start_bit + 1) < tmpcnt ? (bitcnt - ppos_start_bit + 1) : tmpcnt;
//
//			if (gap == gobject_bytes*8 && !flag) {
////				free(out);
////				out = NULL;
//				return 0;
//			} else if (gap == gobject_bytes*8 && flag) {
//				if (ini_0 > 0) {
//					memcpy(&out[idx*GAP_SIZE_BYTES+1], &ini_0, GAP_SIZE_BYTES);
//					idx++;
//				}
//				memcpy(&out[idx*GAP_SIZE_BYTES+1], &gap, GAP_SIZE_BYTES);
//				idx++;
//				if (later_0 > 0) {
//					memcpy(&out[idx*GAP_SIZE_BYTES+1], &later_0, GAP_SIZE_BYTES);
//					idx++;
//				}
//				return (idx*GAP_SIZE_BYTES+1);
//			}
//
//			if (idx == 0) {
//				if (!flag) {
//					ini_0 += gap;
//					memcpy(&out[idx*GAP_SIZE_BYTES+1], &ini_0, GAP_SIZE_BYTES);
//				} else {
//					if (ini_0 > 0) {
//						memcpy(&out[idx*GAP_SIZE_BYTES+1], &ini_0, GAP_SIZE_BYTES);
//						idx++;
//					}
//					memcpy(&out[idx*GAP_SIZE_BYTES+1], &gap, GAP_SIZE_BYTES);
//				}
//				idx++;
//				continue;
//			} else {
//				if (bitcnt == ppos_end_bit) {
//					//this predicate ends here.. so no more processing
//					inside = false;
//					if (!flag) {
//						later_0 += gap;
//					} else {
//						memcpy(&out[idx*GAP_SIZE_BYTES+1], &gap, GAP_SIZE_BYTES);
//						idx++;
//					}
//					if (later_0 > 0) {
//						memcpy(&out[idx*GAP_SIZE_BYTES+1], &later_0, GAP_SIZE_BYTES);
//						idx++;
//					}
//					return (idx*GAP_SIZE_BYTES+1);
//				}
//				memcpy(&out[idx*GAP_SIZE_BYTES+1], &gap, GAP_SIZE_BYTES);
//				idx++;
//			}
//
//		} else if (inside) {
//			//take care of the last gap
//			inside = false;
//			gap = ppos_end_bit - (bitcnt - tmpcnt);
//
//			if (!flag) {
//				later_0 += gap;
//			} else {
//				memcpy(&out[idx*GAP_SIZE_BYTES+1], &gap, GAP_SIZE_BYTES);
//				idx++;
//			}
//			if (later_0 > 0) {
//				memcpy(&out[idx*GAP_SIZE_BYTES+1], &later_0, GAP_SIZE_BYTES);
//				idx++;
//			}
//			return (idx*GAP_SIZE_BYTES+1);
//
//		} else if (bitcnt - tmpcnt < ppos_start_bit && bitcnt > ppos_end_bit) {
//			if (!flag) {
//				//all 0s
////				unsigned int t = num_preds * gobject_bytes *8;
////				memcpy(&out[idx*4+1], &t, 4);
////				idx++;
////				free(out);
////				out = NULL;
//				return 0;
//			} else {
//				//ini_0s then gobject_bytes *8 1s and then later_0s
//				if (ini_0 > 0) {
//					memcpy(&out[idx*GAP_SIZE_BYTES+1], &ini_0, GAP_SIZE_BYTES);
//					idx++;
//				}
//				memcpy(&out[idx*GAP_SIZE_BYTES+1], &t, GAP_SIZE_BYTES);
//				idx++;
//				if (later_0 > 0) {
//					memcpy(&out[idx*GAP_SIZE_BYTES+1], &later_0, GAP_SIZE_BYTES);
//					idx++;
//				}
//				return (idx*GAP_SIZE_BYTES+1);
//			}
//		}
//	}
//
//	return (idx*GAP_SIZE_BYTES+1);
//
//}

////////////////////////////////////////////////////////////

// has got tested
//unsigned int var_p_fixed_o(unsigned char *in, unsigned int size, unsigned char *out,
//							unsigned int opos, unsigned char *and_array, unsigned int and_array_size)
//{
//
//
//   ////////////////////////////////////
//   // DEBUG
//   ////////////////////////////////////
////	unsigned int cnt = 0;
////	unsigned int tmpcnt, bitcnt;
////
////	bitcnt = 0;
////	printf("gobject_bytes * 8 = %u\n", gobject_bytes*8);
////	printf("[%u] ", and_array[0]);
////	while (cnt < (idx-1)/4) {
////		memcpy(&tmpcnt, &and_array[cnt*4+1], 4);
////		bitcnt += tmpcnt;
////		printf("%u ", tmpcnt);
////		cnt++;
////	}
////	printf("\nbitcnt %u %u\n", bitcnt, num_preds*gobject_bytes*8);
//
//	unsigned int out_size = dgap_AND(in, size, and_array, and_array_size, out);
//
////	printf("exiting var_p_fixed_o out_size %u\n", out_size);
//
//	return out_size;
//
//}

////////////////////////////////////////////////////////////

unsigned long count_triples_in_bitmat(BitMat *bitmat)
{
	unsigned int i;
	unsigned char *resptr;
	unsigned int rowcnt = 0, size = 0;
	unsigned int ressize = 0;

//	rowcnt = bitmat->num_subs;

	bitmat->num_triples = 0;

	if (bitmat->bm.size() == 0 && bitmat->vbm.size() == 0) {
		return 0;
	}

	unsigned char *data = NULL;

	for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); it++) {
		ressize = 0; size = 0;
		data = (*it).data;
		memcpy(&ressize, data, BM_ROW_SIZE);
		size = TOTAL_ROWSIZE(ressize);
//		cout << "Rowid " << (*it).rowid << " size " << ressize << endl;
		bitmat->num_triples += count_triples_in_row(&data[BM_ROW_SIZE], size);
	}

	if (bitmat->num_triples == 0) {
		for (vector<struct row>::iterator it = bitmat->vbm.begin(); it != bitmat->vbm.end(); it++) {
			ressize = 0; size = 0;
			data = (*it).data;
			memcpy(&ressize, data, BM_ROW_SIZE);
			size = TOTAL_ROWSIZE(ressize);
			bitmat->num_triples += count_triples_in_row(&data[BM_ROW_SIZE], size);
		}
	}

//	printf("Total triples %u\n", bitmat->num_triples);

	return bitmat->num_triples;

}

void threaded_simple_col_fold(std::list<struct row>::iterator it_begin,
		std::list<struct row>::iterator it_end, unsigned char *foldarr, unsigned int foldarr_size)
{
	for (std::list<struct row>::iterator it = it_begin; it != it_end; it++) {
		unsigned char *data = (*it).data;
//			unsigned int rowsize = 0;
		unsigned int sz = 0;
		memcpy(&sz, data, BM_ROW_SIZE);
		dgap_uncompress(data + BM_ROW_SIZE, TOTAL_ROWSIZE(sz), foldarr, foldarr_size);
	}
}

void threaded_simple_row_fold(std::list<struct row>::iterator it_begin,
		std::list<struct row>::iterator it_end, unsigned char *foldarr, unsigned int foldarr_size)
{
	for (std::list<struct row>::iterator it = it_begin; it != it_end; it++) {
		unsigned int rowbit = (*it).rowid - 1;
		assert(rowbit/8 < foldarr_size);
		foldarr[rowbit/8] |= (0x80 >> (rowbit%8));
	}
}

////////////////////////////////////////////////////////////
/*
 * Rule
 * folded S is always foldarr[NUM_OF_ALL_SUB/8]
 * folded P is always foldarr[NUM_OF_ALL_PRED/8]
 * folded O is always foldarr[column_bytes]
 * the caller function has to make sure that the foldarr
 * array is initialized properly
 */
void simple_fold(BitMat *bitmat, int ret_dimension, unsigned char *foldarr, unsigned int foldarr_size)
{
//	printf("Inside fold\n");
	memset(foldarr, 0, foldarr_size);

	if (ret_dimension == ROW) {
		/*
		 * Parallelized
		 */
//		for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); it++) {
//			unsigned int rowbit = (*it).rowid - 1;
//			assert(rowbit/8 < foldarr_size);
//			foldarr[rowbit/8] |= (0x80 >> (rowbit%8));
//		}

		thread myThreads[NUMTHREADS];
		std::list<struct row>::iterator it = bitmat->bm.begin(), it_begin, it_end;
		unsigned int bitmat_rows = bitmat->bm.size();
		unsigned int blocksz = bitmat_rows / NUMTHREADS;
		unsigned int remaining = bitmat_rows;
		unsigned int total_thr = 0;

		for (unsigned short i = 0; i < NUMTHREADS && it != bitmat->bm.end(); i++, total_thr++) {
			unsigned int curr_alloc = min(blocksz, remaining);
			it_begin = it;
			advance(it, curr_alloc);
			it_end = it;
			remaining -= curr_alloc;

			while (it_end != bitmat->bm.end() && ((*it_end).rowid - 1) % 8 > 0) {
				it_end++;
				remaining--;
			}

			myThreads[i] = thread(threaded_simple_row_fold, it_begin, it_end, foldarr, foldarr_size);
			it = it_end;
		}
		for (unsigned short i = 0; i < total_thr; i++) {
			myThreads[i].join();
		}

///////////////////////////////////////////////////////////////////////////////////////

	} else if (ret_dimension == COLUMN) {
//		if (bitmat->last_unfold == COLUMN) {
//			assert(foldarr_size == bitmat->column_bytes);
//			memcpy(foldarr, bitmat->colfold, bitmat->column_bytes);
//		} else {
		thread myThreads[NUMTHREADS];
		unsigned char *local_foldarr[NUMTHREADS];
		std::list<struct row>::iterator it = bitmat->bm.begin(), it_begin, it_end;
		unsigned int bitmat_rows = bitmat->bm.size();

		unsigned int blocksz = bitmat_rows % NUMTHREADS > 0 ?
				bitmat_rows / NUMTHREADS + 1 : bitmat_rows / NUMTHREADS;
		unsigned int remaining = bitmat_rows, total_thr = 0;

		for (unsigned short i = 0; i < NUMTHREADS && it != bitmat->bm.end(); i++, total_thr++) {
			unsigned int curr_alloc = min(blocksz, remaining);
			it_begin = it;
//			unsigned int remaining = bitmat_rows - i*blocksz;
			advance(it, curr_alloc);
			it_end = it;
			remaining -= curr_alloc;
			local_foldarr[i] = (unsigned char *) malloc (foldarr_size);
			memset(local_foldarr[i], 0x00, foldarr_size);
			myThreads[i] = thread(threaded_simple_col_fold, it_begin, it_end, local_foldarr[i], foldarr_size);
			it = it_end;
		}

		for (unsigned short i = 0; i < total_thr; i++) {
			myThreads[i].join();
		}

		for (unsigned int j = 0; j < foldarr_size; j++) {
			for (unsigned short i = 0; i < total_thr; i++) {
				foldarr[j] |= local_foldarr[i][j];
			}
		}

		for (unsigned short i = 0; i < total_thr; i++) {
			free (local_foldarr[i]);
		}

//		for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); it++) {
//			unsigned char *data = (*it).data;
////			unsigned int rowsize = 0;
//			unsigned short sz = 0;
//			memcpy(&sz, data, BM_ROW_SIZE);
//
//			dgap_uncompress(data + BM_ROW_SIZE, ((sz-1)*GAP_SIZE_BYTES+1), foldarr, foldarr_size);
//
//		}
//		}
	} else {
		cout << "simple_fold: **** ERROR unknown dimension " << ret_dimension << endl;
		assert(0);
	}
}


unsigned char * rowwise_column_unfold(unsigned char *data, unsigned char *maskarr, unsigned int maskarr_size)
{
//	if (maskarr_bits == bitmat->num_objs) {
//		printf("All objects present.. didn't unfold anything\n");
//		return data;
//	}
	if (maskarr_size == 0) {
		free(data);
		data = NULL;
		return data;
	}
	unsigned int maskarr_bits = count_bits_in_row(maskarr, maskarr_size);
	unsigned char *andres = (unsigned char *) malloc (GAP_SIZE_BYTES * 2 * maskarr_bits + 1 + 1024);
//	unsigned int orig_bm_triples = bitmat->num_triples;
	unsigned int andres_size = 0;
//	unsigned int rowsize = 0;
	unsigned int sz = 0;

	memcpy(&sz, data, BM_ROW_SIZE);
	data += BM_ROW_SIZE;

	/*
	 * prev0_bit: position of the most recent 0 bit
	 * prev1_bit: position of the most recent 1 bit
	 * p0bit_set: if most recent bit seen in prev gap was 0
	 * p1bit_set: if most recent bit seen in prev gap was 1
	 *
	 * The logic is that whenever you are in 1 bit zone
	 * you only check if there was a previously set 1 bit
	 * more than 1 bit away from the current one,
	 * so as to push the intermediate 0s and vice versa.
	 *
	 * p0bit_set and p1bit_set are NEVER set to false once
	 * set to true.
	 */

	unsigned int cnt = 0, total_gaps = (sz-1), tmpcnt = 0, triplecnt = 0,
			prev1_bit = 0, prev0_bit = 0, tmpval = 0;
	bool flag = data[0] & 0x01, begin = true, p1bit_set = false, p0bit_set = false;//, done = false;
	unsigned char format = data[0] & 0x02;
	if (format == 0x02) {
		begin = false;
		andres[0] = 0x02;
		andres_size = 1;
	}

	while (cnt < total_gaps) {
		//tmpcnt is going to be the position of 1s in case of alternate format,
		//positions begin with index 1, and NOT with 0
		memcpy(&tmpcnt, &data[cnt*GAP_SIZE_BYTES+1], GAP_SIZE_BYTES);
//				cout << "triplecnt " << triplecnt << " tmpcnt " << tmpcnt << endl;

		if (format == 0x02) {
			if ((tmpcnt-1)/8 >= maskarr_size) {
			} else if (maskarr[(tmpcnt-1)/8] & (0x80 >> ((tmpcnt-1)%8))) {
				//The bit (position) under consideration is present in maskarr
				memcpy(&andres[andres_size], &tmpcnt, GAP_SIZE_BYTES);
				andres_size += GAP_SIZE_BYTES;
			}
		} else if (format == 0x00) {
			//triplecnt is number of bits already read from gaps.
			//if maskarr size is already overshot by the number of bits read
			//then we don't need to do anything.
			if (flag) {
				if (triplecnt/8 >= maskarr_size) {
					//p0bit_set = true;
					prev0_bit = triplecnt + tmpcnt - 1;
//							cout << "prev0_bit set1 " << prev0_bit << endl;
				} else {
					for (unsigned int i = triplecnt; i < triplecnt + tmpcnt; i++) {
						if (i/8 >= maskarr_size) {
							//An existing bit is unset to due to maskarr
							//so decreasing the triple count
							assert(p0bit_set || p1bit_set);
							if (!p0bit_set && p1bit_set) {
								tmpval = prev1_bit + 1;
								memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
								andres_size += GAP_SIZE_BYTES;
								p0bit_set = true;
							} else if (prev1_bit > prev0_bit) {
								tmpval = prev1_bit - prev0_bit;
								memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
								andres_size += GAP_SIZE_BYTES;
							}
							prev0_bit = i;
//									cout << "prev0_bit set2 " << prev0_bit << endl;
//									if (numpass == 1) {
//										done = true;
//										break;
//									}
						} else if ((maskarr[i/8] & (0x80 >> (i%8))) == 0x00) {
							//An existing bit is unset to due to maskarr
							if (begin) {
								begin = false;
								andres[0] = 0x00;
								andres_size = 1;
								p0bit_set = true;
							} else {
								if (!p0bit_set) {//this means all prev bits are 1s
									//pushing previous # of 1s
									memcpy(&andres[andres_size], &i, GAP_SIZE_BYTES);
									andres_size += GAP_SIZE_BYTES;
									p0bit_set = true;
								} else if (prev0_bit != (i - 1)) {//this means there are one or more 1s between this 0 bit and prev one
									//pushing previous # of 1s
									tmpval = i - prev0_bit -1;
									memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
									andres_size += GAP_SIZE_BYTES;
								}
							}
							prev0_bit = i;
//									cout << "prev0_bit set3 " << prev0_bit << endl;

						} else if (maskarr[i/8] & (0x80 >> (i%8))) {
							if (begin) {
								begin = false;
								andres[0] = 0x01;
								andres_size = 1;
								p1bit_set = true;
							} else {
								if (!p1bit_set) {//means all prev bits are 0s only
									//pushing previous # of 0s
									memcpy(&andres[andres_size], &i, GAP_SIZE_BYTES);
									andres_size += GAP_SIZE_BYTES;
									p1bit_set = true;
								} else if (prev1_bit != (i - 1)) {
									//pushing previous # of 0s
									tmpval = i - prev1_bit -1;
									memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
									andres_size += GAP_SIZE_BYTES;
								}
							}
							prev1_bit = i;
//									cout << "prev1_bit set " << prev1_bit << endl;
						}
					}
				}
			} else {
				if (begin) {
					begin = false;
					andres[0] = 0x00;
					andres_size = 1;
					p0bit_set = true;
				} else {
					if (!p0bit_set) {
						//pushing previous # of 1s
						tmpval = prev1_bit+1;
						memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
						andres_size += GAP_SIZE_BYTES;
						p0bit_set = true;
					} else if (prev0_bit != (triplecnt - 1)) {
						//pushing previous # of 1s
						tmpval = triplecnt - prev0_bit -1;
						memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
						andres_size += GAP_SIZE_BYTES;
					}
				}
				prev0_bit = triplecnt + tmpcnt - 1;
//						cout << "prev0_bit set4 " << prev0_bit << " triplecnt " << triplecnt << " tmpcnt " << tmpcnt << endl;
			}
		} else {
			printf("Unrecognized format %x\n", format);
			assert(0);
		}

		flag = !flag;
		cnt++;
		triplecnt += tmpcnt;
	}
	if (format == 0x00) {
		//COMMENT: no need to push last 0 bits.
		if (prev1_bit > prev0_bit) {
			//COMMENT: Push last 1 bits.
			if (!p0bit_set) {
				tmpval = (prev1_bit + 1);
				memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
			} else {
				tmpval = (prev1_bit - prev0_bit);
				memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
			}
			andres_size += GAP_SIZE_BYTES;
		} else if (prev1_bit == prev0_bit && p1bit_set && p0bit_set) {
			//if both bits have been set before and both are 0, then something is wrong here
//					cout << "prev0_bit " << prev0_bit << " prev1_bit " << prev1_bit << endl;
			assert(0);
		}

		if ((!andres[0]) && andres_size == 1) {
			//empty AND res
			data -= BM_ROW_SIZE;
			free(data);
			data = NULL;
		}
	} else if (andres_size == 1) {
		data -= BM_ROW_SIZE;
		free(data);
		data = NULL;
	}

	if (data != NULL) {
//		tmpval = andres_size;
//		memcpy(andres, &andres_size, ROW_SIZE_BYTES);
		data -= BM_ROW_SIZE;
		free(data);
//		assert(ROW_SIZE_HEADER(andres_size) < 0xffff);
		assert(IS_ROW_SIZE_SHORT(andres_size));
//		assert(((andres_size-1)/GAP_SIZE_BYTES + 1) < 0xffff);
		data = (unsigned char *) calloc (1, andres_size + BM_ROW_SIZE);
		sz = ROW_SIZE_HEADER(andres_size);// (andres_size-1)/GAP_SIZE_BYTES + 1;
		memcpy(data, &sz, BM_ROW_SIZE);
		memcpy(&data[BM_ROW_SIZE], andres, andres_size);
	}

	return data;

}

void threaded_simple_col_unfold(BitMat *bitmat, list<struct row>::iterator it_begin,
		list<struct row>::iterator it_end, unsigned char *maskarr, unsigned int maskarr_size)
{
	for (list<struct row>::iterator it = it_begin; it != it_end; ) {
		(*it).data = rowwise_column_unfold((*it).data, maskarr, maskarr_size);
		if ((*it).data == NULL && it != it_begin) {
			it = bitmat->bm.erase(it);
			continue;
		}
		it++;

	}
}

//////////////////////////////////////////////////////////
/*
 * Modified simple_unfold updates the count of triples in bitmat
 * dynamically as it unfolds on the "column" dimension.
 * For the "row" dimension, it does so only if it is asked
 * to do so through "numpass" variable.
 */
void simple_unfold(BitMat *bitmat, unsigned char *maskarr, unsigned int maskarr_size,
					int dimension)
{
//	cout << "simple_unfold " << numpass << " dimension " << ret_dimension << endl;

	if (dimension == ROW) {
		if (count_bits_in_row(maskarr, maskarr_size) == bitmat->num_rows) {
			//don't need to unfold anything just return
			printf("All subjects present.. didn't unfold anything\n");
			return;
		}
		if (maskarr_size == 0) {
			for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); ) {
				free((*it).data);
				it = bitmat->bm.erase(it);
			}
			bitmat->num_triples = 0;
//			bitmat->last_unfold = ROW;
			return;
		}

//		unsigned int bm_orig_rows = bitmat->bm.size();

		for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); ) {
			unsigned int rowbit = (*it).rowid - 1;
			if (rowbit/8 >= maskarr_size) {
				free((*it).data);
				it = bitmat->bm.erase(it);
				continue;
			}
			if (!( maskarr[rowbit/8] & (0x80 >> (rowbit%8)))) {
				free((*it).data);
				it = bitmat->bm.erase(it);
				continue;
			}
			it++;
		}
//		if (bm_orig_rows - bitmat->bm.size() > 0) {
//			bitmat->last_unfold = ROW;
//			if (numpass > 1) {
//				bitmat->num_triples = count_triples_in_bitmat(bitmat);
//			}
//		}

////////////////////////////////////////////////////////////////////////////

	} else if (dimension == COLUMN) {
		//use concatenation function
//		cout << "*********** simple_unfold OBJ_DIMENSION **********" << endl;
		unsigned int andres_size = 0;
		unsigned int maskarr_bits = count_bits_in_row(maskarr, maskarr_size);

		if (maskarr_bits == bitmat->num_columns) {
			printf("All objects present.. didn't unfold anything\n");
			return;
		}
//		if (count_triples_in_row(maskarr, maskarr_size) == bitmat->num_objs) {
//			//don't need to do anything
//			printf("All objects present.. didn't unfold anything\n");
//			return;
//		}

//		cumulative_dgap(maskarr, maskarr_size, cumu_maskarr);

//		cout << "bitmat->numtriples " << bitmat->num_triples << endl;

		if (maskarr_size == 0) {
			printf("simple_unfold: maskarr_size is 0\n");
			for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); ) {
				unsigned char *data  = (*it).data;
				free(data);
				it = bitmat->bm.erase(it);
			}
//			bitmat->last_unfold = COLUMN;
			bitmat->num_triples = 0;
			return;
		}

		//////////
		//DEBUG
		//TODO: remove later
		///////
//		cout << "---------- simple_unfold: maskarr_bits" << endl;
//		unsigned int maskarr_set_bits =	print_set_bits_in_row(maskarr, maskarr_size);
		///////

		thread myThreads[NUMTHREADS];
		std::list<struct row>::iterator it = bitmat->bm.begin(), it_begin, it_end;
		std::list<struct row>::iterator it_start[NUMTHREADS];

		unsigned int bitmat_rows = bitmat->bm.size();

		unsigned int blocksz = bitmat_rows % NUMTHREADS > 0 ?
				bitmat_rows / NUMTHREADS + 1 : bitmat_rows / NUMTHREADS;

		for (unsigned short i = 0; i < NUMTHREADS && it != bitmat->bm.end(); i++) {
			it_begin = it;
			unsigned int remaining = bitmat_rows - i*blocksz;
			advance(it, min(blocksz, remaining));
			it_end = it;
			it_start[i] = it_begin;

			myThreads[i] = thread(threaded_simple_col_unfold, bitmat, it_begin, it_end, maskarr, maskarr_size);
			it = it_end;
		}

		for (unsigned short i = 0; i < NUMTHREADS; i++) {
			myThreads[i].join();
		}

		/*
		 * This is required after all threads are done,
		 * because we cannot erase start points of threads
		 * from bitmat->bm list in a threaded cond
		 * due to possible race
		 */
		for (unsigned short i = 0; i < NUMTHREADS; i++) {
			if ((*(it_start[i])).data == NULL) {
				bitmat->bm.erase(it_start[i]);
			}
		}

//		for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); ) {
//			if ((*it).data == NULL) {
//				it = bitmat->bm.erase(it);
//				continue;
//			}
//			it++;
//		}
//
		/*
		 * This has been parallelized
		 */

//		unsigned char *andres = (unsigned char *) malloc (GAP_SIZE_BYTES * 2 * maskarr_bits + BM_ROW_SIZE + 1 + 1024);
//
//		unsigned int orig_bm_triples = bitmat->num_triples;
//
//		for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); ) {
//			unsigned char *data  = (*it).data;
////			unsigned int rowsize = 0;
//			unsigned short sz = 0;
//			memcpy(&sz, data, BM_ROW_SIZE);
//			data += BM_ROW_SIZE;
//			unsigned int cnt = 0, total_cnt = sz-1, tmpcnt = 0, triplecnt = 0,
//					prev1_bit = 0, prev0_bit = 0, tmpval = 0;
//			bool flag = data[0] & 0x01, begin = true, p1bit_set = false, p0bit_set = false;//, done = false;
//			unsigned char format = data[0] & 0x02;
//			if (format == 0x02) {
//				begin = false;
//				andres[ROW_SIZE_BYTES] = 0x02;
//				andres_size = ROW_SIZE_BYTES + 1;
//			}
//
////			cout << "cnt " << cnt << " total_cnt " << total_cnt << " rowid " << (*it).rowid << endl;
////			printf("format %x\n", format);
//
//			while (cnt < total_cnt) {
//				//tmpcnt is going to be the position of 1s in case of alternate format,
//				//positions begin with index 1, and NOT with 0
//				memcpy(&tmpcnt, &data[cnt*GAP_SIZE_BYTES+1], GAP_SIZE_BYTES);
////				cout << "triplecnt " << triplecnt << " tmpcnt " << tmpcnt << endl;
//
//				if (format == 0x02) {
//					if ((tmpcnt-1)/8 >= maskarr_size) {
////						if (numpass == 1) {
////							break;
////						}
////						bitmat->num_triples--;
//					} else if (maskarr[(tmpcnt-1)/8] & (0x80 >> ((tmpcnt-1)%8))) {
//						//The bit (position) under consideration is present in maskarr
//						memcpy(&andres[andres_size], &tmpcnt, GAP_SIZE_BYTES);
//						andres_size += GAP_SIZE_BYTES;
//					} else {
////						bitmat->num_triples--;
//					}
//				} else if (format == 0x00) {
//					//triplecnt is number of bits already read from gaps.
//					//if maskarr size is already overshot by the number of bits read
//					//then we don't need to do anything.
////					if (numpass == 1 && triplecnt/8 >= maskarr_size) {
////						break;
////					}
//					if (flag) {
//						if (triplecnt/8 >= maskarr_size) {
////							bitmat->num_triples -= tmpcnt;
//							//p0bit_set = true;
//							prev0_bit = triplecnt + tmpcnt - 1;
////							cout << "prev0_bit set1 " << prev0_bit << endl;
//						} else {
//							for (unsigned int i = triplecnt; i < triplecnt + tmpcnt; i++) {
//								if (i/8 >= maskarr_size) {
//									//An existing bit is unset to due to maskarr
//									//so decreasing the triple count
////									bitmat->num_triples--;
//									assert(p0bit_set || p1bit_set);
//									if (!p0bit_set && p1bit_set) {
//										tmpval = prev1_bit + 1;
//										memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//										andres_size += GAP_SIZE_BYTES;
//										p0bit_set = true;
//									} else if (prev1_bit > prev0_bit) {
//										tmpval = prev1_bit - prev0_bit;
//										memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//										andres_size += GAP_SIZE_BYTES;
//									}
//									prev0_bit = i;
////									cout << "prev0_bit set2 " << prev0_bit << endl;
////									if (numpass == 1) {
////										done = true;
////										break;
////									}
//								} else if ((maskarr[i/8] & (0x80 >> (i%8))) == 0x00) {
//									//An existing bit is unset to due to maskarr
////									bitmat->num_triples--;
//									if (begin) {
//										begin = false;
//										andres[ROW_SIZE_BYTES] = 0x00;
//										andres_size = ROW_SIZE_BYTES +1;
//										p0bit_set = true;
//									} else {
//										if (!p0bit_set) {//this means all prev bits are 1s
//											//pushing previous # of 1s
//											memcpy(&andres[andres_size], &i, GAP_SIZE_BYTES);
//											andres_size += GAP_SIZE_BYTES;
//											p0bit_set = true;
//										} else if (prev0_bit != (i - 1)) {//this means there are one or more 1s between this 0 bit and prev one
//											//pushing previous # of 1s
//											tmpval = i - prev0_bit -1;
//											memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//											andres_size += GAP_SIZE_BYTES;
//										}
//									}
//									prev0_bit = i;
////									cout << "prev0_bit set3 " << prev0_bit << endl;
//
//								} else if (maskarr[i/8] & (0x80 >> (i%8))) {
//									if (begin) {
//										begin = false;
//										andres[ROW_SIZE_BYTES] = 0x01;
//										andres_size = ROW_SIZE_BYTES + 1;
//										p1bit_set = true;
//									} else {
//										if (!p1bit_set) {//means all prev bits are 0s only
//											//pushing previous # of 0s
//											memcpy(&andres[andres_size], &i, GAP_SIZE_BYTES);
//											andres_size += GAP_SIZE_BYTES;
//											p1bit_set = true;
//										} else if (prev1_bit != (i - 1)) {
//											//pushing previous # of 0s
//											tmpval = i - prev1_bit -1;
//											memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//											andres_size += GAP_SIZE_BYTES;
//										}
//									}
//									prev1_bit = i;
////									cout << "prev1_bit set " << prev1_bit << endl;
//								}
//							}
//						}
//
////						if (done) break;
//
//					} else {
//						if (begin) {
//							begin = false;
//							andres[ROW_SIZE_BYTES] = 0x00;
//							andres_size = ROW_SIZE_BYTES +1;
//							p0bit_set = true;
//						} else {
//							if (!p0bit_set) {
//								//pushing previous # of 1s
//								tmpval = prev1_bit+1;
//								memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//								andres_size += GAP_SIZE_BYTES;
//								p0bit_set = true;
//							} else if (prev0_bit != (triplecnt - 1)) {
//								//pushing previous # of 1s
//								tmpval = triplecnt - prev0_bit -1;
//								memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//								andres_size += GAP_SIZE_BYTES;
//							}
//						}
//						prev0_bit = triplecnt + tmpcnt - 1;
////						cout << "prev0_bit set4 " << prev0_bit << " triplecnt " << triplecnt << " tmpcnt " << tmpcnt << endl;
//					}
//				} else {
//					printf("Unrecognized format %x\n", format);
//					assert(0);
//				}
//
//				flag = !flag;
//				cnt++;
//				triplecnt += tmpcnt;
//			}
//			if (format == 0x00) {
//				//COMMENT: no need to push last 0 bits.
//				if (prev1_bit > prev0_bit) {
//					//COMMENT: Push last 1 bits.
//					if (!p0bit_set) {
//						tmpval = (prev1_bit + 1);
//						memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//					} else {
//						tmpval = (prev1_bit - prev0_bit);
//						memcpy(&andres[andres_size], &tmpval, GAP_SIZE_BYTES);
//					}
//					andres_size += GAP_SIZE_BYTES;
//				} else if (prev1_bit == prev0_bit && p1bit_set && p0bit_set) {
//					//if both bits have been set before and both are 0, then something is wrong here
////					cout << "prev0_bit " << prev0_bit << " prev1_bit " << prev1_bit << endl;
//					assert(0);
//				}
//
//				if ((!andres[ROW_SIZE_BYTES]) && andres_size == 1+ROW_SIZE_BYTES) {
//					//empty AND res
////					cout << "removed rowid1 " << (*it).rowid << endl;
//					data -= BM_ROW_SIZE;
//					free(data);
//					it = bitmat->bm.erase(it);
//					continue;
//				}
//			} else {
//				if (andres_size == ROW_SIZE_BYTES+1) {
////					cout << "removed rowid2 " << (*it).rowid << endl;
//					data -= BM_ROW_SIZE;
//					free(data);
//					it = bitmat->bm.erase(it);
//					continue;
//				}
//			}
//
////			cumulative_dgap(data + ROW_SIZE_BYTES, rowsize, data + ROW_SIZE_BYTES);
////			andres_size = dgap_AND(data + ROW_SIZE_BYTES, rowsize, cumu_maskarr, maskarr_size, &andres[ROW_SIZE_BYTES]);
////			de_cumulative_dgap(&andres[ROW_SIZE_BYTES], andres_size, &andres[ROW_SIZE_BYTES]);
////			if (!andres[ROW_SIZE_BYTES] && andres_size == 1+GAP_SIZE_BYTES) {
////				//empty AND res
////				free(data);
////				it = bitmat->bm.erase(it);
////				continue;
////			}
//
////			bool format_changed = false;
//			/*
//			 * To be enabled if memory performance suffers
//			if ((andres[ROW_SIZE_BYTES] & 0x02) == 0x00) {
//				unsigned int total_gaps = count_dgaps_in_row(&andres[ROW_SIZE_BYTES], andres_size - ROW_SIZE_BYTES);
//				unsigned int total_setbits = count_triples_in_row(&andres[ROW_SIZE_BYTES], andres_size - ROW_SIZE_BYTES);
//				if (total_gaps > total_setbits) {
//					unsigned char *tmpdata = (unsigned char *) malloc (1 + total_setbits*GAP_SIZE_BYTES + ROW_SIZE_BYTES);
//					bool flag = andres[ROW_SIZE_BYTES] & 0x01;
//					unsigned int cnt = 0, total_bits = 0, gapcnt = 0;
//					tmpdata[ROW_SIZE_BYTES] = 0x02; //format
//					unsigned int size = 1;
//					while (cnt < total_gaps) {
//						memcpy(&gapcnt, &andres[ROW_SIZE_BYTES + cnt*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
//						if (flag) {
//							for (unsigned int i = total_bits+1; i <= total_bits+gapcnt; i++) {
//								memcpy(&tmpdata[size + ROW_SIZE_BYTES], &i, GAP_SIZE_BYTES);
//								size += GAP_SIZE_BYTES;
//							}
//						}
//						flag = !flag;
//						cnt++;
//						total_bits += gapcnt;
//					}
//					assert(size == (1 + total_setbits * GAP_SIZE_BYTES));
//					memcpy(tmpdata, &size, ROW_SIZE_BYTES);
//					data -= ROW_SIZE_BYTES;
//					free(data);
//					(*it).data = tmpdata;
//					format_changed = true;
//				}
//			}
//	 */
//
////			if (!format_changed) {
//			tmpval = (andres_size - ROW_SIZE_BYTES);
////				memcpy(andres, &tmpval, ROW_SIZE_BYTES);
//			data -= BM_ROW_SIZE;
//			free(data);
//			assert((tmpval-1)/GAP_SIZE_BYTES < 0xffff);
//			sz = (tmpval-1)/GAP_SIZE_BYTES + 1;
//			(*it).data = (unsigned char *) malloc (tmpval + BM_ROW_SIZE);
//			memcpy((*it).data, &sz, BM_ROW_SIZE);
//			memcpy(&(*it).data[BM_ROW_SIZE], &andres[ROW_SIZE_BYTES], tmpval);
////			}
//			it++;
//		}
////		if (orig_bm_triples - bitmat->num_triples > 0) {
////			bitmat->last_unfold = COLUMN;
////		}
//		free(andres);

	} else {
		assert(0);
	}
}

void conditional_fold(BitMat *bitmat, int ret_dimension, unsigned char *foldarr, unsigned int foldarr_size)
{
	if (ret_dimension == ROW) {
		if (bitmat->last_op == ROW_VEC || bitmat->last_op == ROW_UNFOLD || bitmat->last_op == BMUPTODATE) {
			if (foldarr != bitmat->rowfold) {
				memset(foldarr, 0x00, foldarr_size);
				memcpy(foldarr, bitmat->rowfold, bitmat->row_bytes < foldarr_size ?
						bitmat->row_bytes : foldarr_size);
			}
		} else {
			if (bitmat->last_op == COLUMN_VEC) {
				simple_unfold(bitmat, bitmat->colfold, bitmat->column_bytes, COLUMN);
				bitmat->last_op = COLUMN_UNFOLD;
			}
			simple_fold(bitmat, ROW, bitmat->rowfold, bitmat->row_bytes);
			bitmat->last_op = BMUPTODATE;
		}
	} else if (ret_dimension == COLUMN) {
		if (bitmat->last_op == COLUMN_VEC || bitmat->last_op == COLUMN_UNFOLD || bitmat->last_op == BMUPTODATE) {
			if (foldarr != bitmat->colfold) {
				memset(foldarr, 0x00, foldarr_size);
				memcpy(foldarr, bitmat->colfold, bitmat->column_bytes < foldarr_size ?
						bitmat->column_bytes : foldarr_size);
			}
		} else {
			if (bitmat->last_op == ROW_VEC) {
				simple_unfold(bitmat, bitmat->rowfold, bitmat->row_bytes, ROW);
				bitmat->last_op = ROW_UNFOLD;
			}
			simple_fold(bitmat, COLUMN, bitmat->colfold, bitmat->column_bytes);
			bitmat->last_op = BMUPTODATE;
		}
	}
}
////////////////////////////////////////////////////////////
//unsigned int count_size_of_bitmat(unsigned char **bitmat)
//{
//#if USE_MORE_BYTES
//	unsigned long tmpcnt = 0;
//#else
//	unsigned int tmpcnt = 0;
//#endif
//	unsigned long size = 0, i;
//	unsigned int ptrsize = sizeof(unsigned char *);
//
//	if (bitmat == NULL)
//		return 0;
//
//	size += (gnum_subs * ptrsize);
//	for (i = 0; i < gnum_subs; i++) {
//		if (bitmat[i] == NULL)
//			continue;
//		memcpy(&tmpcnt, bitmat[i], ROW_SIZE_BYTES);
//		size += tmpcnt;
//	}
//
//	return size;
//}
////////////////////////////////////////////////////////////
//void explore_other_nodes(struct node *n, int curr_nodenum, struct node *parent, FILE *outfile)
//{
//	//COMMENT: instead of maintaining unvisited list
//	//every time when you don't find an unmapped neighbor
//	//of the node or its parent's neighbor, just go over the mapped
//	//nodes in the q_to_gr_node_map and pick an unmapped neighbor
//	//of the first mapped node
//
//	LIST *ptr = n->nextTP;
//
//	//First insert all the unvisited neighbors in the "unvisited" set
//	//every time you call match_query_graph on an unvisited node
//	//remove it from the set
//	for (; ptr; ) {
//
//		assert(ptr->gnode->isNeeded);
////		if (!ptr->gnode->isNeeded) {
////			ptr = ptr->next;
////			continue;
////		}
//
//		if (q_to_gr_node_map.find((TP *)ptr->gnode->data) == q_to_gr_node_map.end()) {
//
//			unvisited.insert(unvisited.end(), ptr->gnode);
//			
//		}
//		ptr = ptr->next;
//	}
//
//	//3.2: Medha: now call match_query_graph on the first unvisited neighbor
//	ptr = n->nextTP;
//	bool explored = false;
////	bool newly_mapped = true;
//
////	for (int k = 0; k < n->numTPs && newly_mapped; k++) 
//	for (; ptr; ) {
//
//		assert(ptr->gnode->isNeeded);
////		if (!ptr->gnode->isNeeded) {
////			ptr = ptr->next;
////			continue;
////		}
//
//		if (q_to_gr_node_map.find((TP *)ptr->gnode->data) == q_to_gr_node_map.end()) {
//			//this neighbor node is unmapped yet.. so explore that and 
//			
//			//make sure that the edge that you
//			//are exploring is a "tree" edge and not backedge
//			vector<struct node*>::iterator it = find(tree_edges_map[n].begin(), tree_edges_map[n].end(), ptr->gnode);
//			if (it != tree_edges_map[n].end()) {
//				explored = true;
//				//remove this neighbor being explored from the unvisited set
//				//before calling match_query_graph on it
//				unvisited.erase(ptr->gnode);
//				TP *tp = (TP *)ptr->gnode->data;
////				cout << "Exploring neighbor node " << tp->sub << " " << tp->pred << " " << tp->obj  << endl;
//				match_query_graph(ptr->gnode, curr_nodenum+1, n, outfile);
////				newly_mapped = false;
//				return;
//			}
//		}
//		ptr = ptr->next;
//	}
//
//	//3.3: Medha: The very fact that you came here means you don't have a pathlike
//	//spanning tree edges so just start exploring the parent of the current node
//
//	if (!explored && parent != NULL) {
////		for (vector<struct node*>::iterator it = tree_edges_map[parent].begin();
////				it != tree_edges_map[parent].end() && newly_mapped; it++) 
//		for (vector<struct node*>::iterator it = tree_edges_map[parent].begin();
//				it != tree_edges_map[parent].end(); it++) {
//			if (q_to_gr_node_map.find((TP *)(*it)->data) == q_to_gr_node_map.end()) {
//				//cout << "4. calling match_query_graph on " << ((TP *)(*it)->data)->sub 
//				//	<< " " << ((TP *)(*it)->data)->pred << " " << ((TP *)(*it)->data)->obj << endl;
//				explored = true;
//				unvisited.erase((*it));
//				TP *tp = (TP *)(*it)->data;
////				cout << "Exploring2 neighbor node " << tp->sub << " " << tp->pred << " " << tp->obj  << endl;
//				match_query_graph((*it), curr_nodenum+1, parent, outfile);
////				newly_mapped = false;
//				return;
//			}
//
//		}
//	}
//
//	//3.4: Medha: if you come here and you still have "explored" == false, that means you have
//	//exhausted all the unvisited neighbors of your parent too but still you haven't
//	//finished visited all query nodes.. so pick a node from the "unvisited" set
//	//This can happen in case of a graph as shown at the top while defining "unvisited" set
//	if (!explored) {
//		//TODO: recheck
//		if (unvisited.size() == 0) {
//			for (int i = 0; i < graph_tp_nodes; i++) {
//				if (!graph[MAX_JVARS_IN_QUERY + i].isNeeded)
//					continue;
//				TP *tp = (TP *)graph[MAX_JVARS_IN_QUERY + i].data;
//				//Find an already mapped node and then find its unmapped
//				//neighbor
//				if (q_to_gr_node_map.find(tp) != q_to_gr_node_map.end()) {
//					ptr = graph[MAX_JVARS_IN_QUERY + i].nextTP;
//					for (; ptr; ) {
//
//						if (!ptr->gnode->isNeeded) {
//							ptr = ptr->next;
//							continue;
//						}
//
//						if (q_to_gr_node_map.find((TP *)ptr->gnode->data) == q_to_gr_node_map.end()) {
//							unvisited.insert(unvisited.end(), ptr->gnode);
//							break;
//						}
//						ptr = ptr->next;
//					}
//				}
//			}
//		}
//
//		set<struct node *>::iterator it = unvisited.begin();
//		explored = true;
//		unvisited.erase((*it));
//		//parent is set to NULL as we are randomly picking
//		//it up from unvisited set, hence we don't know who
//		//is its parent node.						
//		TP *tp = (TP *)(*it)->data;
////		cout << "Exploring3 neighbor node " << tp->sub << " " << tp->pred << " " << tp->obj  << endl;
//		match_query_graph((*it), curr_nodenum+1, NULL, outfile);
////		newly_mapped = false;
//	}
//
//}
//////////////////////////////////////////////////////////
//void explore_other_nodes(struct node *n, int curr_nodenum, struct node *parent, FILE *outfile)
//{
//
//	//TODO: fix here.. you need 2 vectors -- vec1, vec2, one is drained and
//	//another is populated with the draining values from the other.
//	map<struct node *, vector<struct node *> >::iterator nitr = tree_edges_map.find(n);
//	vector<struct node*> nVec ;
//
//	if (nitr != tree_edges_map.end()) {
//		nVec = tree_edges_map[n];
//	}
//
//	for (vector<struct node*>::iterator it = nVec.begin(); it != nVec.end(); it++) {
//		assert((*it)->isNeeded);
//
//		if (q_to_gr_node_map.find((TP *)(*it)->data) == q_to_gr_node_map.end()) {
//			//this neighbor node is unmapped yet.. so explore that and 
//			unvisited.push_back(*it);
//		}
//	}
//	assert(unvisited.size() != 0);
//	struct node *qfront = unvisited.front();
//	unvisited.pop_front();
//
//	match_query_graph(qfront, curr_nodenum+1, NULL, outfile);
//
//}

//////////////////////////////////////////////////////

//TODO: fix for OPTIONAL queries
//bool check_back_edges(struct node *curr_node, TP curr_matched_triple)
//{
//	struct node *n = curr_node;
//	bool matched = true;
//
//	LIST *ptr = n->nextTP;
//	struct triple t = curr_matched_triple;
//	TP *q_node = (TP *)n->data;
//
//	for(; ptr; ) {
//		// go over all the neighbors.
//
//		assert(ptr->gnode->isNeeded);
////		if (!ptr->gnode->isNeeded) {
////			ptr = ptr->next;
////			continue;
////		}
//
//		TP *q_nt = (TP *)ptr->gnode->data;
//
////		cout << "check_back_edges: Exploring [" << q_nt->sub << " " << q_nt->pred << " " << q_nt->obj << endl; 
//
//		if(q_to_gr_node_map.find(q_nt) != q_to_gr_node_map.end()) {
//			// Found occurrence of mapping for a neighbor in the map.
//			// nt == neighbor triple q_nt == query node associated with that triple
//
//			struct triple nt = q_to_gr_node_map[q_nt];
//
//			// Compare the corresponding entries based on the edge type.
//			int edge_type = ptr->edgetype;
//			switch(edge_type) {
//				case SUB_EDGE:
//					if(t.sub != nt.sub)
//						matched=false;
//					break;
//
//				case PRED_EDGE:
//					if(t.pred != nt.pred)
//						matched=false;
//					break;
//
//				case OBJ_EDGE:
//					if(t.obj != nt.obj)
//						matched=false;
//					break;
//
//				case SO_EDGE:
//					if((q_node->sub < 0) && (q_node->sub == q_nt->obj)) {
//						if(t.sub != nt.obj)
//							matched=false;
//					} else if((q_node->obj < 0) && (q_node->obj == q_nt->sub)) {
//						if(t.obj != nt.sub)
//							matched=false;
//					}
//					break;
//			}
//
//			if(!matched)
//				break;
//				
//		} 
//		ptr=ptr->next;
//
//	}
//
//	return matched;
//
//}
//void set_neighbors_null2(struct node *n, set<string> &null_tmp)
//{
//	LIST *ptr = n->nextTP;
//	struct triple t = {0, 0, 0};
//	TP *q_node = (TP *)n->data;
//
//	for(; ptr; ) {
//		assert(ptr->gnode->isNeeded);
//
//		TP *q_nt = (TP *)ptr->gnode->data;
//		map<struct node*, struct triple>::iterator it = q_to_gr_node_map.find(ptr->gnode);
//
//		if(it != q_to_gr_node_map.end()) {
//			struct triple tn = (*it).second;
//			if (tn.sub != 0 && tn.pred != 0 && tn.obj != 0) {
//				if ((q_nt->isSlave(q_node->strlevel))
//					||
//					(!q_node->isSlave(q_nt->strlevel) && !q_nt->isSlave(q_node->strlevel))
//					) {
//					q_to_gr_node_map[ptr->gnode] = t;
//					null_tmp.insert(q_nt->strlevel);
//					set_neighbors_null2(ptr->gnode, null_tmp);
//				}
//			}
//			//TODO: remove later
//			else if (tn.sub == 0 && tn.pred == 0 && tn.obj == 0) {} else {assert(0);}
//		}
//		ptr=ptr->next;
//	}
//}


//void set_neighbors_null(struct node *n, set<string> &null_tmp)
//{
//	//TODO: you have to do this recursively until entire map is consistent
//	LIST *ptr = n->nextTP;
//	struct triple t = {0, 0, 0};
//	TP *q_node = (TP *)n->data;
//
//	for(; ptr; ) {
//		// go over all the neighbors.
//
//		assert(ptr->gnode->isNeeded);
//
//		TP *q_nt = (TP *)ptr->gnode->data;
//		map<struct node*, struct triple>::iterator it = q_to_gr_node_map.find(ptr->gnode);
//
//		if(it != q_to_gr_node_map.end()) {
//			struct triple tn = (*it).second;
//			if (tn.sub != 0 && tn.pred != 0 && tn.obj != 0) {
//				if (q_node->strlevel.compare(q_nt->strlevel) == 0 ||
//					(!q_node->isSlave(q_nt->strlevel) && !q_nt->isSlave(q_node->strlevel))) {
//					q_to_gr_node_map[ptr->gnode] = t;
//					null_tmp.insert(q_nt->strlevel);
//					set_neighbors_null2(ptr->gnode, null_tmp);
//				}
//
//			}
//		}
//		ptr=ptr->next;
//	}
//}

/////////////////////////////////////////////////////

/*
 * This logic replaces the previous logic which was going
 * over the nullified nodes again and again.
 */
void set_neighbors_null(struct node **bfsarr, map<struct node*, struct triple>&resmap)
{
	deque<struct node *> bfsq;

	for (int i = 0; bfsarr[i]; i++) {
		TP *tp = (TP *)bfsarr[i]->data;
		struct triple t = resmap[bfsarr[i]];
		if (t.sub == 0 || t.pred == 0 || t.obj == 0) {
			bfsq.push_back(bfsarr[i]);
		}
	}

	while (bfsq.size() > 0) {
		struct node *node = bfsq.front(); bfsq.pop_front();
		TP *tp = (TP *)node->data;
		LIST *next = node->nextTP;

		for (; next; next = next->next) {
			TP *nbr = (TP *)next->gnode->data;

			/*
			 * With the changed structure that always expects a query graph
			 * in the well-designed form, you must not check for mutually
			 * excluside neighbors, as they do not exist
			 */
			if (nbr->isSlave(tp->strlevel)) {
				struct triple t = resmap[next->gnode];// it->second;
				if (t.sub != 0 || t.pred != 0 || t.obj != 0) {
					t.sub = 0; t.pred = 0; t.obj = 0;
					resmap[next->gnode] = t;
					bfsq.push_back(next->gnode);
				}
			}
		}
	}
}


void set_neighbors_null2(struct node *gnode, map<struct node*, struct triple>&resmap)
{
	deque<struct node *> bfsq;
	bfsq.push_back(gnode);
	struct node *node = NULL;

	while (bfsq.size() > 0) {
		node = bfsq.front(); bfsq.pop_front();
		TP *tp = (TP *)node->data;
		LIST *next = node->nextTP;

		for (; next; next = next->next) {
			/*
			 * You should not check this because now this
			 * function gets called only while printing the results
			 * when all the required TPs are already mapped to the
			 * respective triples.
			 */
//			map<struct node*, struct triple>::iterator it = resmap.find(next->gnode);
//			assert(it != resmap.end());
//			if (it == resmap.end())
//				continue;

			TP *neigh = (TP *)next->gnode->data;

			/*
			 * With the changed structure that always expects a query graph
			 * in the well-designed form, you must not check for mutually
			 * excluside neighbors, as they do not exist
			 */
//			if (neigh->isSlave(tp->strlevel) ||
//					(!neigh->isSlave(tp->strlevel) && !tp->isSlave(neigh->strlevel)))
			if (neigh->isSlave(tp->strlevel)) {
				struct triple t = resmap[next->gnode];// it->second;
				if (t.sub != 0 || t.pred != 0 || t.obj != 0) {
					t.sub = 0; t.pred = 0; t.obj = 0;
					resmap[next->gnode] = t;
					bfsq.push_back(next->gnode);
				}
			}
		}
	}
}

void print_res(FILE **outfile, FILE **intfile, FILE **intfile_preds, struct node **bfsarr, bool bestm, map<struct node *, struct triple> &q_to_gr_node_map,
		map<struct node *, struct triple> &copy_of_resmap)
{
//	for(map<struct node*, struct triple>::iterator it = copy_of_resmap.begin(); it != copy_of_resmap.end(); it++) {
//		struct triple t = it->second;
//		if (t.sub == 0 || t.pred == 0 || t.obj == 0) {
//			set_neighbors_null(it->first, copy_of_resmap);
//		}
//	}

	if (!bestm) {
		for (int i = (-1)*q_to_gr_node_map.size(); i < 0; i++) {
	//		cout << i << " tp " << ((TP*)bfsarr[i]->data)->toString() << endl;
			struct triple t = q_to_gr_node_map[bfsarr[i]];
			if (((TP*)bfsarr[i]->data)->sub < 0) {
				if (t.sub <= gnum_subs) {
					fprintf(*outfile, "%u ", t.sub);
					fprintf(*intfile, "%u\n", t.sub);
				} else {
					assert(0);
				}
			}
			if (((TP*)bfsarr[i]->data)->pred < 0) {
				fprintf(*outfile, "%u ", t.pred);
				fprintf(*intfile_preds, "%u\n", t.pred);
			}
			if (((TP*)bfsarr[i]->data)->obj < 0) {
				if (t.obj <= gnum_comm_so) {
					fprintf(*outfile, "%u ", t.obj);
					fprintf(*intfile, "%u\n", t.obj);
				} else {
					unsigned int tmp = t.obj - gnum_comm_so + gnum_subs;
					fprintf(*outfile, "%u ", tmp);
					fprintf(*intfile, "%u\n", tmp);
				}
			}
		}
	} else {
		for(map<struct node*, struct triple>::iterator it = q_to_gr_node_map.begin();
				it != q_to_gr_node_map.end(); it++) {
			copy_of_resmap[it->first] = it->second;
		}

		set_neighbors_null(&bfsarr[(-1)*copy_of_resmap.size()], copy_of_resmap);

//		for (int i = (-1)*copy_of_resmap.size(); i < 0; i++) {
//			struct triple t = copy_of_resmap[bfsarr[i]];
//			if (t.sub == 0 || t.pred == 0 || t.obj == 0) {
//				//TODO: rewrite this logic to avoid unnecessary function calls
//				set_neighbors_null(bfsarr[i], copy_of_resmap);
//			}
//		}

		for (int i = (-1)*copy_of_resmap.size(); i < 0; i++) {
	//		cout << i << " tp " << ((TP*)bfsarr[i]->data)->toString() << endl;
			struct triple t = copy_of_resmap[bfsarr[i]];
			if (((TP*)bfsarr[i]->data)->sub < 0) {
				if (t.sub <= gnum_subs) {
					fprintf(*outfile, "%u ", t.sub);
					if (t.sub > 0) fprintf(*intfile, "%u\n", t.sub);
				} else assert(0);
			}
			if (((TP*)bfsarr[i]->data)->pred < 0) {
				fprintf(*outfile, "%u ", t.pred);
				if (t.pred > 0) fprintf(*intfile_preds, "%u\n", t.pred);
			}
			if (((TP*)bfsarr[i]->data)->obj < 0) {
				if (t.obj <= gnum_comm_so) {
					fprintf(*outfile, "%u ", t.obj);
					if (t.obj > 0) fprintf(*intfile, "%u\n", t.obj);
				} else {
					unsigned int tmp = t.obj - gnum_comm_so + gnum_subs;
					fprintf(*outfile, "%u ", tmp);
					fprintf(*intfile, "%u\n", tmp);
				}
			}
		}
	}

//	for(map<struct node*, struct triple>::iterator it = copy_of_resmap.begin(); it != copy_of_resmap.end(); it++) {
//		if (((TP*)(it->first)->data)->sub < 0) {
//			fprintf(outfile, "%u ", (it->second).sub);
//		}
//		if (((TP*)it->first->data)->pred < 0) {
//			fprintf(outfile, "%u ", (it->second).pred);
//		}
//		if (((TP*)it->first->data)->obj < 0) {
//			fprintf(outfile, "%u ", (it->second).obj);
//		}
//	}
//	fprintf(*outfile, "%s", num_null_pads);
//	for (int i = 0; i < num_null_pads; i++) {
//		fprintf(outfile, "0 ");
//	}
	fprintf(*outfile, "\n");
}

/////////////////////////////////////////////////////

bool check_bit(struct row *r, unsigned int bit)
{
	unsigned char *resptr = r->data + BM_ROW_SIZE;
	bool flag = resptr[0] & 0x01;
	unsigned char format = resptr[0] & 0x02;
//	unsigned int ressize = 0;
	unsigned int sz = 0;
	memcpy(&sz, r->data, BM_ROW_SIZE);
	unsigned int total_cnt = sz - 1;
	unsigned int cnt = 0;
	unsigned int bitcnt = 0, tmpcnt = 0;

	while (cnt < total_cnt) {
		memcpy(&tmpcnt, &resptr[cnt*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);

		if (format == 0x02) {
			if (bit == tmpcnt)
				return true;
		} else {
			if (flag && (bitcnt+1 <= bit && bit <= bitcnt+tmpcnt)) {
				return true;
			} else if (!flag && (bitcnt+1 <= bit && bit <= bitcnt+tmpcnt)) {
				return false;
			}
		}

		bitcnt += tmpcnt;
		cnt++;
		flag = !flag;
	}
	return false;
}

FILE * print_triples_to_file(unsigned int rownum, unsigned int bmnum, unsigned int column, int bmdim, FILE *outfile)
{
	switch(bmdim) {
		case(SPO_BITMAT):
			if (column > gnum_comm_so) {
				column = column - gnum_comm_so + gnum_subs;
			}
			fprintf(outfile, "%d:%d:%d\n", rownum, bmnum, column);
			break;
		case(OPS_BITMAT):
			if (rownum > gnum_comm_so) {
				rownum = rownum - gnum_comm_so + gnum_subs;
			}

			fprintf(outfile, "%d:%d:%d\n", column, bmnum, rownum);
			break;
		case(PSO_BITMAT):
			if (column > gnum_comm_so) {
				column = column - gnum_comm_so + gnum_subs;
			}
			fprintf(outfile, "%d:%d:%d\n", bmnum, rownum, column);
			break;
		case(POS_BITMAT):
			if (bmnum > gnum_comm_so) {
				bmnum = bmnum - gnum_comm_so + gnum_subs;
			}
			fprintf(outfile, "%d:%d:%d\n", column, rownum, bmnum);
			break;
	}
	return outfile;
}

/////////////////////////////////////////////////////////
void list_enctrips_bitmat_new(BitMat *bitmat, unsigned int bmnum, vector<twople> &twoplist, FILE *outfile)
{
#if USE_MORE_BYTES
	unsigned long opos = 0;
	unsigned long tmpcnt = 0, bitcnt = 0;
#else
	unsigned int opos = 0;
	unsigned int tmpcnt = 0, bitcnt = 0;
#endif
	unsigned int total_cnt = 0;
	bool flag;
	unsigned int cnt = 0;
//	FILE *fp;
	unsigned int rowcnt;

	if (outfile == NULL)
		twoplist.clear();

	for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); it++) {
		unsigned char *data = (*it).data;
		unsigned int rownum = (*it).rowid;
		flag = data[BM_ROW_SIZE] & 0x01;
		unsigned char format = data[BM_ROW_SIZE] & 0x02;
		unsigned int sz = 0;
		memcpy(&sz, data, BM_ROW_SIZE);
		total_cnt = sz - 1;
		cnt = 0;
		bitcnt = 0;
		data += BM_ROW_SIZE;

		while (cnt < total_cnt) {
			memcpy(&tmpcnt, &data[cnt*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);

			if (format == 0x02) {
//				cout << rownum << ":" << bmnum << ":"<< tmpcnt << endl;
				if (outfile != NULL) {
					outfile = print_triples_to_file(rownum, bmnum, tmpcnt, bitmat->dimension, outfile);
//					fprintf(outfile, "%d:%d:%d\n", rownum, bmnum, tmpcnt);
				} else {
					struct twople t = {rownum, tmpcnt};
					twoplist.push_back(t);
				}
			} else {
				if (flag) {
					for (opos = bitcnt + 1; opos <= bitcnt+tmpcnt; opos++) {
						if (outfile != NULL) {
							outfile = print_triples_to_file(rownum, bmnum, opos, bitmat->dimension, outfile);
//							fprintf(outfile, "%d:%d:%d\n", rownum, bmnum, opos);
						} else {
							struct twople t = {rownum, opos};
							twoplist.push_back(t);
						}
					}
				}
			}
			bitcnt += tmpcnt;
			cnt++;
			flag = !flag;
		}
	}
}
/////////////////////////////
void list_enctrips_bitmat2(BitMat *bitmat, unsigned int bmnum, vector<twople> &twoplist, FILE *outfile)
{
#if USE_MORE_BYTES
	unsigned long opos = 0;
	unsigned long tmpcnt = 0, bitcnt = 0;
#else
	unsigned int opos = 0;
	unsigned int tmpcnt = 0, bitcnt = 0;
#endif
	unsigned int total_cnt, cnt;
	bool flag;
//	FILE *fp;
	unsigned int rowcnt;

	if (outfile == NULL)
		twoplist.clear();

	//Depending on the dimension of the bitmat
	//it's either a SPO or OPS bitmat
	//TODO: this will change for large datasets
	//if you are maining a hash_map instead of char**

	for (std::vector<struct row>::iterator it = bitmat->vbm.begin(); it != bitmat->vbm.end(); it++) {

		unsigned char *data = (*it).data;
		unsigned int rownum = (*it).rowid;
		flag = data[BM_ROW_SIZE] & 0x01;
		unsigned char format = data[BM_ROW_SIZE] & 0x02;
		unsigned int sz = 0;
		memcpy(&sz, data, BM_ROW_SIZE);
		total_cnt = sz - 1;
		cnt = 0;
		bitcnt = 0;
		data += BM_ROW_SIZE;

		while (cnt < total_cnt) {
			memcpy(&tmpcnt, &data[cnt*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
			if (format == 0x02) {
				if (outfile != NULL) {
					outfile = print_triples_to_file(rownum, bmnum, tmpcnt, bitmat->dimension, outfile);
//					fprintf(outfile, "%u:%u:%u\n", rownum, bmnum, tmpcnt);
				} else {
					struct twople t = {rownum, tmpcnt};
					twoplist.push_back(t);
				}
			} else {
				if (flag) {
					for (opos = bitcnt+1; opos <= bitcnt+tmpcnt; opos++) {
						if (outfile != NULL) {
							outfile = print_triples_to_file(rownum, bmnum, opos, bitmat->dimension, outfile);
//							fprintf(outfile, "%u:%u:%u\n", rownum, bmnum, opos);
						} else {
							struct twople t = {rownum, opos};
							twoplist.push_back(t);
						}
					}
				}
			}
			bitcnt += tmpcnt;
			cnt++;
			flag = !flag;
		}

	}

}

unsigned long get_size_of_bitmat(int dimension, unsigned int node)
{
	char dumpfile[1024];
	off_t off1 = 0, off2 = 0;
	switch (dimension) {
	case SPO_BITMAT:
		sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_SPO")].c_str());
		off1 = get_offset(dumpfile, node);
		if (node == gnum_preds) {
			off2 = vectfd[0].second;
		} else {
			off2 = get_offset(dumpfile, node+1);
		}
		return (off2 - off1);
		break;
	case OPS_BITMAT:
		sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_OPS")].c_str());
		off1 = get_offset(dumpfile, node);
		if (node == gnum_preds) {
			off2 = vectfd[2].second;
		} else {
			off2 = get_offset(dumpfile, node+1);
		}
		return (off2 - off1);
		break;
	case PSO_BITMAT:
		sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
		off1 = get_offset(dumpfile, node);
		if (node == gnum_subs) {
			off2 = vectfd[4].second;
		} else {
			off2 = get_offset(dumpfile, node+1);
		}
		return (off2 - off1);
		break;
	case POS_BITMAT:
		sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
		off1 = get_offset(dumpfile, node);
		if (node == gnum_objs) {
			off2 = vectfd[6].second;
		} else {
			off2 = get_offset(dumpfile, node+1);
		}
		return (off2 - off1);
		break;
	default:
		assert(0);
		break;
	}
	return 0;
}

void list_all_data(int dimension, unsigned int node, char *filename)
{
	BitMat bitmat;
	char dumpfile[1024];
	FILE *file = NULL;
	if (filename != NULL) {
		file = fopen(filename, "a");
		assert(file != NULL);
		setvbuf(file, NULL, _IOFBF, 0x8000000);

	}
	switch (dimension) {
		case SPO_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_SPO")].c_str());
			init_bitmat(&bitmat, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
			cout << "Listing triples from SPO bitmats" << endl;
			if (node != 0) {
//				unsigned long offset = get_offset(dumpfile, node);
//				if (offset == 0 && node > 1) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);
//				}
//				if (offset == 0 && node > 1 && offset == get_offset(dumpfile, node+1)) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);
//				}
				load_from_dump_file(dumpfile, node, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				vector<twople> twoples;
				list_enctrips_bitmat_new(&bitmat, node, twoples, file);
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", (*it).sub, node, (*it).obj);
				}
				twoples.clear();
				bitmat.freebm();
//				clear_rows(&bitmat, true, true, true);
				if (file != NULL) {
					fclose(file);
				}

				return;
			}
//			{
//				unsigned long offset1 = get_offset(dumpfile, 1);
//				unsigned long offset2 = get_offset(dumpfile, 2);
//				cout << "Size of bitmat 1 " << (offset2 - offset1) << endl;
//				offset1 = get_offset(dumpfile, 3);
//				cout << "Size of bitmat 2 " << (offset1 - offset2) << endl;
//				offset2 = get_offset(dumpfile, 4);
//				cout << "Size of bitmat 3 " << (offset2 - offset1) << endl;
//				offset1 = get_offset(dumpfile, 8);
//				cout << "Size of bitmat 4 " << (offset1 - offset2) << endl;
//			}
			for (unsigned int i = 1; i <= gnum_preds; i++) {
//				unsigned long offset = get_offset(dumpfile, i);
				//This offset is just a placeholder
//				if (offset == 0 && i > 1) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}
//				if (offset == 0 && i == 1 && offset == get_offset(dumpfile, i+1)) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}

				cout << "Loading bitmat " << i << endl;
				load_from_dump_file(dumpfile, i, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				cout << "Done loading bitmat " << i << endl;
				vector<twople> twoples;
				cout << "Listing triples in bitmat " << i << endl;
				list_enctrips_bitmat_new(&bitmat, i, twoples, file);
				cout << "Done listing triples in bitmat " << i << endl;
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", (*it).sub, i, (*it).obj);
//					cerr << (*it).sub << ":" << i << ":" << (*it).obj << endl;
				}
				twoples.clear();
				cout << "Clearing rows for bitmat " << i << endl;
				bitmat.freebm();
//				clear_rows(&bitmat, true, true, true);
			}
			break;

		case OPS_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_OPS")].c_str());
			init_bitmat(&bitmat, gnum_objs, gnum_preds, gnum_subs, gnum_comm_so, OPS_BITMAT);
			cout << "Listing triples from OPS bitmats" << endl;

			if (node != 0) {
//				unsigned long offset = get_offset(dumpfile, node);
//				cout << "Offset for bmnum " << node << " is "<< offset << endl;
				//This offset is just a placeholder
//				if (offset == 0 && node > 1) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);
//				}
//				if (offset == 0 && node == 1 && offset == get_offset(dumpfile, node+1)) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);
//				}

				load_from_dump_file(dumpfile, node, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				vector<twople> twoples;
				list_enctrips_bitmat_new(&bitmat, node, twoples, file);
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", (*it).obj, node, (*it).sub);
//					cerr << (*it).obj << ":" << node << ":" << (*it).sub << endl;
				}
				twoples.clear();
				bitmat.freebm();
//				clear_rows(&bitmat, true, true, true);
				if (file != NULL) {
					fclose(file);
				}

				return;

			}

			for (unsigned int i = 1; i <= gnum_preds; i++) {
//				unsigned long offset = get_offset(dumpfile, i);
				//This offset is just a placeholder
//				if (offset == 0 && i > 1) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}
//				if (offset == 0 && i == 1 && offset == get_offset(dumpfile, i+1)) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}
//
				load_from_dump_file(dumpfile, i, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				vector<twople> twoples;
				list_enctrips_bitmat_new(&bitmat, i, twoples, file);
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", (*it).obj, i, (*it).sub);
//					cerr << (*it).obj << ":" << i << ":" << (*it).sub << endl;
				}
				twoples.clear();
				bitmat.freebm();
//				clear_rows(&bitmat, true, true, true);
			}
			break;

		case PSO_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
			init_bitmat(&bitmat, gnum_preds, gnum_subs, gnum_objs, gnum_comm_so, PSO_BITMAT);
			cout << "Listing triples from PSO bitmats" << endl;

			if (node != 0) {
//				unsigned long offset = get_offset(dumpfile, node);
				//This offset is just a placeholder
//				if (offset == 0 && node > 1) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);;
//				}
//
//				if (offset == 0 && node == 1 && offset == get_offset(dumpfile, node+1)) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);
//				}
				load_from_dump_file(dumpfile, node, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				vector<twople> twoples;
				list_enctrips_bitmat_new(&bitmat, node, twoples, file);
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", node, (*it).sub, (*it).obj);
//					cerr << node << ":" << (*it).sub << ":" << (*it).obj << endl;
				}
				twoples.clear();
//				clear_rows(&bitmat, true, true, true);
				bitmat.freebm();

				if (file != NULL) {
					fclose(file);
				}
				return;

			}

			for (unsigned int i = 1; i <= gnum_subs; i++) {
//				unsigned long offset = get_offset(dumpfile, i);
				//This offset is just a placeholder
//				if (offset == 0 && i > 1) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}
//
//				if (offset == 0 && i == 1 && offset == get_offset(dumpfile, i+1)) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}

				load_from_dump_file(dumpfile, i, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				vector<twople> twoples;
				list_enctrips_bitmat_new(&bitmat, i, twoples, file);
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", i, (*it).sub, (*it).obj);
//					cerr << i << ":" << (*it).sub << ":" << (*it).obj << endl;
				}
				twoples.clear();
				bitmat.freebm();
//				clear_rows(&bitmat, true, true, true);
			}
			break;

		case POS_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
			init_bitmat(&bitmat, gnum_preds, gnum_objs, gnum_subs, gnum_comm_so, POS_BITMAT);
			cout << "Listing triples from POS bitmats" << endl;
			if (node != 0) {
//				unsigned long offset = get_offset(dumpfile, node);
				//This offset is just a placeholder
//				if (offset == 0 && node > 1) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);;
//				}
//				if (offset == 0 && node == 1 && offset == get_offset(dumpfile, node+1)) {
//					cout << "Non-existing predicate " << node << endl;
//					assert(0);
//				}

				load_from_dump_file(dumpfile, node, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				vector<twople> twoples;
				list_enctrips_bitmat_new(&bitmat, node, twoples, file);
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", (*it).obj, (*it).sub, node);
//					cerr << (*it).obj << ":" << (*it).sub << ":" << node << endl;
				}
				twoples.clear();
//				clear_rows(&bitmat, true, true, false);
				bitmat.freebm();
				if (file != NULL) {
					fclose(file);
				}

				return;
			}

			for (unsigned int i = 1; i <= gnum_objs; i++) {
//				unsigned long offset = get_offset(dumpfile, i);
				//This offset is just a placeholder
//				if (offset == 0 && i > 1) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}
//
//				if (offset == 0 && i == 1 && offset == get_offset(dumpfile, i+1)) {
//					cout << "Non-existing predicate " << i << endl;
//					continue;
//				}

				load_from_dump_file(dumpfile, i, &bitmat, true, true, NULL, 0, 0, NULL, 0, true);
				vector<twople> twoples;
				list_enctrips_bitmat_new(&bitmat, i, twoples, file);
				for (vector<twople>::iterator it = twoples.begin(); it != twoples.end(); it++) {
					fprintf(file, "%u:%u:%u\n", (*it).obj, (*it).sub, i);
//					cerr << (*it).obj << ":" << (*it).sub << ":" << i << endl;
				}
				twoples.clear();
				clear_rows(&bitmat, true, true, false);
			}
			break;
		default:
			break;

	}
	if (file != NULL) {
		fclose(file);
	}

}

////////////////////////////////////////////////////////////

//unsigned long list_enctriples_in_bitmat(unsigned char **bitmat, unsigned int dimension, unsigned int num_triples, char *triplefile)
//{
//#if USE_MORE_BYTES
//	unsigned long ppos_start_bit, ppos_end_bit, opos = 0;
//	unsigned long tmpcnt = 0, bitcnt = 0, gap;
//	unsigned long triplecnt = 0;
//#else
//	unsigned int ppos_start_bit, ppos_end_bit, opos = 0;
//	unsigned int tmpcnt = 0, bitcnt = 0, gap;
//	unsigned int triplecnt = 0;
//#endif
//	unsigned long pred_cnt;
//	unsigned int i, total_cnt;
//	bool flag, inside;
//	unsigned int cnt;
//	FILE *fp;
//	unsigned int ressize = 0, rowcnt, chunk_bytes;
//
////	if (num_triples == 0) {
////		return;
////	}
//
//	fp = fopen(triplefile, "w");
//	if (fp == NULL) {
//		printf("****ERROR opening triplefile %s\n", triplefile);
//		assert(0);
//	}
//
//	//Depending on the dimension of the bitmat
//	//it's either a SPO or OPS bitmat
//	if (dimension == SUB_DIMENSION) {
//		rowcnt = gnum_subs;
//		chunk_bytes = gobject_bytes;
//	} else if (dimension == OBJ_DIMENSION) {
//		rowcnt = gnum_objs;
//		chunk_bytes = gsubject_bytes;
//	}
//	for (i = 0; i < rowcnt; i++) {
//		if (bitmat[i] == NULL)
//			continue;
//		unsigned char *resptr = bitmat[i] + ROW_SIZE_BYTES;
//		flag = resptr[0] & 0x01;
//		ressize = 0;
//		memcpy(&ressize, bitmat[i], ROW_SIZE_BYTES);
//		total_cnt = (ressize - 1)/GAP_SIZE_BYTES;
//		cnt = 0;
//		bitcnt = 0;
//		pred_cnt = 0;
//		ppos_start_bit = pred_cnt * (chunk_bytes << 3) + 1;
//		ppos_end_bit = ppos_start_bit + (chunk_bytes << 3) - 1;
//
//		memcpy(&tmpcnt, &resptr[cnt*GAP_SIZE_BYTES + 1], GAP_SIZE_BYTES);
//		bitcnt += tmpcnt;
//		inside = false;
//
//		while (cnt < total_cnt) {
//			if (bitcnt >= ppos_start_bit && bitcnt <= ppos_end_bit) {
//				inside = true;
//				if (flag) {
//					gap = (bitcnt - ppos_start_bit + 1) < tmpcnt ?
//								(bitcnt - ppos_start_bit + 1)	:  tmpcnt;
//					opos = bitcnt - gap + 1 - (pred_cnt * (chunk_bytes << 3));
//					while (gap > 0) {
////						triples[triplecnt] = (unsigned char *) malloc (TRIPLE_STR_SPACE);
////						memset(triples[triplecnt], 0, TRIPLE_STR_SPACE);
//						if (dimension == SUB_DIMENSION)
//							fprintf(fp, "%u:%u:%u\n", i+1, pred_cnt+1, opos);
//						else if (dimension == OBJ_DIMENSION)
//							fprintf(fp, "%u:%u:%u\n", opos, pred_cnt+1, i+1);
//
//						triplecnt++;
//						gap--;
//						opos++;
//					}
//				}
//				if (bitcnt == ppos_end_bit) {
//					inside = false;
//					pred_cnt++;
//					ppos_start_bit = pred_cnt * (chunk_bytes << 3) + 1;
//					ppos_end_bit = ppos_start_bit + (chunk_bytes << 3) - 1;
//				}
//				cnt++;
//				if (cnt >= total_cnt)
//					break;
//				flag = !flag;
//				memcpy(&tmpcnt, &resptr[cnt*GAP_SIZE_BYTES+1], GAP_SIZE_BYTES);
//				bitcnt += tmpcnt;
//			} else if (inside) {
//				inside = false;
//				if (flag) {
//					gap = ppos_end_bit - (bitcnt - tmpcnt);
//					opos = ppos_end_bit - gap + 1 - (pred_cnt * (chunk_bytes << 3));
//
//					while (gap > 0) {
////						triples[triplecnt] = (unsigned char *) malloc (TRIPLE_STR_SPACE);
////						memset(triples[triplecnt], 0, TRIPLE_STR_SPACE);
//						if (dimension == SUB_DIMENSION)
//							fprintf(fp, "%u:%u:%u\n", i+1, pred_cnt+1, opos);
//						else if (dimension == OBJ_DIMENSION)
//							fprintf(fp, "%u:%u:%u\n", opos, pred_cnt+1, i+1);
//						triplecnt++;
//						gap--;
//						opos++;
//					}
//				}
//				pred_cnt++;
//				ppos_start_bit = pred_cnt * (chunk_bytes << 3) + 1;
//				ppos_end_bit = ppos_start_bit + (chunk_bytes << 3) - 1;
//
//			} else {
//				pred_cnt++;
//				ppos_start_bit = pred_cnt * (chunk_bytes << 3) + 1;
//				ppos_end_bit = ppos_start_bit + (chunk_bytes << 3) - 1;
//			}
//
//		}
//
//	}
//
//	fclose(fp);
//
//	return triplecnt;
//
////	printf("Total triples %u\n", num_triples);
//
//}
////////////////////////////////////////////////////////////

void dump_out_data(FILE *fdump_fp, BitMat *bitmat, char *fname_tmp)
{
//	int fd = 0;
	unsigned int i = 0;
//	unsigned int size = 0;

	if (bitmat->num_triples != 0) {
//		write(fd, &bitmat->num_triples, sizeof(unsigned int));
//		cout << "Num triples - " << bitmat->num_triples << endl;;
		fwrite(&bitmat->num_triples, sizeof(unsigned long), 1, fdump_fp);
		gtotal_size += sizeof(unsigned int);
	}
	if (bitmat->rowfold != NULL && bitmat->colfold != NULL) {
		if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
//			write(fd, comp_rowfold, comp_rowfold_size + ROW_SIZE_BYTES);
			fwrite(comp_rowfold, comp_rowfold_size + ROW_SIZE_BYTES, 1, fdump_fp);
			gtotal_size += comp_rowfold_size + ROW_SIZE_BYTES;
		} else {
			fwrite(bitmat->rowfold, bitmat->row_bytes, 1, fdump_fp);
			gtotal_size += bitmat->row_bytes;
		}

		if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
			if (comp_folded_arr) {
				//compress the array
				unsigned int comp_colfold_size = dgap_compress_new(bitmat->colfold, bitmat->column_bytes, comp_colfold);
				fwrite(&comp_colfold_size, ROW_SIZE_BYTES, 1, fdump_fp);
				fwrite(comp_colfold, comp_colfold_size, 1, fdump_fp);
				gtotal_size += ROW_SIZE_BYTES + comp_colfold_size;
			} else {
				fwrite(bitmat->colfold, bitmat->column_bytes, 1, fdump_fp);
				gtotal_size += bitmat->column_bytes;
			}
		}

	}

	if (0 == bitmat->bm.size()) {
		int tmpdump = open(fname_tmp, O_RDONLY);
		if (tmpdump < 0)
			assert(0);

		unsigned char readbuf[0x1000000];
		ssize_t rbytes = 1;

		while (rbytes > 0) {
			rbytes = read(tmpdump, readbuf, 0x1000000);
			if (rbytes > 0) {
				fwrite(readbuf, rbytes, 1, fdump_fp);
				gtotal_size += rbytes;
			}
		}
		close(tmpdump);

	} else {
		for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); it++) {
			unsigned int sz = 0;
			memcpy(&sz, (*it).data, BM_ROW_SIZE);
			fwrite((*it).data, TOTAL_ROWSIZE(sz) + BM_ROW_SIZE, 1, fdump_fp);
			gtotal_size += (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
		}
	}

//	printf("dump_out_data: Total bytes to file %u #triples %u\n", gtotal_size, bitmat->num_triples);

}

///////////////////////////////////////////////

void load_data_vertically(char *file, vector<struct twople> &twoplelist, BitMat *bitmat, char *fname_dump,
						bool ondisk, bool invert, bool listload, char *tmpfile)
{
	FILE *fp = NULL, *fp_table = NULL;
	FILE *fdump_fp = NULL;
	FILE *tmpdump = NULL;
	unsigned int spos = 0, ppos = 0, opos = 0;
	unsigned int sprev = 0, pprev = 0;
	unsigned char l, m, n;
	bool start = true;
//	unsigned char **table;
	//TODO: change later
//	unsigned char table[bitmat->num_preds][table_col_bytes];
//	unsigned char table[20000000][table_col_bytes];
//	int fdump_fd = -1;

//	cout << "Inside load_data_vertically" << endl;

	if (file != NULL) {
		fp = fopen(file, "r");
		if (fp == NULL) {
			printf("Error opening file %s\n", file);
			return;
		}
	}

//	cout << "load_data_vertically : HERE1" << endl;
	if (grow != NULL) {
		free(grow); grow = NULL;
	}
	grow = (unsigned char *) malloc (GAP_SIZE_BYTES * (bitmat->num_columns / 2));
	memset(grow, 0, GAP_SIZE_BYTES * (bitmat->num_columns / 2));
	grow_size = GAP_SIZE_BYTES * (bitmat->num_columns / 2);
//	cout << "load_data_vertically : HERE2" << endl;

	if (ondisk) {
		//TODO: while loading in chunks no need to alloc num_preds
		char tablefile[1024];
		sprintf(tablefile, "%s_table", fname_dump);

		fp_table = fopen(tablefile, "wb");
		if (fp_table == NULL) {
			cout << "*** ERROR: Could not open tablefile " << tablefile << endl;
			assert(0);
		}

		if (buf_table != NULL) {
			free(buf_table); buf_table = NULL;
		}

		buf_table = (char *) malloc (0xf000000 * sizeof(char));

		if (setvbuf(fp_table, buf_table, _IOFBF, 0xf000000) != 0)
			assert(0);

//		table = (unsigned char **) malloc (bitmat->num_preds * sizeof(unsigned char *));
//
//		for (unsigned int i = 0; i < bitmat->num_preds; i++) {
//			table[i] = (unsigned char *) malloc (table_col_bytes * sizeof (unsigned char));
//			memset(table[i], 0, table_col_bytes * sizeof (unsigned char));
//		}
//		memset(table, 0, bitmat->num_preds * table_col_bytes);
	}

//	cout << "load_data_vertically : HERE3" << endl;
	unsigned int total_triples = 0;
	prev_bitset = 0;
	prev_rowbit = 1, comp_rowfold_size = 0;
	gtotal_size = 0;

	bitmat->bm.clear();

//	cout << "load_data_vertically : HERE4" << endl;
	if (ondisk) {

		if (tmpfile != NULL) {
			tmpdump = fopen(tmpfile, "wb");

			if (buf_tmp != NULL) {
				free(buf_tmp); buf_tmp = NULL;
			}

			buf_tmp = (char *) malloc (0xf000000 * sizeof(char));

			if (setvbuf(tmpdump, buf_tmp, _IOFBF, 0xf000000) != 0)
				assert(0);

			if (tmpdump == NULL) {
				cout << "Cannot open " << config[string("TMP_STORAGE")] << endl;
				assert(0);
			}
		}

		if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
			comp_rowfold = (unsigned char *) malloc (GAP_SIZE_BYTES * bitmat->num_rows);
			memset(comp_rowfold, 0, GAP_SIZE_BYTES * bitmat->num_rows);
			if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
				comp_colfold = (unsigned char *) malloc (GAP_SIZE_BYTES * bitmat->num_columns);
				memset(comp_colfold, 0, GAP_SIZE_BYTES * bitmat->num_columns);
			}
		} else if (bitmat->rowfold == NULL) {
			//TODO:
			bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof (unsigned char));
			memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof (unsigned char));
		}

		if (bitmat->colfold == NULL) {
			bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof (unsigned char));
			memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof (unsigned char));
		}

		fdump_fp = fopen(fname_dump, "wb");
		if (fdump_fp == NULL) {
			cout << "Could not open " << fname_dump << endl;
			assert(0);
		}

		if (buf_dump != NULL) {
			free(buf_dump); buf_dump = NULL;
		}

		buf_dump = (char *) malloc (sizeof(char) * 0xf000000);

		if (setvbuf(fdump_fp, buf_dump, _IOFBF, 0xf000000) != 0)
			assert(0);
	}

//	unsigned int pcnt = 0;
	unsigned int swap = 0;

	if (file == NULL) {
		//Load from the triplelist vector
//		cout << "load_data_vertically : HERE5" << endl;
		bitmat->num_triples = twoplelist.size();
		for (vector<struct twople>::iterator it = twoplelist.begin(); it != twoplelist.end(); it++) {
			spos = (*it).sub;
			opos = (*it).obj;
			if (invert) {
				swap = spos;
				spos = opos;
				opos = swap;
			}
			if (sprev != spos) {
				if (start) {
					map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, false, true, ondisk,
							listload, NULL, false);
					start = false;
				} else {
					map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, true, true, ondisk,
							listload, NULL, false);
				}
				sprev = spos;
			} else {
				map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, false, false, ondisk,
						listload, NULL, false);
			}
		}

		//For the last row entries
		map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, true, false, ondisk, listload, NULL, false);

		if (grow != NULL) {
			free(grow);
			grow = NULL;
		}

		return;

	}

	unsigned int tcount = 0;
	while (!feof(fp)) {
//		cout << "loading line" << endl;
		char line[50];
		char s[10], p[10], o[10];

		memset(line, 0, 50);
		if(fgets(line, sizeof(line), fp) != NULL) {
			char *c = NULL, *c2=NULL;
			c = strchr(line, ':');
			*c = '\0';
			strcpy(s, line);
			c2 = strchr(c+1, ':');
			*c2 = '\0';
			strcpy(p, c+1);
			c = strchr(c2+1, '\n');
			*c = '\0';
			strcpy(o, c2+1);
			tcount++;
			spos = strtoul(s, NULL, 10); ppos = strtoul(p, NULL, 10); opos = strtoul(o, NULL, 10);
			if (invert) {
				swap = spos;
				spos = opos;
				opos = swap;
			}
			if (pprev != ppos || sprev != spos) {
				if (start) {
					tmpdump = map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, false, true,
								ondisk, listload, tmpdump, false);
					start = false;
					if (ppos != 1) {
						unsigned long tmpval = 0;
						for (unsigned int k = 1; k < ppos; k++) {
							fwrite(&tmpval, table_col_bytes, 1, fp_table);
						}
					}
					fwrite(&gtotal_size, table_col_bytes, 1, fp_table);

				} else {
					if (pprev != ppos) {
						if (ondisk) {
							tmpdump = map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, true,
										false, ondisk, listload, tmpdump, true);
//							bitmat->num_triples = count_triples_in_bitmat(bitmat, bitmat->dimension);
							bitmat->num_triples = tcount-1;
							//For new BitMat triple count is reset
							tcount = 1;
//							cout << fname_dump << " " << bitmat->num_triples << endl;
							total_triples += bitmat->num_triples;

							//TODO: remove this later
//							memcpy(table[pprev-1], &gtotal_size, table_col_bytes);
//							memcpy(table[pcnt], &gtotal_size, table_col_bytes);
//							pcnt++;

							if (comp_rowfold_size != 0) {
								//COMMENT: No need to add last zeros
//								if (prev_rowbit != (bitmat->row_bytes * 8)) {
//									//You will have to append later_0s
//									unsigned int later_0 = (bitmat->row_bytes * 8) - prev_rowbit;
//									memcpy(&comp_rowfold[ROW_SIZE_BYTES + comp_rowfold_size], &later_0, GAP_SIZE_BYTES);
//									comp_rowfold_size += GAP_SIZE_BYTES;
//								}
								memcpy(comp_rowfold, &comp_rowfold_size, ROW_SIZE_BYTES);
							}
							if (bitmat->dimension == SPO_BITMAT || bitmat->dimension == OPS_BITMAT) {
								cout << "load_data_vertically: Dumping data for BM " << pprev << endl;
							}

							if (tmpdump) {
								fclose(tmpdump);
								tmpdump = NULL;
							}

							dump_out_data(fdump_fp, bitmat, (char *)config[string("TMP_STORAGE")].c_str());

							if (ppos - 1 != pprev) {
								unsigned long tmpval = 0;
								for (unsigned int k = pprev+1; k < ppos; k++) {
									fwrite(&tmpval, table_col_bytes, 1, fp_table);
								}
							}

							fwrite(&gtotal_size, table_col_bytes, 1, fp_table);

							if (tmpfile != NULL) {
								tmpdump = fopen(tmpfile, "wb");
								if (setvbuf(tmpdump, buf_tmp, _IOFBF, 0xf000000) != 0)
									assert(0);

								if (tmpdump == NULL) {
									cout << "Cannot open " << config[string("TMP_STORAGE")] << endl;
									assert(0);
								}
							}
							//Empty the earlier bitmat
							clear_rows(bitmat, true, true, true);
							comp_rowfold_size = 0;
							prev_rowbit = 1;

							tmpdump = map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, false,
												true, ondisk, listload, tmpdump, false);
						}
//						cout << "Cleared rows after pred " <<  pprev << endl;

					} else {
						tmpdump = map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, true, true,
											ondisk, listload, tmpdump, false);
					}
				}
				sprev = spos;
				pprev = ppos;
			} else {
				tmpdump = map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, false, false, ondisk,
											listload, tmpdump, false);
			}
		}
//		cout << "Done line" << endl;
	}
	//for the last entries
	if (spos != 0 && ppos != 0 && opos != 0) {
		tmpdump = map_to_row_wo_dgap_vertical(bitmat, spos, opos, sprev, true, false, ondisk, listload, tmpdump, false);

		if (ondisk) {
//			bitmat->num_triples = count_triples_in_bitmat(bitmat, bitmat->dimension);
			bitmat->num_triples = tcount;
			tcount = 0;
//			cout << fname_dump << " " << bitmat->num_triples << endl;
			total_triples += bitmat->num_triples;
//			memcpy(table[pprev-1], &gtotal_size, table_col_bytes);
//			memcpy(table[pcnt], &gtotal_size, table_col_bytes);
//			pcnt++;
			//TODO:
			if (comp_rowfold_size != 0) {
//				if (prev_rowbit != (bitmat->row_bytes * 8)) {
//					//You will have to append later_0s
//					unsigned int later_0 = (bitmat->row_bytes * 8) - prev_rowbit;
//					memcpy(&comp_rowfold[ROW_SIZE_BYTES + comp_rowfold_size], &later_0, GAP_SIZE_BYTES);
//					comp_rowfold_size += GAP_SIZE_BYTES;
//				}
				memcpy(comp_rowfold, &comp_rowfold_size, ROW_SIZE_BYTES);
			}
			cout << "load_data_vertically: Dumping data for the last BM " << ppos << " ";
			if (tmpdump) {
				fclose(tmpdump);
				tmpdump = NULL;
			}
			dump_out_data(fdump_fp, bitmat, (char *)config[string("TMP_STORAGE")].c_str());

//			if ( (pcnt % 1048576) == 0)
//				cout << "**** Done with BM num " << pcnt << endl;

			//Empty the earlier bitmat
			clear_rows(bitmat, true, true, true);
			comp_rowfold_size = 0;
			prev_rowbit = 1;
			fclose(fdump_fp);
			fclose(fp_table);
		}

	}

	fclose(fp);

	if (ondisk) {

		cout << "***Total triples in all bitmats " << total_triples <<  endl;
//		cout << "*** Now writing out the table" << endl;
//
//		char tablefile[1024];
//		sprintf(tablefile, "%s_table", fname_dump);
//
//		FILE *fp = fopen(tablefile, "wb");
//		if (fp == NULL) {
//			cout << "*** ERROR: Could not open tablefile " << tablefile << endl;
//			exit (-1);
//		}
//
//		for (unsigned int i = 0; i < bitmat->num_preds; i++) {
//			fwrite(table[i], table_col_bytes, 1, fp);
//			free(table[i]);
//		}

		if (comp_rowfold != NULL)
			free(comp_rowfold);
		comp_rowfold = NULL;
		comp_rowfold_size = 0;
//		fclose(fp);
//		free(table);
		if (bitmat->rowfold != NULL) {
			free(bitmat->rowfold);
			bitmat->rowfold = NULL;
		}
		if (bitmat->colfold != NULL) {
			free(bitmat->colfold);
			bitmat->colfold = NULL;
		}
	}

	if (grow != NULL) {
		free(grow);
		grow = NULL;
	}
}

unsigned int last_set_bit(unsigned char *in, unsigned int size)
{
	if (size == 0)
		return 0;

	unsigned int last_set_bit = size * 8;

//	for (unsigned int i = 0; i < size; i++) {
//		if (in[i] == 0xff) {
//			last_set_bit = (i+1)*8;
//		} else if (in[i] > 0x00) {
//			for (unsigned short j = 0; j < 8; j++) {
//				if (in[i] & (0x80 >> j)) {
//					last_set_bit = (i*8)+j+1;
//				}
//			}
//		}
//	}

	for (unsigned int i = size-1; i >= 0; i--) {
		if (in[i] == 0x00) {
			continue;
		} else if (in[i] > 0x00) {
			for (unsigned short j = 0; j < 8; j++) {
				if (in[i] & (0x01 << j)) {
					last_set_bit = (8 * (i+1)) - j;
					return last_set_bit;
				}
			}
		}
	}

	return 0;

}

unsigned int load_from_dump_file2(char *fname_dump, unsigned int bmnum, BitMat *bitmat, bool readtcnt,
			bool readarray, unsigned char *maskarr, unsigned int maskarr_size, int maskarr_dim,
			char *fpos, int fd, bool loadvbm)
{
//	unsigned int size = 0;

	assert((readtcnt && readarray) || (!readtcnt && !readarray && ((fd == 0) ^ (NULL == fpos))));

#if MMAPFILES
	if (fpos == NULL) {
		fpos = mmap_table[string(fname_dump)];
		fpos += get_offset(fname_dump, bmnum);
	}
#else
	if (fd == 0) {
		fd = open(fname_dump, O_RDONLY);
		if (fd < 0) {
			printf("*** ERROR opening dump file %s\n", fname_dump);
			assert(0);
		}
		unsigned long long offset = get_offset(fname_dump, bmnum);
		///////
		//DEBUG: remove later
		////////
//		cout << "load_from_dump_file: offset is " << offset << endl;
		if (offset > 0) {
			lseek(fd, offset, SEEK_CUR);
		}
	}
#endif

	assert((fd == 0) ^ (NULL == fpos));
//	cout << "load_from_dump_file: offset " << offset << endl;

	if (readtcnt) {
		bitmat->num_triples = 0;
#if MMAPFILES
		memcpy(&bitmat->num_triples, fpos, sizeof(unsigned int));
		fpos += sizeof(unsigned int);
#else
		read(fd, &bitmat->num_triples, sizeof(unsigned int));
#endif
//		total_size += sizeof(unsigned int);
		if (bitmat->num_triples == 0) {
			cout << "load_from_dump_file: 0 results" << endl;
			return bitmat->num_triples;
		}
	}
//	cout << "Num triples1 " << bitmat->num_triples << endl;
	if (readarray) {
		if (bitmat->rowfold == NULL) {
			bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof (unsigned char));
			memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof (unsigned char));

		}
		if (bitmat->colfold == NULL) {
			bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof (unsigned char));
			memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof (unsigned char));
		}

		if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
			//first read size of the compressed array
			unsigned int comp_arr_size = 0;
#if MMAPFILES
			memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
			fpos += ROW_SIZE_BYTES;
#else
			read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//			total_size += ROW_SIZE_BYTES;
			unsigned char *comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
			memcpy(comp_arr, fpos, comp_arr_size);
			fpos += comp_arr_size;
#else
			read(fd, comp_arr, comp_arr_size);
#endif
			dgap_uncompress(comp_arr, comp_arr_size, bitmat->rowfold, bitmat->row_bytes);
//			total_size += comp_arr_size;

			free(comp_arr);

			if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
				comp_arr_size = 0;
#if MMAPFILES
				memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
				fpos += ROW_SIZE_BYTES;
#else
				read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//				total_size += ROW_SIZE_BYTES;

				comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
				memcpy(comp_arr, fpos, comp_arr_size);
				fpos += comp_arr_size;
#else
				read(fd, comp_arr, comp_arr_size);
#endif
				dgap_uncompress(comp_arr, comp_arr_size, bitmat->colfold, bitmat->column_bytes);
	//			total_size += comp_arr_size;

				free(comp_arr);
			}

		} else {
#if MMAPFILES
			memcpy(bitmat->rowfold, fpos, bitmat->row_bytes);
			fpos += bitmat->row_bytes;
			memcpy(bitmat->colfold, fpos, bitmat->column_bytes);
			fpos += bitmat->column_bytes;
#else
			read(fd, bitmat->rowfold, bitmat->row_bytes);
			read(fd, bitmat->colfold, bitmat->column_bytes);
#endif
//			total_size += bitmat->row_bytes;
//			total_size += bitmat->column_bytes;
		}
	}
//	cout << "Num bits in rowfold " << count_bits_in_row(bitmat->rowfold, bitmat->row_bytes) << endl;
//	cout << "Num bits in colfold " << count_bits_in_row(bitmat->colfold, bitmat->column_bytes) << endl;
//	print_set_bits_in_row(bitmat->colfold, bitmat->column_bytes);
	unsigned int limit_bytes = 0;

//	bool skipAtLeastOnce = false, skip = false;
	unsigned int last_bit = last_set_bit(maskarr, maskarr_size);
	unsigned int total_set_bits = count_bits_in_row(maskarr, maskarr_size);
	unsigned int rowfold_bits = count_bits_in_row(bitmat->rowfold, bitmat->row_bytes);
	/////////
	//DEBUG
	//TODO: remove later
//	cout << "--------- load_from_dump_file: maskarr_size " << maskarr_size << " maskarr_dim " << maskarr_dim << " last_setbit " << last_bit
//		<< " total_set_bits " << total_set_bits << " rowfold_bits " << rowfold_bits << endl;
	/////////////////////////////
	maskarr_size = last_bit/8 + ((last_bit%8) > 0 ? 1 : 0);

	if (maskarr != NULL && maskarr_dim == ROW) {
		limit_bytes = bitmat->row_bytes > maskarr_size ? maskarr_size : bitmat->row_bytes;
	} else {
		limit_bytes = bitmat->row_bytes;
	}

	unsigned int rownum = 1;
	//////////
	//DEBUG
//	off_t total_size = 0;

	unsigned int sz = 0;

	if (maskarr == NULL || maskarr_dim == COLUMN) {
		/*
		 * Load normally and don't do any unfolding
		 */
		for (unsigned int i = 0; i < bitmat->row_bytes; i++) {
			if (bitmat->rowfold[i] == 0x00) {
				rownum += 8;
			} else {
				for (short j = 0; j < 8; j++) {
					//Row does not exist
					if (!(bitmat->rowfold[i] & (0x80 >> j))) {
						rownum++;
						continue;
					}
	#if MMAPFILES
					memcpy(&sz, fpos, BM_ROW_SIZE);
					fpos += BM_ROW_SIZE;
	#else
					read(fd, &sz, BM_ROW_SIZE);
	#endif
					unsigned char *data = (unsigned char *) malloc (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
					memcpy(data, &sz, BM_ROW_SIZE);
					////////
					//DEBUG
//					total_size += (8+size+ROW_SIZE_BYTES);
//					cout << "total_size " << total_size << endl;
	#if MMAPFILES
					memcpy(data + BM_ROW_SIZE, fpos, TOTAL_ROWSIZE(sz));
					fpos += TOTAL_ROWSIZE(sz);
	#else
					read(fd, data + BM_ROW_SIZE, TOTAL_ROWSIZE(sz));
	#endif
					struct row r = {rownum, data};
					if (loadvbm) {
						bitmat->vbm.push_back(r);
					} else {
						bitmat->bm.push_back(r);
					}
					rownum++;
				}
			}
		}

		return bitmat->num_triples;
	}
	if (maskarr_dim == ROW) {
		if ((total_set_bits <= ((maskarr_size * 8)/SELECTIVITY_THRESHOLD)) &&
				rowfold_bits > 10*total_set_bits &&
				(bitmat->dimension == SPO_BITMAT || bitmat->dimension == OPS_BITMAT)) {
//			cout << "--------- load_from_dump_file: NOT loading the traditional way" << endl;
			char dumpfile[1024];
			BitMat bm_tmp;
			unsigned int add_row_dim = 0;
			if (bitmat->dimension == SPO_BITMAT) {
				sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
				shallow_init_bitmat(&bm_tmp, gnum_preds, gnum_subs, gnum_objs, gnum_comm_so, PSO_BITMAT);
				add_row_dim = PSO_BITMAT;

			} else if (bitmat->dimension == OPS_BITMAT) {
				sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
				shallow_init_bitmat(&bm_tmp, gnum_preds, gnum_objs, gnum_subs, gnum_comm_so, POS_BITMAT);
				add_row_dim = POS_BITMAT;
			} else assert(0);

			rownum = 1;
			for (unsigned int i = 0; i < limit_bytes; i++) {
				if (maskarr[i] == 0x00) {
					rownum += 8;
				} else {
					for (short j = 0; j < 8; j++) {
						if (!(maskarr[i] & (0x80 >> j))) {
							rownum++;
							continue;
						}
//						load_one_row(&bm_tmp, dumpfile, bmnum, rownum, false);
						add_row(&bm_tmp, dumpfile, rownum, add_row_dim, bmnum, false);
						assert(bm_tmp.bm.size() == 1);
						list<struct row>::iterator bmit = bm_tmp.bm.begin();
//						unsigned int size = 0;
						unsigned int sz = 0;
						memcpy(&sz, (*bmit).data, BM_ROW_SIZE);
						unsigned char *data = (unsigned char *) malloc (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
						memcpy(data, (*bmit).data, TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
						struct row r = {rownum, data};
						if (loadvbm) {
							bitmat->vbm.push_back(r);
						} else {
							bitmat->bm.push_back(r);
						}
//						bitmat->bm.push_back(r);
						bm_tmp.reset();

						rownum++;
					}
				}
			}

			return bitmat->num_triples;

		}

		unsigned int sz = 0;

		for (unsigned int i = 0; i < limit_bytes; i++) {
			if (bitmat->rowfold[i] == 0x00) {
				rownum += 8;
				continue;
			}
			for (short j = 0; j < 8; j++) {
				if (!(bitmat->rowfold[i] & (0x80 >> j))) {
					rownum++;
					continue;
				}
#if MMAPFILES
				memcpy(&sz, fpos, BM_ROW_SIZE);
				fpos += BM_ROW_SIZE;
#else
				read(fd, &sz, BM_ROW_SIZE);
#endif
				if ((maskarr[i] & (0x80 >> j)) == 0x00) {
#if MMAPFILES
					fpos += TOTAL_ROWSIZE(sz);
#else
					lseek(fd, sz, SEEK_CUR);
#endif
				} else {
					unsigned char *data = NULL;
					data = (unsigned char *) malloc (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
					memcpy(data, &sz, BM_ROW_SIZE);
#if MMAPFILES
					memcpy(data + BM_ROW_SIZE, fpos, TOTAL_ROWSIZE(sz));
					fpos += TOTAL_ROWSIZE(sz);
#else
					read(fd, data + BM_ROW_SIZE, TOTAL_ROWSIZE(sz));
#endif
					struct row r = {rownum, data};
					if (loadvbm) {
						bitmat->vbm.push_back(r);
					} else {
						bitmat->bm.push_back(r);
					}
				}

				rownum++;
			}

		}
	}

#if MMAPFILES
#else
	close(fd);
#endif

	return bitmat->num_triples;
}


/*
 * This is for the new delayed loading of BitMat only in the pruning step
 */
unsigned int wrapper_load_from_dump_file2(BitMat *bitmat, unsigned int bmnum)
{
	char dumpfile[254];

	switch (bitmat->dimension) {
		case (SPO_BITMAT):
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_SPO")].c_str());
			init_bitmat(bitmat, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
			break;
		case (OPS_BITMAT):
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_OPS")].c_str());
			init_bitmat(bitmat, gnum_objs, gnum_preds, gnum_subs, gnum_comm_so, OPS_BITMAT);
			break;
		case (PSO_BITMAT):
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
			init_bitmat(bitmat, gnum_preds, gnum_subs, gnum_objs, gnum_comm_so, PSO_BITMAT);
			break;
		case (POS_BITMAT):
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
			init_bitmat(bitmat, gnum_preds, gnum_objs, gnum_subs, gnum_comm_so, POS_BITMAT);
			break;

		default:
			assert(0);
	}

	int fd = 0;
	unsigned int size = 0;
	unsigned int total_size = 0;
	char *fpos = NULL;

#if MMAPFILES
	fpos = mmap_table[string(dumpfile)];
	fpos += get_offset(dumpfile, bmnum);
#else
	fd = open(fname_dump, O_RDONLY);
	if (fd < 0) {
		printf("*** ERROR opening dump file %s\n", fname_dump);
		assert(0);
	}
	unsigned long offset = get_offset(fname_dump, bmnum);
	if (offset > 0) {
		lseek(fd, offset, SEEK_CUR);
	}
#endif

	assert((fd == 0) ^ (NULL == fpos));

#if MMAPFILES
	//		memcpy(&bitmat->num_triples, fpos, sizeof(unsigned int));
	fpos += sizeof(unsigned int);
#else
	read(fd, &bitmat->num_triples, sizeof(unsigned int));
#endif
	if (bitmat->num_triples == 0) {
		cout << "wrapper_load_from_dump_file2: 0 results" << endl;
		return bitmat->num_triples;
	}

	if (bitmat->rowfold == NULL) {
		bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof (unsigned char));
		memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof (unsigned char));

	}
	if (bitmat->colfold == NULL) {
		bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof (unsigned char));
		memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof (unsigned char));
	}

	if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
		//first read size of the compressed array
		unsigned int comp_arr_size = 0;
#if MMAPFILES
		memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
		fpos += ROW_SIZE_BYTES;
#else
		read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//			total_size += ROW_SIZE_BYTES;
		unsigned char *comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
		memcpy(comp_arr, fpos, comp_arr_size);
		fpos += comp_arr_size;
#else
		read(fd, comp_arr, comp_arr_size);
#endif
		dgap_uncompress(comp_arr, comp_arr_size, bitmat->rowfold, bitmat->row_bytes);
//			total_size += comp_arr_size;

		free(comp_arr);

		if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
			comp_arr_size = 0;
#if MMAPFILES
			memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
			fpos += ROW_SIZE_BYTES;
#else
			read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//				total_size += ROW_SIZE_BYTES;

			comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
			memcpy(comp_arr, fpos, comp_arr_size);
			fpos += comp_arr_size;
#else
			read(fd, comp_arr, comp_arr_size);
#endif
			dgap_uncompress(comp_arr, comp_arr_size, bitmat->colfold, bitmat->column_bytes);
//			total_size += comp_arr_size;

			free(comp_arr);
		}

	} else {
#if MMAPFILES
		memcpy(bitmat->rowfold, fpos, bitmat->row_bytes);
		fpos += bitmat->row_bytes;
		memcpy(bitmat->colfold, fpos, bitmat->column_bytes);
		fpos += bitmat->column_bytes;
#else
		read(fd, bitmat->rowfold, bitmat->row_bytes);
		read(fd, bitmat->colfold, bitmat->column_bytes);
#endif
	}

	unsigned char *bmfold = NULL;
	unsigned int limit_bytes = 0, bmfold_size = 0;

//	if (joinres != NULL && !joinres->empty) {
//		assert(joinres->joinres_dim != PRED_DIMENSION);
//
//		if (joinres->joinres_dim == SO_DIMENSION) {
//			 limit_bytes = bitmat->common_so_bytes;
//		}
//		if (bmjoindim == ROW) {
//			if ((bitmat->dimension == SPO_BITMAT && joinres->joinres_dim == SUB_DIMENSION)
//					||
//				(bitmat->dimension == OPS_BITMAT && joinres->joinres_dim == OBJ_DIMENSION)) {
//
//				assert(joinres->joinres_size == bitmat->row_bytes);
//				limit_bytes = bitmat->row_bytes;
//			} else if ((bitmat->dimension == SPO_BITMAT && joinres->joinres_dim == OBJ_DIMENSION)
//					||
//					(bitmat->dimension == OPS_BITMAT && joinres->joinres_dim == SUB_DIMENSION)) {
//
//				limit_bytes = bitmat->common_so_bytes;
//				memset(&joinres->joinres[bitmat->common_so_bytes], 0x00,
//						joinres->joinres_size - bitmat->common_so_bytes);
//
//				joinres->joinres[bitmat->common_so_bytes-1] &= (0xff << (8-(bitmat->num_comm_so%8)));
//				joinres->joinres_size = bitmat->common_so_bytes;
//				joinres->joinres_dim = SO_DIMENSION;
//
//			} else if (bitmat->dimension == POS_BITMAT || bitmat->dimension == PSO_BITMAT) {
//				assert(0);
//			}
//
//			bmfold = bitmat->rowfold;
//			bmfold_size = bitmat->row_bytes;
//
//		} else if (bmjoindim == COLUMN) {
//			if (((bitmat->dimension == SPO_BITMAT || bitmat->dimension == PSO_BITMAT)
//					&& joinres->joinres_dim == OBJ_DIMENSION)
//					||
//					((bitmat->dimension == OPS_BITMAT || bitmat->dimension == POS_BITMAT)
//					&& joinres->joinres_dim == SUB_DIMENSION)) {
//
//				assert(joinres->joinres_size == bitmat->column_bytes);
//				limit_bytes = bitmat->column_bytes;
//			} else if ((bitmat->dimension == SPO_BITMAT && joinres->joinres_dim == SUB_DIMENSION)
//						||
//						(bitmat->dimension == OPS_BITMAT && joinres->joinres_dim == OBJ_DIMENSION)) {
//
//				limit_bytes = bitmat->common_so_bytes;
//				memset(&joinres->joinres[bitmat->common_so_bytes], 0x00,
//						joinres->joinres_size - bitmat->common_so_bytes);
//
//				joinres->joinres[bitmat->common_so_bytes-1] &= (0xff << (8-(bitmat->num_comm_so%8)));
//				joinres->joinres_size = bitmat->common_so_bytes;
//				joinres->joinres_dim = SO_DIMENSION;
//
//			} else {
//				assert(0);
//			}
//
//			bmfold = bitmat->colfold;
//			bmfold_size = bitmat->column_bytes;
//
//		} else {
//			cout << "*** ERROR: joindim " << bmjoindim << " ****";
//			assert(0);
//		}
//
//		for (unsigned int i = 0; i < limit_bytes; i++) {
//			joinres->joinres[i] = bmfold[i] & joinres->joinres[i];
//		}
//	}

	unsigned int ret = load_from_dump_file2(dumpfile, bmnum, bitmat, false, false,
//							 	(joinres == NULL ? NULL : joinres->joinres),
								NULL, limit_bytes,
//								bmjoindim,
								0, fpos, fd, false);

//	if (joinres != NULL && !joinres->empty) {
//		memcpy(bmfold, joinres->joinres, limit_bytes);
//		if (limit_bytes > bmfold_size)
//			memset(&bmfold[limit_bytes], 0, bmfold_size - limit_bytes);
//
//		if (bmjoindim == ROW) {
//			bitmat->last_op = ROW_UNFOLD;
//		} else if (bmjoindim == COLUMN) {
//			bitmat->last_op = COLUMN_VEC;
//		}
//	} else {
//		bitmat->last_op = BMUPTODATE;
//	}

	return ret;
}

///////////////////////////////////////////////////////////

unsigned int wrapper_load_from_dump_file(char *fname_dump, unsigned int bmnum, struct node *gnode,
		bool readtcnt, bool readarray)
{
//	if(tp->bitmat.dimension == POS_BITMAT || tp->bitmat.dimension == PSO_BITMAT) {
//		return load_from_dump_file(fname_dump, offset, &tp->bitmat, readtcnt, readarray, NULL, 0);
//	}

	BitMat *bitmat = &((TP *)gnode->data)->bitmat;
	int fd = 0;
	unsigned int size = 0;
	unsigned int total_size = 0;
	char *fpos = NULL;

#if MMAPFILES
	fpos = mmap_table[string(fname_dump)];
	fpos += get_offset(fname_dump, bmnum);
#else
	fd = open(fname_dump, O_RDONLY);
	if (fd < 0) {
		printf("*** ERROR opening dump file %s\n", fname_dump);
		assert(0);
	}
	unsigned long offset = get_offset(fname_dump, bmnum);
	if (offset > 0) {
		lseek(fd, offset, SEEK_CUR);
	}
#endif

	assert((fd == 0) ^ (NULL == fpos));

	if (readtcnt) {
		bitmat->num_triples = 0;
#if MMAPFILES
		memcpy(&bitmat->num_triples, fpos, sizeof(unsigned int));
		fpos += sizeof(unsigned int);
#else
		read(fd, &bitmat->num_triples, sizeof(unsigned int));
#endif
//		total_size += sizeof(unsigned int);
		if (bitmat->num_triples == 0) {
			cout << "wrapper_load_from_dump_file: 0 results" << endl;
			return bitmat->num_triples;
		}
	}
	if (readarray) {
		if (bitmat->rowfold == NULL) {
			bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof (unsigned char));
			memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof (unsigned char));

		}
		if (bitmat->colfold == NULL) {
			bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof (unsigned char));
			memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof (unsigned char));
		}

		if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
			//first read size of the compressed array
			unsigned int comp_arr_size = 0;
#if MMAPFILES
			memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
			fpos += ROW_SIZE_BYTES;
#else
			read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//			total_size += ROW_SIZE_BYTES;
			unsigned char *comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
			memcpy(comp_arr, fpos, comp_arr_size);
			fpos += comp_arr_size;
#else
			read(fd, comp_arr, comp_arr_size);
#endif
			dgap_uncompress(comp_arr, comp_arr_size, bitmat->rowfold, bitmat->row_bytes);
//			total_size += comp_arr_size;

			free(comp_arr);

			if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
				comp_arr_size = 0;
#if MMAPFILES
				memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
				fpos += ROW_SIZE_BYTES;
#else
				read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//				total_size += ROW_SIZE_BYTES;

				comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
				memcpy(comp_arr, fpos, comp_arr_size);
				fpos += comp_arr_size;
#else
				read(fd, comp_arr, comp_arr_size);
#endif
				dgap_uncompress(comp_arr, comp_arr_size, bitmat->colfold, bitmat->column_bytes);
	//			total_size += comp_arr_size;

				free(comp_arr);
			}

		} else {
#if MMAPFILES
			memcpy(bitmat->rowfold, fpos, bitmat->row_bytes);
			fpos += bitmat->row_bytes;
			memcpy(bitmat->colfold, fpos, bitmat->column_bytes);
			fpos += bitmat->column_bytes;
#else
			read(fd, bitmat->rowfold, bitmat->row_bytes);
			read(fd, bitmat->colfold, bitmat->column_bytes);
#endif
		}

	}

//#if MMAPFILES
//#else
//	close(fd);
//#endif

//	cout << "Num objects here " << count_bits_in_row(bitmat->colfold, bitmat->column_bytes) << " Num subjects here " << count_bits_in_row(bitmat->rowfold, bitmat->row_bytes) << endl;

//	cout << "-----------------" << endl;
//	print_set_bits_in_row(bitmat->colfold, bitmat->column_bytes);
//	cout << "-----------------" << endl;

//	char dumpfile[1024];
//	strcpy(dumpfile, fname_dump);

//	unsigned int rows = count_bits_in_row(bitmat->rowfold, bitmat->row_bytes);
//	unsigned int columns = count_bits_in_row(bitmat->colfold, bitmat->column_bytes);
//	unsigned char *maskarr = NULL;
//	unsigned int maskarr_size = 0;
	TP *tp = (TP *)gnode->data;
	LIST *next = gnode->nextTP;
	unsigned char *rowarr = NULL;
	int rowarr_dim = 0, rowarr_pos = 0, rowarr_jvar = 0;
	unsigned int rowarr_size = 0, limit_bytes = 0;

	//TODO: Needs a fix for dimensions of bitmat, e.g. LUBM query9
	for (; next; next=next->next) {
		TP *other = (TP *)next->gnode->data;
		if (other->bitmat.bm.size() == 0)
			continue;

		if (!tp->isSlave(other->strlevel) && other->isSlave(tp->strlevel))
			continue;

		//other is either the peer or master of tp
		//S-S join
		if (tp->sub < 0 && other->sub == tp->sub) {
			if (other->bitmat.dimension == SPO_BITMAT && tp->bitmat.dimension == SPO_BITMAT) {
//				cout << "------taking bindings from other " << other->toString() << endl;
				if (rowarr == NULL) {
					rowarr = (unsigned char *) malloc (tp->bitmat.row_bytes);
					rowarr_size = tp->bitmat.row_bytes;
					rowarr_dim = ROW;
					rowarr_pos = SUB_DIMENSION;
					rowarr_jvar = tp->sub;
					for (unsigned int i = 0; i < tp->bitmat.row_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.row_bytes ? rowarr_size : other->bitmat.row_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == POS_BITMAT || other->bitmat.dimension == OPS_BITMAT)
							&& tp->bitmat.dimension == SPO_BITMAT) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.row_bytes);
					rowarr_size = tp->bitmat.row_bytes;
					assert(other->bitmat.column_bytes == tp->bitmat.row_bytes);
					rowarr_dim = ROW;
					rowarr_pos = SUB_DIMENSION;
					rowarr_jvar = tp->sub;
					for (unsigned int i = 0; i < tp->bitmat.row_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.column_bytes ? rowarr_size : other->bitmat.column_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			} else if (other->bitmat.dimension == SPO_BITMAT && tp->bitmat.dimension == POS_BITMAT) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.column_bytes);
					rowarr_size = tp->bitmat.column_bytes;
					assert(other->bitmat.row_bytes == tp->bitmat.column_bytes);
					rowarr_dim = COLUMN;
					rowarr_pos = SUB_DIMENSION;
					rowarr_jvar = tp->sub;
					for (unsigned int i = 0; i < tp->bitmat.column_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.row_bytes ? rowarr_size : other->bitmat.row_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == POS_BITMAT || other->bitmat.dimension == OPS_BITMAT)
							&& tp->bitmat.dimension == POS_BITMAT) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.column_bytes);
					rowarr_size = tp->bitmat.column_bytes;
					assert(other->bitmat.column_bytes == tp->bitmat.column_bytes);
					rowarr_dim = COLUMN;
					rowarr_pos = SUB_DIMENSION;
					rowarr_jvar = tp->sub;
					for (unsigned int i = 0; i < tp->bitmat.column_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.column_bytes ? rowarr_size : other->bitmat.column_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			}
		} else if (tp->pred < 0 && other->pred == tp->pred) {
			//P-P join
			if ((other->bitmat.dimension == POS_BITMAT || other->bitmat.dimension == PSO_BITMAT) &&
				(tp->bitmat.dimension == POS_BITMAT || tp->bitmat.dimension == PSO_BITMAT)) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.row_bytes);
					rowarr_size = tp->bitmat.row_bytes;
					rowarr_dim = ROW;
					rowarr_pos = PRED_DIMENSION;
					rowarr_jvar = tp->pred;
					for (unsigned int i = 0; i < tp->bitmat.row_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->pred) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					assert(rowarr_size == other->bitmat.row_bytes);
					for (unsigned int i = 0; i < other->bitmat.row_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else {
				assert(0);
			}
		} else if (tp->obj < 0 && other->obj == tp->obj) {
			//O-O join
			if (other->bitmat.dimension == OPS_BITMAT && tp->bitmat.dimension == OPS_BITMAT) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.row_bytes);
					rowarr_size = tp->bitmat.row_bytes;
					assert(other->bitmat.row_bytes == tp->bitmat.row_bytes);
					rowarr_dim = ROW;
					rowarr_pos = OBJ_DIMENSION;
					rowarr_jvar = tp->obj;
					for (unsigned int i = 0; i < tp->bitmat.row_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.row_bytes ? rowarr_size : other->bitmat.row_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == PSO_BITMAT || other->bitmat.dimension == SPO_BITMAT) && tp->bitmat.dimension == OPS_BITMAT) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.row_bytes);
					rowarr_size = tp->bitmat.row_bytes;
					assert(other->bitmat.column_bytes == tp->bitmat.row_bytes);
					rowarr_dim = ROW;
					rowarr_pos = OBJ_DIMENSION;
					rowarr_jvar = tp->obj;
					for (unsigned int i = 0; i < tp->bitmat.row_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.column_bytes ? rowarr_size : other->bitmat.column_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			} else if (other->bitmat.dimension == OPS_BITMAT && tp->bitmat.dimension == PSO_BITMAT) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.column_bytes);
					rowarr_size = tp->bitmat.column_bytes;
					assert(other->bitmat.row_bytes == tp->bitmat.column_bytes);
					rowarr_dim = COLUMN;
					rowarr_pos = OBJ_DIMENSION;
					rowarr_jvar = tp->obj;
					for (unsigned int i = 0; i < other->bitmat.row_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.row_bytes ? rowarr_size : other->bitmat.row_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == PSO_BITMAT || other->bitmat.dimension == SPO_BITMAT)
							&& tp->bitmat.dimension == PSO_BITMAT) {
				if (rowarr == NULL) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					rowarr = (unsigned char *) malloc (tp->bitmat.column_bytes);
					rowarr_size = tp->bitmat.column_bytes;
					assert(other->bitmat.column_bytes == tp->bitmat.column_bytes);
					rowarr_dim = COLUMN;
					rowarr_pos = OBJ_DIMENSION;
					rowarr_jvar = tp->obj;
					for (unsigned int i = 0; i < tp->bitmat.column_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					limit_bytes = rowarr_size < other->bitmat.column_bytes ? rowarr_size : other->bitmat.column_bytes;
					for (unsigned int i = 0; i < limit_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			}
		} else if (tp->sub < 0 && other->obj == tp->sub) {
			//S-O join1
			if (tp->bitmat.dimension != SPO_BITMAT && tp->bitmat.dimension != POS_BITMAT)
				continue;
			bool rowarr_newly_set = false;
			if (rowarr == NULL) {
				rowarr = (unsigned char *) malloc (tp->bitmat.common_so_bytes);
				rowarr_size = tp->bitmat.common_so_bytes;
				rowarr_pos = SO_DIMENSION;
				rowarr_jvar = tp->sub;
				switch (tp->bitmat.dimension) {
				case (SPO_BITMAT):
					rowarr_dim = ROW;
					break;
				case (POS_BITMAT):
					rowarr_dim = COLUMN;
					break;
				default:
					assert(0);
					break;
				}
				rowarr_newly_set = true;
			} else if (rowarr_jvar == tp->sub) {
				if (rowarr_pos != SO_DIMENSION) {
					assert(rowarr_size >= tp->bitmat.common_so_bytes);
					memset(&rowarr[tp->bitmat.common_so_bytes], 0, rowarr_size - tp->bitmat.common_so_bytes);
					rowarr[tp->bitmat.common_so_bytes-1] &= (0xff << (8-(tp->bitmat.num_comm_so%8)));
					rowarr_pos = SO_DIMENSION;
					rowarr_size = tp->bitmat.common_so_bytes;
				} else {
					assert(rowarr_size == tp->bitmat.common_so_bytes);
				}
			}

			if (other->bitmat.dimension == OPS_BITMAT && tp->bitmat.dimension == SPO_BITMAT) {
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == PSO_BITMAT || other->bitmat.dimension == SPO_BITMAT)
						&& tp->bitmat.dimension == SPO_BITMAT) {
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			} else if (other->bitmat.dimension == OPS_BITMAT && tp->bitmat.dimension == POS_BITMAT) {
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == PSO_BITMAT || other->bitmat.dimension == SPO_BITMAT)
							&& tp->bitmat.dimension == POS_BITMAT) {
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->sub) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			}
			if (rowarr_jvar == tp->sub) {
				rowarr[bitmat->common_so_bytes-1] &= (0xff << (8-(bitmat->num_comm_so%8)));
			}

		} else if (tp->obj < 0 && other->sub == tp->obj) {
			//S-O join2
			if (tp->bitmat.dimension != OPS_BITMAT && tp->bitmat.dimension != PSO_BITMAT)
				continue;
			bool rowarr_newly_set = false;
			if (rowarr == NULL) {
				rowarr = (unsigned char *) malloc (tp->bitmat.common_so_bytes);
				rowarr_size = tp->bitmat.common_so_bytes;
				rowarr_pos = SO_DIMENSION;
				rowarr_jvar = tp->obj;
				switch (tp->bitmat.dimension) {
				case (OPS_BITMAT):
					rowarr_dim = ROW;
					break;
				case (PSO_BITMAT):
					rowarr_dim = COLUMN;
					break;
				default:
//					cout << "tp->bitmat.dim " << tp->bitmat.dimension << endl;
//					continue;
					assert(0);
					break;
				}
				rowarr_newly_set = true;
			} else if (rowarr_jvar == tp->obj) {
				if (rowarr_pos != SO_DIMENSION) {
					assert(rowarr_size >= tp->bitmat.common_so_bytes);
					memset(&rowarr[tp->bitmat.common_so_bytes], 0, rowarr_size - tp->bitmat.common_so_bytes);
					rowarr[tp->bitmat.common_so_bytes-1] &= (0xff << (8-(tp->bitmat.num_comm_so%8)));
					rowarr_pos = SO_DIMENSION;
					rowarr_size = tp->bitmat.common_so_bytes;
				} else {
					assert(rowarr_size == tp->bitmat.common_so_bytes);
				}
			}

			if (other->bitmat.dimension == SPO_BITMAT && tp->bitmat.dimension == OPS_BITMAT) {
//				rowarr_dim = ROW;
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == POS_BITMAT || other->bitmat.dimension == OPS_BITMAT)
						&& tp->bitmat.dimension == OPS_BITMAT) {
//				rowarr_dim = ROW;
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.rowfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			} else if (other->bitmat.dimension == SPO_BITMAT && tp->bitmat.dimension == PSO_BITMAT) {
//				rowarr_dim = COLUMN;
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.rowfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.rowfold[i];
					}
				}
			} else if ((other->bitmat.dimension == POS_BITMAT || other->bitmat.dimension == OPS_BITMAT)
							&& tp->bitmat.dimension == PSO_BITMAT) {
//				rowarr_dim = COLUMN;
				if (rowarr_newly_set) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = tp->bitmat.colfold[i] & other->bitmat.colfold[i];
					}
				} else if (rowarr_jvar == tp->obj) {
//					cout << "------taking bindings from other " << other->toString() << endl;
					for (unsigned int i = 0; i < tp->bitmat.common_so_bytes; i++) {
						rowarr[i] = rowarr[i] & other->bitmat.colfold[i];
					}
				}
			}
			if (rowarr_jvar == tp->obj) {
				rowarr[bitmat->common_so_bytes-1] &= (0xff << (8-(bitmat->num_comm_so%8)));
			}
		}
	}

	///////////
	//DEBUG
	//TODO: remove later
	//////////
//	cout << "--------Wrapper_load_from_dump_file: " << tp->toString() << " bitmat.dim " <<
//		tp->bitmat.dimension << " rowarr is " << (rowarr == NULL ? "NULL " : "NOT NULL ") << " rowarr_size "
//		<< rowarr_size << endl;
//	unsigned int setbits = print_set_bits_in_row(rowarr, rowarr_size);
//	cout << "--------setbits in maskarr/rowarr " << count_bits_in_row(rowarr, rowarr_size) << endl;
	///////////////////////////////////

	unsigned int ret = load_from_dump_file(fname_dump, bmnum, bitmat, false, false, rowarr,
								rowarr_size, rowarr_dim, fpos, fd, true);

//	cout << "--------triples in bitmat " << count_triples_in_bitmat(bitmat) << endl;

	free(rowarr);

	return ret;
}

///////////////////////////////////////////////////////////
unsigned int load_from_dump_file(char *fname_dump, unsigned int bmnum, BitMat *bitmat, bool readtcnt,
			bool readarray, unsigned char *maskarr, unsigned int maskarr_size, int maskarr_dim,
			char *fpos, int fd, bool fold_objdim)
{
//	struct timeval start_time, stop_time;
//	clock_t  t_start, t_end;
//	double curr_time;
//	double st, en;
//	gettimeofday(&start_time, (struct timezone *)0);
//	cout << "Inside load_from_dump_file" << endl;
//	int fd = 0;
//	unsigned int size = 0;
//	unsigned int total_size = 0;
//	char *fpos = NULL;

	assert((readtcnt && readarray) || (!readtcnt && !readarray && ((fd == 0) ^ (NULL == fpos))));

#if MMAPFILES
	if (fpos == NULL) {
		fpos = mmap_table[string(fname_dump)];
		fpos += get_offset(fname_dump, bmnum);
	}
#else
	if (fd == 0) {
		fd = open(fname_dump, O_RDONLY);
		if (fd < 0) {
			printf("*** ERROR opening dump file %s\n", fname_dump);
			assert(0);
		}
		unsigned long offset = get_offset(fname_dump, bmnum);
		///////
		//DEBUG: remove later
		////////
//		cout << "load_from_dump_file: offset is " << offset << endl;
		if (offset > 0) {
			lseek(fd, offset, SEEK_CUR);
		}
	}
#endif

	assert((fd == 0) ^ (NULL == fpos));
//	cout << "load_from_dump_file: offset " << offset << endl;

	if (readtcnt) {
		bitmat->num_triples = 0;
#if MMAPFILES
		memcpy(&bitmat->num_triples, fpos, sizeof(unsigned int));
		fpos += sizeof(unsigned int);
#else
		read(fd, &bitmat->num_triples, sizeof(unsigned int));
#endif
//		total_size += sizeof(unsigned int);
		if (bitmat->num_triples == 0) {
			cout << "load_from_dump_file: 0 results" << endl;
			return bitmat->num_triples;
		}
	}
//	cout << "Num triples1 " << bitmat->num_triples << endl;
	if (readarray) {
		if (bitmat->rowfold == NULL) {
			bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof (unsigned char));
			memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof (unsigned char));

		}
		if (bitmat->colfold == NULL) {
			bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof (unsigned char));
			memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof (unsigned char));
		}

		if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
			//first read size of the compressed array
			unsigned int comp_arr_size = 0;
#if MMAPFILES
			memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
			fpos += ROW_SIZE_BYTES;
#else
			read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//			total_size += ROW_SIZE_BYTES;
			unsigned char *comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
			memcpy(comp_arr, fpos, comp_arr_size);
			fpos += comp_arr_size;
#else
			read(fd, comp_arr, comp_arr_size);
#endif
			dgap_uncompress(comp_arr, comp_arr_size, bitmat->rowfold, bitmat->row_bytes);
//			total_size += comp_arr_size;

			free(comp_arr);

			if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
				comp_arr_size = 0;
#if MMAPFILES
				memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
				fpos += ROW_SIZE_BYTES;
#else
				read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//				total_size += ROW_SIZE_BYTES;

				comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
				memcpy(comp_arr, fpos, comp_arr_size);
				fpos += comp_arr_size;
#else
				read(fd, comp_arr, comp_arr_size);
#endif
				dgap_uncompress(comp_arr, comp_arr_size, bitmat->colfold, bitmat->column_bytes);
	//			total_size += comp_arr_size;

				free(comp_arr);
			}

		} else {
#if MMAPFILES
			memcpy(bitmat->rowfold, fpos, bitmat->row_bytes);
			fpos += bitmat->row_bytes;
			memcpy(bitmat->colfold, fpos, bitmat->column_bytes);
			fpos += bitmat->column_bytes;
#else
			read(fd, bitmat->rowfold, bitmat->row_bytes);
			read(fd, bitmat->colfold, bitmat->column_bytes);
#endif
//			total_size += bitmat->row_bytes;
//			total_size += bitmat->column_bytes;
		}
	}
//	cout << "Num bits in rowfold " << count_bits_in_row(bitmat->rowfold, bitmat->row_bytes) << endl;
//	cout << "Num bits in colfold " << count_bits_in_row(bitmat->colfold, bitmat->column_bytes) << endl;
//	print_set_bits_in_row(bitmat->colfold, bitmat->column_bytes);
	unsigned int limit_bytes = 0;
	if (maskarr_size == 0) {
		assert(maskarr == NULL);
		limit_bytes = bitmat->row_bytes;
	}

	bool skipAtLeastOnce = false, skip = false;
	unsigned int last_bit = last_set_bit(maskarr, maskarr_size);
	unsigned int total_set_bits = count_bits_in_row(maskarr, maskarr_size);
	unsigned int rowfold_bits = count_bits_in_row(bitmat->rowfold, bitmat->row_bytes);
	/////////
	//DEBUG
	//TODO: remove later
//	cout << "--------- load_from_dump_file: maskarr_size " << maskarr_size << " maskarr_dim " << maskarr_dim << " last_setbit " << last_bit
//		<< " total_set_bits " << total_set_bits << " rowfold_bits " << rowfold_bits << endl;
	/////////////////////////////
	maskarr_size = last_bit/8 + ((last_bit%8) > 0 ? 1 : 0);
	if (maskarr_dim == ROW) {
		limit_bytes = bitmat->row_bytes > maskarr_size ? maskarr_size : bitmat->row_bytes;
	} else {
		limit_bytes = bitmat->row_bytes;
	}

	unsigned int rownum = 1;
	//////////
	//DEBUG
//	off_t total_size = 0;

	unsigned int sz = 0;

	if (maskarr == NULL || maskarr_dim == COLUMN) {
		for (unsigned int i = 0; i < bitmat->row_bytes; i++) {
			if (bitmat->rowfold[i] == 0x00) {
				rownum += 8;
			} else {
				for (short j = 0; j < 8; j++) {
					if (!(bitmat->rowfold[i] & (0x80 >> j))) {
						rownum++;
						continue;
					}
	#if MMAPFILES
					memcpy(&sz, fpos, BM_ROW_SIZE);
					fpos += BM_ROW_SIZE;
	#else
					read(fd, &sz, BM_ROW_SIZE);
	#endif
					unsigned char *data = (unsigned char *) malloc (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
					memcpy(data, &sz, BM_ROW_SIZE);
					////////
					//DEBUG
//					total_size += (8+size+ROW_SIZE_BYTES);
//					cout << "total_size " << total_size << endl;
	#if MMAPFILES
					memcpy(data + BM_ROW_SIZE, fpos, TOTAL_ROWSIZE(sz));
					fpos += sz;
	#else
					read(fd, data + BM_ROW_SIZE, TOTAL_ROWSIZE(sz));
	#endif
					struct row r = {rownum, data};
//					if (!vbm_load)
					bitmat->bm.push_back(r);
//					else
//						bitmat->vbm.push_back(r);

					rownum++;
				}
			}
		}

	} else if (maskarr_dim == ROW) {
		if ((total_set_bits <= ((maskarr_size * 8)/SELECTIVITY_THRESHOLD)) &&
				rowfold_bits > 10*total_set_bits &&
				(bitmat->dimension == SPO_BITMAT || bitmat->dimension == OPS_BITMAT)) {
//			cout << "--------- load_from_dump_file: NOT loading the traditional way" << endl;
			char dumpfile[1024];
			BitMat bm_tmp;
			unsigned int add_row_dim = 0;
			if (bitmat->dimension == SPO_BITMAT) {
				sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
				shallow_init_bitmat(&bm_tmp, gnum_preds, gnum_subs, gnum_objs, gnum_comm_so, PSO_BITMAT);
				add_row_dim = PSO_BITMAT;

			} else if (bitmat->dimension == OPS_BITMAT) {
				sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
				shallow_init_bitmat(&bm_tmp, gnum_preds, gnum_objs, gnum_subs, gnum_comm_so, POS_BITMAT);
				add_row_dim = POS_BITMAT;
			}

			rownum = 1;
			for (unsigned int i = 0; i < limit_bytes; i++) {
				if (maskarr[i] == 0x00) {
					rownum += 8;
				} else {
					for (short j = 0; j < 8; j++) {
						if (!(maskarr[i] & (0x80 >> j))) {
							rownum++;
							continue;
						}
						add_row(&bm_tmp, dumpfile, rownum, add_row_dim, bmnum, false);
						assert(bm_tmp.bm.size() == 1);
						list<struct row>::iterator bmit = bm_tmp.bm.begin();
//						unsigned int size = 0;
						unsigned int sz = 0;
						memcpy(&sz, (*bmit).data, BM_ROW_SIZE);
						unsigned char *data = (unsigned char *) malloc (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
						memcpy(data, (*bmit).data, TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
						struct row r = {rownum, data};
//							if (!vbm_load)
						bitmat->bm.push_back(r);
//							else
//								bitmat->vbm.push_back(r);
						bm_tmp.reset();

						rownum++;
					}
				}
			}

			free(bitmat->rowfold);
			bitmat->rowfold = maskarr;
			bitmat->row_bytes = maskarr_size;

/*
			if (bitmat->dimension == SPO_BITMAT) {
				char dumpfile[1024];
				sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());

				BitMat bm_pso;
				shallow_init_bitmat(&bm_pso, gnum_preds, gnum_subs, gnum_objs, gnum_comm_so, PSO_BITMAT);
				rownum = 1;
				for (unsigned int i = 0; i < limit_bytes; i++) {
					if (maskarr[i] == 0x00) {
						rownum += 8;
					} else {
						for (short j = 0; j < 8; j++) {
							if (!(maskarr[i] & (0x80 >> j))) {
								rownum++;
								continue;
							}
							add_row(&bm_pso, dumpfile, rownum, PSO_BITMAT, bmnum, false);
							assert(bm_pso.bm.size() == 1);
							list<struct row>::iterator bmit = bm_pso.bm.begin();
							unsigned int size = 0; memcpy(&size, (*bmit).data, ROW_SIZE_BYTES);
							unsigned char *data = (unsigned char *) malloc (size + ROW_SIZE_BYTES);
							memcpy(data, (*bmit).data, size+ROW_SIZE_BYTES);
							struct row r = {rownum, data};
//							if (!vbm_load)
							bitmat->bm.push_back(r);
//							else
//								bitmat->vbm.push_back(r);
							bm_pso.reset();

							rownum++;
						}
					}
				}

			} else if (bitmat->dimension == OPS_BITMAT) {
				char dumpfile[1024];
				sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());

				BitMat bm_pos;
				shallow_init_bitmat(&bm_pos, gnum_preds, gnum_objs, gnum_subs, gnum_comm_so, POS_BITMAT);
				rownum = 1;
				for (unsigned int i = 0; i < limit_bytes; i++) {
					if (maskarr[i] == 0x00) {
						rownum += 8;
					} else {
						for (short j = 0; j < 8; j++) {
							if (!(maskarr[i] & (0x80 >> j))) {
								rownum++;
								continue;
							}
							add_row(&bm_pos, dumpfile, rownum, POS_BITMAT, bmnum, false);
							assert(bm_pos.bm.size() == 1);
							list<struct row>::iterator bmit = bm_pos.bm.begin();
							unsigned int size = 0; memcpy(&size, (*bmit).data, ROW_SIZE_BYTES);
							unsigned char *data = (unsigned char *) malloc (size + ROW_SIZE_BYTES);
							memcpy(data, (*bmit).data, size+ROW_SIZE_BYTES);
							struct row r = {rownum, data};
//							if (!vbm_load)
							bitmat->bm.push_back(r);
//							else
//								bitmat->vbm.push_back(r);
							bm_pos.reset();

							rownum++;
						}
					}
				}

			}
*/
			simple_fold(bitmat, ROW, bitmat->rowfold, bitmat->row_bytes);
			simple_fold(bitmat, COLUMN, bitmat->colfold, bitmat->column_bytes);
			count_triples_in_bitmat(bitmat);

			return bitmat->num_triples;

		}
		for (unsigned int i = 0; i < limit_bytes; i++) {
			if (bitmat->rowfold[i] == 0x00) {
				rownum += 8;
			} else {
				for (short j = 0; j < 8; j++) {
					if (!(bitmat->rowfold[i] & (0x80 >> j))) {
						rownum++;
						continue;
					}
					skip = false;
	#if MMAPFILES
					memcpy(&sz, fpos, BM_ROW_SIZE);
					fpos += BM_ROW_SIZE;
	#else
					read(fd, &sz, BM_ROW_SIZE);
	#endif
					if ((maskarr[i] & (0x80 >> j)) == 0x00) {
						skip = true;
//						skip_cnt++;
						bitmat->rowfold[i] &= ~(0x80 >> j);
						skipAtLeastOnce = true;
					}
					unsigned char *data = NULL;
					if (!skip) {
						data = (unsigned char *) malloc (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
						memcpy(data, &sz, BM_ROW_SIZE);
					}
	#if MMAPFILES
					if (!skip) {
						memcpy(data + BM_ROW_SIZE, fpos, TOTAL_ROWSIZE(sz));
					}
					fpos += TOTAL_ROWSIZE(sz);
	#else
					if (!skip) {
						read(fd, data + BM_ROW_SIZE, TOTAL_ROWSIZE(sz));
					} else {
						lseek(fd, TOTAL_ROWSIZE(sz), SEEK_CUR);
					}
	#endif
					if (!skip) {
						struct row r = {rownum, data};
//						if (!vbm_load)
						bitmat->bm.push_back(r);
//						else
//							bitmat->vbm.push_back(r);
					}
					rownum++;
				}
			}
		}
		free(bitmat->rowfold);
		bitmat->rowfold = maskarr;
		bitmat->row_bytes = maskarr_size;

//		cout << "load_from_dump_file: skip_cnt " << skip_cnt << endl;
	}

	if (maskarr_dim == COLUMN) {
		//modified simple_unfold dynamically updates triple count
		//in bitmat while unfolding on COLUMN dim.
//		cout << "load_from_dump_file: maskarr_dim COLUMN" << endl;
		simple_unfold(bitmat, maskarr, maskarr_size, COLUMN);
		simple_fold(bitmat, ROW, bitmat->rowfold, bitmat->row_bytes);
//		cout << "bitmat->num_triples " << bitmat->num_triples << endl;
	}

//	gettimeofday(&stop_time, (struct timezone *)0);
//	st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
//	en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
//	curr_time = en-st;
//
//	printf("Time for load_from_dump_file(gettimeofday)1: %f\n", curr_time);

	if (skipAtLeastOnce) {
//		for (unsigned int i = 0; i < limit_bytes; i++) {
//			bitmat->rowfold[i] = bitmat->rowfold[i] & maskarr[i];
//		}
		memset(&bitmat->rowfold[limit_bytes], 0, bitmat->row_bytes - limit_bytes);

		//TODO: Should go away no need with Bitmat->last_op
		if (fold_objdim) {
			simple_fold(bitmat, COLUMN, bitmat->colfold, bitmat->column_bytes);
		}
		count_triples_in_bitmat(bitmat);
	}
//	printf("Total size loaded from file %u\n", total_size);

#if MMAPFILES
#else
	close(fd);
#endif

//	cout << "Exiting load_from_dump_file" << endl;
//	gettimeofday(&stop_time, (struct timezone *)0);
//	st = start_time.tv_sec + (start_time.tv_usec/MICROSEC);
//	en = stop_time.tv_sec + (stop_time.tv_usec/MICROSEC);
//	curr_time = en-st;
//
//	printf("Time for load_from_dump_file(gettimeofday): %f\n", curr_time);

//	cout << "bitmat->num_triples " << bitmat->num_triples << endl;

	return bitmat->num_triples;

//	assert(bitmat->num_triples == testtrip);
//	cout << "Num triples " << bitmat->num_triples << endl;
//	unsigned char *test = (unsigned char *) malloc (10);
//	test[0] = 0xff;
//	test[1] = 0x00;
//	test[2] = 0xff;
//	test[3] = 0x00;
//	test[4] = 0xff;
//	test[5] = 0xff;
//	cout << "Num objects " << count_bits_in_row(test, 10) << endl;
//	print_set_bits_in_row(test, 10);


//	vector<struct triple> triplist;
//	list_enctrips_bitmat_new(bitmat, 44, triplist);
//	for (vector<struct triple>::iterator itr = triplist.begin(); itr != triplist.end(); itr++) {
//		cout << (*itr).sub << ":"<< (*itr).pred << ":" << (*itr).obj << endl;
//	}

//	cout << "222 Num objects here " << count_bits_in_row(bitmat->colfold, bitmat->column_bytes) << " Num subjects here " << count_bits_in_row(bitmat->rowfold, bitmat->row_bytes) << endl;
}

////////////////////////////////////
bool filter_and_load_bitmat(BitMat *bitmat, int fd, char *fpos, unsigned char *and_array,
							unsigned int and_array_size)
{
	assert((fd == 0) ^ (NULL == fpos));

	unsigned int total_size = 0;
	unsigned int i = 0;

//	cout << "Inside filter_and_load_bitmat" << endl;

	if (bitmat->rowfold == NULL) {
		bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes * sizeof (unsigned char));
		memset(bitmat->rowfold, 0, bitmat->row_bytes * sizeof (unsigned char));

	}
	if (bitmat->colfold == NULL) {
		bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes * sizeof (unsigned char));
		memset(bitmat->colfold, 0, bitmat->column_bytes * sizeof (unsigned char));
	}

	if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
		//first read size of the compressed array
		unsigned int comp_arr_size = 0;
#if MMAPFILES
		memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
		fpos += ROW_SIZE_BYTES;
#else
		read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif
//		total_size += ROW_SIZE_BYTES;
		unsigned char *comp_arr = (unsigned char *) malloc (comp_arr_size * sizeof(unsigned char));
#if MMAPFILES
		memcpy(comp_arr, fpos, comp_arr_size);
#else
		read(fd, comp_arr, comp_arr_size);
#endif
		dgap_uncompress(comp_arr, comp_arr_size, bitmat->rowfold, bitmat->row_bytes);
//		total_size += comp_arr_size;

		free(comp_arr);

		if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
			comp_arr_size = 0;
#if MMAPFILES
			memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
			fpos += (ROW_SIZE_BYTES + comp_arr_size);
#else
			read(fd, &comp_arr_size, ROW_SIZE_BYTES);
			lseek(fd, comp_arr_size, SEEK_CUR);
#endif
	//		total_size += ROW_SIZE_BYTES;

			//No need of reading object array as that's going to
			//change after "filtering"
		}
//		total_size += comp_arr_size;

	} else {
#if MMAPFILES
		memcpy(bitmat->rowfold, fpos, bitmat->row_bytes);
		fpos += (bitmat->row_bytes + bitmat->column_bytes);
#else
		read(fd, bitmat->rowfold, bitmat->row_bytes);
		lseek(fd, bitmat->column_bytes, SEEK_CUR);
#endif
//		total_size += bitmat->row_bytes;
//		total_size += bitmat->column_bytes;
	}

//	unsigned char *andres = NULL;
//	unsigned int andres_size = 0;

	unsigned int rownum = 1;

//	cumulative_dgap(and_array, and_array_size, and_array);
//	andres = (unsigned char *) malloc (GAP_SIZE_BYTES * bitmat->num_objs + ROW_SIZE_BYTES);

	unsigned int size = 0;

	for (unsigned int i = 0; i < bitmat->row_bytes; i++) {
		if (bitmat->rowfold[i] == 0x00) {
			rownum += 8;
		} else {
			for (short j = 0; j < 8; j++) {
				if (bitmat->rowfold[i] & (0x80 >> j)) {
#if MMAPFILES
					memcpy(&size, fpos, BM_ROW_SIZE);
					fpos += BM_ROW_SIZE;
#else
					read(fd, &size, BM_ROW_SIZE);
#endif
					unsigned char *data = (unsigned char *) malloc (size + BM_ROW_SIZE);
					memcpy(data, &size, BM_ROW_SIZE);
#if MMAPFILES
					memcpy(data + BM_ROW_SIZE, fpos, size);
					fpos += size;
#else
					read(fd, data + BM_ROW_SIZE, size);
#endif

//					cumulative_dgap(&data[ROW_SIZE_BYTES], size, &data[ROW_SIZE_BYTES]);
//					andres_size = dgap_AND(data + ROW_SIZE_BYTES, size, and_array, and_array_size, &andres[ROW_SIZE_BYTES]);
//					de_cumulative_dgap(&andres[ROW_SIZE_BYTES], andres_size, &andres[ROW_SIZE_BYTES]);

//					if ((andres_size == 1 + GAP_SIZE_BYTES) && !andres[ROW_SIZE_BYTES]) {
//						//AND result is all 0s
//						free(data);
//						//Unset bit in rowfold
//						bitmat->rowfold[i/8] &= ((0x80 >> (i%8)) ^ 0xff);
//						continue;
//					}
//
//					memcpy(andres, &andres_size, ROW_SIZE_BYTES);
//					free(data);
//					data = (unsigned char *) malloc (andres_size + ROW_SIZE_BYTES);
//					memcpy(data, andres, andres_size + ROW_SIZE_BYTES);
					//Populate colfold array too

					unsigned char *res = rowwise_column_unfold(data, and_array, and_array_size);

					if (res != NULL) {
						struct row r = {rownum, res};
						bitmat->bm.push_back(r);
					} else {
						//Unset the bit in rowfold
						bitmat->rowfold[i/8] &= ((0x80 >> j) ^ 0xff);
					}
//					dgap_uncompress(data + ROW_SIZE_BYTES, andres_size, bitmat->colfold, bitmat->column_bytes);
				}
				rownum++;
			}

		}
	}

	//TODO: Optimization required if this method is going to be heavily used.. right now ignoring
//	simple_unfold(bitmat, and_array, and_array_size, COLUMN, 3);
	unsigned int limit_bytes = and_array_size < bitmat->column_bytes ? and_array_size : bitmat->column_bytes;
	for (unsigned int i = 0; i < limit_bytes; i++) {
		bitmat->colfold[i] = bitmat->colfold[i] & and_array[i];
	}
	if (bitmat->column_bytes > limit_bytes)
		memset(&bitmat->colfold[limit_bytes], 0x00, bitmat->column_bytes - limit_bytes);

	bitmat->last_op = BMUPTODATE;

//	simple_fold(bitmat, ROW, bitmat->rowfold, bitmat->row_bytes);
//	simple_fold(bitmat, COLUMN, bitmat->colfold, bitmat->column_bytes);

//	printf("Total size loaded from file %u\n", total_size);

#if MMAPFILES
#else
	close(fd);
#endif

	if (bitmat->num_triples != 0)
		return true;

	return false;

}

////////////////////////////////////////////////////////////
void clear_rows(BitMat *bitmat, bool clearbmrows, bool clearfoldarr, bool optimize)
{
	if (clearbmrows) {
		if (bitmat->bm.size() > 0) {
			for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); ){
				free((*it).data);
				it = bitmat->bm.erase(it);
			}
		}

	}

	if (clearfoldarr) {

		//This was probably done this way for optimizing the load
		//operation, but this is clearly wrong if clear_rows
		//is used independently
		if (optimize) {
			if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
				if (comp_folded_arr == 0 && bitmat->rowfold != NULL)
					memset(bitmat->rowfold, 0, bitmat->row_bytes);
				if (bitmat->colfold != NULL)
					memset(bitmat->colfold, 0, bitmat->column_bytes);
			}
		} else {

			if (bitmat->rowfold != NULL)
				memset(bitmat->rowfold, 0, bitmat->row_bytes);
			if (bitmat->colfold != NULL)
				memset(bitmat->colfold, 0, bitmat->column_bytes);
		}

	}

}
///////////////////////////////////////////////////////////////////////
//TODO: make use of the fact that at any load pt
//if any bitmat is compl empty the result is empty
//since it's a join

bool load_one_row(BitMat *bitmat, char *fname, unsigned int rownum, unsigned int bmnum, bool load_colfold)
{
	int fd = 0;
	unsigned long long total_size = 0;//, size = 0;
	char *fpos = NULL;

	unsigned long offset = get_offset(fname, bmnum);

#if MMAPFILES
	fpos = mmap_table[string(fname)];
#else
	fd = open(fname, O_RDONLY);
	if (fd < 0) {
		cout << "*** ERROR opening dump file " << fname << endl;
		assert(0);
	}

	//Move over number of triples
	lseek(fd, offset + sizeof(unsigned int), SEEK_CUR);
	total_size += sizeof(unsigned int);
#endif

	assert((fpos == NULL) ^ (fd == 0));

	if (bitmat->rowfold == NULL) {
		bitmat->rowfold = (unsigned char *) malloc (bitmat->row_bytes);
		memset(bitmat->rowfold, 0, bitmat->row_bytes);
	}
	if (bitmat->colfold == NULL && load_colfold) {
		bitmat->colfold = (unsigned char *) malloc (bitmat->column_bytes);
		memset(bitmat->colfold, 0, bitmat->column_bytes);
	}

	/*
	 * We don't store compressed colfold (columnfold) for PSO and POS
	 * bitmats.
	 */
	if (bitmat->dimension == PSO_BITMAT || bitmat->dimension == POS_BITMAT || comp_folded_arr) {
		unsigned int comp_arr_size = 0;
#if MMAPFILES
		//sizeof(unsigned int) is to walk over the number of triples bytes
		memcpy(&comp_arr_size, &fpos[offset + sizeof(unsigned int)], ROW_SIZE_BYTES);
		fpos += offset + sizeof(unsigned int) + ROW_SIZE_BYTES;
#else
		read(fd, &comp_arr_size, ROW_SIZE_BYTES);
#endif

		unsigned char *comp_arr =  (unsigned char *) malloc (comp_arr_size);

#if MMAPFILES
		memcpy(comp_arr, fpos, comp_arr_size);
		fpos += comp_arr_size;
#else
		read(fd, comp_arr, comp_arr_size);
#endif
		total_size += ROW_SIZE_BYTES + comp_arr_size;

		dgap_uncompress(comp_arr, comp_arr_size, bitmat->rowfold, bitmat->row_bytes);

		free(comp_arr);

		/*
		 * Only if BitMat is of type SPO or OPS we need to move
		 * over colfold (columnfold) as PSO and POS don't store those
		 */
		if (bitmat->dimension != PSO_BITMAT && bitmat->dimension != POS_BITMAT) {
			comp_arr_size = 0;
#if MMAPFILES
			memcpy(&comp_arr_size, fpos, ROW_SIZE_BYTES);
			fpos += (ROW_SIZE_BYTES + comp_arr_size);
#else
			read(fd, &comp_arr_size, ROW_SIZE_BYTES);
			lseek(fd, comp_arr_size, SEEK_CUR);
#endif
			total_size += ROW_SIZE_BYTES + comp_arr_size;
		}

	} else {
#if MMAPFILES
		memcpy(bitmat->rowfold, fpos, bitmat->row_bytes);
		fpos += (bitmat->row_bytes + bitmat->column_bytes);
#else
		read(fd, bitmat->rowfold, bitmat->row_bytes);
		lseek(fd, bitmat->column_bytes, SEEK_CUR);
#endif

		total_size += (bitmat->row_bytes + bitmat->column_bytes);
	}

	unsigned int sz = 0;
	unsigned int fwd = 0;

	if (bitmat->rowfold[(rownum-1)/8] & (0x80 >> ((rownum-1)%8)) ) {
		//this row exists
		for (unsigned int i = 0; i < rownum - 1; i++) {
			if (bitmat->rowfold[i/8] & (0x80 >> (i%8)) ) {
				//this row exists
#if MMAPFILES
				memcpy(&sz, fpos, BM_ROW_SIZE);
				fwd = (BM_ROW_SIZE + TOTAL_ROWSIZE(sz));
				fpos += fwd;
#else
				read(fd, &sz, BM_ROW_SIZE);
				fwd = TOTAL_ROWSIZE(sz);
				lseek(fd, fwd, SEEK_CUR);
#endif
				total_size += (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
			}
		}
#if MMAPFILES
		memcpy(&sz, fpos, BM_ROW_SIZE);
		fpos += BM_ROW_SIZE;
#else
		read(fd, &sz, BM_ROW_SIZE);
#endif

//		if (bitmat->bm == NULL) {
			//just load one row
//			bitmat->bm = (unsigned char **) malloc (sizeof(unsigned char *));
//		}
		unsigned char *data = (unsigned char *) malloc (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
		memcpy(data, &sz, BM_ROW_SIZE);
#if MMAPFILES
		memcpy(data + BM_ROW_SIZE, fpos, TOTAL_ROWSIZE(sz));
		fpos += (TOTAL_ROWSIZE(sz));
#else
		read(fd, data + BM_ROW_SIZE, TOTAL_ROWSIZE(sz));
#endif
		struct row r = {rownum, data};
		bitmat->bm.push_back(r);
		total_size += (TOTAL_ROWSIZE(sz) + BM_ROW_SIZE);
		//TODO: make rowfold and colfold consistent to this row
		memset(bitmat->rowfold, 0, bitmat->row_bytes);
		bitmat->rowfold[(rownum-1)/8] |= (0x80 >> ((rownum-1)%8));
		bitmat->last_op = ROW_UNFOLD;
		//Populating colfold
		if (load_colfold) {
			dgap_uncompress(data + BM_ROW_SIZE, TOTAL_ROWSIZE(sz),
							bitmat->colfold, bitmat->column_bytes);
			bitmat->last_op = BMUPTODATE;
			count_triples_in_bitmat(bitmat);
			bitmat->last_op = ALLUPTODATE;
		}
//		bitmat->single_row = true;
#if MMAPFILES
#else
		close(fd);
#endif
		return true;
	}

#if MMAPFILES
#else
		close(fd);
#endif

	//since you came here, it means that the
	//row doesn't exist
	bitmat->freebm();
	bitmat->last_op = BMUPTODATE;
	return false;

}
/*
 * Just add one row
 */
bool add_row2(BitMat *bitmat, unsigned int dimension, unsigned int bmnum, unsigned int rownum)
{

	char dumpfile[254];

	switch (dimension) {
		case SPO_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_SPO")].c_str());
			shallow_init_bitmat(bitmat, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
			break;

		case OPS_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_OPS")].c_str());
			shallow_init_bitmat(bitmat, gnum_objs, gnum_preds, gnum_subs, gnum_comm_so, OPS_BITMAT);
			break;
		case PSO_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
			shallow_init_bitmat(bitmat, gnum_preds, gnum_subs, gnum_objs, gnum_comm_so, PSO_BITMAT);
			break;
		case POS_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
			shallow_init_bitmat(bitmat, gnum_preds, gnum_objs, gnum_subs, gnum_comm_so, POS_BITMAT);
			break;
		default:
			assert(0);
			break;
	}
	return load_one_row(bitmat, dumpfile, rownum, bmnum, true);
}

/*
 * Just read the triplecnt
 */
void read_bitmat_triplecnt(BitMat *bitmat, unsigned int dimension, unsigned int bmnum)
{
	char dumpfile[254];

	switch (dimension) {
		case SPO_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_SPO")].c_str());
			break;
		case OPS_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_OPS")].c_str());
			break;
		case PSO_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
			break;
		case POS_BITMAT:
			sprintf(dumpfile, "%s", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
			break;
		default:
			assert(0);
			break;
	}

	unsigned long offset = get_offset(dumpfile, bmnum);

#if MMAPFILES
	char *fpos = mmap_table[string(dumpfile)];
	memcpy(&bitmat->num_triples, &fpos[offset], sizeof(unsigned int));
#else
	int fd = open(dumpfile, O_RDONLY);
	assert(fd != -1);
	if (offset > 0)
		lseek(fd, offset, SEEK_CUR);
	read(fd, &bitmat->num_triples, sizeof(unsigned int));
	close(fd);
#endif

}


///////////////////////////////////////////////////////////
bool add_row(BitMat *bitmat, char *fname, unsigned int bmnum, unsigned int readdim, unsigned int rownum, bool load_colfold)
{

//	cout << "Inside add_row" << endl;

//	cout << "Advanced file pointer" << endl;

	//TODO: this will go away later
	if (readdim == bitmat->dimension) {

//		cout << "Calling load_one_row" << endl;
		return load_one_row(bitmat, fname, rownum, bmnum, load_colfold);
	} else {
		assert(0);
		return false;
	}
/*
	size = 0;

	if (tmpsubf[(rownum-1)/8] & (0x80 >> ((rownum-1)%8)) ) {
		//"bmnum" is the predicate num in the original SPO BM
		if (bitmat->bm == NULL) {
			init_bitmat_rows(bitmat, true, true);
		}
		bitmat->rowfold[(bmnum-1)/8] |= 0x80 >> ((bmnum-1)%8);

		for (unsigned int i = 0; i < rownum - 1; i++) {
			if (tmpsubf[i/8] & (0x80 >> (i%8)) ) {
				//this row exists
				read(fd, &size, ROW_SIZE_BYTES);
				lseek(fd, size, SEEK_CUR);
				total_size += (size + ROW_SIZE_BYTES);
			}
		}
		read(fd, &size, ROW_SIZE_BYTES);
		bitmat->bm[bmnum-1] = (unsigned char *) malloc (size + ROW_SIZE_BYTES);
		memcpy(bitmat->bm[bmnum-1], &size, ROW_SIZE_BYTES);
		read(fd, bitmat->bm[bmnum-1] + ROW_SIZE_BYTES, size);
		total_size += (size + ROW_SIZE_BYTES);
		return true;
	}
	return false;
*/
}
/////////////////////////////////////////////////////
unsigned int get_mask_array(unsigned char *and_array, unsigned int setbit)
{
	unsigned int and_array_size = setbit/8 + (setbit%8 > 0 ? 1:0);
	and_array = (unsigned char *) malloc (and_array_size);
	memset(and_array, 0, and_array_size);
	and_array[(setbit-1)/8] |= (0x80 >> ((setbit-1)%8));

	return and_array_size;
}

//////////////////////////////////////////////////////
unsigned int get_and_array(BitMat *bitmat, unsigned char *and_array, unsigned int bit)
{
	unsigned int t = 1;
	unsigned int and_array_size = 1;
	unsigned later_0 = (bitmat->column_bytes << 3) - bit;
	unsigned ini_0 = bit - 1;

	if (bit == 1) {
		and_array[0] = 0x01;
		memcpy(&and_array[1], &t, GAP_SIZE_BYTES);
		and_array_size += GAP_SIZE_BYTES;
	} else {
		and_array[0] = 0x00;
		memcpy(&and_array[1], &ini_0, GAP_SIZE_BYTES);
		memcpy(&and_array[GAP_SIZE_BYTES+1], &t, GAP_SIZE_BYTES);
		and_array_size += 2*GAP_SIZE_BYTES;
	}
	if (later_0 > 0) {
		memcpy(&and_array[and_array_size], &later_0, GAP_SIZE_BYTES);
		and_array_size += GAP_SIZE_BYTES;
	}

	return and_array_size;

}
////////////////////////////////
unsigned long count_size_of_bitmat(BitMat *bitmat)
{
	unsigned long size = 0;

	if (bitmat->bm.size() > 0) {
		for (std::list<struct row>::iterator it = bitmat->bm.begin(); it != bitmat->bm.end(); it++) {
//			unsigned int rowsize = 0;
			unsigned int rsz = 0;
			memcpy(&rsz, (*it).data, BM_ROW_SIZE);
			size += (TOTAL_ROWSIZE(rsz) + BM_ROW_SIZE) + sizeof((*it).rowid);
		}
	} else if (bitmat->vbm.size() > 0) {
		for (vector<struct row>::iterator it = bitmat->vbm.begin(); it != bitmat->vbm.end(); it++) {
//			unsigned int rowsize = 0;
			unsigned int rsz = 0;
			memcpy(&rsz, (*it).data, BM_ROW_SIZE);
			size += (TOTAL_ROWSIZE(rsz) + BM_ROW_SIZE + sizeof((*it).rowid));
		}

	}
	if (bitmat->rowfold != NULL) {
		size += sizeof(bitmat->row_bytes);
	}
//	if (bitmat->rowfold_prev != NULL) {
//		size += sizeof(bitmat->row_bytes);
//	}
	if (bitmat->colfold != NULL) {
		size += sizeof(bitmat->column_bytes);
	}
//	if (bitmat->colfold_prev != NULL) {
//		size += sizeof(bitmat->column_bytes);
//	}

	return size;
}

//////////////////////////////////////////
unsigned long get_offset(char *dumpfile, unsigned int bmnum)
{
	//Get the offset
	char tablefile[1024];
	sprintf(tablefile, "%s_table", dumpfile);
	unsigned long offset = 0;
#if MMAPFILES
	char *mmapfile = mmap_table[string(tablefile)];
	memcpy(&offset, &mmapfile[(bmnum-1)*table_col_bytes], table_col_bytes);
#else
	int fd = open(tablefile, O_RDONLY);
	if (fd < 0) {
		cout << "*** ERROR opening " << tablefile << endl;
		assert (0);
	}

	lseek(fd, (bmnum-1)*table_col_bytes, SEEK_CUR);
	unsigned char tablerow[table_col_bytes];
	read(fd, tablerow, table_col_bytes);
	memcpy(&offset, tablerow, table_col_bytes);
	close(fd);
#endif
	return offset;

}

bool mmap_all_files()
{
	int fd = 0; off_t size = 0; char *fstream; char tablefile[1024];

	fd = open((char *)config[string("BITMATDUMPFILE_SPO")].c_str(), O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[0] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[config[string("BITMATDUMPFILE_SPO")]] = fstream;
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_SPO")].c_str());
	fd = open(tablefile, O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[1] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[string(tablefile)] = fstream;

	fd = open((char *)config[string("BITMATDUMPFILE_OPS")].c_str(), O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[2] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[config[string("BITMATDUMPFILE_OPS")]] = fstream;
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_OPS")].c_str());
	fd = open(tablefile, O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[3] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[string(tablefile)] = fstream;

	fd = open((char *)config[string("BITMATDUMPFILE_PSO")].c_str(), O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[4] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[config[string("BITMATDUMPFILE_PSO")]] = fstream;
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
	fd = open(tablefile, O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[5] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[string(tablefile)] = fstream;


	fd = open((char *)config[string("BITMATDUMPFILE_POS")].c_str(), O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[6] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[config[string("BITMATDUMPFILE_POS")]] = fstream;
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
	fd = open(tablefile, O_RDONLY);
	assert(fd >= 0);
	size = lseek(fd, 0, SEEK_END);
	vectfd[7] = make_pair(fd, size);
	fstream = (char *) mmap(0, size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(fstream != (void *)-1);
	mmap_table[string(tablefile)] = fstream;

	return true;
}

bool munmap_all_files()
{

	char tablefile[1024];

	munmap(mmap_table[config[string("BITMATDUMPFILE_SPO")]], vectfd[0].second);
	close(vectfd[0].first);
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_SPO")].c_str());
	munmap(mmap_table[string(tablefile)], vectfd[1].second);
	close(vectfd[1].first);

	munmap(mmap_table[config[string("BITMATDUMPFILE_OPS")]], vectfd[2].second);
	close(vectfd[2].first);
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_OPS")].c_str());
	munmap(mmap_table[string(tablefile)], vectfd[3].second);
	close(vectfd[3].first);

	munmap(mmap_table[config[string("BITMATDUMPFILE_PSO")]], vectfd[4].second);
	close(vectfd[4].first);
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_PSO")].c_str());
	munmap(mmap_table[string(tablefile)], vectfd[5].second);
	close(vectfd[5].first);

	munmap(mmap_table[config[string("BITMATDUMPFILE_POS")]], vectfd[6].second);
	close(vectfd[6].first);
	sprintf(tablefile, "%s_table", (char *)config[string("BITMATDUMPFILE_POS")].c_str());
	munmap(mmap_table[string(tablefile)], vectfd[7].second);
	close(vectfd[7].first);

	return true;
}

void print_stats(char *fname_dump, unsigned int numbms, bool compfold, unsigned int numsubs, unsigned int numobjs)
{
	char tablefile[2048];
	sprintf(tablefile, "%s_table", fname_dump);

	int fd_table = open(tablefile, O_RDONLY);

	assert(fd_table != -1);

	unsigned long offset = 0;

	int fd = open(fname_dump, O_RDONLY);

	assert(fd != -1);

	unsigned int numtriples = 0;

	unsigned int row_bytes = (numsubs%8>0 ? numsubs/8+1 : numsubs/8);
	unsigned int column_bytes = (numobjs%8>0 ? numobjs/8+1 : numobjs/8);

	unsigned char rowfold[row_bytes];
	unsigned char colfold[column_bytes];

	memset(rowfold, 0, row_bytes);
	memset(colfold, 0, column_bytes);

	cout << "Subject_bytes " << row_bytes << endl;
	cout << "column_bytes " << column_bytes << endl;

	cout << fname_dump << " " << tablefile << endl;

	unsigned int pred_subs = 0, pred_objs = 0;
	unsigned long prev_offset = 0, size = 0;

	for (unsigned int i = 0; i < numbms; i++) {
//		offset = get_offset(fname_dump, i+1);
//		cout << "Offset is " << offset << " ";
		read(fd_table, &offset, table_col_bytes);
		size = offset - prev_offset;

		if (pred_subs != 0 && pred_objs != 0 && size != 0) {
			cout << i << " #Subjects " << pred_subs << " #Objects " << pred_objs << " #Size " << size << " #Triples " << numtriples << endl;
		}

		assert (lseek(fd, offset, SEEK_SET) != -1);
		//Moving over the numtriples field
		read(fd, &numtriples, sizeof(unsigned int));
//		cout << "***#triples " << numtriples << " ";
//		cout << "***Counting rows for bmnum " << i + 1 << " -- ";
		if (compfold) {
			unsigned int comp_arr_size = 0;
			read(fd, &comp_arr_size, ROW_SIZE_BYTES);
			unsigned char *comp_arr = (unsigned char *) malloc (comp_arr_size);
			dgap_uncompress(comp_arr, comp_arr_size, rowfold, row_bytes);
		} else {
			read(fd, rowfold, row_bytes);
			read(fd, colfold, column_bytes);
		}

		pred_subs = count_bits_in_row(rowfold, row_bytes);
		pred_objs = count_bits_in_row(colfold, column_bytes);

//		cout << i+1 << " #Subjects " << pred_subs << " #Objects " << pred_objs << " #Size " << size << " #Triples " << numtriples << endl;

		prev_offset = offset;
		offset = 0;
		//Now count set bits in rowfold
//		cout << "#Subjects " << count_bits_in_row(rowfold, row_bytes) << " #Objects " << count_bits_in_row(colfold, column_bytes) << endl;
	}

	offset = lseek(fd, 0L, SEEK_END);
	size = offset - prev_offset;

	cout << numbms << " #Subjects " << pred_subs << " #Objects " << pred_objs << " #Size " << size << " #Triples " << numtriples << endl;
}

unsigned int count_intersect(unsigned int bm1, unsigned int bm2)
{
	BitMat bitmat_spo1;
	init_bitmat(&bitmat_spo1, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
	BitMat bitmat_spo2;
	init_bitmat(&bitmat_spo2, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
	char *fpos1 = NULL; char *fpos2 = NULL;
	unsigned long offset1 = get_offset((char *)config[string("BITMATDUMPFILE_SPO")].c_str(), bm1);
	unsigned long offset2 = get_offset((char *)config[string("BITMATDUMPFILE_SPO")].c_str(), bm2);

#if MMAPFILES
	fpos1 = mmap_table[string((char *)config[string("BITMATDUMPFILE_SPO")].c_str())];
	fpos2 = fpos1;
	fpos1 += offset1;
	fpos2 += offset2;
#endif

#if MMAPFILES
	memcpy(&bitmat_spo1.num_triples, fpos1, sizeof(unsigned int));
	fpos1 += sizeof(unsigned int);
	memcpy(&bitmat_spo2.num_triples, fpos2, sizeof(unsigned int));
	fpos2 += sizeof(unsigned int);
#endif

#if MMAPFILES
	memcpy(bitmat_spo1.rowfold, fpos1, bitmat_spo1.row_bytes);
	fpos1 += bitmat_spo1.row_bytes;
	memcpy(bitmat_spo1.colfold, fpos1, bitmat_spo1.column_bytes);
	fpos1 += bitmat_spo1.column_bytes;

	memcpy(bitmat_spo2.rowfold, fpos2, bitmat_spo2.row_bytes);
	fpos2 += bitmat_spo2.row_bytes;
	memcpy(bitmat_spo2.colfold, fpos2, bitmat_spo2.column_bytes);
	fpos2 += bitmat_spo2.column_bytes;
#endif

	unsigned char res[bitmat_spo1.common_so_bytes];

	for (unsigned int i=0; i < bitmat_spo1.common_so_bytes; i++) {
		res[i] = bitmat_spo1.colfold[i] & bitmat_spo2.rowfold[i];
	}

	res[bitmat_spo1.common_so_bytes-1] &= (0xff << (8-(bitmat_spo1.num_comm_so%8)));

	return count_bits_in_row(res, bitmat_spo1.common_so_bytes);
}

unsigned int count_intersect_ss(unsigned int bm1, unsigned int bm2)
{
	BitMat bitmat_spo1;
	init_bitmat(&bitmat_spo1, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
	BitMat bitmat_spo2;
	init_bitmat(&bitmat_spo2, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
	char *fpos1 = NULL; char *fpos2 = NULL;
	unsigned long offset1 = get_offset((char *)config[string("BITMATDUMPFILE_SPO")].c_str(), bm1);
	unsigned long offset2 = get_offset((char *)config[string("BITMATDUMPFILE_SPO")].c_str(), bm2);

#if MMAPFILES
	fpos1 = mmap_table[string((char *)config[string("BITMATDUMPFILE_SPO")].c_str())];
	fpos2 = fpos1;
	fpos1 += offset1;
	fpos2 += offset2;
#endif

#if MMAPFILES
	memcpy(&bitmat_spo1.num_triples, fpos1, sizeof(unsigned int));
	fpos1 += sizeof(unsigned int);
	memcpy(&bitmat_spo2.num_triples, fpos2, sizeof(unsigned int));
	fpos2 += sizeof(unsigned int);
#endif

#if MMAPFILES
	memcpy(bitmat_spo1.rowfold, fpos1, bitmat_spo1.row_bytes);
	memcpy(bitmat_spo2.rowfold, fpos2, bitmat_spo2.row_bytes);
#endif

	unsigned char res[bitmat_spo1.row_bytes];

	for (unsigned int i=0; i < bitmat_spo1.row_bytes; i++) {
		res[i] = bitmat_spo1.rowfold[i] & bitmat_spo2.rowfold[i];
	}

	return count_bits_in_row(res, bitmat_spo1.row_bytes);
}

unsigned int count_intersect_oo(unsigned int bm1, unsigned int bm2)
{
	BitMat bitmat_spo1;
	init_bitmat(&bitmat_spo1, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
	BitMat bitmat_spo2;
	init_bitmat(&bitmat_spo2, gnum_subs, gnum_preds, gnum_objs, gnum_comm_so, SPO_BITMAT);
	char *fpos1 = NULL; char *fpos2 = NULL;
	unsigned long offset1 = get_offset((char *)config[string("BITMATDUMPFILE_SPO")].c_str(), bm1);
	unsigned long offset2 = get_offset((char *)config[string("BITMATDUMPFILE_SPO")].c_str(), bm2);

#if MMAPFILES
	fpos1 = mmap_table[string((char *)config[string("BITMATDUMPFILE_SPO")].c_str())];
	fpos2 = fpos1;
	fpos1 += offset1;
	fpos2 += offset2;
#endif

#if MMAPFILES
	memcpy(&bitmat_spo1.num_triples, fpos1, sizeof(unsigned int));
	fpos1 += sizeof(unsigned int);
	memcpy(&bitmat_spo2.num_triples, fpos2, sizeof(unsigned int));
	fpos2 += sizeof(unsigned int);
#endif

#if MMAPFILES
//	memcpy(bitmat_spo1.rowfold, fpos1, bitmat_spo1.row_bytes);
	fpos1 += bitmat_spo1.row_bytes;
	memcpy(bitmat_spo1.colfold, fpos1, bitmat_spo1.column_bytes);
	fpos1 += bitmat_spo1.column_bytes;

//	memcpy(bitmat_spo2.rowfold, fpos2, bitmat_spo2.row_bytes);
	fpos2 += bitmat_spo2.row_bytes;
	memcpy(bitmat_spo2.colfold, fpos2, bitmat_spo2.column_bytes);
	fpos2 += bitmat_spo2.column_bytes;
#endif

	unsigned char res[bitmat_spo1.column_bytes];

	for (unsigned int i=0; i < bitmat_spo1.column_bytes; i++) {
		res[i] = bitmat_spo1.colfold[i] & bitmat_spo2.colfold[i];
	}

	return count_bits_in_row(res, bitmat_spo1.column_bytes);
}

