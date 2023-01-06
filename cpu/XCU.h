#ifndef _XCU_H_
#define _XCU_H_
#include <bits/stdc++.h>
#include "Param.h"
#include "CorrectDetector.h"
#include "hash.h"
#include "XSketch.h"

#define inf 2147483647
#define eps 1e-6

// const int max_range[] = {
// 	2, 14, 65534
// };

// const int counter_bit[] = {
// 	2, 4, 16
// };


template<typename ID_TYPE>
class Stage1_CU {
public: 
	Stage1_CU(uint32_t memory) {
		cell_num = memory * 1024 / sizeof(uint32_t) / 3 / S;
		for (int i = 0; i < 3; ++i) {
			size[i] = cell_num * 32 / counter_bit[i];
		}
		TowerSketch = new uint32_t** [S];
		for (int i = 0; i < S; ++i) {
			TowerSketch[i] = new uint32_t* [3];
			for (int j = 0; j < 3; ++j) {
				TowerSketch[i][j] = new uint32_t [cell_num];
				memset(TowerSketch[i][j], 0, cell_num * sizeof(uint32_t));
			}
		}
		std::cout<<"in CU"<<std::endl;
	}
	Stage1_CU(uint32_t memory, int S1Win){
		NumOfWin = S1Win;
		
		cell_num = memory * 1024 / sizeof(uint32_t) / 3 / NumOfWin;
		for (int i = 0; i < 3; ++i) {
			size[i] = cell_num * 32 / counter_bit[i];
		}
		TowerSketch = new uint32_t** [NumOfWin];
		for (int i = 0; i < NumOfWin; ++i) {
			TowerSketch[i] = new uint32_t* [3];
			for (int j = 0; j < 3; ++j) {
				TowerSketch[i][j] = new uint32_t [cell_num];
				memset(TowerSketch[i][j], 0, cell_num * sizeof(uint32_t));
			}
		}
	}
	~Stage1_CU() {
		for (int i = 0; i < S; ++i) {
			for (int j = 0; j < 3; ++j) {
				delete[] TowerSketch[i][j];
			}
			delete[] TowerSketch[i];
		}
		delete[] TowerSketch;
	}
	void add(uint32_t k, uint32_t i, uint32_t j, uint32_t start, uint32_t end){
        uint32_t cur_bit = get(k, i, j, start, end);
		if (cur_bit > max_range[i])
			return;
		cur_bit++;
        TowerSketch[k][i][j] &= (~(((1 << (end - start)) - 1) << start));
		TowerSketch[k][i][j] |= (cur_bit << start);
		assert(cur_bit == get(k, i, j, start, end));
    }

    uint32_t get(uint32_t k, uint32_t i, uint32_t j, uint32_t start, uint32_t end) {
		return (TowerSketch[k][i][j] & (((1 << (end - start)) - 1) << start)) >> start;
	}
	
   void insert(ID_TYPE id, uint32_t win) {
        // std::cout<<"insert in CU"<<std::endl;
        uint32_t freq[_NumOfArray] = {};
        uint32_t array_cell[_NumOfArray];
        uint32_t array_res[_NumOfArray];
		for (int i = 0; i < _NumOfArray; ++i) {
			uint32_t index = hash(id, i) % size[i];
			uint32_t cell = index * counter_bit[i] / 32;
			uint32_t res = index - cell * 32 / counter_bit[i];

            freq[i] =get(win % S, i, cell, res, res + counter_bit[i]);
            if (freq[i] > max_range[i]) 
                freq[i] = UINT32_MAX;
            array_cell[i] = cell;
            array_res[i] = res;
		}
        uint32_t cur_bit = *std::min_element(freq, freq + _NumOfArray);

        for(int i = 0; i < _NumOfArray; i++){
            if(freq[i] == cur_bit)
                add(win % S, i, array_cell[i], array_res[i], array_res[i] + counter_bit[i]);
        }

	}
	int query(ID_TYPE id, uint32_t win, uint32_t* c) {
		win = win % S;
		for (int k = S - 1; k >= 0; --k, win = (win + S - 1) % S) {
			c[k] = inf;
			for (int i = 0; i < 3; ++i) {
				uint32_t index = hash(id, i) % size[i];
				uint32_t cell = index * counter_bit[i] / 32;
				uint32_t res = index - cell * 32 / counter_bit[i];
				uint32_t temp = get(win, i, cell, res, res + counter_bit[i]);
				if (temp <= max_range[i]) 
					c[k] = MIN(c[k], temp);
			}
			assert(c[k] <= window_size);
			if (c[k] == 0)
				return 0;
		}
		return 1;
	}
	void clear(uint32_t win) {
		for (int i = 0; i < 3; ++i) {
			memset(TowerSketch[win % S][i], 0, cell_num * sizeof(uint32_t));
		}
	}
protected: 
	uint32_t*** TowerSketch;
	uint32_t size[3];
	uint32_t cell_num;
	int _NumOfArray = 3;
	int NumOfWin = S;
};



template<typename ID_TYPE>
class XSketch_CU {
public: 
	XSketch_CU() {}
	XSketch_CU(uint32_t memory) : win_cnt(0), last_timestamp(0) {

		double stage1_mem = memory * stage_ratio;
		// double stage1_mem = 100000;
		double stage2_mem = memory * (1 - stage_ratio);
		// double stage2_mem = 200;
		stage1 = new Stage1_CU<ID_TYPE>(stage1_mem);
		stage2 = new Stage2<ID_TYPE>(stage2_mem);
	}
	XSketch_CU(uint32_t memory, double var_input, double error_input, double ratio_input, int CellNum, int S1Wins, double potential):win_cnt(0),last_timestamp(0){
		stage_ratio = ratio_input;
		var_thers = var_input;
		error_thres = error_input;
		bucket_size = CellNum;
		double stage1_mem = memory * stage_ratio;
		// double stage1_mem = 100000 ;
		double stage2_mem = memory * (1 - stage_ratio);
		// double stage2_mem = 200;
		stage1 = new Stage1_CU<ID_TYPE>(stage1_mem, S1Wins);
		stage2 = new Stage2<ID_TYPE>(stage2_mem, var_input,error_input,CellNum, potential);

	}
	~XSketch_CU() {
		delete stage1;
		delete stage2;
	}
	void insert(ID_TYPE id, uint32_t timestamp) {
		if (last_timestamp + window_size < timestamp) 
			transition();
		if (stage2->query(id) >= 0) {
			stage2->insert(id, win_cnt);
			return;
		}
		stage1->insert(id, win_cnt);
		if (win_cnt < S - 1) 
		 	return;
		uint32_t c[S] = {};
		if (stage1->query(id, win_cnt, c))
		// stage1->query(id, win_cnt, c);
			stage2->push(id, win_cnt, c);
	}
	void transition() {
		// if(win_cnt == 0){
		// 	std::cout<<"stage ratio: "<<stage_ratio<<std::endl;
		// }
		stage2->check(win_cnt);
		win_cnt++;
		stage1->clear(win_cnt);
		stage2->clear(win_cnt);
		last_timestamp += window_size;
	}
	std::vector<std::pair<ID_TYPE, uint32_t>> query() {
		transition();
		stage2->get_all(win_cnt);
		return stage2->result;
	}
	std::vector<Report_Slot<ID_TYPE>> report() {
		std::sort(stage2->report_top_k.begin(), stage2->report_top_k.end());
		reverse(stage2->report_top_k.begin(), stage2->report_top_k.end());
		return stage2->report_top_k;
	}
private: 
	uint32_t win_cnt;
	uint32_t last_timestamp;
	Stage1_CU<ID_TYPE>* stage1;
	Stage2<ID_TYPE>* stage2;
	double error_thres = error_thres_p;
	double var_thers = var_thres_p;
	double stage_ratio = stage_ratio_p;
	int bucket_size = bucket_size_p;
	int S = S_p;
};



#endif