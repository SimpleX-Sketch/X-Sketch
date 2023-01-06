#ifndef _X_H_
#define _X_H_
#include <bits/stdc++.h>
#include "Param.h"
#include "CorrectDetector.h"
#include "Common/hash.h"
#define inf 2147483647
#define eps 1e-6

const int max_range[] = {
	2, 14, 65534
};

const int counter_bit[] = {
	2, 4, 16
};

template<typename ID_TYPE>
class Slot {
public:
	Slot() : id(0) {
		memset(counter, 0, P * sizeof(uint32_t));
	}
	Slot(ID_TYPE _id) : id(_id) {
		memset(counter, 0, P * sizeof(uint32_t));
	}
	~Slot() {}
	void clear() {
		id = 0;
		start_window = 0;
		memset(counter, 0, P * sizeof(uint32_t));
	}
	void add(uint32_t win) {
		counter[win % P]++;
	}
	void insert(ID_TYPE ID, uint32_t win, uint32_t* c) {
		// assert(id == 0);
		id = ID;
		for (int k = S - 1, w = win; k >= 0; --k, --w) {
			counter[w % P] = c[k];
		}
		// assert(counter[win % P]);
		start_window = -1;
		for (int k = S - 1, w = win; k >= 0; --k, --w) {
			if (!c[k]) {
				start_window = w + 1;
				break;
			}
		}
		if (start_window == -1) 
			start_window = win - S + 1;
	}
	int query(uint32_t win, uint32_t* c) {
		for (int k = P - 1, w = win; k >= 0; --k, --w) {
			c[k] = counter[w % P];
			if (c[k] == 0) 
				return 0;
		}
		return 1;
	}
	ID_TYPE id;
	uint32_t counter[P];
	// int S = S_p;
	uint32_t start_window = 0;
};

template<typename ID_TYPE>
class Stage1 {
public: 
	Stage1(uint32_t memory) {
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

	Stage1(uint32_t memory, int S1Win){
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

	~Stage1() {
		for (int i = 0; i < NumOfWin; ++i) {
			for (int j = 0; j < 3; ++j) {
				delete[] TowerSketch[i][j];
			}
			delete[] TowerSketch[i];
		}
		delete[] TowerSketch;
	}
	void add(uint32_t k, uint32_t i, uint32_t j, uint32_t start, uint32_t end) {
		// increment bit [start, end) in j-th counter in i-th level of k-th TowerSketch
		uint32_t cur_bit = get(k, i, j, start, end);
		if (cur_bit > max_range[i])
			return;
		cur_bit++;
		TowerSketch[k][i][j] &= (~(((1 << (end - start)) - 1) << start));
		TowerSketch[k][i][j] |= (cur_bit << start);
		// assert(cur_bit == get(k, i, j, start, end));
	}
	uint32_t get(uint32_t k, uint32_t i, uint32_t j, uint32_t start, uint32_t end) {
		return (TowerSketch[k][i][j] & (((1 << (end - start)) - 1) << start)) >> start;
	}
	void insert(ID_TYPE id, uint32_t win) {
		// if(win == 0)
			// std::cout<<"S: "<<NumOfWin<<std::endl;

		for (int i = 0; i < 3; ++i) {
			uint32_t index = hash(id, i) % size[i];
			uint32_t cell = index * counter_bit[i] / 32;
			uint32_t res = index - cell * 32 / counter_bit[i];
			add(win % NumOfWin, i, cell, res, res + counter_bit[i]);
		}
	}
	int query(ID_TYPE id, uint32_t win, uint32_t* c) {
		win = win % NumOfWin;
		for (int k = NumOfWin - 1; k >= 0; --k, win = (win + NumOfWin - 1) % NumOfWin) {
			c[k] = inf;
			for (int i = 0; i < 3; ++i) {
				uint32_t index = hash(id, i) % size[i];
				uint32_t cell = index * counter_bit[i] / 32;
				uint32_t res = index - cell * 32 / counter_bit[i];
				uint32_t temp = get(win, i, cell, res, res + counter_bit[i]);
				if (temp <= max_range[i]) 
					c[k] = MIN(c[k], temp);
				
			}
			if (c[k] == 0)
				return 0;
			// assert(c[k] <= window_size);
		}
		return 1;
	}
	void clear(uint32_t win) {
		for (int i = 0; i < 3; ++i) {
			memset(TowerSketch[win % NumOfWin][i], 0, cell_num * sizeof(uint32_t));
		}
	}
private: 
	uint32_t*** TowerSketch;
	uint32_t size[3];
	uint32_t cell_num;
	int NumOfWin = S;
};

template<typename ID_TYPE>
class Stage2 {
public: 
	Stage2() {}
	Stage2(uint32_t memory) {
		size = memory * 1024 / sizeof(Slot<ID_TYPE>) / bucket_size;
		slot = new Slot<ID_TYPE>* [size];
		for (int i = 0; i < size; ++i) {
			slot[i] = new Slot<ID_TYPE> [bucket_size];
			// memset(slot[i], 0, bucket_size * sizeof(Slot<ID_TYPE>));
		}
		oFile.open("./ttmp2.csv", std::ios::app);
		oFile<<"X"<<std::endl;

	}
	Stage2(uint32_t memory, double var_input, double error_input, int CellNum, double potential) {
		var_thres = var_input;
		error_thres = error_input;
		bucket_size = CellNum;
		potential_thres = potential;

		size = memory * 1024 / sizeof(Slot<ID_TYPE>) / bucket_size;
		std::cout<<"size: "<<size<<std::endl;
		slot = new Slot<ID_TYPE>* [size];
		for (int i = 0; i < size; ++i) {
			slot[i] = new Slot<ID_TYPE> [bucket_size];
			// memset(slot[i], 0, bucket_size * sizeof(Slot<ID_TYPE>));
		}
		oFile.open("./ttmp2.csv", std::ios::app);
		oFile<<"X"<<std::endl;
		
	}

	~Stage2() {
		for (int i = 0; i < size; ++i) 
			delete slot[i];
		delete slot;
		oFile.close();
	}
	int query(ID_TYPE id) {
		// return which cell id is in, -1 if does not exist
		uint32_t index = hash(id, 33) % size;
		for (int i = 0; i < bucket_size; ++i) {
			if (slot[index][i].id == id)
				return i;
		}
		return -1;
	}
	void insert(ID_TYPE id, uint32_t win) {
		uint32_t index = hash(id, 33) % size;
		int i = query(id);
		// assert(i >= 0);
		slot[index][i].add(win);
	}
	void clear(uint32_t win) {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < bucket_size; ++j) {
				slot[i][j].counter[win % P] = 0;
			}
		}
	}
	void push(ID_TYPE id, uint32_t win, uint32_t* c) {

		// assert(win >= S - 1);
		uint32_t index = hash(id, 33) % size;
		// assert(query(id) < 0);
		int min_counter = inf, min_index = -1;
		for (int i = 0; i < bucket_size; ++i) {
			if (slot[index][i].id == 0) {
				min_index = i;
				break;
			}
		}
		// empty cell
		if (min_index >= 0) {
			slot[index][min_index].insert(id, win, c);
			return;
		}
		double l[K + 1] = {}, error = 0;
		double temp[S] = {};
		double min_weight = inf, id_weight;
		for (int j = 0; j < S; ++j)
			temp[j] = c[j];
		linear_regressing_try(temp, l);
		for (int j = 0; j < S; ++j) {
			double t = temp[j];
			for (int k = 0; k <= K; ++k)
				t -= l[k] * pow(j, k);
			error += t * t;
		}
		id_weight = abs(l[K]) / ((error + eps)*(error + eps));
		

		if (id_weight < potential_thres)
			return;
		int min_start_window = inf;
		min_index = -1;
		for (int i = 0; i < bucket_size; ++i) {
			if (slot[index][i].start_window < min_start_window) {
				min_start_window = slot[index][i].start_window;
				min_index = i;
			}
		}
		// assert(min_index >= 0);
		if (win != min_start_window && rand() % (win - min_start_window) == 0) {
			if (win >= min_start_window + P) {
				report_top_k.push_back(Report_Slot(slot[index][min_index].id, slot[index][min_index].start_window, win - 1));
			}
			slot[index][min_index].clear();
			slot[index][min_index].insert(id, win, c);
		}
	}
	void check(uint32_t win) {
		/*****/
		if(win == 0) {
			std::cout<<"var_thres: "<<var_thres<<std::endl;
		}
		/******/
	#ifdef PREDICT_MODE
		std::map<ID_TYPE, uint32_t> predict_result;
	#endif
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < bucket_size; ++j) {
				if (slot[i][j].id == 0)
					continue;
				if (slot[i][j].counter[win % P] == 0) {
					if (win >= P + slot[i][j].start_window) {
						assert(std::find(result.begin(), result.end(), std::make_pair(slot[i][j].id, win - 1)) != result.end());
						report_top_k.push_back(Report_Slot(slot[i][j].id, slot[i][j].start_window, win - 1));
					}
					slot[i][j].clear();
					continue;
				}
				if (win < P - 1)
					continue;
				uint32_t c[P] = {};
				double y[P] = {}, z[P] = {}, b[K + 1] = {};
				if (slot[i][j].query(win, c)) {
					// linear regressing
					for (int k = 0; k < P; ++k) 
						y[k] = c[k];
					linear_regressing(y, b);
					if (abs(b[K]) < var_thres) {
						if (win >= P + slot[i][j].start_window) {
							// slot[i][j].print(win - 1);
							report_top_k.push_back(Report_Slot(slot[i][j].id, slot[i][j].start_window, win - 1));
						}
						slot[i][j].start_window = win - P + 2;
						continue;
					}
					double error = 0;
					for (int k = 0; k < P; ++k) {
						z[k] = y[k];
						for (int l = 0; l <= K; ++l) 
							z[k] -= b[l] * pow(k, l);
						error += z[k] * z[k];
					}
					if (error / P <= error_thres) {
						result.emplace_back(std::make_pair(slot[i][j].id, win));
						#ifdef PREDICT_MODE
							double pred_res = 0.0;
							for (int k = 0; k <= K; ++k) {
								pred_res += b[k] * pow(P, k);
							} 
							if (pred_res < 0) 
								predict_result[slot[i][j].id] = 0;
							else 
								predict_result[slot[i][j].id] = (int)(pred_res + 0.5);
						#endif
					}
					else {
						if (win >= P + slot[i][j].start_window) {
							report_top_k.push_back(Report_Slot(slot[i][j].id, slot[i][j].start_window, win - 1));
						}
						slot[i][j].start_window = win - P + 2;
					}
				}
			}
		}
		#ifdef PREDICT_MODE
			predict[win + 1] = predict_result;
		#endif
	}
	void get_all(uint32_t win) {
		for (int i = 0 ; i < size; ++i) {
			for (int j = 0; j < bucket_size; ++j) {
				if (slot[i][j].id == 0)
					continue;
				if (win - slot[i][j].start_window >= P)
					report_top_k.push_back(Report_Slot(slot[i][j].id, slot[i][j].start_window, win - 1));
			}
		}
	}
private: 
	Slot<ID_TYPE>** slot;
	uint32_t size;
	int bucket_size = bucket_size_p;
	double error_thres = error_thres_p;
	double var_thres = var_thres_p;
	// int S = S_p;
	double potential_thres = potential_thres_p;
	std::ofstream oFile;
public:
	std::vector<std::pair<ID_TYPE, uint32_t>> result;
	std::vector<Report_Slot<ID_TYPE>> report_top_k;
#ifdef PREDICT_MODE
	std::map<uint32_t, std::map<ID_TYPE, uint32_t>> predict;
#endif
};

template<typename ID_TYPE>
class XSketch {
public: 
	XSketch() {}
	XSketch(uint32_t memory) : win_cnt(0), last_timestamp(0) {

		double stage1_mem = memory * stage_ratio;
		// double stage1_mem = memory * S;
		double stage2_mem = memory * (1 - stage_ratio);
		// double stage2_mem = 100;
		stage1 = new Stage1<ID_TYPE>(stage1_mem);
		stage2 = new Stage2<ID_TYPE>(stage2_mem);
	}
	XSketch(uint32_t memory, double var_input, double error_input, double ratio_input, int CellNum, int S1Wins, double potential):win_cnt(0),last_timestamp(0){
		stage_ratio = ratio_input;
		var_thers = var_input;
		error_thres = error_input;
		bucket_size = CellNum;
		double stage1_mem = memory * stage_ratio;
		// double stage1_mem = memory * S;
		double stage2_mem = memory * (1 - stage_ratio);
		// double stage2_mem = 100;
		stage1 = new Stage1<ID_TYPE>(stage1_mem, S1Wins);
		stage2 = new Stage2<ID_TYPE>(stage2_mem, var_input,error_input,CellNum, potential);

	}
	~XSketch() {
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
		// if(win_cnt == 0)
			// std::cout<<"stage ratio: "<<stage_ratio<<std::endl;
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
#ifdef PREDICT_MODE
	std::map<uint32_t, std::map<ID_TYPE, uint32_t>> predict() {
		return stage2->predict;
	}
#endif
private: 
	uint32_t win_cnt;
	uint32_t last_timestamp;
	Stage1<ID_TYPE>* stage1;
	Stage2<ID_TYPE>* stage2;
	double error_thres = error_thres_p;
	double var_thers = var_thres_p;
	double stage_ratio = stage_ratio_p;
	int bucket_size = bucket_size_p;
	int S = S_p;
};


#endif