#ifndef STRAWMAN_H
#define STRAWMAN_H

#include <bits/stdc++.h>
#include "Common/hash.h"
#include "Param.h"


template<typename ID_TYPE>
class CountMinSketch{
private:
	int _NumOfArray = 3;
	uint32_t **CMSketch = {};
	int HashSeed = 100;
	uint32_t CounterOfArray = 0;
	

public:
	CountMinSketch(int memory){
		CounterOfArray = memory * 1024 / (_NumOfArray * sizeof(uint32_t));
		// std::cout<<"counter of each array: "<<CounterOfArray<<std::endl;
		CMSketch = new uint32_t*[_NumOfArray];
		for(int i = 0; i < _NumOfArray; ++i){
			CMSketch[i] = new uint32_t[CounterOfArray]();
		}

	}
	~CountMinSketch(){
		for(int i= 0; i < _NumOfArray; ++i){
			delete[] CMSketch[i];
		}
		delete[] CMSketch;
	}

	int insert(ID_TYPE ItemID){
		for(int i = 0; i < _NumOfArray; ++i){
			uint32_t HashValue = hash(ItemID,i);
			// uint32_t HashValue = MurmurHash32((const void*)ItemID.data(),KEY_LEN, 50 + i);
			uint32_t loc = HashValue % CounterOfArray;
			// std::cout<<loc<<","<<CounterOfArray<<std::endl;
			CMSketch[i][loc] ++;
		}
		return 1;
	}

	uint32_t query(ID_TYPE ItemID){
		uint32_t freq_min = UINT32_MAX;
		for(int i = 0; i < _NumOfArray; ++i){
			uint32_t loc = hash(ItemID, i) % CounterOfArray;
			// uint32_t loc = MurmurHash32((const void*)ItemID.data(), KEY_LEN, 50 + i) % CounterOfArray;
			uint32_t freq_loc = CMSketch[i][loc];
			freq_min = std::min(freq_min,freq_loc);
		}
		return freq_min;
	}

	void clear(){
		for(int i = 0; i < _NumOfArray; ++i){
			memset(CMSketch[i], 0, CounterOfArray * sizeof(uint32_t));
		}
	}

};

template<typename ID_TYPE, typename RESULT_ID>
class Baseline{
private:
	CountMinSketch<ID_TYPE>* Basic[P];
	std::set<ID_TYPE> X_Queue = {};
	uint32_t last_timestamp;
	uint32_t win_cnt;
	std::vector<std::pair<RESULT_ID, uint32_t>> result;
	std::map<ID_TYPE, uint32_t> now_window;
	std::map<ID_TYPE, uint32_t> last_window;
	std::vector<Report_Slot<ID_TYPE>> report_top_k;

	uint64_t debugID;
	uint32_t debugWin;

	double var_thres = var_thres_p;
	double error_thres = error_thres_p;

public:
	Baseline(int memory):win_cnt(0), last_timestamp(0){
		int MemPerSketch = memory / P;
		for(int i = 0; i < P; ++i)
		 	Basic[i] = new CountMinSketch<ID_TYPE>(MemPerSketch);

	}
	Baseline(int memory, double var_input, double error_input):win_cnt(0), last_timestamp(0){
		var_thres = var_input;
		error_thres = error_input;
		
		int MemPerSketch = memory / P;
		for(int i = 0; i < P; ++i)
		 	Basic[i] = new CountMinSketch<ID_TYPE>(MemPerSketch);

	}
	~Baseline(){
		for(int i = 0; i < P; ++i){
			delete[] Basic[i];
		}
	}

	void insert(ID_TYPE id, uint32_t timestamp){
		if(last_timestamp + window_size < timestamp){//时间切换
			transition();
		}
		Basic[win_cnt % P]->insert(id);
		X_Queue.insert(id);
	}

	void traverse_query(){
		last_window = now_window;
		now_window.clear();
		uint32_t freq[P] = {};
		for(auto i = X_Queue.begin(); i != X_Queue.end(); i++){
			double y[P] = {}, b[K + 1] = {}, z[P] = {};
			int flag = 0;
			ID_TYPE ItemID = *i;
			for(int j = P - 1, k = win_cnt; j >= 0; j--, k--){
				// freq[j] = Basic[j]->query(&ItemID);
				y[j] = Basic[(P + k) % P]->query(ItemID);
				z[j] = y[j];
				if(y[j] == 0){
					flag = 1;
				}
			}

			//debug
			// if(debugID == ItemID && debugWin == win_cnt){
			// 	std::cout<<"Basic Solution : find debug id in time:"<<win_cnt<<std::endl;
			// 	for(int t = 0; t < P; t++){
			// 		std::cout<<int(y[t])<<" ";
			// 	}
			// 	std::cout<<std::endl;
			// }

			if(flag){
				continue;
			}
			linear_regressing(z,b);

			// if(debugID == ItemID && debugWin == win_cnt){
			// 	std::cout<<"Basic: b[k] is:"<<abs(b[K])<<" var_thres: "<<var_thres <<" id is: " << ItemID<<" in time :"<<win_cnt<<std::endl;
			// }

			if(abs(b[K]) < var_thres){
				continue;
			}

			double error = 0;
			for(int j = 0; j < P; ++j){
				z[j] = y[j];
				for(int k = 0; k <= K; ++k){
					z[j] -= pow(j,k) * b[k];
				}
				error += z[j] * z[j];
			}

			// if(debugID == ItemID && debugWin == win_cnt){
			// 	std::cout<<"Basic: error is:"<<error <<" error_thres: "<<error_thres<<std::endl;
			// }

			if (error / P <= error_thres){
				result.emplace_back(std::make_pair(ItemID, win_cnt));
				if(last_window.find(ItemID) != last_window.end()){
					now_window[ItemID] = last_window[ItemID];
					last_window.erase(ItemID);
				}
				else{
					now_window[ItemID] = win_cnt - P + 1;
				}
				// if(debugID == ItemID && debugWin == win_cnt){
				// 	std::cout<<"Basic Solution: we have insert the item into result, id is:"<< ItemID<<" in time :"<<win_cnt<<std::endl;
				// }
			}
		}
		for(auto i : last_window){
			report_top_k.push_back(Report_Slot(i.first, i.second, win_cnt - 1));
		}
	}

	void transition(){
		traverse_query();
		win_cnt ++;
		Basic[win_cnt % P]->clear();
		X_Queue.clear();
		last_timestamp += window_size;
	}

	std::vector<std::pair<RESULT_ID, uint32_t>> query(){
		traverse_query();
		return result;
	}

	std::vector<Report_Slot<ID_TYPE>> report(){
		std::sort(report_top_k.begin(), report_top_k.end());
		reverse(report_top_k.begin(), report_top_k.end());
		return report_top_k;
	}
};

#endif