#ifndef _CORRECT_H_
#define _CORRECT_H_
#include <bits/stdc++.h>
#include "Param.h"
#define NEXT(x) (((x) + 1) % P)
#define PREV(x) (((x) + P - 1) % P)
#define NOW(x) ((x) % P)

template<typename ID_TYPE>
class CorrectItem {
public: 
	CorrectItem() {
		id = 0;
		memset(counter, 0, P * sizeof(uint32_t));
	}
	CorrectItem(ID_TYPE _id) : id(_id) {
		memset(counter, 0, P * sizeof(uint32_t));
	}
	void add(uint32_t win) {
		counter[NOW(win)]++;
		for (int i = 0; i < P; ++i)
			assert(counter[i] <= 40000);
	} 
	void remove(uint32_t win) {
		counter[NOW(win)] = 0;
	}
	void set(ID_TYPE ID, uint32_t win, uint32_t* c) {
		id = ID;
		for (int w = win, k = K; k >= 0 && w >= 0; --w, --k) {
			counter[NOW(w)] = c[k];
		}
		for (int i = 0; i < P; ++i)
			assert(counter[i] <= 40000);
	}
	void clear() {
		id = 0;
		memset(counter, 0, P * sizeof(uint32_t));
	}
	
	ID_TYPE id;
	uint32_t counter[P];
};

template<typename ID_TYPE>
class CorrectDetector {
public: 
	CorrectDetector() : win_cnt(0), last_timestamp(0) {}
	CorrectDetector(double error_input, double var_input) : win_cnt(0), last_timestamp(0) {
		error_thres = error_input;
		var_thres = var_input;
	}
	~CorrectDetector() {}
	void insert(ID_TYPE id, uint32_t timestamp) {
		if (timestamp > last_timestamp + window_size) {
			transition();
		}
		if (id_map.find(id) == id_map.end()) {
			id_map[id] = detect.size();
			detect.emplace_back(CorrectItem(id));
		}
		assert(id == detect[id_map[id]].id);
		detect[id_map[id]].add(win_cnt);
	#ifdef PREDICT_MODE
		history[id][win_cnt]++;
	#endif
	}
	void transition() {
		last_window = now_window;
		now_window.clear();
		for (auto &i : detect) {
			if (win_cnt < P - 1) {
				i.remove(win_cnt + 1);
				continue;
			}
			double y[P] = {}, b[K + 1] = {}, z[P] = {};
			int flag = 0;
			for (int j = P - 1, k = win_cnt; j >= 0; j--, k--) {
				y[j] = i.counter[NOW(k)];
				z[j] = y[j];
				if (y[j] == 0) {
					flag = 1;
				}
			}
			if (flag) {
				i.remove(win_cnt + 1);
				continue;
			}
			linear_regressing(z, b);
			if (abs(b[K]) < var_thres) {
				i.remove(win_cnt + 1);
				continue;
			}
			double error = 0;
			for (int j = 0; j < P; ++j) {
				z[j] = y[j];
				for (int k = 0; k <= K; ++k) {
					z[j] -= pow(j, k) * b[k];
				}
				error += z[j] * z[j];
			}
			if (error / P <= error_thres) {
				result.emplace_back(std::make_pair(i.id, win_cnt));
				if (last_window.find(i.id) != last_window.end()) {
					now_window[i.id] = last_window[i.id];
					last_window.erase(i.id);
				}
				else {
					now_window[i.id] = win_cnt - P + 1;
				}
			}
			i.remove(win_cnt + 1);
		}
		for (auto i : last_window) {
			report_top_k.push_back(Report_Slot(i.first, i.second, win_cnt - 1));
		}
		win_cnt++;
		last_timestamp += window_size;
	}
	std::vector<std::pair<ID_TYPE, uint32_t>> query() {
		transition();
		return result;
	}
	std::vector<Report_Slot<ID_TYPE>> report() {
		for (auto i : now_window) {
			report_top_k.push_back(Report_Slot(i.first, i.second, win_cnt - 1));
		}
		std::sort(report_top_k.begin(), report_top_k.end());
		reverse(report_top_k.begin(), report_top_k.end());
		return report_top_k;
	}
#ifdef PREDICT_MODE
	std::map<ID_TYPE, std::map<uint32_t, uint32_t>> get_history() {
		return history;
	}
#endif
private:
	uint32_t win_cnt;
	uint32_t last_timestamp;
	std::map<ID_TYPE, uint32_t> id_map;
	std::vector<CorrectItem<ID_TYPE>> detect;
	std::vector<std::pair<ID_TYPE, uint32_t>> result;
	std::vector<Report_Slot<ID_TYPE>> report_top_k;
	std::map<ID_TYPE, uint32_t> now_window;
	std::map<ID_TYPE, uint32_t> last_window; 
#ifdef PREDICT_MODE
	std::map<ID_TYPE, std::map<uint32_t, uint32_t>> history;
#endif
	double var_thres = var_thres_p;
	double error_thres = error_thres_p;
};


#endif