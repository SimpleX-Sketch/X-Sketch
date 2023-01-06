#ifndef _PARAM_H_
#define _PARAM_H_
#include <bits/stdc++.h>
#include "Common/Matrix.h"
#define PREDICT_MODE

const uint32_t window_size = 10000;

// memory
const int memory_p = 300;

// memory ratio of stage 1 to stage 2
const double stage_ratio_p = 0.8;

// the minimum value for a_K
const double var_thres_p = 0.1;

// number of consecutive windows
const int P = 7;

// degree of polynomial
const int K = 0;

// number of recorded windows in stage 1
const int S = 4;

const int S_p = S;

// threshold for mean square error
const double error_thres_p = 0.5;

// number of cells in each bucket
const int bucket_size_p = 4;

// threshold for potential
const double potential_thres_p = 0.25;

double calcu_matrix[K + 1][P] = {};
double calcu_matrix_try[K + 1][S_p] = {};

template<typename ID_TYPE>
class Report_Slot {
public:
	ID_TYPE id;
	uint32_t start_window;
	uint32_t end_window;
	Report_Slot() {}
	Report_Slot(ID_TYPE _id, uint32_t st, uint32_t et): id(_id), start_window(st), end_window(et) {}
	~Report_Slot() {}
	bool operator < (const Report_Slot &r) {
		if (end_window - start_window == r.end_window - r.start_window) {
			if (start_window == r.start_window)
				return id < r.id;
			return start_window < r.start_window;
		}
		return end_window - start_window < r.end_window - r.start_window;
	}
};

void init_matrix() {
	// by normal equation, the LSE satisfies X'Xb = X'Y, hence b = (X'X)^{-1}X'Y
	// to get b for an arbitrary Y, we only need to find (X'X)^{-1}X'
	double X[P][K + 1] = {}, XX[K + 1][K + 1] = {}, temp[(K + 1) * (K + 1)];
	for (int i = 0; i <= P - 1; ++i) {
		for (int j = 0; j <= K; ++j) {
			X[i][j] = pow(i, j);
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			for (int k = 0; k <= P - 1; ++k) {
				XX[i][j] += X[k][i] * X[k][j];
			}
		}
	}

	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			temp[i * (K + 1) + j] = XX[i][j];
		}
	}
	// find the inverse matrix of X'X
	inverse(K + 1, temp);

	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			XX[i][j] = temp[i * (K + 1) + j];
		}
	}

	// get (X'X)^{-1}X'
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= P - 1; ++j) {
			for (int k = 0; k <= K; ++k) {
				calcu_matrix[i][j] += XX[i][k] * X[j][k];
			}
		}
	}


	double Z[S_p][K + 1] = {};
	for (int i = 0; i <= S_p - 1; ++i) {
		for (int j = 0; j <= K; ++j) {
			X[i][j] = pow(i, j);
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			for (int k = 0; k <= S_p - 1; ++k) {
				XX[i][j] += X[k][i] * X[k][j];
			}
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			temp[i * (K + 1) + j] = XX[i][j];
		}
	}
	inverse(K + 1, temp);
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			XX[i][j] = temp[i * (K + 1) + j];
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			XX[i][j] = temp[i * (K + 1) + j];
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= S_p - 1; ++j) {
			for (int k = 0; k <= K; ++k) {
				calcu_matrix_try[i][j] += XX[i][k] * Z[j][k];
			}
		}
	}
}

void linear_regressing(double* y, double* b) {
	memset(b, 0, (K + 1) * sizeof(double));
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= P - 1; ++j) {
			b[i] += calcu_matrix[i][j] * y[j];
		}
	}
}

void calcu_variation(uint32_t c[], int d[], uint32_t size, uint32_t k) {
	int e[size];
	memcpy(e, c, size * sizeof(uint32_t));
	while (k--) {
		for (int i = 0; i < size - 1; ++i)
			e[i] = e[i + 1] - e[i];
		--size;
	}
	memcpy(d, e, size * sizeof(int));
}

void linear_regressing_try(double* y, double* b) {
	// only use S number of windows for linear regressing
	memset(b, 0, (K + 1) * sizeof(double));
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= S_p - 1; ++j) {
			b[i] += calcu_matrix[i][j] * y[j];
		}
	}
}

#endif
