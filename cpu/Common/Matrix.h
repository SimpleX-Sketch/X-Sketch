#ifndef _MATRIX_H_
#define _MATRIX_H_

double det(int n, double *aa) {
	if (n == 1)
		return aa[0];
	double *bb = new double[(n - 1) * (n - 1)];
	int mov = 0;
	double sum = 0.0;
	for (int arow = 0; arow < n; arow++) {
		for (int brow = 0; brow < n - 1; brow++) {
			mov = arow > brow ? 0 : 1;
			for (int j = 0; j < n - 1; j++) {
				bb[brow * (n - 1) + j] = aa[(brow + mov) * n + j + 1];
			}
		}
		int flag = (arow % 2 == 0 ? 1: -1);
		sum += flag * aa[arow*n] * det(n - 1, bb);
	}
	delete[] bb;
	return sum;
}
 
void inverse(int n, double* aa) {
	if (n == 1) {
		aa[0] = 1 / aa[0];
		return;
	}
	double det_aa = det(n, aa);
	assert(det_aa != 0);
	double *adjoint = new double[n * n];
	double *bb = new double[(n - 1) * (n - 1)];
	int pi, pj, q;
	for (int ai = 0; ai < n; ai++) {
		for (int aj = 0; aj < n; aj++) {
			for (int bi = 0; bi < n - 1; bi++) {
				for (int bj = 0; bj < n - 1; bj++) {
					if (ai > bi)
						pi = 0;
					else
						pi = 1; 
					if (aj > bj)
						pj = 0;
					else
						pj = 1;
					bb[bi * (n - 1) + bj] = aa[(bi + pi) * n + bj + pj];
				}
			}
			if ((ai + aj) % 2 == 0) 
				q = 1;
			else 
				q = -1;
			adjoint[ai * n + aj] = q * det(n - 1, bb);
		}
	}
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < i; j++) {
            double tem = adjoint[i * n + j];
            adjoint[i * n + j] =  adjoint[j * n + i];
            adjoint[j * n + i] =  tem;
        }
    }
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			aa[i * n + j] = adjoint[i * n + j] / det_aa;
		}
	}

	delete[] adjoint;
	delete[] bb;
}

#endif


