#include <stdio.h>
#include <time.h>
#include <malloc.h> 
#include <mcheck.h>

#include "cublas.h"

#include "mdgpu.h"
#include "util.h"
#include "cutil.h"
#include "cs.h"
#include "etree.h"

#define blasint int
#define _Complex 

#ifndef _WIN32
	extern "C" {
#endif
		#include "cblas.h"
#ifndef _WIN32
	}
#endif

#define MatrizGetFk MatrizGetFk_3
#define FactorAux_CPU FactorAux_CPU_3
#define FactorAux_GPU FactorAux_GPU_2

bool GPU = true;

double flops = 0;
clock_t ticksMalloc = 0;
clock_t ticksMallocGPU = 0;
clock_t ticksFree = 0;
clock_t ticksFreeGPU = 0;
clock_t ticksFactorAux = 0;
clock_t ticksFactorAux1 = 0;
clock_t ticksFactorAux2 = 0;
clock_t ticksFactorAux3 = 0;
clock_t ticksMemcpy = 0;
clock_t ticksMemcpy2 = 0;
clock_t ticksMemcpy21 = 0;
clock_t ticksMemcpyX = 0;
clock_t ticksMerge = 0;
clock_t ticksExtendAdd = 0;
clock_t ticksSymbolic = 0;
clock_t ticksGetFk = 0;
clock_t ticksMemcpyHost = 0;

clock_t ticksTRSM_GPU = 0;
clock_t ticksGEMM_GPU = 0;

long countExtendAdd = 0;
int countGetFk = 0;

long bytesMemcpy2 = 0;
long bytesMemcpy21 = 0;

int block_size = -1;
bool pad = false;

FILE* logger;

Matriz* MatrizGetFk_1(cs* A, int k, int cols) {
    clock_t t = tic();

    Matriz* m = (Matriz*)my_malloc(sizeof(Matriz));
    m->k = k;
	m->offset = 0;
	m->n = 0;
    for (int i = k+cols; i < A->n; i++) {
		bool b = false;
		for (int h = 0; h < cols; h++) {
			if (MatrizGet(A, i, k+h) != 0) {
				b = true;
				break;
			}
		}
        if (b) {
            m->n++;
        }
    }
	m->n += cols;
    MatrizAlloc(m);
    m->q = (int*)my_malloc((m->n)*sizeof(int));
	for (int i = 0; i < cols; i++) {
		for (int h = 0; h < cols; h++) {
			MatrizSet(m, i, h, MatrizGet(A, k+i, k+h));
		}
        m->q[i] = k+i;
	}
    int j = 0;
	for (int i = k+cols; i < A->n; i++) {
		bool b = false;
		for (int h = 0; h < cols; h++) {
			if (MatrizGet(A, i, k+h) != 0) {
				b = true;
				break;
			}
		}
        if (b) {
			for (int h = 0; h < cols; h++) {
				MatrizSet(m, j+cols, h, MatrizGet(A, i, k+h));
			}
            m->q[j+cols] = i;
            j++;
        }
    }
    ticksGetFk += toc(t);
    return m;
}

Matriz* MatrizGetFk_2(cs* A, int k, int cols) {
    clock_t t = tic();

	int* pos = (int*) my_malloc(cols*sizeof(int));
	memset(pos, 0, cols*sizeof(int));

    Matriz* m = (Matriz*)my_malloc(sizeof(Matriz));
    m->k = k;
	m->offset = 0;
	m->n = 0;

	for (int h = 0; h < cols; h++) {
		while (A->i[A->p[k+h] + pos[h]] < k+h) {
			pos[h]++;
		}
	}

    for (int i = k+cols; i < A->n; i++) {
		bool b = false;
		for (int h = 0; h < cols; h++) {
			while (pos[h] < A->p[k+h+1] - A->p[k+h] && A->i[A->p[k+h] + pos[h]] < i) {
				pos[h]++;
			}
			if (pos[h] < A->p[k+h+1] - A->p[k+h] && A->i[A->p[k+h] + pos[h]] == i) {
				b = true;
			}
		}
        if (b) {
            m->n++;
		}
    }
	m->n += cols;
    MatrizAlloc(m);
    m->q = (int*)my_malloc((m->n)*sizeof(int));
	for (int i = 0; i < cols; i++) {
		for (int h = 0; h < cols; h++) {
			MatrizSet(m, i, h, MatrizGet(A, k+i, k+h));
		}
        m->q[i] = k+i;
	}
	
	memset(pos, 0, cols*sizeof(int));

    int j = 0;
	for (int i = k+cols; i < A->n; i++) {
		bool b = false;
		for (int h = 0; h < cols; h++) {
			while (pos[h] < A->p[k+h+1] - A->p[k+h] && A->i[A->p[k+h] + pos[h]] < i) {
				pos[h]++;
			}
			if (pos[h] < A->p[k+h+1] - A->p[k+h] && A->i[A->p[k+h] + pos[h]] == i) {
				MatrizSet(m, j+cols, h, /*MatrizGet(A, i, k+h)*/ A->x[A->p[k+h] + pos[h]]);
				b = true;
			}
		}
        if (b) {
            m->q[j+cols] = i;
            j++;
        }
    }

	my_free(pos);

    ticksGetFk += toc(t);
    return m;
}

Matriz* MatrizGetFk_3(cs* C, int k, int c) {

	const int MAX_INT = 0x7fffffff;

	Matriz* M = (Matriz*)my_malloc(sizeof(Matriz));
	M->k = k;
	M->offset = 0;
	M->n = 0;

	//int* pos = (int*)malloc(c*sizeof(int));
	//int* min_i = (int*)malloc(c*sizeof(int));
	//int* min_i2 = (int*)malloc(c*sizeof(int));

	int pos[c];
	int min_i[c];
	int min_i2[c];

	int minr;
	int* C_p = &C->p[k];
	int* C_i = C->i;
	FLOTANTE* C_x = C->x;

	static int* merge = NULL;
	if (merge == NULL) {	
		merge = (int*) malloc(C->n*sizeof(int));
	}

	clock_t t = tic();

	for (int i = 0; i < c; i++) {
		merge[M->n] = k+i;
		M->n++;
	}

	for (int i = 0; i < c; i++) {
		int pos_i = C_p[i];
		while (C_i[pos_i] < k+c && pos_i < C_p[i+1]) {
			pos_i++;
		}
		pos[i] = pos_i;
	}

	int h = 0;

	while (true) {
		minr = MAX_INT;

		if (h == 0) {

			for (int i = 0; i < c; i++) {
				if (pos[i]+1 < C_p[i+1]) {
					min_i[i] = C_i[pos[i]];
					min_i2[i] = C_i[pos[i]+1];
				} else 
				if (pos[i] < C_p[i+1]) {
					min_i[i] = C_i[pos[i]];
					min_i2[i] = MAX_INT;
				} else {
					min_i[i] = MAX_INT;
					min_i2[i] = MAX_INT;
				}	
			}

		}

		for (int i = 0; i < c; i++) {
			if (min_i[i] < minr) {
				minr = min_i[i];
			}
		}

		if (minr == MAX_INT) {
			break;
		}

		merge[M->n] = minr;
		M->n++;

		for (int i = 0; i < c; i++) {
			if (min_i[i] == minr) {
				pos[i]++;
				min_i[i] = min_i2[i];
			}
		}

		h = h == 0 ? 1 : 0;

	}

	MatrizAlloc(M);
	M->q = (int*)my_malloc(M->n*sizeof(int));

	memcpy(M->q, merge, M->n*sizeof(int));

	for (int i = 0; i < c; i++) {
		int p = 0;
		for (int j = C_p[i]; j < C_p[i+1]; j++) {
			while (merge[p] < C_i[j]) {
				p++;
			}
			_MatrizSet(M, p, i, C_x[j]);
		}
	}

	ticksGetFk += toc(t);

	//my_free(pos);
	//my_free(min_i);	
	//my_free(min_i2);

	return M;
}

inline int log2(int x) {
	int square = 2;
	int count = 1;
	while (square < x) {
		square *= 2;
		++count;
	}
	return count;
}	

inline int pow2(int x) {
	int count = 1;
	for (int i = 0; i < x; i++) {
		count *= 2;
	}
	return count;
}

inline static int sizeNivel(int x) {
	return pow2(x);
}    

inline static int comienzoNivel(int x) {
	return pow2(x)-1;
}

typedef struct {
	int x;
	int ptr;
} KMergeTree;

const int INF = 0x7fffffff;

#define kmerge_min(x, y) (x <= y)

Matriz* MatrizGetFk_4(cs* C, int k, int c) {

	clock_t t;

	Matriz* M = (Matriz*)my_malloc(sizeof(Matriz));
	M->k = k;
	M->offset = 0;
	M->n = 0;

	int* C_p = &C->p[k];
	int* C_i = C->i;
	int i, x1, x2, n1, n2;

	int niv = log2(c);
	int n_tree = sizeNivel(niv+1)-1;

	int pos[c];
	KMergeTree tree[n_tree];
	KMergeTree treeLoser[n_tree];
	int sizeNiveles[niv+1];
	int comienzoNiveles[niv+1];

	/*int* pos = (int*)malloc(c*sizeof(int));
	KMergeTree* tree = (KMergeTree*) my_malloc(n_tree*sizeof(KMergeTree));
	KMergeTree* treeLoser = (KMergeTree*) my_malloc(n_tree*sizeof(KMergeTree));
	int* sizeNiveles = (int*) my_malloc((niv+1)*sizeof(int));
	int* comienzoNiveles = (int*) my_malloc((niv+1)*sizeof(int));*/
	
	for (i = 0; i <= niv; i++) {
		sizeNiveles[i] = sizeNivel(i);
		comienzoNiveles[i] = comienzoNivel(i);
	}

	t = tic();

	M->n = c;

	for (i = 0; i < c; i++) {
		pos[i] = C_p[i];
		while (pos[i] < C_p[i+1] && C_i[pos[i]] < k+c) {
			pos[i]++;
		}
	}

	for (i = 0; i < n_tree; i++) {
		tree[i].x = INF;
		tree[i].ptr = INF;
		treeLoser[i].x = INF;
		treeLoser[i].ptr = INF;
	}

	for (i = 0; i < c; i++) {
		tree[comienzoNiveles[niv]+i].x = treeLoser[comienzoNiveles[niv]+i].x = pos[i] < C_p[i+1] ? C_i[pos[i]] : INF;
		tree[comienzoNiveles[niv]+i].ptr = treeLoser[comienzoNiveles[niv]+i].ptr = i;
	}

	for (int n = niv-1; n >= 0; n--) {
		for (i = 0; i < sizeNiveles[n]; i++) {
			x1 = tree[comienzoNiveles[n+1]+i*2].x;
			x2 = tree[comienzoNiveles[n+1]+i*2+1].x;

			n1 = tree[comienzoNiveles[n+1]+i*2].ptr;
			n2 = tree[comienzoNiveles[n+1]+i*2+1].ptr;

			tree[comienzoNiveles[n]+i].x = kmerge_min(x1,x2) ? x1 : x2;
			tree[comienzoNiveles[n]+i].ptr = kmerge_min(x1,x2) ? n1 : n2;

			treeLoser[comienzoNiveles[n]+i].x = kmerge_min(x1,x2) ? x2 : x1;
			treeLoser[comienzoNiveles[n]+i].ptr = kmerge_min(x1,x2) ? n2 : n1;
		}
	}

	KMergeTree winner = tree[0];

	int ant;
	
	ant = INF;

	while (winner.x != INF) {

		if (ant == INF || ant != winner.x) {
			ant = winner.x;
			M->n++;
		}

		int h = winner.ptr;
		int m = comienzoNiveles[niv]+h;
		if (pos[h]+1 < C_p[h+1]) {
			pos[h]++;
			treeLoser[m].x = C_i[pos[h]];
		} else {
			treeLoser[m].x = INF;
		}
		
		winner = treeLoser[m];		

		for (i = niv-1; i >= 0; i--) {

			//countGetFk++;

			int p = (m-1)/2;	

			if (!kmerge_min(winner.x, treeLoser[p].x)) {
				KMergeTree tmp;

				tmp = winner;
				winner = treeLoser[p];
				treeLoser[p] = tmp;

			}	
			m = p;
		}
	}

	MatrizAlloc(M);
	M->q = (int*)my_malloc(M->n*sizeof(int));

	for (i = 0; i < c; i++) {
		pos[i] = C_p[i];
	}

	for (i = 0; i < n_tree; i++) {
		tree[i].x = INF;
		tree[i].ptr = INF;

		treeLoser[i].x = INF;
		treeLoser[i].ptr = INF;
	}

	for (i = 0; i < c; i++) {
		tree[comienzoNiveles[niv]+i].x = treeLoser[comienzoNiveles[niv]+i].x = C_i[pos[i]];
		tree[comienzoNiveles[niv]+i].ptr = treeLoser[comienzoNiveles[niv]+i].ptr = i;
	}

	for (int n = niv-1; n >= 0; n--) {
		for (i = 0; i < sizeNiveles[n]; i++) {
			x1 = tree[comienzoNiveles[n+1]+i*2].x;
			x2 = tree[comienzoNiveles[n+1]+i*2+1].x;

			n1 = tree[comienzoNiveles[n+1]+i*2].ptr;
			n2 = tree[comienzoNiveles[n+1]+i*2+1].ptr;

			tree[comienzoNiveles[n]+i].x = kmerge_min(x1,x2) ? x1 : x2;
			tree[comienzoNiveles[n]+i].ptr = kmerge_min(x1,x2) ? n1 : n2;

			treeLoser[comienzoNiveles[n]+i].x = kmerge_min(x1,x2) ? x2 : x1;
			treeLoser[comienzoNiveles[n]+i].ptr = kmerge_min(x1,x2) ? n2 : n1;
		}
	}

	winner = tree[0];

	int row = 0;
	ant = INF;

	KMergeTree tmp;

	while (winner.x != INF) {

		countGetFk++;

		if (ant == INF || ant != winner.x) {
			ant = winner.x;
			row++;
		}

		int h = winner.ptr;
		int m = comienzoNiveles[niv]+h;

		if (pos[h] < C_p[h+1]) {
			MatrizSet(M, row-1, h, C->x[pos[h]]);
			M->q[row-1] = C_i[pos[h]];
		}

		if (pos[h]+1 < C_p[h+1]) {
			pos[h]++;
			treeLoser[m].x = C_i[pos[h]];
		} else {
			treeLoser[m].x = INF;
		}

		winner = treeLoser[m];

		for (i = niv-1; i >= 0; i--) {

			int p = (m-1)/2;	

			if (!kmerge_min(winner.x, treeLoser[p].x)) {
				tmp = winner;
				winner = treeLoser[p];
				treeLoser[p] = tmp;
			}	
			m = p;
		}
	}

	ticksGetFk += toc(t);

	/*my_free(comienzoNiveles);
	my_free(sizeNiveles);
	my_free(treeLoser);
	my_free(tree);
	my_free(pos);*/

	return M;
}

void VectorMerge(int* q1, int n1, int* q2, int n2, /* out */ int** q, /* out */ int* n) {
	
	int* qq = (int*) my_malloc(sizeof(int)*(n1+n2));
	*q = qq;

	clock_t tick = tic();

	int nn = 0;
	int k1 = 0, k2 = 0;
	while (k1 < n1 || k2 < n2) {

		if (k1 < n1 && k2 < n2) {
			if (q1[k1] < q2[k2]) {
				if (nn == 0 || qq[nn-1] != q1[k1]) {
					qq[nn] = q1[k1];
					nn++;
					k1++;
				} else {
					k1++;
				}
			} else {
				if (nn == 0 || qq[nn-1] != q2[k2]) {
					qq[nn] = q2[k2];
					nn++;
					k2++;
				} else {
					k2++;
				}
			}
		} else if (k1 < n1) {
			if (qq[nn-1] != q1[k1]) {
				qq[nn] = q1[k1];
				nn++;
			}
			k1++;
		} else {
			if (qq[nn-1] != q2[k2]) {
				qq[nn] = q2[k2];
				nn++;
			}
			k2++;
		}

	}

	*n = nn;

	ticksMerge += toc(tick);
}

void VectorMergeIndices(int* Q, int N, int* q, int n, /* out */ int** ql) {
	*ql = (int*) my_malloc(sizeof(int)*n);
	int i = 0, j = 0;
	while (i < n) {
		if (q[i] == Q[j]) {
			(*ql)[i] = j;
			i++;
		}
		j++;
	}
}

void AllocFrenteGPU(Frente** F, int nF) {

	clock_t tick = tic();

	int size_x = 0;
	int size_q = 0;
	for (int i = 0; i < nF; i++) {
		int n = F[i]->h_frente->w;
		size_x += n*n*sizeof(FLOTANTE);
		size_q += n*sizeof(int);
	}

	FLOTANTE* x;
	cutilSafeCall( cudaMalloc((void**)&x, size_x) );
	cutilSafeCall( cudaMemset(x, 0, size_x) );

	/*int* q;	
	cutilSafeCall( cudaMalloc((void**)&q, size_q) );*/

	//Matriz* f;
	///*cutilSafeCall(*/ cudaMalloc((void**)&f, size_f) /*)*/;

	size_x = 0;
	size_q = 0;
	for (int i = 0; i < nF; i++) {

		int n = F[i]->h_frente->w;

		F[i]->hd_frente = (Matriz*) my_malloc(sizeof(Matriz));
		F[i]->hd_frente->w = F[i]->h_frente->w;
		F[i]->hd_frente->n = F[i]->h_frente->n;
		F[i]->hd_frente->k = F[i]->h_frente->k;
		F[i]->hd_frente->offset = F[i]->h_frente->offset;

		F[i]->hd_frente->x = &x[size_x];
		size_x += n*n;

		/*F[i]->hd_frente->q = &q[size_q];
		size_q += n;*/		

	}

	ticksMallocGPU += toc(tick);
}

void AllocFrenteGPU(Frente* F) {

	AllocFrenteGPU(&F, 1);

}

void FreeFrenteGPU(Frente** F, int n) {

	clock_t tick = tic();

	for (int i = 0; i < n; i++) {

		if (F[i]->GPU) {

			my_free(F[i]->h_frente->ql);			
			my_free(F[i]->h_frente->q);
			my_free(F[i]->h_frente->x);
			my_free(F[i]->h_frente);

			if (i == 0) {
				cutilSafeCall( cudaFree(F[i]->hd_frente->x) );
				//cutilSafeCall( cudaFree(F[i]->hd_frente->q) );
			}			
			cutilSafeCall( cudaFree(F[i]->hd_frente->ql) );
			my_free(F[i]->hd_frente);
			my_free(F[i]);
		} else {
			my_free(F[i]->h_frente->x);
			my_free(F[i]->h_frente);
			my_free(F[i]);
		}

	}

	ticksFreeGPU += toc(tick);
}

void FreeFrenteGPU(Frente* F) {
	FreeFrenteGPU(&F, 1);
}

void MoverFrenteAGPU(Frente* F/*, cudaStream_t stream*/) {

	int n = F->h_frente->w;

	/*static FLOTANTE* buffer = NULL;
	static int sizeBuffer = 0;		
	if (buffer == NULL || n*F->h_frente->offset*sizeof(FLOTANTE) > sizeBuffer) {
		if (buffer != NULL) {
			cudaFreeHost(buffer);
		}
		sizeBuffer = n*F->h_frente->offset*sizeof(FLOTANTE);
		cudaMallocHost((void**)&buffer, sizeBuffer);
	}

	clock_t tick;
	//tick = tic();
	memcpy(buffer, F->h_frente->x, n*F->h_frente->offset*sizeof(FLOTANTE));
	//ticksMemcpyHost += toc(tick);*/

	clock_t tick;
	
	tick = tic();
		cutilSafeCall( cudaMemcpy(F->hd_frente->x, F->h_frente->x, n*F->h_frente->offset*sizeof(FLOTANTE), cudaMemcpyHostToDevice/*, stream*/) );
	ticksMemcpyX += toc(tick);

}

void MoverFrenteQlAGPU(Frente* F) {

	int n = F->h_frente->n /*- F->h_frente->offset*/;

	clock_t tick;

	tick = tic();
	cutilSafeCall( cudaMalloc((void**)&F->hd_frente->ql, n*sizeof(int)) );
	ticksMallocGPU += toc(tick);

	tick = tic();
	cutilSafeCall( cudaMemcpy(F->hd_frente->ql, F->h_frente->ql, n*sizeof(int), cudaMemcpyHostToDevice) );
	ticksMemcpy += toc(tick);
}

/*void MoverFrenteXACPU(Frente* F) {

	int n = F->hd_frente->w;

	clock_t tick = tic();
	cutilSafeCall( cudaMemcpy(F->h_frente->x, F->hd_frente->x, n*n*sizeof(FLOTANTE), cudaMemcpyDeviceToHost) );
	ticksMemcpy += toc(tick);
}*/

void ExtendAdd_CPU_old(Matriz* src, Matriz* dst) {

	int qi;
	int qj;
	for (int j = 0; j < src->n-src->offset; j++) {
		qj = src->ql[j];
		for (int i = 0; i < src->n-src->offset; i++) {
			qi = src->ql[i];
			MatrizSet(dst, qi, qj, MatrizGet(src, i+src->offset, j+src->offset)+MatrizGet(dst, qi, qj));
		}
	}

}

void ExtendAdd_CPU(Matriz* src, Matriz* dst) {

	int qi;
	int qj;
	int i, j, cd, cs;
	int n = src->n-src->offset;
	FLOTANTE* xd = dst->x;
	FLOTANTE* xs = src->x;
	for (j = 0; j < n; j++) {
		qj = src->ql[j];
		cd = qj*dst->w;
		cs = (j+src->offset)*src->w+src->offset;
		for (i = 0; i < n; i++) {
			qi = src->ql[i];
			xd[cd+qi] += xs[cs+i];
			countExtendAdd++;
		}
	}

}

__global__ void ExtendAdd_GPU_device(int src_offset, int src_n, int dst_n, 
									 FLOTANTE* src_x, FLOTANTE* dst_x, int* src_ql, int src_n2) {
	int i = blockIdx.x * 16 + threadIdx.x;
	int j = blockIdx.y * 16 + threadIdx.y;

	__shared__ int p1[16];
	__shared__ int p2[16];

	if (threadIdx.y == 0) {
		p1[threadIdx.x] = src_ql[i];
		p2[threadIdx.x] = src_ql[blockIdx.y * 16 + threadIdx.x] * dst_n;
		//printf("p1[%i] = %i; p2[%i] = %i\n", threadIdx.x, p1[threadIdx.x], threadIdx.x, p2[threadIdx.x]);
	}

	__syncthreads();

	if (i < src_n2 && j < src_n2) {
		dst_x[p2[threadIdx.y]+p1[threadIdx.x]] += src_x[j*src_n+i];
	}
}

void ExtendAdd_GPU(Frente* src, Frente* dst) {
	int n = src->h_frente->n;

	dim3 block(16,16);
	dim3 grid(n / 16 + 1, n / 16 + 1);
	ExtendAdd_GPU_device<<<grid,block>>>(src->h_frente->offset, src->h_frente->w, dst->h_frente->w,
		&src->hd_frente->x[src->h_frente->offset*src->h_frente->w+src->h_frente->offset], 
		dst->hd_frente->x, src->hd_frente->ql, src->h_frente->n-src->h_frente->offset);
	//cudaThreadSynchronize();
	cutilCheckMsg("ExtendAdd_GPU");
}

void ExtendAdd(Frente* src, Frente* dst, bool b) {

	clock_t tick = tic();

	if (GPU && !b) {
		ExtendAdd_GPU(src, dst);
	} else {
		ExtendAdd_CPU(src->h_frente, dst->h_frente);
	}

	ticksExtendAdd += toc(tick);
}

__global__ void FactorAux1_GPU_1_device(int c, Matriz* Uk, int* spL_i, FLOTANTE* spL_x) {
	
	int i = blockIdx.x * 256 + threadIdx.x;
	FLOTANTE raiz;
	int n = Uk->n;
	int w = Uk->w;

	raiz = sqrt(Uk->x[c*w+c]);
	if (i >= c && i < n) {
		spL_i[i-c] = Uk->q[i];
		spL_x[i-c] = Uk->x[c*w+i]/raiz;
	}

}

/*void FactorAux_GPU_1(Frente* Uk, cs* spL) {

	//clock_t tick;

	int cols = Uk->h_frente->offset;
	int n = Uk->hd_frente->n;
	int w = Uk->hd_frente->w;
	int k = Uk->hd_frente->k;
	FLOTANTE* x = Uk->hd_frente->x;
	FLOTANTE f;

	//printf("\n");
	//MatrizPrint(Uk->d_frente, "%f ");

	for (int c = 0; c < cols; c++) {

		//tick = tic();
		//FactorAux1_GPU(c, Uk, spL);
		//ticksFactorAux1 += toc(tick);

		FLOTANTE* d_spL_x;
		cutilSafeCall( cudaMalloc((void**)&d_spL_x, (n-c)*sizeof(FLOTANTE)) );

		int* d_spL_i;
		cutilSafeCall( cudaMalloc((void**)&d_spL_i, (n-c)*sizeof(int)) );

		dim3 block(256);
		int s = n / 256 + 1;
		dim3 grid(s);
		FactorAux1_GPU_1_device<<<grid,block>>>(c, Uk->d_frente, d_spL_i, d_spL_x);
		cutilCheckMsg("FactorAux1_GPU");

		cutilSafeCall( cudaMemcpy(&spL->i[spL->p[k+c]], d_spL_i, (n-c)*sizeof(int), cudaMemcpyDeviceToHost) );
		cutilSafeCall( cudaMemcpy(&spL->x[spL->p[k+c]], d_spL_x, (n-c)*sizeof(FLOTANTE), cudaMemcpyDeviceToHost) );

		cudaFree(d_spL_i);
		cudaFree(d_spL_x);

		//tick = tic();
		//FactorAux2_GPU_CUBLAS(c, Uk);
		//ticksFactorAux2 += toc(tick);

		cutilSafeCall( cudaMemcpy((void**)&f, x+(c*w+c), sizeof(FLOTANTE), cudaMemcpyDeviceToHost) );

		mdgpu_cublasXsyr('L', n-c-1, (FLOTANTE)-1.0/f, &x[(c+0)*w+(c+1)], 1, &x[(c+1)*w+(c+1)], w);

		cublasStatus status = cublasGetError();
		if (status != CUBLAS_STATUS_SUCCESS) {
			printf("Error en cublaSgemm()");
			exit(1);
		}
		//cudaThreadSynchronize();

	}

	//printf("\n");
	//MatrizPrint(Uk->d_frente, "%f ");
}*/

__global__ void FactorAux1_GPU_2_device_1(FLOTANTE* x) {
	
	x[0] = sqrt(x[0]);
}

__global__ void FactorAux1_GPU_2_device_2(int n, FLOTANTE* x) {
	
	int i = blockIdx.x * 256 + threadIdx.x;
	FLOTANTE f = x[0];
	if (i > 0 && i < n) {
		x[i] /= f;
	}
}

void FactorAux_GPU_2(Frente** F, int nF, cs* spL) {

	clock_t tick;
	clock_t tick2;

	cublasStatus status;

	if (block_size == -1) {
		printf("block_size == -1\n");
		exit(1);
	}

	#define cols(i) (F[i]->h_frente->offset)
	int b = block_size;

	int max_bloques = 16;

	static FLOTANTE* bloqueP = NULL;
	static FLOTANTE** bloques = NULL;
	static FLOTANTE* d_bloqueP = NULL;
	static FLOTANTE** d_bloques = NULL;

	if (bloqueP == NULL) {
		cutilSafeCall( cudaMallocHost((void**)&bloqueP, max_bloques*b*b*sizeof(FLOTANTE)) );
		//bloqueP = (FLOTANTE*) my_malloc(max_bloques*b*b*sizeof(FLOTANTE));
		
		bloques = (FLOTANTE**)	my_malloc(max_bloques*sizeof(FLOTANTE*));	
		for (int i = 0; i < max_bloques; i++) {
			bloques[i] = &bloqueP[i*b*b];
		}

		cutilSafeCall( cudaMalloc((void**)&d_bloqueP, max_bloques*b*b*sizeof(FLOTANTE)) );
		d_bloques = (FLOTANTE**) my_malloc(max_bloques*sizeof(FLOTANTE*));	
		for (int i = 0; i < max_bloques; i++) {
			d_bloques[i] = &d_bloqueP[i*b*b];
		}		
	}

	int max_cols = 0;
	for (int i = 0; i < nF; i++) {
		if (cols(i) > max_cols) {
			max_cols = cols(i);
		}
	}

	for (int j = 0; j < nF; j++) {
		F[j]->tiempoFact = 0;
	}
	
	for (int i = 0; i < max_cols; i += b) {

		for (int j = 0; j < nF; j++) {
		
			clock_t tini = tic();
		
			if (i < cols(j)) {
				FLOTANTE* x = F[j]->hd_frente->x;
				int w = F[j]->hd_frente->w;
				int nb = min(i+b,cols(j));
				int b2 = nb-i;
				FLOTANTE* bloque = bloques[j];
				
				bytesMemcpy2 += b2*b2*sizeof(FLOTANTE);
				
				if (j == 0) {
					
					//cutilSafeCall(
						cudaMemcpy2D(bloque, b*sizeof(FLOTANTE), &x[i*w+i], w*sizeof(FLOTANTE), b2*sizeof(FLOTANTE), b2, cudaMemcpyDeviceToHost);
					//);					
					
				} 
				
			}

			if (i < cols(j)) {

				if (j < nF-1) {
				
					int j2 = j+1;
					FLOTANTE* x = F[j2]->hd_frente->x;
					int w = F[j2]->hd_frente->w;
					int nb = min(i+b,cols(j2));
					int b2 = nb-i;
					FLOTANTE* bloque = bloques[j2];
					
					bytesMemcpy2 += b2*b2*sizeof(FLOTANTE);
					
					tick = tic();
					//cutilSafeCall(
						cudaMemcpy2DAsync(bloque, b*sizeof(FLOTANTE), &x[i*w+i], w*sizeof(FLOTANTE), b2*sizeof(FLOTANTE), b2, cudaMemcpyDeviceToHost);
					//);
					//cudaThreadSynchronize();
					ticksMemcpy21 += toc(tick);
				}
			
				int nb = min(i+b,cols(j));
				int b2 = nb-i;

				FLOTANTE* bloque = bloques[j];
				
				tick = tic();
				
				int cb = 0;
				for (int c = 0; c < b2; c++) {
					mdgpu_cblasXscal(b2-c, 1.0/sqrt(bloque[cb+c]), &bloque[cb+c], 1);
					mdgpu_cblasXsyr(CblasColMajor, CblasLower, b2-c-1, -1.0, &bloque[cb+c+1], 1, &bloque[cb+b+c+1], b);
					cb += b;
				}
				
				ticksFactorAux1 += toc(tick);
				
				cudaThreadSynchronize();
				
			}			
		
			if (i < cols(j)) {

				FLOTANTE* x = F[j]->hd_frente->x;
				int w = F[j]->hd_frente->w;

				int nb = min(i+b,cols(j));
				int b2 = nb-i;

				FLOTANTE* bloque = bloques[j];
				
				tick = tic();
				//cutilSafeCall(
					cudaMemcpy2D(&x[i*w+i], w*sizeof(FLOTANTE), bloque, b*sizeof(FLOTANTE), b2*sizeof(FLOTANTE), b2, cudaMemcpyHostToDevice);
				//);
				ticksMemcpy2 += toc(tick);
			}		
		
			if (i < cols(j)) {

				FLOTANTE* x = F[j]->hd_frente->x;
				int w = F[j]->hd_frente->w;
				int n = F[j]->hd_frente->n;

				int nb = min(i+b,cols(j));
				int b2 = nb-i;

				if (n-nb > 0) {

					tick = tic();
					mdgpu_cublasXtrsm('R', 'L', 'T', 'N', n-nb, nb-i, 1.0f, &x[i*w+i], w, &x[i*w+i+b2], w);
					cudaThreadSynchronize();
					status = cublasGetError();
					ticksTRSM_GPU += toc(tick);
					
					if (status != CUBLAS_STATUS_SUCCESS) {
						printf("Error en cublaXtrsm()");
						exit(1);
					}

					tick = tic();
					mdgpu_cublasXgemm('N', 'T', n-nb, n-nb, b2, -1.0f, &x[i*w+i+b2], w, &x[i*w+i+b2], w, 1.0f, &x[(i+b2)*w+i+b2], w);
					cudaThreadSynchronize();
					status = cublasGetError();
					ticksGEMM_GPU += toc(tick);
					
					if (status != CUBLAS_STATUS_SUCCESS) {
						printf("Error en cublaXgemm()");
						exit(1);
					}

				}
			}
		
			cudaThreadSynchronize();
			
			F[j]->tiempoFact += toc(tini);
		}
	
		/*tick = tic();
		
		for (int j = 0; j < nF; j++) {

			if (i < cols(j)) {
				FLOTANTE* x = F[j]->hd_frente->x;
				int w = F[j]->hd_frente->w;
				int nb = min(i+b,cols(j));
				int b2 = nb-i;
				FLOTANTE* bloque = bloques[j];
				
				bytesMemcpy2 += b2*b2*sizeof(FLOTANTE);
				//cutilSafeCall(
					cudaMemcpy2D(bloque, b*sizeof(FLOTANTE), &x[i*w+i], w*sizeof(FLOTANTE), b2*sizeof(FLOTANTE), b2, cudaMemcpyDeviceToHost); 
				//);
			}
		}

		cudaThreadSynchronize();
		ticksMemcpy21 += toc(tick);
		tick = tic();

		for (int j = 0; j < nF; j++) {

			if (i < cols(j)) {

				int nb = min(i+b,cols(j));
				int b2 = nb-i;

				FLOTANTE* bloque = bloques[j];

				
				int cb = 0;
				for (int c = 0; c < b2; c++) {
					mdgpu_cblasXscal(b2-c, 1.0/sqrt(bloque[cb+c]), &bloque[cb+c], 1);
					mdgpu_cblasXsyr(CblasColMajor, CblasLower, b2-c-1, -1.0, &bloque[cb+c+1], 1, &bloque[cb+b+c+1], b);
					cb += b;
				}
			}
		}

		ticksFactorAux1 += toc(tick);
		tick = tic();

		for (int j = 0; j < nF; j++) {

			if (i < cols(j)) {

				FLOTANTE* x = F[j]->hd_frente->x;
				int w = F[j]->hd_frente->w;

				int nb = min(i+b,cols(j));
				int b2 = nb-i;

				FLOTANTE* bloque = bloques[j];
				
				//cutilSafeCall(
					cudaMemcpy2D(&x[i*w+i], w*sizeof(FLOTANTE), bloque, b*sizeof(FLOTANTE), b2*sizeof(FLOTANTE), b2, cudaMemcpyHostToDevice);
				//);
			}
		}
		
		cudaThreadSynchronize();
		ticksMemcpy2 += toc(tick);
		tick2 = tic();

		for (int j = 0; j < nF; j++) {

			if (i < cols(j)) {

				FLOTANTE* x = F[j]->hd_frente->x;
				int w = F[j]->hd_frente->w;
				int n = F[j]->hd_frente->n;

				int nb = min(i+b,cols(j));
				int b2 = nb-i;

				if (n-nb > 0) {

					tick = tic();
					mdgpu_cublasXtrsm('R', 'L', 'T', 'N', n-nb, nb-i, 1.0f, &x[i*w+i], w, &x[i*w+i+b2], w);
					cudaThreadSynchronize();
					status = cublasGetError();
					ticksTRSM_GPU += toc(tick);
					
					if (status != CUBLAS_STATUS_SUCCESS) {
						printf("Error en cublaXtrsm()");
						exit(1);
					}

					tick = tic();
					mdgpu_cublasXgemm('N', 'T', n-nb, n-nb, b2, -1.0f, &x[i*w+i+b2], w, &x[i*w+i+b2], w, 1.0f, &x[(i+b2)*w+i+b2], w);
					cudaThreadSynchronize();
					status = cublasGetError();
					ticksGEMM_GPU += toc(tick);
					
					if (status != CUBLAS_STATUS_SUCCESS) {
						printf("Error en cublaXgemm()");
						exit(1);
					}

				}
			}
		}
		
		cudaThreadSynchronize();
		ticksFactorAux2 += toc(tick2);*/
	}

	printf("\n");
	for (int j = 0; j < nF; j++) {
		printf("Procesando frente %i, size = %i, cols = %i, tiempo = %.10f\n", F[j]->h_frente->k, F[j]->h_frente->n, F[j]->h_frente->offset, ticks2seg(F[j]->tiempoFact));
	}	
	
	tick = tic();

	int* spL_i;

	for (int j = 0; j < nF; j++) {

		FLOTANTE* x = F[j]->hd_frente->x;
		int w = F[j]->hd_frente->w;
		int n = F[j]->hd_frente->n;
		int k = F[j]->hd_frente->k;

		for (int c = 0; c < cols(j); c++) {
			spL_i = &spL->i[spL->p[k+c]-c];		
			memcpy(&spL_i[c], &F[j]->h_frente->q[c], (n-c)*sizeof(int));

			//cutilSafeCall( cudaMemcpy(&spL_i[c], &F[j]->h_frente->q[c], (n-c)*sizeof(int), cudaMemcpyHostToHost) );

			// TODO: traer L en bloques mas grandes
			//cutilSafeCall( cudaMemcpy(&spL->x[spL->p[k+c]], &x[c*w+c], (n-c)*sizeof(float), cudaMemcpyDeviceToHost) );
		}
	}

	ticksFactorAux3 += toc(tick);

	#undef cols

}

void FactorAux_CPU_1(Frente* Uk, cs* spL) {

	//clock_t tick;

	int cols = Uk->h_frente->offset;
	FLOTANTE raiz;
	int x;
	int n = Uk->h_frente->n;
	int k = Uk->h_frente->k;

	for (int c = 0; c < cols; c++) {

		//FactorAux1_CPU(c, Uk->h_frente, spL);

		flops += (n-c) + 1;
		raiz = sqrt(MatrizGet(Uk->h_frente, c, c));
		for (int i = c; i < n; i++) {
			x = spL->p[k+c] + i - c;
			spL->i[x] = Uk->h_frente->q[i];
			spL->x[x] = MatrizGet(Uk->h_frente, i, c)/raiz;

			//printf("%i %i %lf %lf\n", i, c, raiz, MatrizGet(Uk, i, c));
			if (MatrizGet(Uk->h_frente, c, c) < 0) {
				putchar('X');
				fflush(stdout);
				//printf("La matriz no es definida positiva!!!\n");
				//exit(1);
			}
		}

		//FactorAux2_CPU(c, Uk->h_frente);

		for (int j = c+1; j < n; j++) {
			for (int i = j; i < n; i++) {
				MatrizAdd(Uk->h_frente, i, j, - MatrizGet(Uk->h_frente, i, c) * MatrizGet(Uk->h_frente, j, c) / MatrizGet(Uk->h_frente, c, c));
			}
		}

	}
}

void FactorAux_CPU_2(Frente* Uk, cs* spL) {

	//clock_t tick;

	int cols = Uk->h_frente->offset;
	FLOTANTE raiz;
	int x;
	int n = Uk->h_frente->n;
	int k = Uk->h_frente->k;

	for (int c = 0; c < cols; c++) {

		//FactorAux1_CPU(c, Uk->h_frente, spL);

		flops += (n-c) + 1;
		raiz = sqrt(MatrizGet(Uk->h_frente, c, c));
		for (int i = c; i < n; i++) {
			x = spL->p[k+c] + i - c;
			spL->i[x] = Uk->h_frente->q[i];
			spL->x[x] = MatrizGet(Uk->h_frente, i, c)/raiz;

			//printf("%i %i %lf %lf\n", i, c, raiz, MatrizGet(Uk, i, c));
			if (MatrizGet(Uk->h_frente, c, c) < 0) {
				putchar('X');
				fflush(stdout);
				//printf("La matriz no es definida positiva!!!\n");
				//exit(1);
			}
		}

		//FactorAux2_CPU(c, Uk->h_frente);

		FLOTANTE f = (FLOTANTE)-1.0 / MatrizGet(Uk->h_frente, c, c);
		flops += (n-c)*(n-c-1);
		mdgpu_cblasXsyr(CblasColMajor, CblasLower, n-c-1, f, &Uk->h_frente->x[(c+0)*n+(c+1)], 1, &Uk->h_frente->x[(c+1)*n+(c+1)], n);

	}
}

void FactorAux_CPU_3(Frente* Uk, cs* spL) {

	//clock_t tick;

	int cols = Uk->h_frente->offset;
	FLOTANTE raiz;
	int x;
	int n = Uk->h_frente->n;
	int w = Uk->h_frente->w;
	int k = Uk->h_frente->k;

	//printf("\nk = %i\ncols = %i\n", k, cols);
	//MatrizPrint(Uk->h_frente, "%f ");

	int b = block_size == -1 ? n : block_size;

	for (int i = 0; i < cols; i += b) {

		int nb = min(i+b,cols);
		int b2 = nb-i;
		for (int c = i; c < nb; c++) {
		
			raiz = sqrt(MatrizGet(Uk->h_frente, c, c));
			for (int j = c; j < nb; j++) {
				MatrizSet(Uk->h_frente, j, c, MatrizGet(Uk->h_frente, j, c)/raiz);
			}

			mdgpu_cblasXsyr(CblasColMajor, CblasLower, nb-c-1, -1.0, Uk->h_frente->x+((c+0)*w+(c+1)), 1, Uk->h_frente->x+((c+1)*w+(c+1)), w);
		}

		if (n-nb > 0) {

			mdgpu_cblasXtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, 
				n-nb, nb-i, 1.0f, Uk->h_frente->x+(i*w+i), w, Uk->h_frente->x+(i*w+i+b2), w);

			//printf("\n");
			//MatrizPrint(Uk->h_frente, "%f ");

			mdgpu_cblasXgemm(CblasColMajor, CblasNoTrans, CblasTrans, 
				n-nb, n-nb, b2, -1.0f, &Uk->h_frente->x[i*w+i+b2], w, &Uk->h_frente->x[i*w+i+b2], w, 
					1.0f, &Uk->h_frente->x[(i+b2)*w+i+b2], w);

		}
	}

	//printf("\n");
	//MatrizPrint(Uk->h_frente, "%f ");

	for (int c = 0; c < cols; c++) {
		for (int j = c; j < n; j++) {
			x = spL->p[k+c] + j - c;
			spL->i[x] = Uk->h_frente->q[j];
			spL->x[x] = MatrizGet(Uk->h_frente, j, c);
		}
	}

	for (int c = 0; c < cols; c++) {
		
		flops += (n-c) + 1;
		flops += (n-c)*(n-c-1);

	}

}

void FactorAux(Frente* Uk, cs* spL) {

	clock_t tick = tic();

	if (GPU) {
		FactorAux_GPU(&Uk, 1, spL);
	} else {
		FactorAux_CPU(Uk, spL);
	}

	ticksFactorAux += toc(tick);
}

void FactorAux(Frente** U, int n_U, cs* spL) {

	clock_t tick = tic();

	if (GPU) {
		/*for (int i = 0; i < n_U; i++) {
			FactorAux_GPU(&U[i], 1, spL);
		}*/
		FactorAux_GPU(U, n_U, spL);
	} else {
		for (int i = 0; i < n_U; i++) {
			FactorAux_CPU(U[i], spL);
		}
	}

	ticksFactorAux += toc(tick);
}

void Factor(cs* A, cs* spL, struct ETree** listaNodos, int cantNodos) {

	Frente** stack = (Frente**) my_malloc(cantNodos*sizeof(Frente*));
	int n_stack = 0;

	for (int c = 0; c < cantNodos; c++) {
		struct ETree* tree = listaNodos[c];

		if (c % 20 == 0) {
			//putchar('.');
			//fflush(stdout);
		}

		if (tree->nHijos == 0) {
			continue;
		}

		Frente** F = (Frente**) my_malloc(tree->nHijos*sizeof(Frente*));

		for (int i = 0; i < tree->nHijos; i++) {
			F[i] = (Frente*) my_malloc(sizeof(Frente));
			F[i]->h_frente = MatrizGetFk(A, tree->hijos[i]->nodo, tree->hijos[i]->cols);			
			F[i]->GPU = false;
			F[i]->h_frente->offset = tree->hijos[i]->cols;

		}

		if (GPU) {
			AllocFrenteGPU(F, tree->nHijos);
		}

		for (int i = 0; i < tree->nHijos; i++) {
			if (GPU) {
				MoverFrenteAGPU(F[i]);
				F[i]->GPU = true;
			}
		}

		Frente** U = (Frente**) my_malloc(sizeof(Frente*)*tree->nHijos);

		for (int i = 0; i < tree->nHijos; i++) {
			if (tree->hijos[tree->nHijos-1-i]->nHijos == 0) {
				U[tree->nHijos-1-i] = NULL;
			} else {
				U[tree->nHijos-1-i] = stack[--n_stack];
			}
		}

		Frente** U2 = (Frente**) my_malloc(sizeof(Frente*)*tree->nHijos);

		for (int i = 0; i < tree->nHijos; i++) {

			int* q;
			int n;

			if (U[i] == NULL) {

				U2[i] = F[i];
				//U2[i]->h_frente->offset = tree->hijos[i]->cols;

			} else {

				VectorMerge(F[i]->h_frente->q, F[i]->h_frente->n, U[i]->h_frente->q + U[i]->h_frente->offset, U[i]->h_frente->n - U[i]->h_frente->offset, &q, &n);
				VectorMergeIndices(q, n, F[i]->h_frente->q, F[i]->h_frente->n, &F[i]->h_frente->ql);
				VectorMergeIndices(q, n, U[i]->h_frente->q, U[i]->h_frente->n, &U[i]->h_frente->ql);
				
				U2[i] = (Frente*) my_malloc(sizeof(Frente));
				U2[i]->GPU = false;
				U2[i]->h_frente = (Matriz*) my_malloc(sizeof(Matriz));
				U2[i]->h_frente->n = n;
				U2[i]->h_frente->q = q;
				U2[i]->h_frente->k = tree->hijos[i]->nodo;
				U2[i]->h_frente->offset = tree->hijos[i]->cols;
				MatrizAlloc(U2[i]->h_frente);

				if (GPU) {
					AllocFrenteGPU(U2[i]);
					U2[i]->GPU = true;
				}

				if (GPU) {
					MoverFrenteQlAGPU(F[i]);
					MoverFrenteQlAGPU(U[i]);
				}

				F[i]->h_frente->offset = 0;

				ExtendAdd(F[i], U2[i], false);
				ExtendAdd(U[i], U2[i], false);

			}
		}

		FactorAux(U2, tree->nHijos, spL);

		int n = U2[0]->h_frente->n - U2[0]->h_frente->offset;
		int* q = (int*) my_malloc(n*sizeof(int));
		memcpy(q, U2[0]->h_frente->q + U2[0]->h_frente->offset, n*sizeof(int));

		for (int i = 1; i < tree->nHijos; i++) {

			int n2;
			int* q2;

			VectorMerge(q, n, U2[i]->h_frente->q + U2[i]->h_frente->offset, U2[i]->h_frente->n - U2[i]->h_frente->offset, &q2, &n2);
			
			my_free(q);
	
			n = n2;
			q = q2;
		}

		for (int i = 0; i < tree->nHijos; i++) {
			VectorMergeIndices(q, n, U2[i]->h_frente->q + U2[i]->h_frente->offset, U2[i]->h_frente->n - U2[i]->h_frente->offset, &U2[i]->h_frente->ql);
		}

		Frente* Uk = (Frente*) my_malloc(sizeof(Frente));
		Uk->GPU = false;
		Uk->h_frente = (Matriz*) my_malloc(sizeof(Matriz));
		Uk->h_frente->n = n;
		Uk->h_frente->q = q;
		Uk->h_frente->k = tree->nodo;
		Uk->h_frente->offset = 0;
		MatrizAlloc(Uk->h_frente);

		if (GPU) {
			AllocFrenteGPU(Uk);
			Uk->GPU = true;
		}

		for (int i = 0; i < tree->nHijos; i++) {

			if (GPU) {
				MoverFrenteQlAGPU(U2[i]);
			}

			ExtendAdd(U2[i], Uk, false);

		}

		FreeFrenteGPU(F, tree->nHijos);

		for (int i = 0; i < tree->nHijos; i++) {
			if (U[i] != NULL) {
				FreeFrenteGPU(U[i]);
				FreeFrenteGPU(U2[i]);
			}
		}
		my_free(U);
		my_free(U2);
		my_free(F);

		//fprintf(logger, "%i ", Uk->h_frente->n);

		stack[n_stack++] = Uk;
	}

	Frente* U = stack[--n_stack];

	Frente* F = (Frente*) my_malloc(sizeof(Frente));
	F->GPU = false;
	F->h_frente = MatrizGetFk(A, listaNodos[cantNodos-1]->nodo, listaNodos[cantNodos-1]->cols);
	F->h_frente->offset = /*listaNodos[cantNodos-1]->cols*/ 0;

	if (GPU) {
		AllocFrenteGPU(F);
		MoverFrenteAGPU(F/*, 0*/);
		F->GPU = true;
	}

	int n;
	int* q;
	VectorMerge(U->h_frente->q, U->h_frente->n, F->h_frente->q, F->h_frente->n, &q, &n);

	VectorMergeIndices(q, n, U->h_frente->q, U->h_frente->n, &U->h_frente->ql);	
	VectorMergeIndices(q, n, F->h_frente->q, F->h_frente->n, &F->h_frente->ql);

	Frente* U2 = (Frente*) my_malloc(sizeof(Frente));
	U2->GPU = false;
	U2->h_frente = (Matriz*) my_malloc(sizeof(Matriz));
	U2->h_frente->n = n;
	U2->h_frente->q = q;
	U2->h_frente->k = listaNodos[cantNodos-1]->nodo;
	U2->h_frente->offset = listaNodos[cantNodos-1]->cols;
	MatrizAlloc(U2->h_frente);

	if (GPU) {
		AllocFrenteGPU(U2);
		U2->GPU = true;
	}

	if (GPU) {
		MoverFrenteQlAGPU(F);
		MoverFrenteQlAGPU(U);
	}

	ExtendAdd(F, U2, false);
	ExtendAdd(U, U2, false);

	my_free(F->h_frente->x);
	my_free(F->h_frente->q);
	my_free(F->h_frente);
	my_free(F);

	my_free(U->h_frente->x);
	my_free(U->h_frente);
	my_free(U);

	FactorAux(U2, spL);

	my_free(U2->h_frente->x);
	my_free(U2->h_frente->q);
	my_free(U2->h_frente);

}

void testVectorMerge() {
	int q1[] = {2, 5};
	int q2[] = {0, 1, 2, 3, 4, 5};

	int n;
	int* q;

	VectorMerge(q1, 2, q2, 6, &q, &n);

	VectorPrint(q, n);

	int* ql1;
	VectorMergeIndices(q, n, q1, 2, &ql1);
	VectorPrint(ql1, 2);
}

/* true for off-diagonal entries */
static int dropdiag (int i, int j, FLOTANTE aij, void *other) { return (i != j) ;}

/* C = A + triu(A,1)' */
static cs *make_sym (cs *A)
{
    cs *AT, *C ;
    AT = cs_transpose (A, 1) ;          /* AT = A' */
    cs_fkeep (AT, &dropdiag, NULL) ;    /* drop diagonal entries from AT */
    C = cs_add (A, AT, 1, 1) ;          /* C = A+AT */
    cs_spfree (AT) ;
    return (C) ;
}

/* infinity-norm of x */
static FLOTANTE norm (FLOTANTE *x, int n)
{
    int i ;
    FLOTANTE normx = 0 ;
    for (i = 0 ; i < n ; i++) normx = CS_MAX (normx, fabs (x [i])) ;
    return (normx) ;
}

void help() {
	printf("mdgpu [-gpu|-cpu] [-amd] [-sn|-rn c] [-b bs] [-pad] matriz\n");
	exit(1);
}

void test() {
	FILE* f = fopen("../../mat/nos4.mtx", "r");

	cs* T = cs_load(f);
	cs* C = cs_compress(T);
	cs_spfree(T);

	Matriz* m = MatrizGetFk_3(C, 0, 5);
	MatrizPrint(m, "%3.4f ");

}

int main(int argc, char** argv) {

	mallopt(M_MMAP_MAX, 0);

	//test();
	//exit(1);

	//printf("\nMemory: %i KB\n", memory_usage() );

	cublasStatus status;

	status = cublasInit();
	if (status != CUBLAS_STATUS_SUCCESS) {
		printf("Error en cublasInit(%i)", status);
		exit(1);
	}

	logger = fopen("salida.txt", "w");

	clock_t mainTick = 0;
	clock_t factorTick = 0;
	clock_t loadTick = 0;

	mainTick = tic();

	GPU = false;
	bool amd = false;
	bool supernodos = false;
	bool relaxednodos = false;
	int relaxednodosMaxCeros = 0;

	for (int i = 1; i < argc-1; i++) {
		//printf("%i -> %s\n", i, argv[i]);
		if (strcmp(argv[i], "-gpu") == 0) {
			//printf("gpu\n");
			GPU = true;
		} else if (strcmp(argv[i], "-cpu") == 0) {
			//printf("cpu\n");
			GPU = false;
		} else if (strcmp(argv[i], "-amd") == 0) {
			//printf("amd\n");
			amd = true;
		} else if (strcmp(argv[i], "-sn") == 0) {
			//printf("sn\n");
			supernodos = true;
		} else if (strcmp(argv[i], "-rn") == 0) {
			//printf("sn\n");
			relaxednodos = true;
			relaxednodosMaxCeros = atoi(argv[i+1]);
		} else if (strcmp(argv[i], "-b") == 0) {
			block_size = atoi(argv[i+1]);
		} else if (strcmp(argv[i], "-pad") == 0) {
			pad = true;
		} else if (strcmp(argv[i], "-help") == 0) {
			help();
		}
	}

	loadTick = tic();
	FILE* f = fopen(argv[argc-1], "r");

	cs* T = cs_load(f);
	cs* C_low = cs_compress(T);
	cs_spfree(T);

	loadTick = toc(loadTick);

	//printf("\nMemory: %i KB\n", memory_usage() );

	ticksSymbolic = tic();

	cs* A_low = NULL;
	cs* A = NULL;

	//MatrizSpy(C_low);
	
	int* amd_perm;
	int* amd_pinv;

	if (amd) {

		cs* C_up = cs_transpose(C_low, 1);
		//cs_spfree(C_low);
		
		//printf("\nC = \n");
		//MatrizPrint(C, "%10.4lf ");

		amd_perm = cs_amd(1, C_up);
		amd_pinv = cs_pinv(amd_perm, C_up->n);

		cs* A_amd_U = cs_symperm(C_up, amd_pinv, 1);
		cs_spfree(C_up);
		
		A_low = cs_transpose(A_amd_U, 1);
		cs_spfree(A_amd_U);

		//printf("\nA_amd_U = \n");
		//MatrizPrint(A_amd_U, "%10.4lf ");

		A = make_sym(A_low);

		//printf("\nA = \n");
		//MatrizPrint(A, "%10.4lf ");

	} else {
		A_low = C_low;
		A = make_sym(A_low);
	}

	/*if (amd) {
		printf("amd_perm: ");
		VectorPrint(amd_perm, A->n);

		printf("amd_pinv: ");
		VectorPrint(amd_pinv, A->n);
	}*/

	int* parent = cs_etree(A, 0);
	int* post = cs_post(parent, A->n);
	int* count = cs_counts(A, parent, post, 0);

	//printf("counts: ");
	//VectorPrint(count, A->n);

	//printf("sum counts: ");
	//VectorPrint(spL->p, A->n);

	//printf("etree: ");
	//VectorPrint(parent, A->n);	

	//printf("post: ");
	//VectorPrint(post, A->n);	

	struct ETree** listaNodos;

	//printf("\nMemory: %i KB\n", memory_usage() );
	int cantNodos;
	struct ETree* tree = GetETree(A, parent, post, count, &listaNodos, &cantNodos);
	//printf("Memory ETree: %i KB\n", memory_usage() );
	//PrintETree(tree, 0);

	//printf("\n---------------------\n\n");

	//PrintETree(listaNodos[cantNodos-1], 0);

	printf("\nCant Nodos = %i\n", cantNodos);

	struct ETree** listaNodosSuper;
	int cantNodosSuper;
	if (supernodos) {
		GetSuperETree(listaNodos, cantNodos, &listaNodosSuper, &cantNodosSuper);

		printf("\nCant SuperNodos = %i\n", cantNodosSuper);
	} else {
		listaNodosSuper = listaNodos;
		cantNodosSuper = cantNodos;
	}
	
	//PrintETree(listaNodosSuper[cantNodosSuper-1], 0);

	if (relaxednodos) {

		int* permRelaxedETree;

		GetRelaxedETree(listaNodos, cantNodos, relaxednodosMaxCeros, count, &permRelaxedETree, &listaNodosSuper, &cantNodosSuper);

		my_free(listaNodos);

		printf("\nCant RelaxedNodos = %i\n", cantNodosSuper);

		int* ipermRelaxedETree = cs_pinv(permRelaxedETree, cantNodos);

		cs* tmp1 = cs_transpose(A_low, 1);
		cs_spfree(A_low);

		cs* tmp2 = cs_symperm(tmp1, ipermRelaxedETree, 1);
		cs_spfree(tmp1);

		A_low = cs_transpose(tmp2, 1);
		cs_spfree(tmp2);

		int* count2 = (int*)malloc(cantNodos*sizeof(int));
		cs_pvec(permRelaxedETree, count, count2, cantNodos);
		count = count2;
		int* amd_perm2 = (int*)malloc(cantNodos*sizeof(int));
		cs_pvec(permRelaxedETree, amd_perm, amd_perm2, cantNodos);
		amd_perm = amd_perm2;
	}

	//MatrizPrint(A, "%f ");

	//printf("\nMemory: %i KB\n", memory_usage() );

	cs* spL = (cs*) my_malloc(sizeof(cs));
	spL->n = A->n;
	spL->m = A->m;
	spL->p = (int*) my_malloc((spL->n+1)*sizeof(int));
	int nz = 0;
	for (int i = 0; i < spL->n; i++) {
		spL->p[i] = nz;
		nz += count[i];
	}
	spL->p[spL->n] = nz;
	spL->nzmax = nz;
	spL->nz = -1;
	spL->i = (int*) my_malloc(nz*sizeof(int));
	memset(spL->i, 0, nz*sizeof(int));
	//spL->x = (FLOTANTE*) my_malloc(nz*sizeof(FLOTANTE));

	printf("L.nz = %i\n", nz);

	//printf("\nMemory: %i KB\n", memory_usage() );

	my_free(parent);
	my_free(post);
	my_free(count);

	ticksSymbolic = toc(ticksSymbolic);

	cs_spfree(A);

	//mtrace();
	factorTick = tic();
	Factor(A_low, spL, listaNodosSuper, cantNodosSuper);
	factorTick = toc(factorTick);
	//muntrace();

	printf("\n");
	//MatrizSpy(spL);
	//cs_print(spL, 0);
	//MatrizPrint(logger, spL);
	//MatrizPrint(spL, "%f ");

	mainTick = toc(mainTick);

	printf("\n");
	printf("MFlops     = %.2f\n", flops/1000000.0);
	printf("MFlops/s   = %.3f\n", flops/1000000.0/ticks2seg(ticksFactorAux));
	printf("\n");
	printf("Total          = %.6f\n\n", ticks2seg(mainTick));
	printf("Factor         = %.6f Factorizacion numerica total\n", ticks2seg(factorTick));
	printf("FactorAux      = %.6f Factorizacion numerica solo calculo\n", ticks2seg(ticksFactorAux));
	printf("  BLAS CPU     = %.6f\n", ticks2seg(ticksFactorAux1));
	printf("  BLAS GPU     = %.6f\n", ticks2seg(ticksFactorAux2));
	printf("    TRSM GPU   = %.6f\n", ticks2seg(ticksTRSM_GPU));
	printf("    GEMM GPU   = %.6f\n", ticks2seg(ticksGEMM_GPU));
	printf("  Otros        = %.6f\n", ticks2seg(ticksFactorAux3));
	printf("Memcpy         = %.6f\n", ticks2seg(ticksMemcpy)+ticks2seg(ticksMemcpy2)+ticks2seg(ticksMemcpy21)+ticks2seg(ticksMemcpyX));
	printf("  Memcpy2      = %.6f Copia paneles CPU<->GPU\n", ticks2seg(ticksMemcpy2+ticksMemcpy21));
	//printf("  Memcpy21     = %.6f\n", ticks2seg(ticksMemcpy21));
	printf("  MemcpyX      = %.6f Copia inicial del frente CPU->GPU\n", ticks2seg(ticksMemcpyX));
	printf("  Otros        = %.6f\n", ticks2seg(ticksMemcpy));
	//("MemcpyHost     = %.6f\n", ticks2seg(ticksMemcpyHost));

	printf("Merge          = %.6f\n", ticks2seg(ticksMerge));
	printf("GetFk          = %.6f Formar el frente\n", ticks2seg(ticksGetFk));
	printf("ExtendAdd      = %.6f Suma extendida\n", ticks2seg(ticksExtendAdd));

	printf("Malloc CPU     = %.6f\n", ticks2seg(ticksMalloc));
	printf("Free CPU       = %.6f\n", ticks2seg(ticksFree));
	printf("Malloc GPU     = %.6f\n", ticks2seg(ticksMallocGPU));
	printf("Free GPU       = %.6f\n\n", ticks2seg(ticksFreeGPU));
	
	printf("Load           = %.6f\n", ticks2seg(loadTick));
	printf("Symbolic       = %.6f\n", ticks2seg(ticksSymbolic));

	//printf("\ncountExtendAdd   = %li\n", countExtendAdd);
	//printf("countGetFk         = %i\n", countGetFk);
	//printf("bytesMemcpy2       = %li\n", bytesMemcpy2);
	//printf("bytesMemcpy21      = %li\n", bytesMemcpy21);

	fclose(logger);

	exit(1);

	A = make_sym(A_low);

	FLOTANTE* b = (FLOTANTE*) my_malloc(A->n*sizeof(FLOTANTE));
	FLOTANTE* x = (FLOTANTE*) my_malloc(A->n*sizeof(FLOTANTE));
	FLOTANTE* x2 = (FLOTANTE*) my_malloc(A->n*sizeof(FLOTANTE));

	for (int i = 0; i < A->n; i++) {
		b[i] = 1 + ((FLOTANTE) i) / A->n;
		x[i] = b[i];
	}

	if (amd) {
		cs_pvec(amd_perm, b, x, A->n);
	}

	cs_lsolve(spL, x);
	cs_ltsolve(spL, x);

	if (amd) {
		cs_ipvec(amd_perm, x, x2, A->n);
	}	

	FLOTANTE* resid = (FLOTANTE*) my_malloc(A->n*sizeof(FLOTANTE));
	for (int i = 0; i < A->n; i++) {
		resid[i] = -b[i];
	}

 	cs* C = make_sym(C_low);

	cs_gaxpy(C, x2, resid);

	printf ("resid: %8.2e\n", norm (resid,A->n) / (cs_norm (A) * norm (x,A->n) + norm (b,A->n))) ;
	printf ("resid: %8.2e\n", norm (resid,A->n)) ;

	cudaThreadExit();

	cublasShutdown();

	printf("FIN\n");
	//getchar();
	exit(0);
}
