#ifndef _MDGPU_H
#define _MDGPU_H

#include "types.h"
#include "cs.h"

typedef struct {
	int w;
	int n;
	FLOTANTE* x;
	int* q;
	int* ql;
	int k;
	int offset;
} Matriz;

typedef struct {

	bool GPU;

	Matriz* h_frente;
	Matriz* hd_frente;
	//Matriz* d_frente;
	
	clock_t tiempoFact;
	
} Frente;

void MatrizAlloc(Matriz* m);

inline FLOTANTE MatrizGet(Matriz* m, int i, int j) {
	/*if (i < 0 || j < 0 || i >= m->n || j >= m->n) {
		printf("MatrizSet fuera de rango k=%i i=%i j=%i\n", m->k, i, j);
		exit(0);
	}*/
	return m->x[j*(m->w)+i];
}

#define _MatrizSet(m, i, j, v) {m->x[j*(m->w)+i] = v;}

inline void MatrizSet(Matriz* m, int i, int j, FLOTANTE v) {
	m->x[j*(m->w)+i] = v;
}

inline void MatrizAdd(Matriz* m, int i, int j, FLOTANTE v) {
	/*if (i < 0 || j < 0 || i >= m->n || j >= m->n) {
		printf("MatrizSet fuera de rango k=%i i=%i j=%i\n", m->k, i, j);
		exit(0);
	}*/
	m->x[j*(m->w)+i] += v;
}

inline void MatrizSet(cs* M, int i, int j, FLOTANTE v) {
    int p1 = M->p[j];
    int p2 = M->p[j+1];
    for (int k = p1; k < p2; k++) {
        if (M->i[k] == i) {
            M->x[k] = v;
        }
    }
}

inline FLOTANTE MatrizGet(cs* M, int i, int j) {
    int p1 = M->p[j];
    int p2 = M->p[j+1];
    for (int k = p1; k < p2; k++) {
        if (M->i[k] == i) {
            return M->x[k];
        }
    }
    return 0.0;
}

void VectorPrint(FILE* f, int* v, int n);
void VectorPrint(int* v, int n);
void MatrizPrint(FILE* f, Matriz* m, const char* format);
void MatrizPrint(Matriz* m, const char* format);
void MatrizPrint(FILE* f, cs* m, const char* format);
void MatrizPrint(cs* m, const char* format);
void MatrizSpy(Matriz* m);
void MatrizSpy(cs* m);
void MatrizSpy2(Matriz* A, Matriz* L);
int cs_print (FILE* f, const cs *A, int brief);

#endif
