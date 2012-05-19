#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#include <string.h>
#endif

#include "mdgpu.h"

extern bool pad;

//#define my_malloc malloc
//#define my_free free

void* my_malloc(size_t size);
void my_free(void* ptr);

void MatrizAlloc(Matriz* m) {
	if (pad) {
		m->w = (m->n/16 + 1)*16;
	} else {
		m->w = m->n;
	}
	m->x = (FLOTANTE*)my_malloc(m->w*m->w*sizeof(FLOTANTE));
	memset(m->x, 0, m->w*m->w*sizeof(FLOTANTE));
}

void VectorPrint(FILE* f, int* v, int n) {
    for (int i = 0; i < n; i++) {
        fprintf(f, "%i ", v[i]);
    }
    fprintf(f, "\n");
}

void VectorPrint(int* v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%i ", v[i]);
    }
    printf("\n");
}

void MatrizPrint(FILE* f, Matriz* m, const char* format) {
    for (int i = 0; i < m->n; i++) {
        for (int j = 0; j < m->n; j++) {
            FLOTANTE v = MatrizGet(m, i, j);
            fprintf(f, format, v);
        }
        fprintf(f, "\n");
    }
	fflush(f);
}

void MatrizPrint(Matriz* m, const char* format = "%10.4lf ") {
    MatrizPrint(stdout, m, format);
}

void MatrizPrint(FILE* f, cs* m, const char* format = "%10.4lf ") {
    for (int i = 0; i < m->n; i++) {
        for (int j = 0; j < m->n; j++) {
			FLOTANTE v = MatrizGet(m, i, j);
            fprintf(f, format, v);
        }
        fprintf(f, "\n");
    }

}

void MatrizPrint(cs* m, const char* format = "%10.4lf ") {
    MatrizPrint(stdout, m, format);
}

void MatrizSpy(Matriz* m) {
    for (int i = 0; i < m->n; i++) {
        for (int j = 0; j < m->n; j++) {
            if (MatrizGet(m, i, j) == 0) {
				printf(".");          
			} else {
				printf("x");
			} 
        }
        printf("\n");
    }
}

void MatrizSpy(cs* m) {
	for (int j = 0; j < m->n; j++) {
		printf("%i", j/10);
	}
	printf("\n");
	for (int j = 0; j < m->n; j++) {
		printf("%i", j%10);
	}
	printf("\n");
    for (int i = 0; i < m->n; i++) {
        for (int j = 0; j < m->n; j++) {
            if (MatrizGet(m, i, j) == 0) {
				printf(".");          
			} else {
				printf("x");
			} 
        }
        printf("\n");
    }
}

void MatrizSpy2(Matriz* A, Matriz* L) {
    for (int i = 0; i < L->n; i++) {
        for (int j = 0; j < L->n; j++) {
            if (MatrizGet(L, i, j) == 0) {
                  printf(". ");          
            } else if (MatrizGet(A, i, j) == 0) {
                  printf("o ");
            } else {
                  printf("x ");                   
            }
		}
        printf("\n");
    }
}

/* print a sparse matrix */
int cs_print (FILE* f, const cs *A, int brief)
{
    int p, j, m, n, nzmax, nz, *Ap, *Ai ;
    FLOTANTE *Ax ;
    if (!A) { fprintf (f, "(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    fprintf (f, "CSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
        CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
        fprintf (f, "%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n", m, n, nzmax,
                Ap [n], cs_norm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            fprintf (f, "    col %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                fprintf (f, "      %d : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
                if (brief && p > 20) { fprintf (f, "  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        fprintf (f, "triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
        for (p = 0 ; p < nz ; p++)
        {
            fprintf (f, "    %d %d : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
            if (brief && p > 20) { fprintf (f, "  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}
