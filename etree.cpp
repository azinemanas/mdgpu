#include "etree.h"

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#endif

struct ETree* GetETree(cs* A, int* parent, int* post, int* count, struct ETree*** listaNodos, int* cantNodos) {

    struct ETree* nodos = (struct ETree*)malloc(A->n*sizeof(struct ETree));
    for (int i = 0; i < A->n; i++) {
        nodos[i].nodo = i;
        nodos[i].nHijos = 0;
        nodos[i].hijos = NULL;
		nodos[i].rows = count[i];
		nodos[i].cols = 1;
		nodos[i].ceros = 0;
		nodos[i].nodos = (int*)malloc(1*sizeof(int));
		nodos[i].nodos[0] = i;
    }
    struct ETree* root = NULL;
    for (int i = 0; i < A->n; i++) {
        if (parent[i] == -1) {
            if (root != NULL) {
                    printf("Error: el arbol de eliminacion no es un arbol\n");
                    exit(1);
            } else {
                    root = &nodos[i];
            }
        } else {
            nodos[parent[i]].nHijos++;
            nodos[parent[i]].hijos = (struct ETree**)realloc(nodos[parent[i]].hijos, nodos[parent[i]].nHijos*sizeof(struct ETree*));
            nodos[parent[i]].hijos[nodos[parent[i]].nHijos-1] = &nodos[i];
        }
    }

	*listaNodos = (struct ETree**)malloc(A->n*sizeof(struct ETree*));
	for (int i = 0; i < A->n; i++) {
		(*listaNodos)[i] = &nodos[post[i]];
	}
	*cantNodos = A->n;

	return root;
}

void GetSuperETree(struct ETree** listaNodos, int cantNodos, struct ETree*** listaNodosSuper_out, int* cantNodosSuper_out) {

	struct ETree** stack = (struct ETree**) malloc(cantNodos*sizeof(struct ETree*));
	int n_stack = 0;

	struct ETree** listaNodosSuper = (struct ETree**) malloc(cantNodos*sizeof(struct ETree*));
	int cantNodosSuper = 0;

	for (int i = 0; i < cantNodos; i++) {
		struct ETree* tree = listaNodos[i];
		struct ETree* root = NULL;

		bool merge = false;
		if (tree->nHijos >= 1 && tree->hijos[tree->nHijos-1]->rows == tree->rows+1 && tree->hijos[tree->nHijos-1]->nodo == tree->nodo-1) {
			merge = true;
		}

		//printf("nodo = %i, nHijos = %i, rows = %i", tree->nodo, tree->nHijos, tree->rows);

		root = (struct ETree*)malloc(sizeof(struct ETree));
		root->nodos = NULL;

		struct ETree* hijoMerge = NULL;
		if (merge) {
			hijoMerge = stack[--n_stack];
		}

		if (!merge) {
			root->nodo = tree->nodo;
			root->nHijos = tree->nHijos;
			root->cols = tree->cols;
			root->rows = tree->rows;
		} else {
			root->nodo = hijoMerge->nodo;
			root->nHijos = tree->nHijos-1+hijoMerge->nHijos;
			root->cols = hijoMerge->cols+1;
			root->rows = hijoMerge->rows-1;
		}

		root->hijos = (struct ETree**)malloc(root->nHijos*sizeof(struct ETree*));
		
		if (!merge) {
			if (tree->nHijos >= 1) {
				root->hijos[tree->nHijos-1] = stack[--n_stack];
			}
		} else {
			for (int i = 0; i < hijoMerge->nHijos; i++) {
				root->hijos[tree->nHijos-1+i] = hijoMerge->hijos[i];
			}
		}

		for (int i = 0; i < tree->nHijos-1; i++) {
			root->hijos[i] = stack[--n_stack];
		}

		if (!merge) {
			listaNodosSuper[cantNodosSuper] = root;	
		} else {
			cantNodosSuper--;
			listaNodosSuper[cantNodosSuper] = root;	
		}
		cantNodosSuper++;

		stack[n_stack++] = root;
	}

	*listaNodosSuper_out = listaNodosSuper;
	*cantNodosSuper_out = cantNodosSuper;
}

void GetSuperETree_old(struct ETree** listaNodos, int cantNodos, struct ETree*** listaNodosSuper_out, int* cantNodosSuper_out) {

	struct ETree** stack = (struct ETree**) malloc(cantNodos*sizeof(struct ETree*));
	int n_stack = 0;

	struct ETree** listaNodosSuper = (struct ETree**) malloc(cantNodos*sizeof(struct ETree*));
	int cantNodosSuper = 0;

	for (int i = 0; i < cantNodos; i++) {
		struct ETree* tree = listaNodos[i];
		struct ETree* root = NULL;

		bool merge = false;
		if (tree->nHijos == 1 && tree->hijos[tree->nHijos-1]->rows == tree->rows+1 && tree->hijos[tree->nHijos-1]->nodo == tree->nodo-1) {
			merge = true;
		}

		if (merge) {

			root = stack[--n_stack];
			root->cols++;

		} else {

			root = (struct ETree*)malloc(sizeof(struct ETree));
			root->nHijos = tree->nHijos;
			root->hijos = (struct ETree**)malloc(root->nHijos*sizeof(struct ETree*));
			root->nodo = tree->nodo;
			root->cols = tree->cols;
			root->rows = tree->rows;

			for (int h = 0; h < root->nHijos; h++) {
				root->hijos[h] = stack[--n_stack];
			}

			listaNodosSuper[cantNodosSuper++] = root;
		}
		
		stack[n_stack++] = root;
	}

	*listaNodosSuper_out = listaNodosSuper;
	*cantNodosSuper_out = cantNodosSuper;
}

static void dfs(struct ETree* root, struct ETree** listaNodos, int* cantNodos, int* c, int* perm) {
	for (int i = 0; i < root->nHijos; i++) {
		dfs(root->hijos[i], listaNodos, cantNodos, c, perm);
	}

	listaNodos[(*cantNodos)++] = root;

	root->nodo = *c;
	for (int i = root->cols-1; i >= 0; i--) {
		perm[*c] = root->nodos[i];
		(*c)++;
	}

}

static void sort(int* a, int* c, int n) {
	for (int i = 0; i < n; i++) {
		int x = -1;
		for (int j = i; j < n; j++) {
			if (x == -1 || c[a[j]] < c[a[x]]) {
				x = j;
			}
		}
		int tmp = a[x];
		a[x] = a[i];
		a[i] = tmp;
	}
}

void GetRelaxedETree(struct ETree** listaNodos, int cantNodos, int max_ceros, int* count, int** perm, struct ETree*** listaNodosRelaxed_out, int* cantNodosRelaxed_out) {
	
	struct ETree** stack = (struct ETree**)malloc(cantNodos*sizeof(struct ETree*));
	int n_stack = 0;

	struct ETree** listaNodosRelaxed = (struct ETree**)malloc(cantNodos*sizeof(struct ETree*));
	int cantNodosRelaxed = 0;

	for (int i = 0; i < cantNodos; i++) {

		struct ETree* tree = listaNodos[i];

		struct ETree** hijosMerge = (struct ETree**)malloc(tree->nHijos*sizeof(struct ETree*));
		for (int h = 0; h < tree->nHijos; h++) {
			hijosMerge[h] = NULL;
		}

		struct ETree** hijos = &stack[n_stack-tree->nHijos];
		
		int ceros = tree->ceros;
				
		struct ETree* root = (struct ETree*)malloc(sizeof(struct ETree));
		root->nodo = tree->nodo;
		root->cols = tree->cols;
		root->rows = tree->rows;
		root->ceros = 0;
		root->nodos = (int*)malloc(root->cols*sizeof(int));
		root->nodos[0] = root->nodo;

		//printf("\n%i (%i)\n", tree->nodo, tree->rows);

		int cols = 1;

		while (true) {

			bool merge = false;
			int min_h = -1;
			int min_ceros = -1;
			struct ETree* hijoMerge = NULL;
			for (int h = 0; h < tree->nHijos; h++) {
				struct ETree* hijo = hijos[h];
				//int c = tree->rows+1-hijo->rows+ceros;
				int c = (tree->rows + cols - hijo->rows)*hijo->cols + hijo->ceros;
				if (hijosMerge[h] == NULL && c <= max_ceros && (min_ceros == -1 || c < min_ceros) /*tree->nHijos == 1 && hijo->rows == tree->rows+1*/ 
						/*&& hijo->nodo+hijo->cols-1 == tree->nodo-1*/) {
					min_ceros = c;
					hijoMerge = hijo;
					merge = true;
					min_h = h;
				}
			}

			if (!merge) {
				break;
			}

			hijosMerge[min_h] = hijoMerge;
			ceros += min_ceros;
			cols += hijoMerge->cols;
			
			/*printf("  <- ");
			for (int h = 0; h < hijoMerge->cols; h++) {
				printf("%i, ", hijoMerge->nodos[h]);
			}
			printf("\n");*/

			for (int h = 0; h < hijoMerge->cols; h++) {
				//printf("  ceros add %i to col %i\n", (min_ceros-hijoMerge->ceros)/(hijoMerge->cols), hijoMerge->nodos[h]);
				count[hijoMerge->nodos[h]] += (min_ceros-hijoMerge->ceros)/(hijoMerge->cols);
			}
			//printf("\n");

			/*printf("  ceros = %i\n", ceros);
			printf("  min_ceros = %i\n", min_ceros);
			printf("  cols = %i\n", cols);
			printf("\n");*/

		}

		root->ceros = ceros;
		root->nHijos = tree->nHijos;
		for (int h = 0; h < tree->nHijos; h++) {
			if (hijosMerge[h] != NULL) {
				root->nHijos--;
				root->nHijos += hijosMerge[h]->nHijos;
			}
		}

		root->hijos = (struct ETree**)malloc(root->nHijos*sizeof(struct ETree*));

		int p = 0;
		for (int h = 0; h < tree->nHijos; h++) {
			if (hijosMerge[h] == NULL) {
				root->hijos[p++] = hijos[h];
			} else {
				for (int c = 0; c < hijosMerge[h]->cols; c++) {
					root->nodos = (int*)realloc(root->nodos, (root->cols+1)*sizeof(int));
					root->nodos[root->cols] = hijosMerge[h]->nodos[c];
					root->cols++;
				}
				for (int hh = 0; hh < hijosMerge[h]->nHijos; hh++) {
					root->hijos[p++] = hijosMerge[h]->hijos[hh];
				}
			}
		}

		for (int c = 0; c < root->cols; c++) {
			if (root->nodos[c] < root->nodo) {
				root->nodo = root->nodos[c];
			}
		}

		sort(root->nodos, count, root->cols);

		/*printf("  nodos = ");
		for (int h = 0; h < root->cols; h++) {
			printf("%i, ", root->nodos[h]);
		}
		printf("\n");
		printf("  rows = %i\n", root->rows);
		printf("  ceros = %i\n", root->ceros);*/

		/*if (min_ceros > 0) {
			for (int h = 1; h < root->cols; h++) {
				printf("  ceros add %i to col %i\n", min_ceros, root->nodos[h]);
				count[root->nodos[h]] += min_ceros;
			}
		}*/

		n_stack -= tree->nHijos;
		stack[n_stack++] = root;

		free(hijosMerge);
	}

	*perm = (int*)malloc(cantNodos*sizeof(int));
	int c = 0;
	dfs(stack[n_stack-1], listaNodosRelaxed, &cantNodosRelaxed, &c, *perm);

	/*printf("\nPerm:\n");
	for (int i = 0; i < cantNodos; i++) {
		printf("%i ", (*perm)[i]);
	}
	printf("\n");*/
	
	*listaNodosRelaxed_out = listaNodosRelaxed;
	*cantNodosRelaxed_out = cantNodosRelaxed;

}

void PrintETree(struct ETree* root, int nivel) {
	for (int i = 0; i < nivel; i++) {
		printf("-");
	}
	if (root == NULL) {
		printf("> NULL\n");
	} else {
		if (root->cols == 1) {
			printf("> %i\n", root->nodo);
		} else {
			printf("> ");
			if (root->nodos != NULL) {
				for (int c = 0; c < root->cols; c++) {
					printf("%i, ", root->nodos[c]);
				}
			}
			printf("(%i-%i, %i, %i)\n", root->nodo, root->nodo+root->cols-1, root->rows, root->ceros);
		}
		for (int i = 0; i < root->nHijos; i++) {
			PrintETree(root->hijos[i], nivel+1);
		}
	}
}