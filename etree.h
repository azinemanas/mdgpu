#ifndef _ETREE_H
#define _ETREE_H

#include "cs.h"

struct ETree {

	int nodo;
	
	int nHijos;
	struct ETree** hijos;
	
	int rows;
	int cols;

	int* nodos;

	int ceros;
};

struct ETree* GetETree(cs* m, int* parent, int* post, int* count, struct ETree*** listaNodos, int* cantNodos);
void GetSuperETree(struct ETree** listaNodos, int cantNodos, struct ETree*** listaNodosSuper, int* cantNodosSuper);
void GetRelaxedETree(struct ETree** listaNodosSuper, int cantNodosSuper, int max_ceros, int* count, int** perm, struct ETree*** listaNodosRelaxed_out, int* cantNodosRelaxed_out);
void PrintETree(struct ETree* root, int nivel);

#endif