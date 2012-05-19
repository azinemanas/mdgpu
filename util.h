#ifndef _UTIL_H
#define _UTIL_H

#include "types.h"

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

//#define my_malloc malloc
//#define my_free free

extern clock_t ticksMalloc;
extern clock_t ticksFree;

void* my_malloc(size_t size);
void my_free(void* ptr);

clock_t z_clock();
void z_sleep(int seg);
int memory_usage();
clock_t tic(void);
clock_t toc(clock_t t);
double ticks2seg(clock_t t);

#endif
