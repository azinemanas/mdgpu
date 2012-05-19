#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32

#else
#include <malloc.h>
#endif

#include "mdgpu.h"
#include "util.h"


/*int mem = 0;

void* my_malloc(size_t size) {
        mem += size;
        printf("mem += %d = %d\n", size, mem);
        void* ptr = malloc(size+sizeof(int));
        *((int*)ptr) = size;
        return ((int*)ptr)+1;
}

void my_free(void* ptr) {
        void* ptr2 = ((int*)ptr)-1;
        size_t size = *(int*)ptr2;
        mem -= size;
        printf("mem -= %d = %d\n", size, mem);
        free(ptr2);
}*/

void* my_malloc(size_t size) {
        clock_t t = tic();
        void* res = malloc(size);
        ticksMalloc += toc(t);
	if (res == NULL) {
		printf("memory allocation error (size = %i bytes)\n", size);
		exit(0);
	}
        return res;
}

void my_free(void* ptr) {
        clock_t t = tic();
        free(ptr);
        ticksFree += toc(t);
}

#ifdef _WIN32

clock_t z_clock() {

    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;
	long sec, usec;

	GetSystemTimeAsFileTime(&ft);
	li.LowPart  = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;
	t  = li.QuadPart;
	t /= 10;
	sec  = (long)(t / 1000000);
	usec = (long)(t % 1000000);

	return sec*1000000+usec;

}

#else

clock_t z_clock() {

	struct timeval tv;

	gettimeofday(&tv, NULL);
	return tv.tv_sec*1000000+tv.tv_usec;

	/*struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);

	return (ru.ru_utime.tv_sec+ru.ru_stime.tv_sec)*1000000+ru.ru_utime.tv_usec+ru.ru_stime.tv_usec;*/

}

#endif

#ifdef _WIN32

void z_sleep(int seg) {
	Sleep(seg*1000);
}

#else

void z_sleep(int seg) {
	sleep(seg);
}

#endif

#ifdef _WIN32

int memory_usage() {
	PROCESS_MEMORY_COUNTERS pmc;
	if (GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc))) {
		return pmc.WorkingSetSize;
	}
	return -1;
}

#else

int memory_usage() {
	struct mallinfo info = mallinfo();
	return info.arena/1024;
}

#endif

clock_t tic(void) { return (z_clock()); };
clock_t toc(clock_t t) { clock_t s = tic(); return (s-t); }
double ticks2seg(clock_t t) { return (t/1000000.0); }

