/***************************************************/
/*                                                 */
/*      MEMORY ALLOCATOR                           */
/*                                                 */
/*      JG LEE                                     */
/*      24/NOV/2015                                */
/***************************************************/

#ifndef _MEMALLOC_H
#define _MEMALLOC_H

#include "vega.h"
#include "svc.h"

// ALLOCATE MEMORY
template<class T>
T *alloc_vector(int nl,int nh){
	int i;
    T *v;

    v = new T [nh-nl+1];
    if(!v)  Errmsg("allocation failure in vector()", -1);

	return v-nl;
}

template<class T>    
T **alloc_matrix(int nrl,int nrh,int ncl,int nch){
    int i, j;
    T **m = new T*[nrh - nrl + 1];
	m -= nrl;

	m[nrl] = new T[(nrh - nrl + 1) * (nch - ncl + 1)];
	for (i = 1; i < (nrh - nrl + 1); i++) {
		m[nrl + i] = m[nrl] + i * (nch - ncl + 1);
		m[nrl + i] -= ncl;
	}

    /*m = new T* [nrh-nrl+1];
    if(!m) Errmsg("allocation failure 1 in matrix()", -1);
    m -= nrl;

    for(i=nrl;i<=nrh;i++){
        m[i]=new T [nch-ncl+1];
        if(!m[i]) Errmsg("allocation failure 2 in matrix()", -1);
        m[i] -= ncl;
    }*/

    return m;
}

template<class T>
T ***alloc_tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
	int i, j;
	T ***m;

	m = new T**[nrh - nrl + 1];
	if (!m) Errmsg("allocation failure 1 in tensor()\n", -1);
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = new T*[nch - ncl + 1];
		if (!m[i]) Errmsg("allocation failure 2 in tensor()\n", -1);
		m[i] -= ncl;
	}

	for (i = nrl; i <= nrh; i++) for (j = ncl; j <= nch; j++) {
		m[i][j] = new T[ndh - ndl + 1];
		if (!m[i][j]) Errmsg("allocation failure 3 in tensor()\n", -1);
		m[i][j] -= ndl;
	}

	return m;
}

// FREE the MOMORY
template<class T>
void free_vector(T *v,int nl,int nh){
	delete [] (v+nl);
}

template<class T>
void free_matrix(T **m,int nrl,int nrh,int ncl,int nch){
    //int i;
  
	if (nrh - nrl + 1) {
		delete[](m[-nrl]);
	}
	delete[] m;

    /*for(i=nrh;i>=nrl;i--) delete [] (m[i]+ncl);
    delete [] (m+nrl);*/
}

// ALLOCATING FUNCTIONS
void memalloc_geo(GEODATA& geo);

void memalloc_fld(FIELDDATA& fld, int rN, int zN);

void memalloc_mes(MESHDATA& mes);

void memalloc_sys(SYSALL& sys);



#endif

