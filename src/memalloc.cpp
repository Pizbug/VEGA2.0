/***************************************************/
/*                                                 */
/*      MEMORY ALLOCATOR                           */
/*                                                 */
/*      JG LEE                                     */
/*      24/NOV/2015                                */
/***************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "parameter.h"
#include "vega.h"
#include "memalloc.h"


// ALLOCATING FUNCTIONS
void memalloc_geo(GEODATA& geo) {
	geo.r	= alloc_vector<double>(0, geo.rN-1);
	geo.z	= alloc_vector<double>(0, geo.zN-1);
	geo.Psi = alloc_matrix<double>(0, geo.rN-1, 0, geo.zN-1);

	geo.Norpsi	= alloc_vector<double>(0, geo.rN-1);
	geo.Bt		= alloc_vector<double>(0, geo.rN-1);
	geo.Bt0		= alloc_vector<double>(0, geo.rN - 1);
	geo.Lmt_r	= alloc_vector<double>(0, geo.lnum-1);
	geo.Lmt_z	= alloc_vector<double>(0, geo.lnum-1);
	geo.Pfc_r	= alloc_vector<double>(0, 1);   //  # of PFC is 2
	geo.Pfc_z	= alloc_vector<double>(0, 1);
	geo.Divt_r	= alloc_matrix<double>(0, geo.ndiv-1, 0, geo.dnum-1);
	geo.Divt_z	= alloc_matrix<double>(0, geo.ndiv-1, 0, geo.dnum-1);
}

void memalloc_fld(FIELDDATA& fld, int rN, int zN) {
	fld.Br = alloc_matrix<double>(0, rN-1, 0, zN-1);
	fld.Bz = alloc_matrix<double>(0, rN-1, 0, zN-1);
	fld.Bp = alloc_matrix<double>(0, rN-1, 0, zN-1);

	fld.Psi00 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Psi10 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Psi01 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Psi11 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);

	fld.Br00 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Br10 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Br01 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Br11 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);

	fld.Bz00 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Bz10 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Bz01 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);
	fld.Bz11 = alloc_matrix<double>(0, rN - 2, 0, zN - 2);

	fld.dXp = alloc_vector<TESTPOINT>(0, 1);
}

void memalloc_mes(MESHDATA& mes) {
	//	Sr: allocated at FRAME in vega.cpp
	//	Sz: allocated at FRAME in vega.cpp

	mes.Psi  = alloc_vector<double>(0, MGD);
	mes.nD	 = alloc_vector<int>(0, MGD);
	mes.dP   = alloc_vector<double>(0, MGD);
	mes.EDS  = alloc_matrix<double>(0, MGD, 0, 2);

	mes.Mr = alloc_matrix<double>(0, 2 * HALF - 1, 0, KMAX - 1);
	mes.Mz = alloc_matrix<double>(0, 2 * HALF - 1, 0, KMAX - 1);

	mes.Brm = alloc_matrix<double>(0, 2 * HALF - 1, 0, KMAX - 1);
	mes.Bzm = alloc_matrix<double>(0, 2 * HALF - 1, 0, KMAX - 1);
	mes.Btm = alloc_matrix<double>(0, 2 * HALF - 1, 0, KMAX - 1);
}

void memalloc_sys(SYSALL& sys) {
	memalloc_geo(sys.geo);
	memalloc_fld(sys.fld, sys.geo.rN, sys.geo.zN);
	memalloc_mes(sys.mes);

}
