/***************************************************/
/*												   */
/*		POST PROCESS HEADER						   */
/*					for C2 Code					   */
/*		JG LEE									   */
/*		20/JUL/2016								   */
/***************************************************/

#ifndef _POSTPROCESS_H
#define _POSTPROCESS_H

#include "vega.h"

void MESHPART(SYSALL &sys);
void ASSEMBLER(SYSALL &sys, int id, int nS,	int nW, int nH, int &k, int dW, int dH, int jump);
void PRINTMESH(SYSALL &sys, int id, int dom, int k, int nH, int nW, ivector b, ivector t);
void ONEDPRINT(SYSALL &sys, int id, int nH, int nIS, int nOS, int Nedge);
void POSTP(SYSALL &sys);
void POSTP_SNOW(SYSALL &sys);
void POSTP_DDNU(SYSALL &sys);
void POSTP_DDNL(SYSALL &sys);

#endif

