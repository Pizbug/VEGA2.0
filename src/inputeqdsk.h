/***************************************************/
/*												   */
/*		INPUT EQDSK HEADER						   */
/*												   */
/*		JG LEE									   */
/*		09/FEB/2017								   */
/***************************************************/

#ifndef _INPUTEQDSK_H
#define _INPUTEQDSK_H

#include "vega.h"
#include "svc.h"
#include "memalloc.h"

void get_id(FILE* in, char id);							
void read_par(SYSALL& sys);
void read_input(SYSALL& sys);


#endif

