/***************************************************/
/*												   */
/*		VEGA HEADER								   */
/*												   */
/*		JG LEE									   */
/*		28/DEC/2015								   */
/***************************************************/

#ifndef	_VEGA_H
#define	_VEGA_H

typedef double*		dvector; // 1-D pointer variable
typedef double**	dmatrix; // 2-D pointer variable
typedef double***	dtensor; // 3-D pointer variable
typedef int*		ivector;
typedef int**		imatrix;


typedef struct TESTPOINT{
	int 	id 		;
	double 	r  		;
	double	z  		;
	double	v  		;	//	Value
	int 	i  		;
	int 	j  		;
	
	TESTPOINT(){
		id = 0		;
		r  = 0.0	;
		z  = 0.0	;
		v  = 0.0	;
		i  = 0		;
		j  = 0		;
	};
}TESTPOINT;


struct GEODATA{
	int 	rN, zN		;
	dvector	r,z 		;
	dmatrix	Psi			;

	dvector	Norpsi		;	//	Normalized Psi
	dvector	Bt			;	//	Toroidal B field
	dvector Bt0			;	//	Toroidal vacuum B field

	int 	lnum		;	//	# of limiter
	dvector	Lmt_r		;	//	Limiter
	dvector	Lmt_z		;
	dvector	Pfc_r		;	//	Plasma Face Component
	dvector	Pfc_z		;

	int		ndiv		;	//	# of Divertor
	int		dnum		;	//	It should be # of NODES of EACH divertor
	dmatrix	Divt_r		;	//	matrix -> Non-singular divertors
	dmatrix	Divt_z		;
};

struct FIELDDATA{
	dmatrix	Br			;	//	B = Partial psi
	dmatrix	Bz			;
	dmatrix	Bp			;

	dmatrix	Psi00		;	//	Bilinear Coefficients
	dmatrix	Psi10		;
	dmatrix	Psi01		;	
	dmatrix	Psi11		;

	dmatrix	Br00		;
	dmatrix	Br10		;	
	dmatrix	Br01		;
	dmatrix	Br11		;

	dmatrix	Bz00		;
	dmatrix	Bz10		;
	dmatrix	Bz01		;	
	dmatrix	Bz11		;

	int 	xpn			;
	int 	opn			;
	TESTPOINT*	Xp	 	;	//	#, x, y, v, i, j
	TESTPOINT*	Op		;
	TESTPOINT*	dXp	;	//	Using for snowflake divertor
	int 	ro			;	//	Direction of plasma current (1 : counter-clockwise, -1 : clockwise)
};

struct MESHDATA{
	int		Sn			;
	dmatrix	Sr			;	//	Separatrix x
	dmatrix	Sz 			;
	dmatrix	MaxGr		;	//	Maximum Gradient points near X point
	dmatrix	MaxGz		;

	TESTPOINT	Type	;	//	type
							//	id -> 0 : SN, 1 : CDN, 2 : DDN
							//	r, z, i, j -> most far away X point
	dvector Psi			;
	ivector	nD			;	//	Number of points in each separatrix, Max 6 array
	dvector dP			;	//	Psi difference

	dmatrix EDS			;	//	EDS Factors; SIZE = [MGD][3]

	dvector GS			;	//	Grid Scope

	dmatrix	Mr			;
	dmatrix Mz			;

	dmatrix	Brm			;	//	Interpolated Br in mesh point
	dmatrix Bzm			;	//	
	dmatrix Btm			;

	dtensor SMr			;
	dtensor SMz			;
	dtensor SM			;
};

struct SYSALL{
	GEODATA		geo 	;
	FIELDDATA	fld 	;
	MESHDATA	mes 	;
};

/*
#define usealias(g,f,m,sys) \ // s, m, g = sys sub_structure address
	GEODATA&	g 	= sys.geo;	\
	FIELDDATA&	f 	= sys.fld;	\
	MESHDATA&	m 	= sys.mes;	\
	int rN = sys.geo.rN; 		\
	int zN = sys.geo.zN;
	*/

#define usealias_geo(g,sys) 	\
	GEODATA&	g = sys.geo; 	\
	int rN = sys.geo.rN; 		\
	int zN = sys.geo.zN;

#define usealias_fld(f,sys)		\
	FIELDDATA&	f	= sys.fld;

#define usealias_mes(m,sys)		\
	MESHDATA&	m	= sys.mes;

enum {
	Co, 	//	Core ; number, difference
	IBdy,	//	Inner Body ; number
	OBdy,	//	Outer Body
	IDt,	//	Inner Divertor ; number
	ODt,	//	Outer Divertor
	SOL,	//	SOL ; number, psi, difference
	AS, 	//	Inactive SOL; 
	Pr,	//	Private ; number, psi, difference
	APr,
	MGD		//	Max number of Grid distribution factor
};

enum {
	SN,		//	0	Single Null
	CDN,	//	1	Connected Double Null
	DDN,	//	2	Disconnected Double Null
	LMT,	//	3	Limiter plasma
};

double DERIV_R(int rN, dvector r, dvector z, dmatrix V, int i, int j);
double DERIV_Z(int zN, dvector r, dvector z, dmatrix V, int i, int j);
void DERIV(int rN, int zN, dvector r, dvector z, dmatrix V, dmatrix dVr, dmatrix dVz);

TESTPOINT CROSSFUNC(double r1, double z1, double r2, double z2,
	double br1, double bz1, double br2, double bz2);
TESTPOINT GETPOINT(SYSALL & sys, double r1, double z1,
	double r2, double z2, double psi, double ratio);

double BiliInterp(double **f, double *x, double *y,
	double xx, double yy, int i, int j);

#endif
