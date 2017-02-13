/***************************************************/
/*                                                 */
/*      SERVICE ROUTINE                    		   */
/*                                                 */
/*      JG LEE                                     */
/*      29/DEC/2015                                */
/***************************************************/

#ifndef	_SVC_H
#define	_SVC_H


#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b)) 
#define sqr(x) ((x)*(x))
#define pow2(x) ((x)*(x))
#define div_n(i,a,b,n) ((a)+((b)-(a))*((i)-1.0)/((double)(n)-1.0))  
#define del_n(a,b,n) (((b)-(a))/((double)(n)-1.0))                
#define distance(x1,y1,x2,y2) (sqrt(((x2)-(x1))*((x2)-(x1))+((y2)-(y1))*((y2)-(y1))))
#define inorout(ax,ay,bx,by,x,y)	(((bx-ax)*(y-ay)-(by-ay)*(x-ax))*y)
				//	possitive -> inner, negative -> outer

#define DoAllDomain(i,n)	for(i=0;i<n;i++)	//	due to C language array ordering
#define	DoDomain(i,a,n)		for(i=a;i<n;i++)
#define PI      (3.14159265358979324)
#define	BTWN(v, a, b)		((v) >= (a) && (v) < (b)) || ((v) < (a) && (v) >= (b))
#define	LINEAR1D(r1,r2,f1,f2,f)		(r1) + ((f)-(f1))*((r2)-(r1))/((f2)-(f1))

#define Errmsg(X,N) {printf("%s\n",X);exit(N);}
#define MAXSW 100
static time_t start_time[MAXSW],end_time[MAXSW];
static double etime[MAXSW];
#define stopwatch(i) etime[i] = clock()
#define howlong(i) etime[i] = (double)(clock()-etime[i])/(double)CLOCKS_PER_SEC

inline double LINEAR(double* f,double x,double* xgrd,int n){
    int i;
    double dx=xgrd[1]-xgrd[0];
    
    if(x<=xgrd[0]  ) return f[0  ];
    if(x>=xgrd[n-1]) return f[n-1];

    i=(int)((x-xgrd[0])/dx);
    x-=xgrd[i];
    
    return (dx-x)/dx*f[i]+x/dx*f[i+1];
}    

//	Hold
inline double SPLINE1D(double* r, int fi, double ff,
					   double f0, double f1, double f2, double f3) {
	double	m[16], inv[16], I[16];
	double	det, f[4];
	double	a = 0, b = 0, c = 0, d = 0;
	double	x0, x1;
	int		i;

	DoAllDomain(i, 4) {
		m[4 * i + 0] = r[fi - 2 + i] * r[fi - 2 + i] * r[fi - 2 + i];
		m[4 * i + 1] = r[fi - 2 + i] * r[fi - 2 + i];
		m[4 * i + 2] = r[fi - 2 + i];
		m[4 * i + 3] = 1;
	}

	inv[0]	= m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15]
			+ m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

	inv[4]	=-m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15]
			- m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

	inv[8]	= m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15]
			+ m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

	inv[12] =-m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14]
			- m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

	inv[1]	=-m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15]
			- m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];

	inv[5]	= m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15]
			+ m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

	inv[9]	=-m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15]
			- m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

	inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14]
			+ m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

	inv[2]	= m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15]
			+ m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];

	inv[6]	=-m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15]
			- m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

	inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15]
			+ m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

	inv[14] =-m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14]
			- m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

	inv[3]	=-m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11]
			- m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];

	inv[7]	= m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11]
			+ m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];

	inv[11] =-m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11]
			- m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10]
			+ m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0)
		Errmsg("No Exists Inverse matrix", -1);

	det = 1.0 / det;

	for (i = 0; i < 16; i++)
		I[i] = inv[i] * det;

	f[0] = f0; f[1] = f1; f[2] = f2; f[3] = f3;

	DoAllDomain(i, 4) {
		a += I[i] * f[i];
		b += I[4 + i] * f[i];
		c += I[8 + i] * f[i];
		d += I[12 + i] * f[i];
	}
	
	x0 = LINEAR1D(r[fi],r[fi+1],f1,f2,ff);

	DoAllDomain(i, 10) {
		x1 = x0 - (a*x0*x0*x0 + b*x0*x0 + c*x0 + d-ff) / (3.0*a*x0*x0 + 2 * b*x0 + c);
		x0 = x1;
	}

	return x1;
}

#endif