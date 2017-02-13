/***************************************************/
/*												   */
/*		SERVICE ROUTINE							   */
/*												   */
/*		JG LEE									   */
/*		15/JAN/2016								   */
/***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "mpi.h"


#define MAXCOLUMN 256

void c2error(char *msg)
{
	printf("%s",msg);
	printf("\n");
	exit(-1);
}

void c2error_mpi(char *msg)
{
	printf("%s",msg);
	printf("\n");
	MPI_Finalize();
}

void get_fillnum(int i,int n,char* c)
{
	int k,kk,ii;
	
	ii=i;
	for(k=1;k<=n;k++) if( (ii/=10)==0 ) break;
	if(ii!=0) c2error("C2 : error in get_fillnum");
	for(kk=0;kk<n-k;kk++) c[kk]='0';
	sprintf(c+n-k,"%d",i);
}

void remove_space(char *str,int n)
{
	char tmp[MAXCOLUMN];
	char *p1,*p2;
	int  len,k;
	
	p1=str;
	p2=tmp;
	len=0;
	while(*p1!='\0') {
		len++;
		if(*p1==' ') p1++;
		else *p2++=*p1++;
	}
	
	for(k=0;k<len;k++) str[k]=tmp[k];
	str[len]='\0';
}	

void upper(char *name)
{
	char *p;
	p=name;
	while(*p!='\0')
	{
		*p=(char)toupper((int)*p);
		p++;
	}
}

void lower(char *name)
{
	char *p;
	p=name;
	while(*p!='\0')
	{
		*p=(char)tolower((int)*p);
		p++;
	}
}

void get_id(FILE *fp,char *id) //get_id(in,"SovlerType")
{
	char str[MAXCOLUMN];

	rewind(fp); //파일의 현재 위치를 시작 위치로 한번에 이동!
 	do {
		fgets(str,256,fp); //char fgets(char* s,int n, FILE* stream) 문자열 입력 함수
		if(str[0]=='@'&&strncmp(str+1,id,strlen(id))==0) break; //strncmp : 문자열 비교 //strlen : 문자열 길이 반환
	}
 	while(!feof(fp)); //feof : 파일의 끝에 도달했는지 검사
 
 	if(feof(fp)) {
		printf("@%s : no item in input file\n",id);
		exit(-1);
	}
}

