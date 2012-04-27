/*
 * =====================================================================================
 *
 *       Filename:  read_dynamic_file_fast.c
 *
 *    Description:  read correlation data from autosave files created by ALV Light Scattering Instrument 
 *
 *        Version:  1.0
 *        Created:  22.12.2011 16:07:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Daniel Soraruf (), daniel.soraruf@gmail.com
 *        Company:  
 *
 * =====================================================================================
 */



#include	<stdio.h>
#include 	<stdlib.h>
#include 	<string.h>

/* check if calling from matlab and include necessary mex.h*/
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/*  define max_length of correlation data */
#define MAX_CORR_VECTOR_LENGTH 1000
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mexFunction 
 *  Description:  comunicate between matlab and c program
 * =====================================================================================
 */
int read_data(double *t, double *gt, double *dgt,double *temp,double *angle,  char *path);
#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, 
		mxArray *plhs[], 
		int nrhs, 
		const mxArray *prhs[])
{

	char *path;
	int buflen;
	double *gt, *plhs_gt;
	double *dgt, *plhs_dgt;
	double *t, *plhs_t;
	double angle, *plhs_angle, temperature, *plhs_temperature;
	int buf_out_len;
	int status;
	int i ;
	
 	/*  get length of input string */
	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
	/*  allocate memory for path */
	path   = mxMalloc(buflen * sizeof(char));
	/*  get path string from input prhs[0]  */
	status = mxGetString(prhs[0], path, buflen); 
	if (status != 0)
		mexWarnMsgTxt("Not enough space. String is truncated."); 

	/*  allocate memory */
	t   = calloc(MAX_CORR_VECTOR_LENGTH, sizeof(double));
	gt  = calloc(MAX_CORR_VECTOR_LENGTH, sizeof(double));
	dgt = calloc(MAX_CORR_VECTOR_LENGTH, sizeof(double));
	/*  read data from file: */
	/*  buf_out_len = length of correlation data vectors */
	buf_out_len = read_data(t,gt,dgt,&temperature,&angle, path);
	
	/*  allocate Matlab memory: 3 vectors t,gt, dgt*/
	plhs[0] = mxCreateDoubleMatrix(buf_out_len, 1 , mxREAL);
	plhs[1] = mxCreateDoubleMatrix(buf_out_len, 1 , mxREAL);
	plhs[2] = mxCreateDoubleMatrix(buf_out_len, 1 , mxREAL);
	/*  allocate Matlab memory: 2 variables angle,temperature */
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	/*  get pointer to Matlab angle and temparature variables */
	plhs_angle = mxGetPr(plhs[3]);
	plhs_temperature = mxGetPr(plhs[4]);
	/*  populate Matlab angle and temperature variables */
	*plhs_angle = angle;
	*plhs_temperature = temperature;
	/*  get pointer to 3 Matlab vectors t,gt, dgt */
	plhs_t = mxGetPr(plhs[0]);
	plhs_gt = mxGetPr(plhs[1]);
	plhs_dgt = mxGetPr(plhs[2]);
	/*  populate  3 Matlab vectors t,gt, dgt*/
	for ( i = 0 ; i < buf_out_len; i++)
	{
		plhs_t[i] = t[i];
		plhs_gt[i] = gt[i];
		plhs_dgt[i] = dgt[i];
	}
}				/* ----------  end of function mexFunction  ---------- */
#endif
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_data 
 *  Description:  read dynamic data from DLS instrument ALV autosave
 * =====================================================================================
 */
int read_data(double *t, double *gt, double *dgt,double *temp,double *angle,  char *path)
{
	int read_check;
	FILE* file_pointer;
	char *str = (char*) malloc(1000 * sizeof(char));
	float tmp_float;
	int number_of_columns = 5;
	/*  ASSUMPTION : col_number_correlation != 0  !!!*/
	const int col_number_time = 0;
	const int col_number_correlation = 1;
	int loop_col_number_time = col_number_time;
	/* ASSUMPTION std_dev in separate list in 2 columns 0 = time, 1 = std_dev */
	//const int col_number_std_dev = 1;
	int time_index = 0;
	int correlation_index = 0;
	int std_dev_index = 0;
	int i;
	file_pointer = fopen(path, "r+");
	/*find temperature*/
	while( strcmp(str, "Temperature") != 0)
	{
		read_check=fscanf(file_pointer, "%s", str);
	}
	read_check=fscanf(file_pointer, "%s", str);
	read_check=fscanf(file_pointer, "%s", str);
	read_check=fscanf(file_pointer, "%s", str);
	/*  save temperature */
	*temp = atof(str);
//	printf("%lf\n", *temp);
	/*  find angle */
	while( strcmp(str, "Angle") != 0)
	{
		read_check=fscanf(file_pointer, "%s", str);
	}
	read_check=fscanf(file_pointer, "%s", str);
	read_check=fscanf(file_pointer, "%s", str);
	read_check=fscanf(file_pointer, "%s", str);
	/*  save angle   */
	*angle = atof(str);
//	printf("%lf\n", *angle);
	while( strcmp(str, "\"Correlation\"") != 0 )
	{
		read_check=fscanf(file_pointer, "%s", str);
//		printf("%s\n", str);
	}
	/*  read first column of first row (for speed up, less strcmp) */
	read_check=fscanf(file_pointer,"%s", str);
	if (col_number_time == 0)
	{
		t[time_index++] = atof(str); 
		loop_col_number_time = number_of_columns;
	}
	/* first column already read */
	while( strcmp(str, "\"Count") != 0 )
	{
		for ( i = 0 ; i < number_of_columns; i++)
		{
			if ( i == col_number_correlation - 1)
			{
			 	/*  -1 necessary since i = col_number + 1 for speed up */
				read_check=fscanf(file_pointer,"%lf", & gt[correlation_index++]);
			}
			else
			{
				if (i == loop_col_number_time - 1)
				{
				/*  read_check=fscanf at position  col_number_time (if col_number_time ==0 -> next column)  */
					read_check=fscanf(file_pointer,"%s", str);
					t[time_index++] = atof(str);
				}
				else
				{
					read_check=fscanf(file_pointer, "%f", & tmp_float);
				}
			}
		}
	}
	
	while( strcmp(str, "\"StandardDeviation\"") != 0 && !feof(file_pointer) )
	{
		read_check=fscanf(file_pointer, "%s", str);
//		printf("%s\n", str);
	}
	while(!feof(file_pointer))
	{
		read_check=fscanf(file_pointer, "%f", & tmp_float);
		read_check=fscanf(file_pointer, "%lf",& dgt[std_dev_index++]);
	}
	fclose(file_pointer);
	time_index --;
	/*  check whether std_dev in file, if not fill with ones */
	if (! std_dev_index ) 
	{
		for ( i = 0 ; i < time_index ; i++)
		{
			dgt[i] = 1;
		}
	}
	else
	{
		std_dev_index --;
	}
	return time_index;
}				/* ----------  end of function read_data  ---------- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  For test / launch purposes
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	
	char *path = "/Users/daniel/Documents/tesi/data/data_raw/LS/2012_03_02/BSA_5gl_Nosalt0000_0001.ASC";
	double *t = (double * ) malloc(MAX_CORR_VECTOR_LENGTH * sizeof(double));
	double *gt = (double * ) malloc(MAX_CORR_VECTOR_LENGTH * sizeof(double));
	double *dgt = (double * ) malloc(MAX_CORR_VECTOR_LENGTH * sizeof(double));
	double angle, temperature;
	int len, i;
	i = 1;
	len = read_data(t,gt,dgt,&temperature, &angle,path);
	printf("\n\nangle : %lf \ntemperature: %lf\nlen : %d\n", angle, temperature, len);
	/*  head   */
	for ( i = 0 ; i < 5 ; i++) 
		printf("\n %lf %lf %lf" , t[i], gt[i], dgt[0]);
	printf("\n");
	/*  tail */
	for ( i = len-5 ; i < len ; i++) 
		printf("\n %lf %lf %lf" , t[i], gt[i], dgt[i]);
	printf("\n");
	return 0;
}				/* ----------  end of function main  ---------- */
