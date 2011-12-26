/*
 * =====================================================================================
 *
 *       Filename:  read_static_from_autosave.c
 *
 *    Description:  read static count rate and monitor diode intensity from autosave file of ALV light scattering instrument.
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

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mexFunction 
 *  Description:  comunicate between matlab and c program
 * =====================================================================================
 */
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

	
//	t   = mxGetPr(plhs[0]);
//	gt  = mxGetPr(plhs[1]);
//	dgt = mxGetPr(plhs[2]);
	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
	path = mxMalloc(buflen * sizeof(char));
	status = mxGetString(prhs[0], path, buflen); 
	if (status != 0)
		mexWarnMsgTxt("Not enough space. String is truncated."); 
//	gt  = mxMalloc(1000 * sizeof(double));
//	dgt = mxMalloc(1000 * sizeof(double));
//	t   = mxMalloc(1000 * sizeof(double));
	t = calloc(1000, sizeof(double));
	gt = calloc(1000, sizeof(double));
	dgt = calloc(1000, sizeof(double));
	buf_out_len = read_data(t,gt,dgt,&temperature,&angle, path);

	plhs[0] = mxCreateDoubleMatrix(buf_out_len, 1 , mxREAL);
	plhs[1] = mxCreateDoubleMatrix(buf_out_len, 1 , mxREAL);
	plhs[2] = mxCreateDoubleMatrix(buf_out_len, 1 , mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	
	plhs_angle = mxGetPr(plhs[3]);
	plhs_temperature = mxGetPr(plhs[4]);

	*plhs_angle = angle;
	*plhs_temperature = temperature;
	
	plhs_t = mxGetPr(plhs[0]);
	plhs_gt = mxGetPr(plhs[1]);
	plhs_dgt = mxGetPr(plhs[2]);
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
	FILE* file_pointer;
	char *str = (char*) malloc(1000 * sizeof(char));
	float tmp_float;
	int number_of_columns = 5;
	int col_number_time = 0;
	int col_number_correloation = 1;
	int col_number_std_dev = 1;
	int time_index = 0;
	int correlation_index = 0;
	int std_dev_index = 0;
	int i;
	file_pointer = fopen(path, "r+");
	/*find temperature*/
	while( strcmp(str, "Temperature") != 0)
	{
		fscanf(file_pointer, "%s", str);
	}
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	/*  save temperature */
	*temp = atof(str);
//	printf("%lf\n", *temp);
	/*  find angle */
	while( strcmp(str, "Angle") != 0)
	{
		fscanf(file_pointer, "%s", str);
	}
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	/*  save angle   */
	*angle = atof(str);
//	printf("%lf\n", *angle);
	while( strcmp(str, "\"Correlation\"") != 0 )
	{
		fscanf(file_pointer, "%s", str);
//		printf("%s\n", str);
	}
	fscanf(file_pointer,"%s", str);
	t[time_index++] = atof(str); 
	while( strcmp(str, "\"Count") != 0 )
	{
		for ( i = 0 ; i < number_of_columns; i++)
		{
			if ( i == 0)
			{
				fscanf(file_pointer,"%lf", & gt[correlation_index++]);
//				printf(" gt[%d]: %lf\n",correlation_index -1, \
//					gt[correlation_index -1 ]);
			}
			else
			{
				if (i == number_of_columns - 1)
				{
					fscanf(file_pointer,"%s", str);
					t[time_index++] = atof(str);
//					printf(" t[%d]: %lf \t ", time_index-1,\
//							t[time_index -1]);
				}
				else
				{
					fscanf(file_pointer, "%f", & tmp_float);
				}
			}
		}
	}
	
	while( strcmp(str, "\"StandardDeviation\"") != 0 )
	{
		fscanf(file_pointer, "%s", str);
//		printf("%s\n", str);
	}
	while(!feof(file_pointer))
	{
		fscanf(file_pointer, "%f", & tmp_float);
		fscanf(file_pointer, "%lf",& dgt[std_dev_index++]);
//		printf("t[%d] = %lf ; dgt[%d] = %lf\n", std_dev_index -1, t[std_dev_index-1],std_dev_index -1,  dgt[std_dev_index-1]);
	}
	fclose(file_pointer);
	time_index --;
	std_dev_index --;
//	printf("%d %d %d", std_dev_index, correlation_index, time_index);
	return time_index;
}				/* ----------  end of function read_data  ---------- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	char *path = "/Users/daniel/Documents/tesi/data/data_raw/LS/2011_11_04/BSA_1gl_NaCl_200mM0003.ASC";
	double *t = (double * ) malloc(1000 * sizeof(double));
	double *gt = (double * ) malloc(1000 * sizeof(double));
	double *dgt = (double * ) malloc(1000 * sizeof(double));
	double angle, temperature;
	int len;
	len = read_data(t,gt,dgt,&temperature, &angle,path);
	printf("\n\nangle : %lf \ntemperature: %lf\nlen : %d\n", angle, temperature, len);
	return 0;
}				/* ----------  end of function main  ---------- */
