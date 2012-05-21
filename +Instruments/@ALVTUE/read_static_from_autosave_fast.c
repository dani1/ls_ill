/*
 * =====================================================================================
 *
 *       Filename:  read_static_from_autosave_fast.c
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
int read_data(double *cr1, double *cr2, double *imon,double *temp,double *angle,  char *path);
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
	double angle, *plhs_angle, temperature, *plhs_temperature;
	double count_rate1, *plhs_cr1, count_rate2,*plhs_cr2 , monitor_intensity, *plhs_imon;
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

	/*  read data from file: */
	i = read_data(&count_rate1, &count_rate2,&monitor_intensity,&temperature,&angle, path);
	if (i == 0)
	mexWarnMsgTxt("File not existent / errors during evaluation of function read_data");
	
	/*  allocate Matlab memory: 3 vectors t,gt, dgt*/
	plhs[0] = mxCreateDoubleMatrix(1, 1 , mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1 , mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1 , mxREAL);
	/*  allocate Matlab memory: 2 variables angle,temperature */
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	/*  get pointer to Matlab angle and temparature variables */
	plhs_angle       = mxGetPr(plhs[3]);
	plhs_temperature = mxGetPr(plhs[4]);
	/*  populate Matlab angle and temperature variables */
	*plhs_angle       = angle;
	*plhs_temperature = temperature;
	/*  get pointer to 3 Matlab vectors t,gt, dgt */
	plhs_cr1  = mxGetPr(plhs[0]);
	plhs_cr2  = mxGetPr(plhs[1]);
	plhs_imon = mxGetPr(plhs[2]);
	/*  populate  3 Matlab count rates and monitor intensity*/
	*plhs_cr1  = count_rate1;
	*plhs_cr2  = count_rate2;
	*plhs_imon = monitor_intensity;

	
}				/* ----------  end of function mexFunction  ---------- */
#endif
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  read_data 
 *  Description:  read dynamic data from DLS instrument ALV autosave
 * =====================================================================================
 */
int read_data(double *cr1, double *cr2, double *imon,double *temp,double *angle,  char *path)
{
	FILE* file_pointer;
	char *str = (char*) malloc(1000 * sizeof(char));
	
	if((file_pointer = fopen(path, "r+")) == NULL)
	{
		return 0;
	}
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
	/*  find count_rate 1 */
	while( strcmp(str, "MeanCR0") != 0)
	{
		fscanf(file_pointer, "%s", str);
	}
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	/*  save count_rate 1   */
	*cr1 = atof(str);
	/*  find count_rate 2 */
	while( strcmp(str, "MeanCR1") != 0)
	{
		fscanf(file_pointer, "%s", str);
	}
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	/*  save count_rate 2   */
	*cr2 = atof(str);
	/*  find monitor diode */
	while( strcmp(str, "Monitor") != 0)
	{
		fscanf(file_pointer, "%s", str);
	}
	fscanf(file_pointer, "%s", str);
	fscanf(file_pointer, "%s", str);
	/*  save monitor diode intensity   */
	*imon = atof(str);

	fclose(file_pointer);
	return 1;
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
	
	char *path = "/Users/daniel/Documents/tesi/data/data_raw/LS/2011_11_04/BSA_1gl_NaCl_200mM0003.ASC";
	double cr1;
	double cr2;
	double imon; 
	double angle, temperature;
	int len;
	len = read_data(&cr1,&cr2,&imon,&temperature, &angle,path);
	printf("results: %f %f %f %f %f", cr1, cr2, imon, temperature, angle);
	return 0;
}				/* ----------  end of function main  ---------- */
