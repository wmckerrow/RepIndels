"""
Copyright Wilson McKerrow, 2018
This file is part of RepIndel.
RepIndel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
RepIndel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with RepIndel.  If not, see <http://www.gnu.org/licenses/>.
"""

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>


static PyObject* PEhmm2_forward_sum (PyObject *self, PyObject *args)
{
	/* convert Python arguments */
	PyObject *sequence_arg=NULL,*genome_profile_arg=NULL,*del_open_arg=NULL,*del_extend_arg=NULL,*ins_open_arg=NULL,*ins_extend_arg=NULL,*fm_init_arg=NULL,*fi_init_arg=NULL;
	PyObject *sequence_arr=NULL,*genome_profile_arr=NULL,*del_open_arr=NULL,*del_extend_arr=NULL,*ins_open_arr=NULL,*ins_extend_arr=NULL,*fm_init_arr=NULL,*fi_init_arr=NULL;
	
	PyObject *fm_arg=NULL,*fd_arg=NULL,*fi_arg=NULL;
	PyObject *fm_arr=NULL,*fd_arr=NULL,*fi_arr=NULL;
	
	if (!PyArg_ParseTuple(args, "OOOOOOOOO!O!O!", &sequence_arg, &genome_profile_arg, &del_open_arg, &del_extend_arg, &ins_open_arg, &ins_extend_arg, &fm_init_arg, &fi_init_arg,
		&PyArray_Type, &fm_arg,&PyArray_Type, &fd_arg,&PyArray_Type, &fi_arg)) return NULL;
	
	sequence_arr = PyArray_FROM_OTF(sequence_arg, NPY_INT, NPY_IN_ARRAY);
	genome_profile_arr = PyArray_FROM_OTF(genome_profile_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_open_arr = PyArray_FROM_OTF(del_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_extend_arr = PyArray_FROM_OTF(del_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_open_arr = PyArray_FROM_OTF(ins_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_extend_arr = PyArray_FROM_OTF(ins_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	fm_init_arr = PyArray_FROM_OTF(fm_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	fi_init_arr = PyArray_FROM_OTF(fi_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	fm_arr = PyArray_FROM_OTF(fm_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	fd_arr = PyArray_FROM_OTF(fd_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	fi_arr = PyArray_FROM_OTF(fi_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
		
	if (sequence_arr == NULL) return Py_BuildValue("i", 1);
	if (genome_profile_arr == NULL) return Py_BuildValue("i", 2);
	if (del_open_arr == NULL) return Py_BuildValue("i", 3);
	if (del_extend_arr == NULL) return Py_BuildValue("i", 4);
	if (ins_open_arr == NULL) return Py_BuildValue("i", 5);
	if (ins_extend_arr == NULL) return Py_BuildValue("i", 6);
	if (fm_init_arr == NULL) return Py_BuildValue("i", 7);
	if (fi_init_arr == NULL) return Py_BuildValue("i", 8);
	if (fm_arr == NULL) return Py_BuildValue("i", 9);
	if (fd_arr == NULL) return Py_BuildValue("i", 10);
	if (fi_arr == NULL) return Py_BuildValue("i", 11);
	
	/* Get c pointers and length */
	int *sequence;
	double *genome_profile, *del_open, *del_extend, *ins_open, *ins_extend, *fm_init, *fi_init;
	sequence = (int *)PyArray_DATA(sequence_arr);
	genome_profile = (double *)PyArray_DATA(genome_profile_arr);
	del_open = (double *)PyArray_DATA(del_open_arr);
	del_extend = (double *)PyArray_DATA(del_extend_arr);
	ins_open = (double *)PyArray_DATA(ins_open_arr);
	ins_extend = (double *)PyArray_DATA(ins_extend_arr);
	fm_init = (double *)PyArray_DATA(fm_init_arr);
	fi_init = (double *)PyArray_DATA(fi_init_arr);
	
	double *fm, *fd, *fi;
	fm = (double *)PyArray_DATA(fm_arr);
	fd = (double *)PyArray_DATA(fd_arr);
	fi = (double *)PyArray_DATA(fi_arr);
	
	int seqlength = PyArray_DIMS(sequence_arr)[0];
	int prolength = PyArray_DIMS(genome_profile_arr)[0];
	if (prolength != PyArray_DIMS(del_open_arr)[0]) return Py_BuildValue("i", 12);
	if (prolength != PyArray_DIMS(del_extend_arr)[0]) return Py_BuildValue("i", 13);
	if (prolength != PyArray_DIMS(ins_open_arr)[0]) return Py_BuildValue("i", 14);
	if (prolength != PyArray_DIMS(ins_extend_arr)[0]) return Py_BuildValue("i", 15);
	if (prolength != PyArray_DIMS(fm_init_arr)[0]) return Py_BuildValue("i", 16);
	if (prolength != PyArray_DIMS(fi_init_arr)[0]) return Py_BuildValue("i", 17);
	if (prolength != PyArray_DIMS(fm_arr)[0]) return Py_BuildValue("i", 18);
	if (seqlength != PyArray_DIMS(fm_arr)[1]) return Py_BuildValue("i", 19);
	
	/* do function */	
	int i,j;
	for (i = 0; i<prolength; i++)
	{
		fm[i*seqlength+0] = fm_init[i]*genome_profile[i*4 + sequence[0]];
		fd[i*seqlength+0] = 0.0;
		fi[i*seqlength+0] = fi_init[i]*0.25;
	}
	for (j=1; j<seqlength; j++)
	{
		fi[0*seqlength+j] = 0.25*(fm[0*seqlength+j-1]*ins_open[0] + fi[0*seqlength+j-1]*ins_extend[0]);
	}
	for (i=1; i<prolength; i++)
	{
		fd[i*seqlength+0] = fm[(i-1)*seqlength+0]*del_open[i] + fd[(i-1)*seqlength+0]*del_extend[i];
		for (j=1; j<seqlength; j++) {
			fm[i*seqlength+j] = ( fm[(i-1)*seqlength+j-1]*(1-del_open[i]-ins_open[i-1]) + fd[(i-1)*seqlength+j-1]*(1.0-del_extend[i]) + fi[(i-1)*seqlength+j-1]*(1.0-ins_extend[i-1]) )* genome_profile[i*4 + sequence[j]];
			fd[i*seqlength+j] = fm[(i-1)*seqlength+j]*del_open[i] + fd[(i-1)*seqlength+j]*del_extend[i];
			fi[i*seqlength+j] = 0.25*(fm[i*seqlength+j-1]*ins_open[i] + fi[i*seqlength+j-1]*ins_extend[i]);
		}
	}
		
	/* return something */
	Py_DECREF(sequence_arr);
	Py_DECREF(genome_profile_arr);
	Py_DECREF(del_open_arr);
	Py_DECREF(del_extend_arr);
	Py_DECREF(ins_open_arr);
	Py_DECREF(ins_extend_arr);
	Py_DECREF(fm_init_arr);
	Py_DECREF(fi_init_arr);
	Py_DECREF(fm_arr);
	Py_DECREF(fd_arr);
	Py_DECREF(fi_arr);
	Py_INCREF(Py_None);
	return Py_None;
	
}

static PyObject* PEhmm2_forward_sum_interior (PyObject *self, PyObject *args)
{	
	PyObject *fm_init_arg=NULL, *fi_init_arg=NULL, *del_open_arg=NULL, *del_extend_arg=NULL, *ins_open_arg=NULL, *ins_extend_arg=NULL;
	PyObject *fm_init_arr=NULL, *fi_init_arr=NULL, *del_open_arr=NULL, *del_extend_arr=NULL, *ins_open_arr=NULL, *ins_extend_arr=NULL;
	
	int size;
	double p;
	
	PyObject *fm_arg=NULL,*fd_arg=NULL,*fi_arg=NULL;
	PyObject *fm_arr=NULL,*fd_arr=NULL,*fi_arr=NULL;
	
	if (!PyArg_ParseTuple(args, "OOOOOOidO!O!O!", &fm_init_arg, &fi_init_arg, &del_open_arg, &del_extend_arg, &ins_open_arg, &ins_extend_arg, &size, &p,
		&PyArray_Type, &fm_arg,&PyArray_Type, &fd_arg,&PyArray_Type, &fi_arg)) return NULL;
	
	fm_init_arr = PyArray_FROM_OTF(fm_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	fi_init_arr = PyArray_FROM_OTF(fi_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_open_arr = PyArray_FROM_OTF(del_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_extend_arr = PyArray_FROM_OTF(del_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_open_arr = PyArray_FROM_OTF(ins_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_extend_arr = PyArray_FROM_OTF(ins_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	fm_arr = PyArray_FROM_OTF(fm_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	fd_arr = PyArray_FROM_OTF(fd_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	fi_arr = PyArray_FROM_OTF(fi_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	
	if (fm_init_arr == NULL) return Py_BuildValue("i", 1);
	if (fi_init_arr == NULL) return Py_BuildValue("i", 2);
	if (del_open_arr == NULL) return Py_BuildValue("i", 3);
	if (del_extend_arr == NULL) return Py_BuildValue("i", 4);
	if (ins_open_arr == NULL) return Py_BuildValue("i", 5);
	if (ins_extend_arr == NULL) return Py_BuildValue("i", 6);
	if (fm_arr == NULL) return Py_BuildValue("i", 9);
	if (fd_arr == NULL) return Py_BuildValue("i", 10);
	if (fi_arr == NULL) return Py_BuildValue("i", 11);	
	
	/* Get c pointers and length */
	double *fm_init, *fi_init, *del_open, *del_extend, *ins_open, *ins_extend;
	fm_init = (double *)PyArray_DATA(fm_init_arr);
	fi_init = (double *)PyArray_DATA(fi_init_arr);
	del_open = (double *)PyArray_DATA(del_open_arr);
	del_extend = (double *)PyArray_DATA(del_extend_arr);
	ins_open = (double *)PyArray_DATA(ins_open_arr);
	ins_extend = (double *)PyArray_DATA(ins_extend_arr);
	
	double *fm, *fd, *fi;
	fm = (double *)PyArray_DATA(fm_arr);
	fd = (double *)PyArray_DATA(fd_arr);
	fi = (double *)PyArray_DATA(fi_arr);
	
	int windowlength = PyArray_DIMS(fm_init_arr)[0];
	if (windowlength != PyArray_DIMS(fm_init_arr)[0]) return Py_BuildValue("i", 12);
	if (windowlength != PyArray_DIMS(fi_init_arr)[0]) return Py_BuildValue("i", 13);
	if (windowlength != PyArray_DIMS(del_open_arr)[0]) return Py_BuildValue("i", 14);
	if (windowlength != PyArray_DIMS(del_extend_arr)[0]) return Py_BuildValue("i", 15);
	if (windowlength != PyArray_DIMS(ins_open_arr)[0]) return Py_BuildValue("i", 16);
	if (windowlength != PyArray_DIMS(ins_extend_arr)[0]) return Py_BuildValue("i", 17);
	if (windowlength != PyArray_DIMS(fm_arr)[0]) return Py_BuildValue("i", 18);
	if (windowlength != PyArray_DIMS(fd_arr)[0]) return Py_BuildValue("i", 19);
	if (windowlength != PyArray_DIMS(fi_arr)[0]) return Py_BuildValue("i", 20);
	if (size != PyArray_DIMS(fm_arr)[1]) return Py_BuildValue("i", 21);
	if (size != PyArray_DIMS(fd_arr)[1]) return Py_BuildValue("i", 22);
	if (size != PyArray_DIMS(fi_arr)[1]) return Py_BuildValue("i", 23);
	
	/* Do function */
	int i,k;
	
	fi[0*size+0] = fi_init[0]*ins_extend[0];
	for (i=1; i<windowlength; i++)
	{
		fm[i*size+0] = (fm[(i-1)*size+0]*(1-p)+fm_init[i-1])*(1-del_open[i]-ins_open[i-1]) + fd[(i-1)*size+0]*(1-p)*(1-del_extend[i]) + (fi[(i-1)*size+0]*(1-p)+fi_init[i-1])*(1-ins_extend[i-1]);
		fd[i*size+0] = fm[(i-1)*size+0]*del_open[i] + fd[(i-1)*size+0]*del_extend[i];
		fi[i*size+0] = (fm[i*size+0]*(1-p)+fm_init[i])*ins_open[i] + ins_extend[i]*fi_init[i];
	}
	
	for (k=1; k<size; k++)
	{
		fi[0*size+k] = (fm[0*size+k]*(1-p)+fm[0*size+k-1]*p)*ins_open[0] + ins_extend[0]*fi[0*size+(k-1)];
		for (i=1; i<windowlength; i++)
		{
			fm[i*size+k] = (fm[(i-1)*size+k]*(1-p)+fm[(i-1)*size+(k-1)]*p)*(1-del_open[i]-ins_open[i-1]) + (fd[(i-1)*size+k]*(1-p)+fd[(i-1)*size+k-1]*p)*(1-del_extend[i]) + (fi[(i-1)*size+k]*(1-p)+fi[(i-1)*size+k-1]*p)*(1-ins_extend[i-1]);
			fd[i*size+k] = fm[(i-1)*size+k]*del_open[i] + fd[(i-1)*size+k]*del_extend[i];
			fi[i*size+k] = (fm[i*size+k]*(1-p)+fm[i*size+k-1]*p)*ins_open[i] + ins_extend[i]*fi[i*size+(k-1)];
		}
	}
	
	/* Cleanup */
	Py_DECREF(fm_init_arr);
	Py_DECREF(fi_init_arr);
	Py_DECREF(del_open_arr);
	Py_DECREF(del_extend_arr);
	Py_DECREF(ins_open_arr);
	Py_DECREF(ins_extend_arr);
	Py_DECREF(fm_arr);
	Py_DECREF(fd_arr);
	Py_DECREF(fi_arr);
	Py_INCREF(Py_None);
	
	return Py_None;
	
}

static PyObject* PEhmm2_backward_sum (PyObject *self, PyObject *args)
{	
	/* convert Python arguments */
	PyObject *sequence_arg=NULL,*genome_profile_arg=NULL,*del_open_arg=NULL,*del_extend_arg=NULL,*ins_open_arg=NULL,*ins_extend_arg=NULL,*bm_init_arg=NULL,*bi_init_arg=NULL;
	PyObject *sequence_arr=NULL,*genome_profile_arr=NULL,*del_open_arr=NULL,*del_extend_arr=NULL,*ins_open_arr=NULL,*ins_extend_arr=NULL,*bm_init_arr=NULL,*bi_init_arr=NULL;
	
	PyObject *bm_arg=NULL,*bd_arg=NULL,*bi_arg=NULL;
	PyObject *bm_arr=NULL,*bd_arr=NULL,*bi_arr=NULL;
	
	if (!PyArg_ParseTuple(args, "OOOOOOOOO!O!O!", &sequence_arg, &genome_profile_arg, &del_open_arg, &del_extend_arg, &ins_open_arg, &ins_extend_arg, &bm_init_arg, &bi_init_arg,
		&PyArray_Type, &bm_arg,&PyArray_Type, &bd_arg,&PyArray_Type, &bi_arg)) return NULL;
	
	sequence_arr = PyArray_FROM_OTF(sequence_arg, NPY_INT, NPY_IN_ARRAY);
	if (sequence_arr == NULL) return Py_BuildValue("i", 1);
	genome_profile_arr = PyArray_FROM_OTF(genome_profile_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	if (genome_profile_arr == NULL) return Py_BuildValue("i", 2);
	del_open_arr = PyArray_FROM_OTF(del_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	if (del_open_arr == NULL) return Py_BuildValue("i", 3);
	del_extend_arr = PyArray_FROM_OTF(del_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	if (del_extend_arr == NULL) return Py_BuildValue("i", 4);
	ins_open_arr = PyArray_FROM_OTF(ins_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	if (ins_open_arr == NULL) return Py_BuildValue("i", 5);
	ins_extend_arr = PyArray_FROM_OTF(ins_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	if (ins_extend_arr == NULL) return Py_BuildValue("i", 6);
	bm_init_arr = PyArray_FROM_OTF(bm_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	if (bm_init_arr == NULL) return Py_BuildValue("i", 7);
	bi_init_arr = PyArray_FROM_OTF(bi_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	if (bi_init_arr == NULL) return Py_BuildValue("i", 8);
	bm_arr = PyArray_FROM_OTF(bm_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	if (bm_arr == NULL) return Py_BuildValue("i", 9);
	bd_arr = PyArray_FROM_OTF(bd_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	if (bd_arr == NULL) return Py_BuildValue("i", 10);
	bi_arr = PyArray_FROM_OTF(bi_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	if (bi_arr == NULL) return Py_BuildValue("i", 11);
	
	/* Get c pointers and length */
	int *sequence;
	double *genome_profile, *del_open, *del_extend, *ins_open, *ins_extend, *bm_init, *bi_init;
	sequence = (int *)PyArray_DATA(sequence_arr);
	genome_profile = (double *)PyArray_DATA(genome_profile_arr);
	del_open = (double *)PyArray_DATA(del_open_arr);
	del_extend = (double *)PyArray_DATA(del_extend_arr);
	ins_open = (double *)PyArray_DATA(ins_open_arr);
	ins_extend = (double *)PyArray_DATA(ins_extend_arr);
	bm_init = (double *)PyArray_DATA(bm_init_arr);
	bi_init = (double *)PyArray_DATA(bi_init_arr);
		
	double *bm, *bd, *bi;
	bm = (double *)PyArray_DATA(bm_arr);
	bd = (double *)PyArray_DATA(bd_arr);
	bi = (double *)PyArray_DATA(bi_arr);
	
	int seqlength = PyArray_DIMS(sequence_arr)[0];
	int prolength = PyArray_DIMS(genome_profile_arr)[0];
	if (prolength != PyArray_DIMS(del_open_arr)[0]) return Py_BuildValue("i", 12);
	if (prolength != PyArray_DIMS(del_extend_arr)[0]) return Py_BuildValue("i", 12);
	if (prolength != PyArray_DIMS(ins_open_arr)[0]) return Py_BuildValue("i", 13);
	if (prolength != PyArray_DIMS(ins_extend_arr)[0]) return Py_BuildValue("i", 14);
	if (prolength != PyArray_DIMS(bm_arr)[0]) return Py_BuildValue("i", 15);
	if (seqlength != PyArray_DIMS(bm_arr)[1]) return Py_BuildValue("i", 16);
	
	/* do function */
	int i,j;
	for (i=0;i<prolength;i++)
	{
		bm[i*seqlength+seqlength-1] = bm_init[i];
		bd[i*seqlength+seqlength-1] = 0.0;
		bi[i*seqlength+seqlength-1] = bi_init[i];
	}
	
	for (j=seqlength-2;j>=0;j--)
	{
		bm[(prolength-1)*seqlength+j] = bi[(prolength-1)*seqlength+j+1]*ins_open[prolength-1]*0.25;
		bi[(prolength-1)*seqlength+j] = bi[(prolength-1)*seqlength+j+1]*ins_extend[prolength-1]*0.25;
	}
	for (i=prolength-2;i>=0;i--)
	{
		for (j=seqlength-2;j>=0;j--)
		{
			bm[i*seqlength+j] = bm[(i+1)*seqlength+j+1]*(1-del_open[i+1]-ins_open[i])*genome_profile[(i+1)*4+sequence[j+1]] + bd[(i+1)*seqlength+j]*del_open[i+1] + bi[i*seqlength+j+1]*ins_open[i]*0.25;
			bd[i*seqlength+j] = bm[(i+1)*seqlength+j+1]*(1.0-del_extend[i+1])*genome_profile[(i+1)*4+sequence[j+1]] + bd[(i+1)*seqlength+j]*del_extend[i+1];
			bi[i*seqlength+j] = bm[(i+1)*seqlength+j+1]*(1.0-ins_extend[i])*genome_profile[(i+1)*4+sequence[j+1]] + bi[i*seqlength+j+1]*ins_extend[i]*0.25;
		}
	}
		
	/* return something */
	Py_DECREF(sequence_arr);
	Py_DECREF(genome_profile_arr);
	Py_DECREF(del_open_arr);
	Py_DECREF(del_extend_arr);
	Py_DECREF(ins_open_arr);
	Py_DECREF(ins_extend_arr);
	Py_DECREF(bm_init_arr);
	Py_DECREF(bi_init_arr);
	Py_DECREF(bm_arr);
	Py_DECREF(bd_arr);
	Py_DECREF(bi_arr);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* PEhmm2_backward_sum_interior (PyObject *self, PyObject *args)
{	
	PyObject *bm_init_arg=NULL, *bi_init_arg=NULL, *del_open_arg=NULL, *del_extend_arg=NULL, *ins_open_arg=NULL, *ins_extend_arg=NULL;
	PyObject *bm_init_arr=NULL, *bi_init_arr=NULL, *del_open_arr=NULL, *del_extend_arr=NULL, *ins_open_arr=NULL, *ins_extend_arr=NULL;
	
	int size;
	double p;
	
	PyObject *bm_arg=NULL,*bd_arg=NULL,*bi_arg=NULL;
	PyObject *bm_arr=NULL,*bd_arr=NULL,*bi_arr=NULL;
	
	if (!PyArg_ParseTuple(args, "OOOOOOidO!O!O!", &bm_init_arg, &bi_init_arg, &del_open_arg, &del_extend_arg, &ins_open_arg, &ins_extend_arg, &size, &p,
		&PyArray_Type, &bm_arg,&PyArray_Type, &bd_arg,&PyArray_Type, &bi_arg)) return NULL;
	
	bm_init_arr = PyArray_FROM_OTF(bm_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	bi_init_arr = PyArray_FROM_OTF(bi_init_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_open_arr = PyArray_FROM_OTF(del_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_extend_arr = PyArray_FROM_OTF(del_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_open_arr = PyArray_FROM_OTF(ins_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_extend_arr = PyArray_FROM_OTF(ins_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	bm_arr = PyArray_FROM_OTF(bm_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	bd_arr = PyArray_FROM_OTF(bd_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	bi_arr = PyArray_FROM_OTF(bi_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	
	if (bm_init_arr == NULL) return Py_BuildValue("i", 1);
	if (bi_init_arr == NULL) return Py_BuildValue("i", 2);
	if (del_open_arr == NULL) return Py_BuildValue("i", 3);
	if (del_extend_arr == NULL) return Py_BuildValue("i", 4);
	if (ins_open_arr == NULL) return Py_BuildValue("i", 5);
	if (ins_extend_arr == NULL) return Py_BuildValue("i", 6);
	if (bm_arr == NULL) return Py_BuildValue("i", 9);
	if (bd_arr == NULL) return Py_BuildValue("i", 10);
	if (bi_arr == NULL) return Py_BuildValue("i", 11);	
	
	/* Get c pointers and length */
	double *bm_init, *bi_init, *del_open, *del_extend, *ins_open, *ins_extend;
	bm_init = (double *)PyArray_DATA(bm_init_arr);
	bi_init = (double *)PyArray_DATA(bi_init_arr);
	del_open = (double *)PyArray_DATA(del_open_arr);
	del_extend = (double *)PyArray_DATA(del_extend_arr);
	ins_open = (double *)PyArray_DATA(ins_open_arr);
	ins_extend = (double *)PyArray_DATA(ins_extend_arr);
	
	double *bm, *bd, *bi;
	bm = (double *)PyArray_DATA(bm_arr);
	bd = (double *)PyArray_DATA(bd_arr);
	bi = (double *)PyArray_DATA(bi_arr);
	
	int windowlength = PyArray_DIMS(bm_init_arr)[0];
	if (windowlength == 0) return Py_BuildValue("i", 12);
	if (windowlength != PyArray_DIMS(bi_init_arr)[0]) return Py_BuildValue("i", 13);
	if (windowlength != PyArray_DIMS(del_open_arr)[0]) return Py_BuildValue("i", 14);
	if (windowlength != PyArray_DIMS(del_extend_arr)[0]) return Py_BuildValue("i", 15);
	if (windowlength != PyArray_DIMS(ins_open_arr)[0]) return Py_BuildValue("i", 16);
	if (windowlength != PyArray_DIMS(ins_extend_arr)[0]) return Py_BuildValue("i", 17);
	if (windowlength != PyArray_DIMS(bm_arr)[0]) return Py_BuildValue("i", 18);
	if (windowlength != PyArray_DIMS(bd_arr)[0]) return Py_BuildValue("i", 19);
	if (windowlength != PyArray_DIMS(bi_arr)[0]) return Py_BuildValue("i", 20);
	if (size != PyArray_DIMS(bm_arr)[1]) return Py_BuildValue("i", 21);
	if (size != PyArray_DIMS(bd_arr)[1]) return Py_BuildValue("i", 22);
	if (size != PyArray_DIMS(bi_arr)[1]) return Py_BuildValue("i", 23);
	
	/* Do function */
	int i,k;
	
	bi[(windowlength-1)*size+size-1] = ins_extend[windowlength-1]*bi_init[windowlength-1];
	bm[(windowlength-1)*size+size-1] = bi[(windowlength-1)*size+size-1]*(1-p)*ins_open[windowlength-1];
	for (i=windowlength-2;i>=0;i--)
	{
		bi[i*size+size-1] = (bm[(i+1)*size+size-1]*(1-p)+bm_init[i+1]*p)*(1.0-ins_extend[i]) + ins_extend[i]*bi_init[i];
		bm[i*size+size-1] = (bm[(i+1)*size+size-1]*(1-p)+bm_init[i+1]*p)*(1-del_open[i+1]-ins_open[i]) + bd[(i+1)*size+size-1]*del_open[i+1] + (bi[i*size+size-1]*(1-p)+bi_init[i]*p)*ins_open[i];
		bd[i*size+size-1] = (bm[(i+1)*size+size-1]*(1-p)+bm_init[i+1]*p)*(1.0-del_extend[i+1]) + bd[(i+1)*size+size-1]*del_extend[i+1];
	}
	
	for (k=size-2;k>=0;k--)
	{
		bi[(windowlength-1)*size+k] = ins_extend[windowlength-1]*bi[(windowlength-1)*size+(k+1)];
		bm[(windowlength-1)*size+k] = (bi[(windowlength-1)*size+k]*(1-p)+bi[(windowlength-1)*size+k+1]*p)*ins_open[windowlength-1];
		for (i=windowlength-2;i>=0;i--)
		{
			bi[i*size+k] = (bm[(i+1)*size+k]*(1-p)+bm[(i+1)*size+k+1]*p)*(1.0-ins_extend[i]) + ins_extend[i]*bi[i*size+(k+1)];
			bm[i*size+k] = (bm[(i+1)*size+k]*(1-p)+bm[(i+1)*size+k+1]*p)*(1-del_open[i+1]-ins_open[i]) + bd[(i+1)*size+k]*del_open[i+1] + (bi[i*size+k]*(1-p)+bi[i*size+k+1]*p)*ins_open[i];
			bd[i*size+k] = (bm[(i+1)*size+k]*(1-p)+bm[(i+1)*size+k+1]*p)*(1.0-del_extend[i+1]) + bd[(i+1)*size+k]*del_extend[i+1];
		}
	}
	
	/* Cleanup */
	Py_DECREF(bm_init_arr);
	Py_DECREF(bi_init_arr);
	Py_DECREF(del_open_arr);
	Py_DECREF(del_extend_arr);
	Py_DECREF(ins_open_arr);
	Py_DECREF(ins_extend_arr);
	Py_DECREF(bm_arr);
	Py_DECREF(bd_arr);
	Py_DECREF(bi_arr);
	Py_INCREF(Py_None);
		
	return Py_None;
	
}

static PyObject* PEhmm2_add_counts_interior (PyObject *self, PyObject *args)
{
	PyObject *fm_arg=NULL,*fd_arg=NULL,*fi_arg=NULL,*bm_arg=NULL,*bd_arg=NULL,*bi_arg=NULL,*del_open_arg=NULL,*del_extend_arg=NULL,*ins_open_arg=NULL,*ins_extend_arg=NULL;
	PyObject *fm_arr=NULL,*fd_arr=NULL,*fi_arr=NULL,*bm_arr=NULL,*bd_arr=NULL,*bi_arr=NULL,*del_open_arr=NULL,*del_extend_arr=NULL,*ins_open_arr=NULL,*ins_extend_arr=NULL;
	
	int size;
	double p;
	double pr;
	
	PyObject *match_count_arg=NULL,*del_count_arg=NULL,*ins_count_arg=NULL,*del_open_count_arg=NULL,*del_extend_count_arg=NULL,*ins_open_count_arg=NULL,*ins_extend_count_arg=NULL;
	PyObject *match_count_arr=NULL,*del_count_arr=NULL,*ins_count_arr=NULL,*del_open_count_arr=NULL,*del_extend_count_arr=NULL,*ins_open_count_arr=NULL,*ins_extend_count_arr=NULL;
	
	if (!PyArg_ParseTuple(args, "OOOOOOOOOOiddO!O!O!O!O!O!O!", &fm_arg,&fd_arg,&fi_arg,&bm_arg,&bd_arg,&bi_arg,&del_open_arg,&del_extend_arg,&ins_open_arg,&ins_extend_arg, &size, &p, &pr,
		&PyArray_Type, &match_count_arg,&PyArray_Type, &del_count_arg,&PyArray_Type, &ins_count_arg,&PyArray_Type, &del_open_count_arg,&PyArray_Type, &del_extend_count_arg,&PyArray_Type, &ins_open_count_arg,&PyArray_Type, &ins_extend_count_arg)) return NULL;
	
	fm_arr = PyArray_FROM_OTF(fm_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	fd_arr = PyArray_FROM_OTF(fd_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	fi_arr = PyArray_FROM_OTF(fi_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	bm_arr = PyArray_FROM_OTF(bm_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	bd_arr = PyArray_FROM_OTF(bd_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	bi_arr = PyArray_FROM_OTF(bi_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_open_arr = PyArray_FROM_OTF(del_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	del_extend_arr = PyArray_FROM_OTF(del_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_open_arr = PyArray_FROM_OTF(ins_open_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	ins_extend_arr = PyArray_FROM_OTF(ins_extend_arg, NPY_DOUBLE, NPY_IN_ARRAY);
	
	match_count_arr = PyArray_FROM_OTF(match_count_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	del_count_arr = PyArray_FROM_OTF(del_count_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	ins_count_arr = PyArray_FROM_OTF(ins_count_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	del_open_count_arr = PyArray_FROM_OTF(del_open_count_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	del_extend_count_arr = PyArray_FROM_OTF(del_extend_count_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	ins_open_count_arr = PyArray_FROM_OTF(ins_open_count_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	ins_extend_count_arr = PyArray_FROM_OTF(ins_extend_count_arg, NPY_DOUBLE, NPY_INOUT_ARRAY);
	
	if (fm_arr == NULL) return Py_BuildValue("i", 1);
	if (fd_arr == NULL) return Py_BuildValue("i", 2);
	if (fi_arr == NULL) return Py_BuildValue("i", 3);
	if (bm_arr == NULL) return Py_BuildValue("i", 4);
	if (bd_arr == NULL) return Py_BuildValue("i", 5);
	if (bi_arr == NULL) return Py_BuildValue("i", 6);
	if (del_open_arr == NULL) return Py_BuildValue("i", 7);
	if (del_extend_arr == NULL) return Py_BuildValue("i", 8);
	if (ins_open_arr == NULL) return Py_BuildValue("i", 9);
	if (ins_extend_arr == NULL) return Py_BuildValue("i", 10);
	if (match_count_arr == NULL) return Py_BuildValue("i", 11);
	if (del_count_arr == NULL) return Py_BuildValue("i", 12);
	if (ins_count_arr == NULL) return Py_BuildValue("i", 13);
	if (del_open_count_arr == NULL) return Py_BuildValue("i", 14);
	if (del_extend_count_arr == NULL) return Py_BuildValue("i", 15);
	if (ins_open_count_arr == NULL) return Py_BuildValue("i", 16);
	if (ins_extend_count_arr == NULL) return Py_BuildValue("i", 17);
	
	/* Get c pointers and length */
	double *fm, *fd, *fi, *bm, *bd, *bi, *del_open, *del_extend, *ins_open, *ins_extend;
	fm = (double *)PyArray_DATA(fm_arr);
	fd = (double *)PyArray_DATA(fd_arr);
	fi = (double *)PyArray_DATA(fi_arr);
	bm = (double *)PyArray_DATA(bm_arr);
	bd = (double *)PyArray_DATA(bd_arr);
	bi = (double *)PyArray_DATA(bi_arr);
	del_open = (double *)PyArray_DATA(del_open_arr);
	del_extend = (double *)PyArray_DATA(del_extend_arr);
	ins_open = (double *)PyArray_DATA(ins_open_arr);
	ins_extend = (double *)PyArray_DATA(ins_extend_arr);
	
	double *match_count, *del_count, *ins_count, *del_open_count, *del_extend_count, *ins_open_count, *ins_extend_count;
	match_count = (double *)PyArray_DATA(match_count_arr);
	del_count = (double *)PyArray_DATA(del_count_arr);
	ins_count = (double *)PyArray_DATA(ins_count_arr);
	del_open_count = (double *)PyArray_DATA(del_open_count_arr);
	del_extend_count = (double *)PyArray_DATA(del_extend_count_arr);
	ins_open_count = (double *)PyArray_DATA(ins_open_count_arr);
	ins_extend_count = (double *)PyArray_DATA(ins_extend_count_arr);
	
	int seqlength = PyArray_DIMS(fm_arr)[0];
	if (seqlength != PyArray_DIMS(fd_arr)[0]) return Py_BuildValue("i", 18);
	if (seqlength != PyArray_DIMS(fi_arr)[0]) return Py_BuildValue("i", 19);
	if (seqlength != PyArray_DIMS(bm_arr)[0]) return Py_BuildValue("i", 20);
	if (seqlength != PyArray_DIMS(bd_arr)[0]) return Py_BuildValue("i", 21);
	if (seqlength != PyArray_DIMS(bi_arr)[0]) return Py_BuildValue("i", 22);
	if (seqlength != PyArray_DIMS(del_open_arr)[0]) return Py_BuildValue("i", 23);
	if (seqlength != PyArray_DIMS(del_extend_arr)[0]) return Py_BuildValue("i", 24);
	if (seqlength != PyArray_DIMS(ins_open_arr)[0]) return Py_BuildValue("i", 25);
	if (seqlength != PyArray_DIMS(ins_extend_arr)[0]) return Py_BuildValue("i", 26);
	if (seqlength != PyArray_DIMS(match_count_arr)[0]) return Py_BuildValue("i", 27);
	if (seqlength != PyArray_DIMS(del_count_arr)[0]) return Py_BuildValue("i", 28);
	if (seqlength != PyArray_DIMS(ins_count_arr)[0]) return Py_BuildValue("i", 29);
	if (seqlength != PyArray_DIMS(del_open_count_arr)[0]) return Py_BuildValue("i", 30);
	if (seqlength != PyArray_DIMS(del_extend_count_arr)[0]) return Py_BuildValue("i", 31);
	if (seqlength != PyArray_DIMS(ins_open_count_arr)[0]) return Py_BuildValue("i", 32);
	if (seqlength != PyArray_DIMS(ins_extend_count_arr)[0]) return Py_BuildValue("i", 33);
	
	if (size != PyArray_DIMS(fm_arr)[1]) return Py_BuildValue("i", 34);
	if (size != PyArray_DIMS(fd_arr)[1]) return Py_BuildValue("i", 35);
	if (size != PyArray_DIMS(fi_arr)[1]) return Py_BuildValue("i", 36);
	if (size != PyArray_DIMS(bm_arr)[1]) return Py_BuildValue("i", 37);
	if (size != PyArray_DIMS(bd_arr)[1]) return Py_BuildValue("i", 38);
	if (size != PyArray_DIMS(bi_arr)[1]) return Py_BuildValue("i", 39);
	
	/* Do function */
	int i,k;
	
	for (i=0; i<seqlength-1; i++)
	{
		for (k=0; k<size-1; k++)
		{
			match_count[i] += fm[i*size+k]*bm[i*size+k]/pr;
			del_open_count[i+1] += fm[i*size+k]*del_open[i+1]*bd[(i+1)*size+k]/pr;
			ins_open_count[i] += fm[i*size+k]*ins_open[i]*(1-p)*bi[i*size+k]/pr + fm[i*size+k]*ins_open[i]*p*bi[i*size+k+1]/pr;
			del_count[i] += fd[i*size+k]*bd[i*size+k]/pr;
			del_extend_count[i+1] += fd[i*size+k]*del_extend[i+1]*bd[(i+1)*size+k]/pr;
			ins_count[i] += fi[i*size+k]*bi[i*size+k]/pr;
			ins_extend_count[i] += fi[i*size+k]*ins_extend[i]*bi[i*size+k+1]/pr;
		}
	}
	
//	return Py_BuildValue("i", 107);
	
	/* Cleanup */
	Py_DECREF(fm_arr);
	Py_DECREF(fd_arr);
	Py_DECREF(fi_arr);
	Py_DECREF(bm_arr);
	Py_DECREF(bd_arr);
	Py_DECREF(bi_arr);
	Py_DECREF(del_open_arr);
	Py_DECREF(del_extend_arr);
	Py_DECREF(ins_open_arr);
	Py_DECREF(ins_extend_arr);
	Py_DECREF(match_count_arr);
	Py_DECREF(del_count_arr);
	Py_DECREF(ins_count_arr);
	Py_DECREF(del_open_count_arr);
	Py_DECREF(del_extend_count_arr);
	Py_DECREF(ins_open_count_arr);
	Py_DECREF(ins_extend_count_arr);
	Py_INCREF(Py_None);
	
	return Py_None;
}

static PyMethodDef PEhmm2_funcs[] = {
	{ "forward_sum", PEhmm2_forward_sum,METH_VARARGS,"Perform the forward sum." },
	{ "backward_sum", PEhmm2_backward_sum,METH_VARARGS,"Perform the backward sum." },
	{ "forward_sum_interior", PEhmm2_forward_sum_interior,METH_VARARGS,"Perform the forward sum between read ends."},
	{ "backward_sum_interior", PEhmm2_backward_sum_interior,METH_VARARGS,"Perform the backward sum between read ends."},
	{ "add_counts_interior", PEhmm2_add_counts_interior,METH_VARARGS,"Add expected match/insert/delete counts."},
	{NULL, NULL, 0, NULL}
};

void initPEhmm2(void)
{
	(void)Py_InitModule("PEhmm2", PEhmm2_funcs);
	import_array();
}
