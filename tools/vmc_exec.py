#!/bin/env python
from __future__ import unicode_literals,print_function
import six
import subprocess
import sys
import re
import h5py
import numpy as np
import vmc

def vmc_exec(**kwargs):
    """
    arguments:
    - nprocs (default 4)`
    - hosts (default localhost)
    - Nsamp (default 1)
    plus all those given by runing vmc -h.
    returns a dictionnary with all measurements as
    asked using the kwargs.
    """
    opts=get_vmc_args()
    opts+=['nprocs','hosts','Nsamp']
    kwargs.setdefault('nprocs',4)
    kwargs.setdefault('hosts','localhost')
    kwargs.setdefault('dir','.')
    kwargs.setdefault('Nsamp',1)
    meas_trans={'meas_magnetization':'Magnetization',\
                'meas_projheis':'ProjHeis',\
                'meas_stagmagn':'StagMagn',\
                'meas_statspinstruct':'StatSpinStruct'}
    if kwargs.setdefault('Sztot_conserved',False):
        meas_trans['meas_stagmagn']='StagMagnZ'
    vmcargs=[]
    for key,value in kwargs.items():
        if key!='nprocs' and key!='hosts' and key!='Nsamp':
            if not key in opts:
                raise RuntimeError('Unrecognized option \'{}\'.'.format(key))
            if type(value)==bool:
                if value:
                    vmcargs+=['--{}'.format(key)]
            else:
                vmcargs+=['--{key}={val}'.format(key=key,val=value)]
    vmcproc=subprocess.Popen(['mpiexec','-np',str(kwargs['nprocs']),'--host',kwargs['hosts'],'vmc']+vmcargs,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=vmcproc.communicate()
    if six.PY3:
        stdout=str(stdout,encoding='utf-8')
        stderr=str(stderr,encoding='utf-8')
    else:
        stdout=unicode(stdout,encoding='utf-8')
        stderr=unicode(stderr,encoding='utf-8')
    print(stderr)
    prefix=re.findall(r'output prefix=([0-9]+)',stdout)

    outq=dict()
    for meas in meas_trans.keys():
        if kwargs.setdefault(meas,False):
            outq[meas_trans[meas]]=get_quantity(kwargs['dir']+'/{prefix}-{measname}.h5'.format(prefix=prefix[0],measname=meas_trans[meas]),kwargs['Nsamp'])
    return outq

def get_quantity(filepath,Nsamp):
    print(filepath)
    dpath,args=vmc.GetStat([filepath],Nsamp)
    hfile=h5py.File(dpath[args[0][0]][0],'r')
    dshape=list(hfile[dpath[args[0][0]][1]].shape)
    hfile.close()
    dshape[1]/=2 #because complex numbers a contiguous along a line
    Q=np.zeros([Nsamp]+list(dshape),dtype=complex)
    for si,bi in enumerate(args):
        for d in bi:
            hfile=h5py.File(dpath[d][0],'r')
            dat=hfile[dpath[d][1]]
            Q[si,:]+=dat[:,::2]+1j*dat[:,1::2]
            hfile.close()
        Q[si,:]/=len(bi)
    return Q

def get_vmc_args():
    """
    get all arguments from the vmc command
    """
    vmcquery=subprocess.Popen(['vmc'],stdout=subprocess.PIPE)
    vmcargs,_=vmcquery.communicate()
    if six.PY3:
        vmcargs=str(vmcargs,encoding='utf-8')
    else:
        vmcargs=unicode(vmcargs,encoding='utf-8')
    opts=[]
    for l in vmcargs.split('\n'):
        opts+=re.findall(r'--([a-zA-Z0-9_]*)=?.*',l)
    return opts
