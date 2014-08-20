#!/bin/env python
from __future__ import unicode_literals,print_function
import six
import subprocess
import sys
import re
import h5py
import time
import numpy as np
from vmc_postproc import load

def vmc_exec(**kwargs):
    """
    arguments:
    - nprocs (default 4)`
    - hosts (default localhost)
    - slurmqueue (default False)
    - partition (default qwfall)
    - rundir (default .)
    plus all those given by runing vmc -h.
    returns a dictionnary with all measurements as
    asked using the kwargs.
    """
    opts=get_vmc_args()
    opts+=['nprocs','hosts']
    kwargs.setdefault('nprocs',4)
    kwargs.setdefault('hosts','localhost')
    kwargs.setdefault('dir','.')
    kwargs.setdefault('partition','qwfall')
    kwargs.setdefault('rundir','.')
    kwargs.setdefault('slurmqueue',False)
    add_keys=['nprocs','hosts','slurmqueue','partition','rundir']
    meas_trans={'meas_magnetization':'Magnetization',\
                'meas_projheis':'ProjHeis',\
                'meas_stagmagn':'StagMagn',\
                'meas_statspinstruct':'StatSpinStruct'}
    if kwargs.setdefault('stagflux_wav',False):
        meas_trans['meas_stagmagn']='StagMagnZ'
    vmcargs=[]
    for key,value in kwargs.items():
        if not key in add_keys:
            if not key in opts:
                raise RuntimeError('Unrecognized option \'{}\'.'.format(key))
            if type(value)==bool:
                if value:
                    vmcargs+=['--{}'.format(key)]
            else:
                vmcargs+=['--{key}={val}'.format(key=key,val=value)]
    if kwargs['slurmqueue']:
        vmcargsinline=''
        for a in vmcargs:
            vmcargsinline+=' '+a
        batchscript="""#!/bin/bash
#SBATCH -n {nprocs}
#SBATCH -c 1
#SBATCH -J vmc_exec
#SBATCH -p {partition}
#SBATCH --error=vmc_exec-%j.stderr
#SBATCH --output=vmc_exec-%j.stdout

cd {rundir}

mpirun ./vmc --prefix=$SLURM_JOB_ID""".format(nprocs=kwargs['nprocs'],partition=kwargs['partition'],rundir=kwargs['rundir'])+vmcargsinline
        batch_file=open('vmc_exec_batch.sh','w')
        batch_file.write(batchscript)
        batch_file.close()
        slurmproc=subprocess.Popen(['sbatch','vmc_exec_batch.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=slurmproc.communicate()
        if six.PY3:
            stdout=str(stdout,encoding='utf-8')
            stderr=str(stderr,encoding='utf-8')
        else:
            stdout=unicode(stdout,encoding='utf-8')
            stderr=unicode(stderr,encoding='utf-8')
        print(stderr)
        print(stdout)
        prefix=re.findall(r'Submitted batch job ([0-9]+)',stdout)
        done=False
        while not done:
            squeueproc=subprocess.Popen(['squeue','-h','-j',str(prefix)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            stdout,stderr=squeueproc.communicate()
            done=(len(stdout)==0)
            time.sleep(10)
        print('Finished calculation '+str(prefix))
    else:
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
            outq[meas_trans[meas]]=load.get_quantity(kwargs['dir']+'/{prefix}-{measname}.h5'.format(prefix=prefix[0],measname=meas_trans[meas]),kwargs['samples'])
    return outq

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
