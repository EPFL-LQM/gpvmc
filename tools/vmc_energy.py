#!/bin/env python

import subprocess
import sys
import re
import h5py
import numpy as np

def vmc_energy(**kwargs):
    """
    arguments:
    - nprocs
    - hosts
    plus all those given by runing vmc -h.
    """
    # get all arguments from the vmc command
    vmcquery=subprocess.Popen(['vmc'],stdout=subprocess.PIPE)
    vmcargs,_=vmcquery.communicate()
    opts=[]
    for l in str(vmcargs,encoding='utf-8').split('\n'):
        opts+=re.findall(r'--([a-zA-Z0-9_]*)=?.*',l)
    opts+=['nprocs','hosts']
    kwargs.setdefault('nprocs',4)
    kwargs.setdefault('hosts','localhost')
    kwargs.setdefault('dir','.')
    vmcargs=[]
    for key,value in kwargs.items():
        if key!='nprocs' and key!='hosts':
            if not key in opts:
                raise RuntimeError('Unrecognized option \'{}\'.'.format(key))
            if type(value)==bool:
                if value:
                    vmcargs+=['--{}'.format(key)]
            else:
                vmcargs+=['--{key}={val}'.format(key=key,val=value)]
    vmcproc=subprocess.Popen(['mpiexec','-np',str(kwargs['nprocs']),'--host',kwargs['hosts'],'vmc']+vmcargs,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=vmcproc.communicate()
    print(str(stderr,encoding='utf-8'))
    stdout=str(stdout,encoding='utf-8')
    prefix=re.findall(r'output prefix=([0-9]+)',stdout)
    enfile=kwargs['dir']+'/{}-ProjHeis.h5'.format(prefix[0])
    enfile=h5py.File(enfile,'r')
    energy=[]
    for r in enfile.keys():
        for d in enfile[r].keys():
            energy+=[enfile[r][d][0,0]]
    energy=np.array(energy)
    em=np.mean(energy)
    es=np.std(energy)/np.sqrt(len(energy))
    return em,es,energy
