import numpy as np
from histogram import histo
import math as math

def flat(min,max):
    return np.random.random_sample(1)[0]*(max-min)+min

def norm(mu=0,sigma=1):
    return np.random.normal(mu,sigma,1)[0]

def secret():
    x=flat(0,1)
    y=flat(0,1)
    a=0.1
    b=0.9
    if(y<a*x+b): return x
    else: return secret();

def linbg(xmin,xmax,a,b):
    x=flat(xmin,xmax)
    y=flat(0,max(xmin*a+b,xmax*a+b))
    if(y<a*x+b): return x
    else: return linbg(xmin,xmax,a,b);

mass={}
mass["e-"]=5.109989e-4
mass["e+"]=5.109989e-4
mass["p+"]=0.938272046
mass["pbar-"]=0.938272046
mass["pi+"]=0.139570
mass["pi-"]=0.139570
mass["K+"]=0.493677
mass["K-"]=0.493677
mass["mu-"]=0.1056584
mass["mu+"]=0.1056584
