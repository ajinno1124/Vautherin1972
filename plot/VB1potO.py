#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath


#plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = False
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

d1_pot=pd.read_csv('../data/Z8N8L0_VB2NaN/potential.csv',comment='#')

#print(d1_pot.head)


fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

ax.plot(d1_pot["r(fm)"],d1_pot["VNp(MeV)"],label=r"$U_p$",linewidth=2,color="r")
ax.plot(d1_pot["r(fm)"],d1_pot["VNn(MeV)"],label=r"$U_n$",linewidth=2,color="b")

#ax.text(1,-4,r'SLy4',{'color':'k','fontsize':14})
#ax.text(1,-7,r'HP$\Lambda$2',{'color':'k','fontsize':14})

ax.legend(loc='lower right',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(0,10)
ax.set_ylim(-90,0)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_ylabel('Potential $U$ (MeV)',fontsize=16)
ax.set_xlabel(r'$r$ (fm)',fontsize=16)
ax.tick_params(axis='x', which='both',direction='in',labelsize=14)
ax.tick_params(axis='y', which='both',left='true',right='true', direction='in',labelsize=14)
#plt.tight_layout()

#ax.set_xticks([0,1,2,3,4,5])
#ax.set_yticks([-50,0,50,100,200,300])
#ax.tick_params(labelsize=12)

plt.savefig("VB2potO.pdf",dpi=300)
plt.savefig("VB2potO.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
