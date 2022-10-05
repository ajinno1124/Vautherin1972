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
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

d1_dense=pd.read_csv('../data/Z82N126L0_VB2NaN/density.csv',comment='#')
#print(d1_dense.head)

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

ax.plot(d1_dense["r(fm)"],d1_dense["Rhop"],label=r"$\rho_p$",linewidth=2,color="darkgreen")
ax.plot(d1_dense["r(fm)"],d1_dense["dRhop"],label=r"$d\rho_p$/$dr$",linewidth=2,color="b")
ax.plot(d1_dense["r(fm)"],d1_dense["LapRhop"],label=r"$\Delta \rho_p$",linewidth=2,color="r")
ax.plot(d1_dense["r(fm)"],d1_dense["Jp"],label=r"$J_p$",linewidth=2,color="m")
ax.plot(d1_dense["r(fm)"],d1_dense["DivJp"],label=r"$\nabla \cdot J_p$",linewidth=2,color="k")

ax.text(3,-0.05,r'VB2',{'color':'k','fontsize':14})
#ax.text(3,-0.14,r'HP$\Lambda$2',{'color':'k','fontsize':14})
ax.text(8,0.05,r'$^{208}$Pb',{'color':'k','fontsize':14})

ax.legend(loc='lower right',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(0,10)
#ax.set_ylim(-35,0)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_ylabel('Density',fontsize=16)
ax.set_xlabel(r'$r$ (fm)',fontsize=16)
ax.tick_params(axis='x', which='both', direction='in',labelsize=14)
ax.tick_params(axis='y', which='both', direction='in',labelsize=14)
plt.tight_layout()

#ax.set_xticks([0,1,2,3,4,5])
#ax.set_yticks([-50,0,50,100,200,300])
#ax.tick_params(labelsize=12)

#plt.axis([0.0,1.0, 0.0,0.2])
#plot.ylabel(r'$v_2$',fontsize=20)
#plt.text(1.25,9,'Au+Au b<3.4 fm',{'color':'b','fontsize':20})
#plt.text(0.17,11.5,'$\mu=0.0$',{'color':'b','fontsize':20})
#plt.text(0.17,10,r'hadrons=$(n,\Delta,\pi,\rho$)',{'color':'b','fontsize':20})

#plt.savefig("v2ch.eps",format='eps',dpi=1000)
#plt.savefig("pxe895qmdmsv.eps",format='eps',bbox__inches='tight')
plt.savefig("density.pdf",dpi=300)
plt.savefig("density.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
