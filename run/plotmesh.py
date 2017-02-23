#!/task/imd/lib64/bin/python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri

xdat=[]
ydat=[]
for i in range(7):
	ip1=i+1
	filename='geo%1d.dat'%(ip1)
	f=open(filename,'r')
	line=f.readline()
	linearr=line.split()
	nw=int(linearr[3][0:len(linearr[3])-1])
	nh=int(linearr[5][0:len(linearr[5])-1])
	
	xdat.append(np.zeros((nw,nh)))
	ydat.append(np.zeros((nw,nh)))
	for h in range(nh):
		for w in range(nw):
			line=f.readline()
			linearr=line.split()
			xdat[i][w,h]=float(linearr[0])
			ydat[i][w,h]=float(linearr[1])
	f.close()
	
	if i==0:
		color='k'
	elif i==1:
		color='r'
	elif i==2:
		color='y'
	elif i==3:
		color='b'
	elif i==4:
		color='g'
	elif i==5:
		color='c'
	elif i==6:
		color='m'

	for w in range(nw):
		plt.plot(xdat[i][w,:],ydat[i][w,:],color)
	for h in range(nh):
		plt.plot(xdat[i][:,h],ydat[i][:,h],color)

f=open('input/DDN_DIIID/limiter.in','r')
line=f.readline()
linearr=line.split()
nlim=int(linearr[2])
f.readline()
lim=np.zeros((nlim,2))
for i in range(nlim):
	line=f.readline()
	linearr=line.split()
	lim[i,0]=float(linearr[0])
	lim[i,1]=float(linearr[1])
f.close()
plt.plot(lim[:,0],lim[:,1],'k')
plt.axis('scaled')
plt.show()
