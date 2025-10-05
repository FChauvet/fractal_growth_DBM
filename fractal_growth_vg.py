# -----------------------------------------------------------------------------
# Project: Fractal Growth DBM
# File: fractal_growth_vg.py
# Author: Fabien Chauvet
# License: MIT (see LICENSE file for details)
# Copyright (c) 2025 Fabien Chauvet
# -----------------------------------------------------------------------------



import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import h5py
import time
from solve_adv_diff import*


#-----------Constants---------------

Ld = 20                       #diffusion lenght
m = 200                       #width of the computation domain
hmax = 200
domain_height = 10            #domain height*Ld = size of the moving computation domain above the front
filename = "test.hdf5"

ec = 1e-04/Ld                 #convergence criterion 
rc = ec*10                    #lower limit for mass flux considered





#-----------Functions---------------

def showim(A):

	vminn=A.min(0).min(0)
	vmaxx=A.max(0).max(0)

	plt.imshow(A,vmin=vminn,vmax=vmaxx,interpolation='none',cmap='jet')
	plt.ion()
	plt.show()
	plt.pause(0.001)
	return 


def pt_contour(A):

	[n,m]=A.shape
	k=0
	index=np.array([(0)],dtype='i')

	# computation of the front location (from the matrix top)
	K=A.mean(1).nonzero()
	y_f=K[0].min()

	for i in range(y_f-1,n-1):
		for j in range(1,m-1):
			if A[i-1,j]+A[i+1,j]+A[i,j-1]+A[i,j+1]>0 and A[i,j]<1: 
				k=k+1
				index=np.concatenate((index,np.array([(np.ravel_multi_index([i,j],A.shape,order='F'))])))

	j=0
	for i in range(y_f-1,n-1):
		#if A[i+1,j]+A[i-1,j]+A[i,j+1]>0 and A[i+1,j]+A[i-1,j]+A[i,j+1]<3 and A[i,j]<1:
		if A[i+1,j]+A[i-1,j]+A[i,j+1]+A[i,m-1]>0 and A[i,j]<1:
			k=k+1
			index=np.concatenate((index,np.array([(np.ravel_multi_index([i,j],A.shape,order='F'))])))
	j=m-1
	for i in range(y_f-1,n-1):
		#if A[i+1,j]+A[i-1,j]+A[i,j-1]>0 and A[i+1,j]+A[i-1,j]+A[i,j-1]<3 and A[i,j]<1:
		if A[i+1,j]+A[i-1,j]+A[i,j-1]+A[i,0]>0 and A[i,j]<1:
			k=k+1
			index=np.concatenate((index,np.array([(np.ravel_multi_index([i,j],A.shape,order='F'))])))

	i=n-1
	for j in range(1,m-1):
		if A[i-1,j]+A[i,j-1]+A[i,j+1]>0 and A[i,j]<1:
			k=k+1
			index=np.concatenate((index,np.array([(np.ravel_multi_index([i,j],A.shape,order='F'))])))


	i=n-1
	j=0
	#if A[i-1,j]+A[i,j+1]>0 and A[i,j]<1:
	if A[i-1,j]+A[i,j+1]+A[i,m-1]>0 and A[i,j]<1:
		k=k+1
		index=np.concatenate((index,np.array([(np.ravel_multi_index([i,j],A.shape,order='F'))])))

	i=n-1
	j=m-1
	#if A[i-1,j]+A[i,j-1]>0 and A[i,j]<1:
	if A[i-1,j]+A[i,j-1]+A[i,0]>0 and A[i,j]<1:
		k=k+1
		index=np.concatenate((index,np.array([(np.ravel_multi_index([i,j],A.shape,order='F'))])))

	indexx=index[1:]

	return indexx

def lapl_growth_vg_global_rc(hmax,m,Ld,domain_height,filename,rc,erreur_conv):

	#record results every kkstep
	kkstep=int(0.2*m*Ld/100.)
	#starting
	kkstepp=0

	f=h5py.File(filename,"w")

	n=hmax+domain_height*Ld

	#initialize the lattice (0: empty site, 1: occuped site)
	A=np.zeros((n,m),order='F')
	pl=np.zeros((n,m),order='F')
	A[-1,:]=1
	#Af=A

	c=np.zeros((n,m),order='F')
	y_f=n-1 #initialize deposit front location (distance from the top of the matrix)

	#initialize the concentration matrix with the 1D "unperturbed" solution
	for i in range(0,y_f+1):
		c[i,:]=1.0-np.exp(-(np.float64(y_f-i))/np.float64(Ld))

	#initialize deposit maximal heigth
	h=0
	#initialize the iteration counter
	kk=0

	#initialize point list of the boundary
	indexx=pt_contour(A)


	#main while loop

	print ('start')
	t0=time.time()

	while h<hmax:	
	
		kk=kk+1
		AA=1-A

		# initialize the iteration counter
		k=0
		# initialize the residu
		residu=1.0

		#locations of the occupied sites
		II=(A>0).nonzero()

		# computation of the front location (from the matrix top)
		K=A.mean(1).nonzero()
		y_f=K[0].min()

		# limits of the interior of the domain
		i0=y_f-domain_height*Ld+1
		iff=n-2
		j0=1
		jff=m-2

		# 1D & unperturbed solution above the top boundary condition optionnal
		for i in range(0,y_f-domain_height*Ld-1+1):
			c[i,:]=1.0-np.exp(-(np.float64(y_f-i))/np.float64(Ld)) #1.0-np.exp(-(n-(i+1.0))/Ld)

		# boundary condition at the top of the computation domain (1D solution
		# far from the front
		c[y_f-domain_height*Ld,:]=1.0 - np.exp(-np.float64(domain_height))

		# boundary condition at the front, c=0
		c[II]=0

		#point list of the boundary
		pl=pl*0
		pl[np.unravel_index(indexx,(n,m),order='F')]=1.0

		#computation of the concentration field
		c=solve_adv_diff(Ld,AA,c,pl,n,m,i0,iff,j0,jff,erreur_conv)

		# distribution of flux density of the deposit boundary
		indexx=pt_contour(A)
		EE=0
		liste_sum=np.array([(0)],dtype='float64')
		liste_E0=np.array([(0)],dtype='float64')
		indexxx=np.array([(0)],dtype='i')
		for i in range(0,indexx.size):
			E0=c[np.unravel_index(indexx[i],(n,m), order='F')]
			if E0>rc:
				EE=EE+E0
				liste_sum=np.concatenate((liste_sum,np.array([(EE)],dtype='float64')))
				liste_E0=np.concatenate((liste_E0,np.array([(E0)],dtype='float64')))
				indexxx=np.concatenate((indexxx,np.array([(indexx[i])])))

		liste_sum=liste_sum[1:]
		max_sum=liste_sum.max().copy()
		liste_sum=liste_sum/max_sum
		liste_E0=liste_E0[1:]
		liste_E0_n=liste_E0/max_sum
		indexxx=indexxx[1:]

		if liste_E0_n.min()<2**(-53):
			tmp=input("too small values in liste_E0_n")

		exe_time = time.time()-t0

		# recording of the data
		if kk>=kkstepp:
			print('save')
			kkstepp=kkstepp+kkstep

			#showim(A+c)
			#showim(pl)

			grp=str(kk)
			grp = f.create_group(str(kk))
			grp.create_dataset("data_Ld", data=Ld)
			grp.create_dataset("data_m", data=m)
			grp.create_dataset("data_hmax", data=hmax)
			grp.create_dataset("data_domain_height", data=domain_height)
			grp.create_dataset("data_rc", data=rc)
			grp.create_dataset("data_erreur_conv", data=erreur_conv)
			grp.create_dataset("data_kkstep", data=kkstep)
			grp.create_dataset("data_A", data=A)
			grp.create_dataset("data_c", data=c)
			grp.create_dataset("data_exe_time", data=exe_time)
			grp.create_dataset("data_kk", data=kk)


		#random selection of a point of the deposit boundary following the
		#distribution
		ix=np.random.uniform(0, 1)
		V=liste_sum-ix		
		#ii=(V>0).nonzero()
		#iii=ii[0]
		iiii= np.amin(np.nonzero(V>0)) #iii[0]
		[I,J]=np.unravel_index(indexxx[iiii], (n,m), order='F')
		    
		# filling of the selected lattice site
		A[I,J]=1
		    
		K=A.mean(1).nonzero()
		y_f=K[0].min()
		h=n-1-y_f

		print ('Ld=', Ld, 'm=', m, 'kk=', kk, 'h=', h, '/', hmax)


		  
	#f.close()
	return

#-----------------End functions--------------------


lapl_growth_vg_global_rc(hmax,m,Ld,domain_height,filename,rc,ec)
