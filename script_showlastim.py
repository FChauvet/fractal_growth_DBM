import numpy as np
import matplotlib.pyplot as plt
import h5py


def extract_select(filename,i):

	#kk=100

	f = h5py.File(filename, 'r')
	
	N=len(f.keys())
	last_k=i
	kkstepa=f["/"+str(last_k)+"/data_kkstep/"]
	kkstep=kkstepa[()]
	Nn=N*kkstep-kkstep
	Lda= f["/"+str(last_k)+"/data_Ld/"]
	ma= f["/"+str(last_k)+"/data_m/"]
	hmaxa=f["/"+str(last_k)+"/data_hmax/"]
	domain_heighta=f["/"+str(last_k)+"/data_domain_height/"]
	rca=f["/"+str(last_k)+"/data_rc/"]
	erreur_conva=f["/"+str(last_k)+"/data_erreur_conv/"]
	Aa=f["/"+str(last_k)+"/data_A/"]
	ca=f["/"+str(last_k)+"/data_c/"]
	exe_timea=f["/"+str(last_k)+"/data_exe_time/"]
	kka=f["/"+str(last_k)+"/data_kk/"]

	
	Ld=Lda[()]
	m=ma[()]
	hmax=hmaxa[()]
	domain_height=domain_heighta[()]
	rc=rca[()]
	erreur_conv=erreur_conva[()]
	A=Aa[()]
	c=ca[()]
	exe_time = exe_timea[()]
	kk=kka[()]

	f.close()
	return A,c,Ld,domain_height,Nn, hmax



filename="test.hdf5"

[A,c,Ld,domain_height,Nn,hmax]=extract_select(filename,1)
#Nn=400

[A,c,Ld,domain_height,Nn,hmax]=extract_select(filename,Nn)
[n,m]=A.shape
K=A.mean(1).nonzero()
y_f=K[0].min()
h=n-1-y_f
print (Nn, 'h=', h, '/', hmax)



plt.figure()
plt.imshow((A+c),vmin=0,vmax=1,interpolation='none',cmap='jet')
#plt.imshow(np.concatenate((A+c,A+c),axis=1),vmin=0,vmax=1,interpolation='none',cmap='jet')
plt.ion()
plt.colorbar
plt.show()
plt.pause(0.001)

input()



