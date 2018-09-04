#!/usr/bin/python
import sys
import math, numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import linalg as LA
import numpy as np

#uncomment below if you want to invert the DESC SRD Y10 covariance
# infile =['cov/Y10_3x2pt_clusterN_clusterWL_cov']
# data= ['datav/3x2pt_clusterN_clusterWL_Y10_fid']
# outname=['Y10']

# nggl = 25 	# number of ggl power spectra
# ngcl = 11 	# number of cluster-source galaxy power spectra
# nlens = 10 	# number of lens bins 
# nlenscl= 4 	# number of cluster redshift bins 
# nshear = 15 # number of shear tomographic power spectra
# ncl=20 		# number of ell-bins
# nclgcl=5	# number of cluster ell-bins
# nrich=5 	# number of richness bins

#uncomment below if you want to invert the DESC SRD Y1 covariance
infile =['cov/Y1_3x2pt_clusterN_clusterWL_cov']
data= ['datav/3x2pt_clusterN_clusterWL_Y1_fid']
outname=['Y1']

nggl = 7 	# number of ggl power spectra
ngcl = 6 	# number of cluster-source galaxy power spectra
nlens = 5 	# number of lens bins 
nlenscl= 3 	# number of cluster redshift bins 
nshear = 15 # number of shear tomographic power spectra
ncl=20 		# number of ell-bins
nclgcl=5	# number of cluster ell-bins
nrich=5 	# number of richness bins


ndata = (nshear+nggl+nlens)*ncl+nlenscl*nrich+nrich*ngcl*nclgcl 
n2pt = (nshear+nggl+nlens)*ncl 
ncluster = nlenscl*nrich 
n2ptcl=n2pt+ncluster
nclusterN_WL=ncluster+nrich*ngcl*nclgcl

for k in range(0,1):
  	datafile= np.genfromtxt(data[k])
  	mask = np.zeros(ndata)
	for i in range(0,datafile.shape[0]):
		if (datafile[i,1] >1.0e-15): 
			mask[i]=1.0

  	
  	covfile = np.genfromtxt(infile[k])
	cov = np.zeros((ndata,ndata))

	print ndata,n2pt,int(np.max(covfile[:,0])+1)

	for i in range(0,covfile.shape[0]):
	  	cov[int(covfile[i,0]),int(covfile[i,1])] = covfile[i,8]+covfile[i,9]
	  	cov[int(covfile[i,1]),int(covfile[i,0])] = covfile[i,8]+covfile[i,9]
	 

	cor = np.zeros((ndata,ndata))
	for i in range(0,ndata):
	    for j in range(0,ndata):
	       if (cov[i,i]*cov[j,j] >0):
	         cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])


	a = np.sort(LA.eigvals(cor[:,:]))
	print "min+max eigenvalues full cor:"
	print np.min(a), np.max(a)
	print "neg eigenvalues full cor:"
	for i in range(0,a.shape[0]):
		if (a[i]< 0.0): print a[i]


	# ############### invert shear covariance #################
	inv = LA.inv(cov[0:nshear*ncl,0:nshear*ncl])
	a = np.sort(LA.eigvals(cov[0:nshear*ncl,0:nshear*ncl]))
	print "min+max eigenvalues shear cov:"
	print np.min(a), np.max(a)
	outfile = "cov/"+outname[k]+"_shear_shear_inv"
	f = open(outfile, "w")
	for i in range(0,nshear*ncl):
		inv[i,i]=inv[i,i]*mask[i]
	  	for j in range(0,nshear*ncl):
	  		f.write("%d %d %e\n" %(i,j, inv[i,j]))
	f.close()

	
	# ############### invert clustering covariance #################
	inv = LA.inv(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl])
	a = np.sort(LA.eigvals(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl]))
	print "min+max eigenvalues clustering cov:"
	print np.min(a), np.max(a)
	outfile = "cov/"+outname[k]+"_pos_pos_inv"
	f = open(outfile, "w")
	for i in range(0,nlens*ncl):
		inv[i,i]=inv[i,i]*mask[(nshear+nggl)*ncl+i]
		for j in range(0,nlens*ncl):
	  		f.write("%d %d %e\n" %(i,j, inv[i,j]))
	f.close()

	# ############### invert 2pt covariance #################
	a = np.sort(LA.eigvals(cov[0:n2pt,0:n2pt]))
	print "min+max eigenvalues 2pt cov:"
	print np.min(a), np.max(a)
	inv = LA.inv(cov[0:n2pt,0:n2pt])
	outfile = "cov/"+outname[k]+"_3x2pt_inv" 
	f = open(outfile, "w")
	for i in range(0,n2pt):
		inv[i,i]=inv[i,i]*mask[i]
	  	for j in range(0,n2pt):
	  		f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()



	# # ############### invert full2pt+clusterN+clusterWL covariance #################
	precond = 1.e-7
	for i in range(0,ncluster):
	  cov[n2pt+i,:]*= precond
	  cov[:,n2pt+i]*= precond
	inv = LA.inv(cov)
	a = np.sort(LA.eigvals(cov))
	print "min+max eigenvalues of full 2ptclusterN+clusterWL pre-conditioned matrix:"
	print np.min(a), np.max(a)
	if (np.min(a)<0):
	  print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])
	for i in range(0,ncluster):
	  inv[n2pt+i,:]*= precond
	  inv[:,n2pt+i]*= precond

	outfile = "cov/"+outname[k]+"_3x2pt_clusterN_clusterWL_inv"
	f = open(outfile, "w")
	for i in range(0,ndata):
	  inv[i,i]=inv[i,i]*mask[i]
	  for j in range(0,ndata):
	    f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()



	# # ############### invert clusterN+clusterWL covariance #################
	inv = LA.inv(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL])
	a = np.sort(LA.eigvals(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL]))
	print "min+max eigenvalues of clusterN_WL pre-conditioned matrix:"
	print np.min(a), np.max(a)
	if (np.min(a)<0):
	  print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile[k])

	for i in range(0,ncluster):
	  inv[i,:]*= precond
	  inv[:,i]*= precond

	outfile = "cov/"+outname[k]+"_clusterN_clusterWL_inv"
	f = open(outfile, "w")
	for i in range(0,nclusterN_WL):
	  	inv[i,i]=inv[i,i]*mask[n2pt+i]
	  	for j in range(0,nclusterN_WL):
	  		f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()

	

	cor = np.zeros((ndata,ndata))
	for i in range(0,ndata):
		for j in range(0,ndata):
			if (cov[i,i]*cov[j,j] >0):
				cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

	plt.figure()
	#plt.imshow(np.log10(cov[0:1500,2000:]), interpolation="nearest",vmin=-25, vmax=-10)
	plt.imshow(np.log10(cov[:,:]), interpolation="nearest",vmin=-25, vmax=-10)
	#plt.imshow(cor[n2ptcl:n2ptcl+200,300:nshear*ncl], interpolation="nearest",vmax=0.5)
	plt.colorbar()
	plt.show()
