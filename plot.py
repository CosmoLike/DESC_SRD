import numpy as np
from chainconsumer import ChainConsumer

def sevenchain_multi_plot(filename, out, chainnames, paranames,plotrange):
	fid = np.array([-1.0,0.0])

	d1read = np.genfromtxt(filename[0])
	d2read = np.genfromtxt(filename[1])
	d3read = np.genfromtxt(filename[2])
	d4read = np.genfromtxt(filename[3])
	d5read = np.genfromtxt(filename[4])	
	d6read = np.genfromtxt(filename[5])
	d7read = np.genfromtxt(filename[6])	
	
	d1=d1read[:,(3,4)]-fid
	d2=d2read[:,(3,4)]-fid
	d3=d3read[:,(3,4)]-fid
	d4=d4read[:,(3,4)]-fid
	d5=d5read[:,(3,4)]-fid
	d6=d6read[:,(3,4)]-fid
	d7=d7read[:,(3,4)]-fid
	
	c = ChainConsumer()
	c.add_chain(d1,parameters=paranames, name =chainnames[0])
	c.add_chain(d2,name =chainnames[1])
	c.add_chain(d3,name =chainnames[2])
	c.add_chain(d4,name =chainnames[3])
	c.add_chain(d5,name =chainnames[4])
	c.add_chain(d6,name =chainnames[5])
	c.add_chain(d7,name =chainnames[6])
	#c.configure(kde=[2,2],colors=['c','y'], sigmas=[1,2],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	c.configure(colors=['c','y','r','g','b','k','brown'], sigmas=[1],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	fig = c.plotter.plot(figsize=2.0,extents=plotrange,filename="plots/"+out,truth=[0.0,0.0])


def sixchain_multi_plot(filename, out, chainnames, paranames,plotrange):
	fid = np.array([-1.0,0.0])

	d1read = np.genfromtxt(filename[0])
	d2read = np.genfromtxt(filename[1])
	d3read = np.genfromtxt(filename[2])
	d4read = np.genfromtxt(filename[3])
	d5read = np.genfromtxt(filename[4])	
	d6read = np.genfromtxt(filename[5])	
	
	d1=d1read[:,(3,4)]-fid
	d2=d2read[:,(3,4)]-fid
	d3=d3read[:,(3,4)]-fid
	d4=d4read[:,(3,4)]-fid
	d5=d5read[:,(3,4)]-fid
	d6=d6read[:,(3,4)]-fid
	
	c = ChainConsumer()
	c.add_chain(d1,parameters=paranames, name =chainnames[0])
	c.add_chain(d2,name =chainnames[1])
	c.add_chain(d3,name =chainnames[2])
	c.add_chain(d4,name =chainnames[3])
	c.add_chain(d5,name =chainnames[4])
	c.add_chain(d6,name =chainnames[5])
	#c.configure(kde=[2,2],colors=['c','y'], sigmas=[1,2],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	c.configure(colors=['c','y','r','g','b','k'], sigmas=[1],shade=True,shade_alpha=0.2,shade_gradient=0.0)
	fig = c.plotter.plot(figsize=2.0,extents=plotrange,filename="plots/"+out,truth=[0.0,0.0])
	



# filename=["/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_CL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_SL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_stage3_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_SN_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_3x2_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_ALL_Aug_13"]
# chainnames=[r"Clusters Y10",r"SL Y10",r"Stage III",r"SN Y10",r"3x2pt Y10",r"LSST all+Stage III"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.3,0.5),(-1.5,1.5)]
# sixchain_multi_plot(filename,"LSST_Y10_Aug13_zoomin.pdf",chainnames,paranames,plotrange)


filename=["/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_CL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_SL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_stage3_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_SN_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_3x2_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_ALL_Aug_13"]
chainnames=[r"Clusters Y1",r"SL Y1",r"Stage III",r"SN Y1",r"3x2pt Y1",r"LSST all+Stage III"]
paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
plotrange=[(-0.8,0.8),(-2.5,2.5)]
sixchain_multi_plot(filename,"LSST_Y1_Aug13.pdf",chainnames,paranames,plotrange)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_CL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_SL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_stage3_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_SN_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_3x2_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_ALL_Aug_13"]
# chainnames=[r"Clusters Y10",r"SL Y10",r"Stage III",r"SN Y10",r"3x2pt Y10",r"LSST all+Stage III"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.8,0.8),(-2.5,2.5)]
# sixchain_multi_plot(filename,"LSST_Y10_Aug13.pdf",chainnames,paranames,plotrange)

# filename=["/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_CL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_SL_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_stage3_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_SN_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_3x2_Aug_13","/Users/timeifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_ALL_Aug_13"]
# chainnames=[r"Clusters Y10",r"SL Y10",r"Stage III",r"SN Y10",r"3x2pt Y10",r"LSST all+Stage III"]
# paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
# plotrange=[(-0.3,0.5),(-1.5,1.5)]
# sixchain_multi_plot(filename,"LSST_Y1_Aug13_DESI_1-sigma.pdf",chainnames,paranames,plotrange)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/Tim_April19/lsst_srd_y10_astrosys+stat_noScatter_new_eff.txt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_SN","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_SN+stage3"]
# chainnames=[r"SN Y10 (MCMC)", r"SN Y10 (Fisher)",r"SN Y10+stage3 (Fisher)"]
# paranames=[r"$\Delta w_0$", r"$ \Delta w_a$"]
# plotrange=[(-0.5,0.5),(-1.5,1.5)]
# fivechain_multi_plot(filename,"LSST_Y10_April22_SN_test.pdf",chainnames,paranames,plotrange)


# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_astrosys_shear","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_astrosys_shear_sph005_lph003","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_astrosys_shear_sph005_lph003_photomarg"]
# chainnames=[r"shear shear", r"shear shear contaminated (lph005)",r"shear shear marginalized photo-z"]
# paranames=[r"$\Delta w_0$", r"$ \Delta w_a$"]
# plotrange=[(-1.5,1.5),(-2.5,2.5)]
# fivechain_multi_plot(filename,"LSST_Y1_shear_chain.pdf",chainnames,paranames,plotrange)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_astrosys_clusterN_clusterWL","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_clusters"]
# chainnames=[r"3x2pt", r"3x2pt (Fisher)"]
# paranames=[r"$\Delta w_0$", r"$ \Delta w_a$"]
# plotrange=[(-1.5,1.5),(-2.5,2.5)]
# fivechain_multi_plot(filename,"LSST_Y1_clusters_chainvsFisher.pdf",chainnames,paranames,plotrange)

















# def SNcheck_plot(filename, out, chainnames, paranames):
# 	fid = np.array([-1.0,0.0])
# 	fid2 = np.array([-0.99,0.0])
# 	fid3 = np.array([-1.,0.3])

# 	d1read = np.genfromtxt(filename[0])
# 	d2read = np.genfromtxt(filename[1])
# 	d3read = np.genfromtxt(filename[2])
# 	d4read = np.genfromtxt(filename[3])
	
# 	d1=d1read[:,(3,4)]-fid3
# 	d2=d2read[:,(2,3)]-fid
# 	d3=d3read[:,(3,4)]-fid
# 	d4=d4read[:,(2,3)]-fid2

# 	c = ChainConsumer()
# 	c.add_chain(d1,parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2,name =chainnames[1])
# 	c.add_chain(d3,name =chainnames[2])
# 	c.add_chain(d4,name =chainnames[3])
# 	#c.configure(shade=True, shade_alpha=0.2, bar_shade=True)
# 	c.configure(shade=[False,True,False,True],shade_alpha=[0.2,0.2,0.2,0.2],linestyles=["-", "--", "-", "--"],linewidths=[1.5,1.5,1.5,1.5],colors=["r","r","b","b"])	
# 	fig = c.plotter.plot(figsize=2.5,filename="/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/plots/"+out)


# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_SNY1","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/ReneeSNChains/lsst_y1_jan29_nostarts_varyM_oldomb_rand.txt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_SNY10","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/ReneeSNChains/lsst_y10_jan29_nostarts_varyM_oldomb_rand.txt"]
# chainnames=[r"SN (Fisher) Y1", r"SN (Renee) Y1", r"SN (Fisher) Y10",  r"SN (Renee) Y10"]
# paranames=[r"$w_0$", r"$w_a$"]
# SNcheck_plot(filename,"SN_Y1_Y10_check.pdf",chainnames,paranames)



# def twochain_multi_plot(filename, out, chainnames, paranames):
# 	fid = np.array([-1.0,0.0])
# 	fid = np.array([-1.0,0.0])
# 	d1read = np.genfromtxt(filename[0])
# 	d2read = np.genfromtxt(filename[1])
	
# 	d1=d1read[30000:,(3,4)]-fid
# 	d2=d2read[30000:,(3,4)]-fid
# 	# d3=d3read[:,:7]-fid2
# 	c = ChainConsumer()
# 	c.add_chain(d1,parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2,name =chainnames[1])
# 	c.configure(shade=True, shade_alpha=0.2, bar_shade=True)
# 	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
# 	fig = c.plotter.plot(figsize=2.0,filename="/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/plots/"+out)
	
# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_DES_Y5_2pt_clusterN_clusterWL_Planck_BAO_SN","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_with_sys_2pt_cluster3"]
# chainnames=[r"Planck BOSS JLA", r"DES Y5", r"Y1 3x2pt with sys"]
# paranames=[r"$w_0$", r"$w_a$"]
# twochain_multi_plot(filename,"LSST_3x2pt_Y10_vs_Y1_with_sys.pdf",chainnames,paranames)

 
# def twochain_multi_plot_SN_SL(filename, out, chainnames, paranames):
# 	fid = np.array([0.3156,-1.0,0.0,0.6727])
# 	d1read = np.genfromtxt(filename[0])
# 	d2read = np.genfromtxt(filename[1])
	
# 	d1=d1read-fid
# 	d2=d2read-fid

# 	c = ChainConsumer()
# 	c.add_chain(d1[50000:,(0,1,2,3)], parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2[50000:,(0,1,2,3)], name =chainnames[1])
# 	c.configure(shade=True, shade_alpha=0.2, bar_shade=True)
# 	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
# 	fig = c.plotter.plot(figsize=2.0,filename="/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/plots/"+out)
	
# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_SN_SL_Y10","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_SN_Y10"]
# chainnames=[r"Y10 SL", r"Y10 SN SL"]
# paranames=[r"$\Omega_m$", r"$w_0$", r"$w_a$", r"$h_0$"]
# twochain_multi_plot_SN_SL(filename,"LSST_3x2pt_Y10_vs_Y1_SN_SL.pdf",chainnames,paranames)



# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_with_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_no_sys_2pt"]
# chainnames=[r"Y1 3x2pt with sys", r"Y1 3x2pt no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$", r"$w_a$", r"$\Omega_b$", r"$h_0$"]
# twochain_single_plot(filename,"LSST_Y1_3x2pt_with_sys.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_with_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt"]
# chainnames=[r"Y10 3x2pt with sys", r"Y10 3x2pt no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$", r"$w_a$", r"$\Omega_b$", r"$h_0$"]
# twochain_single_plot(filename,"LSST_Y10_3x2pt_with_sys.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_with_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_with_sys_2pt"]
# chainnames=[r"Y10 3x2pt with sys", r"Y1 3x2pt with sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$", r"$w_a$", r"$\Omega_b$", r"$h_0$"]
# twochain_single_plot(filename,"LSST_Y10_Y1_3x2pt_with_sys.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt_cluster","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_no_sys_2pt_cluster"]
# chainnames=[r"Y10 3x2pt cluster no sys", r"Y1 3x2pt cluster no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$", r"$w_a$", r"$\Omega_b$", r"$h_0$"]
# twochain_single_plot(filename,"LSST_Y10_vs_Y1_3x2pt_cluster_no_sys.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt_cluster"]
# chainnames=[r"Y10 3x2pt no sys", r"Y10 3x2pt clusters no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$", r"$w_a$", r"$\Omega_b$", r"$h_0$"]
# twochain_single_plot(filename,"LSST_Y10_3x2pt_vs_3x2pt_cluster_no_sys.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt_cluster_SN","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_with_sys_2pt_cluster_SN"]
# chainnames=[r"Y10 3x2pt+cluster+SN no sys", r"Y10 3x2pt+cluster+SN with sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$w_0$", r"$w_a$", r"$\Omega_b$", r"$h_0$"]
# twochain_single_plot(filename,"LSST_Y10_3x2pt_cluster_SN_with_sys_vs_no_sys.pdf",chainnames,paranames)

# def twochain_multi_plot(filename, out, chainnames, paranames):
# 	d1 = np.genfromtxt(filename[0])
# 	d2 = np.genfromtxt(filename[1])
# 	c = ChainConsumer()
# 	c.add_chain(d1[50000:,(0,1)], parameters=paranames, name =chainnames[0])
# 	c.add_chain(d2[50000:,(0,1)], name =chainnames[1])
# 	c.configure(shade=True, shade_alpha=0.2, bar_shade=True)
# 	#c.configure(shade=[True,False,False],shade_alpha=[0.2,0.2,0.2],linestyles=["--", "-", "-."],linewidths=[0.5,1.,1.])	
# 	fig = c.plotter.plot(figsize=2.0,filename="/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/plots/"+out)


# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_with_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt"]
# chainnames=[r"Y10 no sys", r"Y1 no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_single_plot(filename,"LSST_Y10_vs_Y1_study.pdf",chainnames,paranames)


# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_with_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_no_sys_2pt"]
# chainnames=[r"Y1 with sys", r"Y1 no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_single_plot(filename,"LSST_Y1_sys_study.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_with_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt"]
# chainnames=[r"Y10 with sys", r"Y10 no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_single_plot(filename,"LSST_Y10_sys_study.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_no_sys_2pt","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt"]
# chainnames=[r"Y1 no sys", r"Y10 no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_single_plot(filename,"LSST_Y10_vs_Y1_study.pdf",chainnames,paranames)

# filename=["/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y1_no_sys_2pt_cluster","/Users/teifler/Dropbox/cosmolike_store/LSSTawakens/like/like_LSST_SRD_Y10_no_sys_2pt_cluster"]
# chainnames=[r"Y1 2pt cluster no sys", r"Y10 2pt cluster no sys"]
# paranames=[r"$\Omega_m$", r"$\sigma_8$"]
# twochain_single_plot(filename,"LSST_Y10_vs_Y1_cluster_study.pdf",chainnames,paranames)
	
