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
	



filename=["/like/like_LSST_SRD_Y1_CL_Aug_13","/like/like_LSST_SRD_Y1_SL_Aug_13","/like/like_LSST_SRD_stage3_Aug_13","/like/like_LSST_SRD_Y1_SN_Aug_13","/like/like_LSST_SRD_Y1_3x2_Aug_13","/like/like_LSST_SRD_Y1_ALL_Aug_13"]
chainnames=[r"Clusters Y1",r"SL Y1",r"Stage III",r"SN Y1",r"3x2pt Y1",r"LSST all+Stage III"]
paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
plotrange=[(-0.8,0.8),(-2.5,2.5)]
sixchain_multi_plot(filename,"LSST_Y1.pdf",chainnames,paranames,plotrange)

filename=["/like/like_LSST_SRD_Y10_CL_Aug_13","/like/like_LSST_SRD_Y10_SL_Aug_13","/like/like_LSST_SRD_stage3_Aug_13","/like/like_LSST_SRD_Y10_SN_Aug_13","/like/like_LSST_SRD_Y10_3x2_Aug_13","/like/like_LSST_SRD_Y10_ALL_Aug_13"]
chainnames=[r"Clusters Y10",r"SL Y10",r"Stage III",r"SN Y10",r"3x2pt Y10",r"LSST all+Stage III"]
paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
plotrange=[(-0.8,0.8),(-2.5,2.5)]
sixchain_multi_plot(filename,"LSST_Y10.pdf",chainnames,paranames,plotrange)

filename=["/like/like_LSST_SRD_Y10_CL_Aug_13","/like/like_LSST_SRD_Y10_SL_Aug_13","/like/like_LSST_SRD_stage3_Aug_13","/like/like_LSST_SRD_Y10_SN_Aug_13","/like/like_LSST_SRD_Y10_3x2_Aug_13","/like/like_LSST_SRD_Y10_ALL_Aug_13"]
chainnames=[r"Clusters Y10",r"SL Y10",r"Stage III",r"SN Y10",r"3x2pt Y10",r"LSST all+Stage III"]
paranames=[r"$\Delta w_0$", r"$\Delta w_a$"]
plotrange=[(-0.5,0.5),(-1.5,1.5)]
sixchain_multi_plot(filename,"LSST_Y10_zoom.pdf",chainnames,paranames,plotrange)

