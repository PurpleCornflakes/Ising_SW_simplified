import numpy as np 
import matplotlib.pyplot as plt 
from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_mutual_info_score
from scipy.stats.stats import pearsonr 
# ising_3D = np.loadtxt("./slice_Tc_3D", delimiter = ",")
# adpNN_r = np.load("r1_64.npy")

def MI(x, y, bins = [[-1.5, 0., 1.5],[-1.5, 0., 1.5]]):
	def H(P):
		return np.sum(-P * np.log2(P))
	Nxy = np.histogram2d(x, y, bins = bins)[0]
	Nx = np.sum(Nxy, axis=0)
	Ny = np.sum(Nxy, axis=1)
	Pxy = Nxy/np.sum(Nxy)
	Px = Nx/np.sum(Nx)
	Py = Ny/np.sum(Ny)
	mi = H(Px) + H(Py) - H(Pxy)
	return mi #/np.sqrt(H(Px)*H(Py))

def pearson_correlation(x, y):
	corr,_ = pearsonr(x, y)
	if corr > 0:
		return corr
	else:
		return 0

def metric_distance(data, method, num_expt = 4):
	N,dmax = data.shape

	funcs = {"NMI": normalized_mutual_info_score,
		   "AMI": adjusted_mutual_info_score,
		   "MI": MI,
		   "correlation": pearson_correlation}
	func = funcs[method]
	MIs = np.zeros(dmax)
	for n in range(num_expt):
		for d in range(dmax):
			MIs[d] += func(data[(N//num_expt)*n:(N//num_expt)*(n+1),0], data[(N//num_expt)*n:(N//num_expt)*(n+1),d])
	MIs /= num_expt
	return MIs

def plot_log(data1D, ax, marker, label):
	dmax = data1D.shape[0]
	ax.plot(np.log10(range(1,dmax)), np.log10(data1D[1:]), marker, label = label)
	ax.legend()
	




def plot_Tc():
	method = "correlation"
	markers = ["go","ro","bo"]
	labels = ["T = 1.5","T = 2.2692(Tc)","T = 2.3"]
	fnames = ["./slices_L1000_D2_T1.5000_pBC", "./slices_L1000_D2_T2.2692_pBC","./slices_L1000_D2_T2.3000_pBC"]
	fig,axes = plt.subplots(3,1, sharey = True)
	fig.set_size_inches(10,5)

	for i in range(3):
		data = np.loadtxt(fnames[i], delimiter=None)
		print("init data shape: ",data[i].shape)
		# convert 0,1 into -1,1
		data = np.where(data[:, 450:550] < 0.5, -1. ,1.)#.astype('Float64')
		print("data shape: ",data.shape)

		plot_log(data1D = metric_distance(data=data, method=method,num_expt = 2), 
				ax = axes[i], 
				marker = markers[i], 
				label = labels[i])
		
	
		# plt.ylim(-4,-0.0)
	# plt.savefig("Ising2D_MI.png")
	plt.xlabel("d")
	plt.ylabel(method)
	plt.show()

plot_Tc()
