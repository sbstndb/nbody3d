import numpy as np
import pylab as pl


def save_position(name):
	faos = open("results/" + name, "r")
	daos = faos.readlines()
	xaos = np.zeros((len(daos), 3), dtype="float64")
	for i in range(len(daos)):
		l = daos[i].rsplit(" ")
		for j in range(3):
			l[j] = float(l[j])
			xaos[i, j] = l[j]
			
	return xaos
	
	
def difference_by_particule(x1, x2):
	diff = np.zeros((len(x1), 1), dtype="float64")
	for i in range(len(x1)):
		d  = abs(x1[i, 0] - x2[i, 0]) + abs(x1[i, 1] - x2[i, 1]) + abs(x1[i, 2] - x2[i, 2])
		diff[i] = d 
	return diff
	
	
def plot_error(dtotal, ntotal):
	for i in range(len(ntotal)):
		pl.plot(dtotal[:,i], linewidth = 0.1, label=ntotal[i])
	pl.legend()
	pl.grid()
	pl.yscale("log")
	pl.xlabel("particle ID")
	pl.ylabel("L1 error VS basis code")	
	#pl.show()
	print("\nPlease open results/error.png\n")
	pl.savefig("results/error.png")
		
	
		
	
	
xaos1 = save_position("aos1x.txt")
xaos2 = save_position("aos2x.txt")
xsoa1 = save_position("soa1x.txt")
xsoa2 = save_position("soa2x.txt")
xsoa3 = save_position("soa3x.txt")
xsoa4 = save_position("soa4x.txt")
#xsoa5 = save("soa5x.txt")		


daos2 = difference_by_particule(xaos1, xaos2)	
daos2 = difference_by_particule(xaos1, xaos2)	
dsoa1 = difference_by_particule(xaos1, xsoa1)	
dsoa2 = difference_by_particule(xaos1, xsoa2)
dsoa3 = difference_by_particule(xaos1, xsoa3)
dsoa4 = difference_by_particule(xaos1, xsoa4)

dtotal = np.zeros((len(daos2),6), dtype="float64")
dtotal[:,0] = daos2[:,0]
dtotal[:,1] = dsoa1[:,0]
dtotal[:,2] = dsoa2[:,0]
dtotal[:,3] = dsoa3[:,0]
dtotal[:,4] = dsoa4[:,0]
ntotal = ["AOS2", "SOA1", "SOA2", "SOA3", "SOA4"]


for i in range(len(ntotal)):
	print("----------------------")
	print(str(ntotal[i])+" VS AOS1 : \n")
	print("Erreur moyenne : " + str(np.mean(dtotal[:,i])))
	print("Erreur minimale : " + str(np.min(dtotal[:,i])))
	print("Erreur maximale : " + str(np.max(dtotal[:,i])))




plot_error(dtotal, ntotal)
	
