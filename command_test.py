from sys import argv
import glob
filename = glob.glob("l=5_*.dat")
file_numb = len(filename)
fp = open("commands.in","w")
#fp.write("module load foss/2017A\n")
for i in range(0,file_numb):
	fp.write("./run "+filename[i]+" observable_%d.dat"%(i)+" 1.0"+" 1"+"\n")
#	fp.write("./run "+filename[i]+" observable.dat\n")

