#./a.out l\=4_normal_gauss_0.dat haha.dat result.dat 0.2
from sys import argv
import glob
filename = glob.glob("2d_L30_*.dat")
file_numb = len(filename)
fp = open("commands.in","w")
#fp.write("module load foss/2017A\n")
for i in range(0,file_numb):
	#fp.write("./run "+filename[i]+" ob_"+filename[i]+" \n")
	fp.write("./run "+filename[i]+" trash.dat\n")

