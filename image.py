import numpy as np
import scipy.misc

with open("./render/index", "r") as f:
	reso = int(f.readline())

image = np.empty([reso, reso, 3], dtype = np.float32)

for i in range(reso):
	with open("./render/{}".format(i)) as f:
		for j in range(reso):
			line = f.readline()
			image[i, j] = [float(x) for x in line.split()]

scipy.misc.imsave("render.png", image)