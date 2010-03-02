import numpy as np
import distribution
import constants
import itertools as it

class Silicon:
	eV = 1.60217646*10**-19
	bandgap = 1.12 * eV
	free_mass = ( 9.10938188 *10**-31)

	electron_mass = 1.08*free_mass
	hole_mass = .58*free_mass
	dielectric = 11.7
	epsilon = dielectric*8.85418782*10**-12
	doping = 10.**24
	doping3d = 10.**24
	doping2d = 10.**22

	def random_momentum(self,type):
		if type == "ptype":
			m = hole_mass
		else:
			m = electron_mass
		energy = distribution(lambda x:(kbT)*math.exp(-x/(kbT)),
					0,kbT*10)
		velocity = map(lambda x: math.sqrt(x)/(2*m),)
		angles = np.arrange(0,2*np.pi,.001)
		retarray = np.ndarray(len(velocity)*len(angles))
		index = it.count()
		for v in velocity:
			for theta in angles:
				i = index.next()
				retarray[i][0] = v*cos[theta]
				retarray[i][1] = v*sin[theta]
		return retarray
	
