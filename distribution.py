import numpy as np
def distribution(f,min,max,boxes=10000.,ddt=1000000.):
	dt = (max-min)/ddt	
	max_area = 1./boxes
	points = [min]
	values = [f(min)]
	area = 0.
	print "dt,min,max:",dt,min,max,f(0)
	for x in np.arange(min,max,dt):
		area += f(x)*dt
		if area >= max_area:
			values.append((f(x)-f(points[-1])/2))
			points.append(x)
			area = 0
	return np.array(points)

def distribution2(f,min,max):
	dx = (max-min)/10000.
	dy = (max-min)/10000.
	max_area = 1/1000.
	lastx = min
	lasty = min
	values = [f(min)]
	area = 0
	for x in range(min,max,dx):
		for y in range(min,max,dy):
			area += f(x,y)*dt
			if area >= max_area:
				values.append((f(x,y)-f(lastx,lasty)/2))
				lastx = x
				lasty = y
				area = 0
	return np.array(values)

	
