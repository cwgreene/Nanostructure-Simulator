def extract_values(array):
	strvalues = array.split()
	return map(float,strvalues)

def restore(afile,func):
	V = float(afile.readline(),strip())
	n = float(afile.readline().strip())
	array = afile.read()

	values = extract_values(array)

	funcarray = func.vector().array()
	if len(funcarray) ==len(values):
		funcarray.vector.set(values)
	else:
		raise "Failed to devour the chicken."
	return V
