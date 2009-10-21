import os

def read_current_output(afile):
	V = float(afile.readline().strip())
	n = float(afile.readline().strip())
	cur_vals = afile.read().split("\n")[:-1]
	return V,n,cur_vals

def mean(x):
	sum = 0.0
	for a in x:
		sum += float(a)
	return sum/(len(x))

def filter_dir(dir,afilter):
	return filter(lambda x: (x.find(afilter) )>0,os.listdir(dir))

def current_values():
	vals = {}
	for name in filter_dir("results","current"):
		afile = open("results/"+name)
		V,n,cur_vals = read_current_output(afile)
		vals[V]=mean(cur_vals)
	voltages = sorted(vals.keys())
	vlist,ilist = [],[]
	for v in voltages:
		vlist.append(str(v))
		ilist.append(str(vals[v]))
		print "V",v,"I",vals[v]
	print "c("+",".join(vlist)+")"
	print "c("+",".join(ilist)+")"
