import os

def read_current_output(afile):
	V = float(afile.readline().strip())
	misc = afile.readline().strip()
	cur_vals = afile.read().split("\n")[:-1]
	return V,misc,cur_vals

def mean(x):
	sum = 0.0
	for a in x:
		sum += float(a)
	return sum/(len(x))

def filter_dir(dir,afilter):
	return filter(lambda x: (x.find(afilter) )>0,os.listdir(dir))

def current_values():
	vals = {}
	tags = {}
	for name in filter_dir("results","current"):
		try:
			afile = open("results/"+name)
			V,misc,cur_vals = read_current_output(afile)
			tags[misc]=0
			if V in vals.keys():
				vals[(V,misc)]+=[mean(cur_vals)]
			else:
				vals[(V,misc)] = [mean(cur_vals)]
		except:
			print "bad input:",name
	all_voltages = sorted(vals.keys())
	for tag in tags:
		voltages = filter(lambda x:x[1]==tag,all_voltages)
		vlist,ilist = [],[]
		print "\n"+tag
		for v in voltages:
			vlist+= map(str,[v[0]]*len(vals[v]))
			ilist+=map(str,vals[v])
			print "V",v[0],"I",vals[v]
		print "c("+",".join(vlist)+")"
		print "c("+",".join(ilist)+")"

if __name__=="__main__":
	current_values()
