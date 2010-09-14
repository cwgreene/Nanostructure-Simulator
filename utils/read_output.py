import os,sys
import optparse

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

def filter_tags(tags,options):
	if options.tag == "None":
		return tags.keys() #return all
	return filter(lambda x: x.find("tag:"+options.tag)>0,tags.keys())

def current_values(options):
	vals = {} #values
	tags = {} #line tags, needs to be made complicated
	#go through filter directory and take all current files
	for name in filter_dir("results","current"):
		try: #assume well formatted
			afile = open("results/"+name)
			V,misc,cur_vals = read_current_output(afile)
			tags[misc]=0
			if V in vals.keys():
				vals[(V,misc)]+=[mean(cur_vals)]
			else:
				vals[(V,misc)] = [mean(cur_vals)]
		except: #not well formmatted
			if options.verbose == True:
				print "bad input:",name
	all_voltages = sorted(vals.keys())
	selected_tags = filter_tags(tags,options)
	for tag in selected_tags:
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
	parser = optparse.OptionParser()
	parser.add_option("-c",type="int",dest="count",default=-1)
	parser.add_option("--tag",type="string",dest="tag",default="None")
	parser.add_option("-v",dest="verbose",action="store_true")
	options,args = parser.parse_args(sys.argv)
	current_values(options)
