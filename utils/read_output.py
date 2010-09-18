import os,sys
import optparse

def read_current_output(afile):
	V = float(afile.readline().strip())
	misc = afile.readline().strip()
	#cur_vals = afile.read().split("\n")[:-1]
	lines = afile.readlines()
	cur_vals = lines[:-1]
	photocurrent=float(cur_vals[-1])
	return V,misc,cur_vals,photocurrent

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
			V,misc,cur_vals,pc = read_current_output(afile)
			tags[misc]=0
			if V in vals.keys():
				vals[(V,misc)]+=[(mean(cur_vals),pc)]
			else:
				vals[(V,misc)] = [(mean(cur_vals),pc)]
		except Exception,e: #not well formmatted
			print name+":",e
			if options.verbose == True:
				print "bad input:",name
	all_voltages = sorted(vals.keys())
	selected_tags = filter_tags(tags,options)
	for tag in selected_tags:
		offset_current = None
		voltages = filter(lambda x:x[1]==tag,all_voltages)
		vlist,ilist,pcurs,dcurs = [],[],[],[]
		print "\n"+tag
		for v in voltages:
			if v[0] == 0.0:
				offset_current=mean(map(lambda x: x[0],vals[v]))
				print offset_current
		for v in voltages:
			vlist+= map(str,[v[0]]*len(vals[v]))
			print vals[v]
			curs = map(lambda x: x[1]-x[0]+offset_current,vals[v])
			pcurs += map(lambda x: str(x[1]),vals[v])
			dcurs += map(lambda x: str(x[0]),vals[v])
			ilist+=map(str,curs)
			print "V",v[0],"I",vals[v]
		print "c("+",".join(vlist)+")"
		print "c("+",".join(ilist)+")"
		print "photocurrent"
		print "c("+",".join(pcurs)+")"
		print "darkcurrent"
		print "c("+",".join(dcurs)+")"
if __name__=="__main__":
	parser = optparse.OptionParser()
	parser.add_option("-c",type="int",dest="count",default=-1)
	parser.add_option("--tag",type="string",dest="tag",default="None")
	parser.add_option("-v",dest="verbose",action="store_true")
	options,args = parser.parse_args(sys.argv)
	current_values(options)
