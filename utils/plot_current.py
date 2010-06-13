import os
import sys
import optparse
import pylab
import re

def create_options(args):
	parser = optparse.OptionParser()
	parser.add_option("-n",type="int",dest="n",default=200)
	options,args = parser.parse_args(sys.argv)
	return options

def load_current(file):
	file.readline()
	file.readline()
	current = []
	for line in file:
		current.append(float(line))
	print len(current)
	return current

def get_file(n):
	for filename in os.listdir('results'):
		if re.match("^"+str(n)+"current",filename):
			return "results/"+filename
	raise "File not found"

opt = create_options(sys.argv)
filename= get_file(opt.n)
print filename
curfile = open(filename)
pylab.plot(load_current(curfile))
raw_input()

