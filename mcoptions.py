import sys,optparse
def get_options():
	parser = optparse.OptionParser()
	parser.add_option("-c",type="int",dest="num",default=200)
	options,args = parser.parse_args(sys.argv)	
	return options
