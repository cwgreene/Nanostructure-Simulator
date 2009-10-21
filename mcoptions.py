import sys,optparse
def get_options():
	parser = optparse.OptionParser()
	parser.add_option("-c",type="int",dest="num",default=200)
	parser.add_option("-V",type="float",dest="V",default=0)
	parser.add_option("-d",type="string",dest="datadir",default="data")
	parser.add_option("--endV",type="float",dest="endV",default=0)
	parser.add_option("--numsteps",type="int",dest="steps",default=1)
	options,args = parser.parse_args(sys.argv)
	return options
