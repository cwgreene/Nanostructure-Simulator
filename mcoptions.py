import sys,optparse
def get_options():
	parser = optparse.OptionParser()
	parser.add_option("-c",type="int",dest="num",default=200)
	parser.add_option("-V",type="float",dest="V",default=0)
	parser.add_option("-d",type="string",dest="datadir",default="data")
	parser.add_option("--endV",type="float",dest="endV",default=0)
	parser.add_option("--numsteps",type="int",dest="steps",default=1)
	parser.add_option("--scale",type="float",dest="scale",default=1)
	parser.add_option("--size",type="int",dest="size",default=1)
	parser.add_option("--particles",type="int",dest="gen_num",default=10)
	parser.add_option("--length",type="float",dest="length",default=10**-6)
	parser.add_option("--doping",type="float",dest="doping",default=10**-24)
	parser.add_option("--tag",type="string",dest="tag",default="None")
	parser.add_option("--dt",type="float",dest="dt",default=10**-13)
	options,args = parser.parse_args(sys.argv)
	return options
