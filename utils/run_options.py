import optparse
def create_options(args):
	parser = optparse.OptionParser()
	parser.add_option("-c",type="int",dest="c",default=200)
	parser.add_option("--doping",type="float",dest="doping",default=None)
	parser.add_option("--start",type="float",dest="start",default=-1.)
	parser.add_option("--end",type="float",dest="end",default=1.0)
	parser.add_option("--runs",type="int",dest="runs",default=10)
	parser.add_option("--length",type="float",dest="length",default=10**-6)
	parser.add_option("--tag",type="string",dest="tag",default="Default")
	parser.add_option("--time",type="float",dest="time",default=10**-15)
