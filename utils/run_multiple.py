import os,time,sys
import run_options

options = run_options.create_options(sys.argv)
start = -1.0
end   = 1.0
max   = 8
total_time = time.time()
tag = "C++_400_scatter"

for x in range(0,max):
	start_time = time.time()
	if(max >1):
		voltage = start+x*(end-start)/(1.*(max-1))
	else:
		voltage = start
	print "running",voltage,"/",x+1,"out of",max
	os.system("python monte.py "
		+"-c 400"       + " "
		+"--size=1"     + " "
		+"--scale=1.0"  + " "
		+"--particles=100"+ " "
		+"-V " +str(voltage) +" "
		+"--tag="+tag +" "
		+"-d data"+str(x)+" "
		)
	elapse = time.time()-start_time
	print "finished",elapse
print "total:",time.time()-total_time
