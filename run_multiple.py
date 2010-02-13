import os,time,sys
import run_options

options = run_options.create_options(sys.argv)
start = -1
end   = 4.
max   = 15
total_time = time.time()
tag = "Test2"

for x in range(0,max):
	start_time = time.time()
	if(max >1):
		voltage = start+x*(end-start)/(1.*(max-1))
	else:
		voltage = start
	print "running",voltage,"/",x+1,"out of",max
	os.system("python monte.py "
		+"-c 200"       + " "
		+"--size=1"     + " "
		+"--scale=1.0"  + " "
		+"--particles=1"+ " "
		+"-V " +str(voltage) +" "
		+"--tag="+tag +" "
		)
	elapse = time.time()-start_time
	print "finished",elapse
print "total:",time.time()-total_time
