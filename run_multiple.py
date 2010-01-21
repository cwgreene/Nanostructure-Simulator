import os,time

start = -1.0
end   = 1.
max   = 10
total_time = time.time()
tag = "CellsAccounted"
for x in range(0,max):
	start_time = time.time()
	if(max >1):
		voltage = start+x*(end-start)/(1.*(max-1))
	else:
		voltage = start
	print "running",voltage,"/",x+1,"out of",max
	os.system("python monte.py "
		+"-c 199"       + " "
		+"--size=1"     + " "
		+"--scale=1.0"  + " "
		+"--particles=1"+ " "
		+"-V " +str(voltage) +" "
		+"--tag="+tag +" "
		)
	elapse = time.time()-start_time
	print "finished",elapse
print "total:",time.time()-total_time
