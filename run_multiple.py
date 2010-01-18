import os,time

start = 0.0
end   = 100000
max   = 10
total_time = time.time()
for x in range(0,max):
	start_time = time.time()
	if(max >1):
		voltage = start+x*(end-start)/(1.*(max-1))
	else:
		voltage = start
	print "running",voltage,"/",x+1,"out of",max
	os.system("python monte2.py "
		+"-c 196 "
		+"--size=1 "
		+"--scale=1.0 "
		+"--particles=5 "
		+"-V " +str(voltage) +" "
		)
	print "finished",time.time()-start_time
print "total:",time.time()-total_time
