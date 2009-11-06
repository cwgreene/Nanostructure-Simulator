import os,time

start = -.001
end   = .001
max   = 4
for x in range(0,max):
	start_time = time.time()
	voltage = start+x*(end-start)/(1.*(max-1))
	print "running",voltage,"/",x+1,"out of",max
	os.system("python monte.py -c 400 --scale=1.2 -V " +str(voltage)) 
	print "finished",time.time()-start_time
