import os,time

start = -.002
end   = .001
max   = 20
for x in range(1,max):
	start_time = time.time()
	voltage = start+(end-start)/(1.*max)
	print "running",voltage,"/",x,"out of",max
	os.system("python monte.py -c 400 -V " +str(voltage)) 
	print "finished",time.time()-start_time
