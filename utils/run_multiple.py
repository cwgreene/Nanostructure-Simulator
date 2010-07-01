import os,time,sys
import run_options

options = run_options.create_options(sys.argv)
start = 0.31
end   = 1.0
max   = 7
total_time = time.time()
tag = "C++_400_triangle_3"

for x in range(0,max):
	start_time = time.time()
	if(max >1):
		voltage = start+x*(end-start)/(1.*(max-1))
	else:
		voltage = start
	print "running",voltage,"/",x+1,"out of",max
	data_dir = str(voltage)[:5].replace(".","_").replace("-","_n")
	os.system("python monte.py "
		+"-c 40"       + " "
		+"--size=1"     + " "
		+"--scale=1.0"  + " "
		+"--particles=100"+ " "
		+"-V " +str(voltage) +" "
		+"--tag="+tag +" "
		+"-d data"+data_dir+" "
		)
	elapse = time.time()-start_time
	print "finished",elapse
print "total:",time.time()-total_time
