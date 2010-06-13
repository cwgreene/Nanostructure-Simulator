import os
import sys
import re
import pylab

def parse_trajectory_line(line):
	trajectory = []
	for x,y in re.findall("\(([0-9.]+), ([0-9.]+)\)",line):
		trajectory.append((float(x),float(y)))
	return trajectory

def generate_trajectories(file):
	#get rid fo two first lines
	file.readline()
	file.readline()
	#parse each line
	for line in file:
		yield parse_trajectory_line(line)

def open_trajectory_file(n):
	for filename in os.listdir("results"):
		if re.match(str(n)+"traj",filename):
			return open("results/"+filename)
	raise "File not found"
	
def display_trajectories(n):
	input =""
	file = open_trajectory_file(n)
	trajectory_gen = generate_trajectories(file)
	trajectory = trajectory_gen.next()
	interactive = True
	i = 0
	while input != 'q':
		first = map(lambda x: x[0],trajectory)
		second = map(lambda x: x[1],trajectory)
		pylab.plot(first,second)
		if interactive:
			input = raw_input()
		if input == "go":
			i += 1
			interactive=False
		if i %100 == 0:
			print i
			raw_input()
		try:
			trajectory=trajectory_gen.next()
		except:
			print "Done"
			break

if __name__=="__main__":
	display_trajectories(sys.argv[1])
