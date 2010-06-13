class Average():
	def __init__(self,value=0,count=0):
		self.count = count
		self.value = value
	def __add__(self,x):
		return Average(self.value+x,self.count+1)
	def __str__(self):
		if self.count == 0:
			return "No count"
		return str(self.value/self.count)
avg_momentum = Average()
avg_force = Average()
avg_dx = Average()
avg_charge = Average()

reap_time = Average()
replenish_time=0
mesh_lookup_time =0 

def print_stats(current):
	print "Current:",current
	
	print "Reaper:",reap_time
	print "Replenish took:",replenish_time
	print "Mesh lookup time:",mesh_lookup_time
	print "Average Momentum:",avg_momentum
	print "Average Force:",avg_force
	print "Average dx:",avg_dx
	print "Average charge:",avg_charge

