import ctypes
import materials
lib = ctypes.cdll.LoadLibrary("c_optimized/move_particles.so")

lib.new_polygon.restype = ctypes.POINTER(ctypes.c_int)
lib.init_particle_list.restype = ctypes.POINTER(ctypes.c_int)
lib.init_dead_list.restype = ctypes.POINTER(ctypes.c_int)
lib.update_density.restype = ctypes.c_double

def print_list(string,list):
	print string,list,len(list)

def move_particles(particles_pos,p_mass,p_charge,p_id,p_live,
					 c_efield,nextDensity,
					 dt,length_scale):
	print_list("Next Density",nextDensity)
	print_list("Particles Pos",particles_pos)
	print_list("P id",p_id)
	print_list("P charge",p_charge)
	print_list("efield",c_efield)
	print "NextDensity",nextDensity,len(nextDensity)
	lib.move_particles(particles_pos.ctypes.data,
			   p_mass.ctypes.data,
			   p_charge.ctypes.data,
			   p_id.ctypes.data,
			   p_live,
			   c_efield.ctypes.data,
			   nextDensity.ctypes.data,
			   ctypes.c_double(dt),
			   ctypes.c_double(length_scale))

def update_density(particles_pos,p_id,nextDensity,p_charge,
		   p_live,p_dead,boundary,kdt):
	return lib.update_density(particles_pos.ctypes.data,
			   p_id.ctypes.data,
			   nextDensity.ctypes.data,
			   p_charge.ctypes.data,
			   p_live,
			   p_dead,
			   boundary,
			   kdt)
def new_polygon(points):
	return lib.new_polygon(points.ctypes.data,len(points))
def init_particle_list(n):
	return lib.init_particle_list(n)
def init_dead_list():
	return lib.init_dead_list()

class System():
	pass

class Particles():
	def __init__(self,nparticles):#particles,p_mass,p_charge,p_id,p_live):
		self.nparticles = nparticles
		self.p_id = zeros((nparticles,1),'int')
		self.p_charge = zeros((nparticles,1),'int')
		self.p_mass = zeros((nparticles,1))
		self.p_live = move_particles_c.init_particle_list(len(mesh.particles))
		self.p_dead = move_particles_c.init_dead_list()


def init_system(mesh):
	system = System()
	ntypepy = materials.Silicon()
	ptypepy = materials.Silicon()

	ntype = lib.new_material(ntypepy.electron_mass,
				  ntypepy.random_momentum("ntype"))
	ptype = lib.new_material(ptypepy.hole_mass,
					ntypepy.random_momentum("ntype"))
	materials = lib.material_pair(ptype,ntype)
	system.materials = materials
	system.
