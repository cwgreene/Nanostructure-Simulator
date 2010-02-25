import numpy as np
from dolfin import *
import traceback
from random import randint as rndi
import itertools as it

def boundary_dict(mesh):
	cells = ArrayUInt()
	bmesh = BoundaryMesh(mesh)
	boundary_coordinates = []
	for x in bmesh.coordinates():
		boundary_coordinates.append(np.array(x))
	return boundary_coordinates

def boundary_id_dict(mesh,boundary):
	cells = ArrayUInt()
	bmesh = BoundaryMesh(mesh)
	boundary_coordinates = {}
	bcprime = {}
	for x in boundary:
		bcprime[tuple(x)] = True
	for x,i in it.izip(mesh.coordinates(),it.count()):
		if tuple(x) in bcprime:
			boundary_coordinates[i] = True
	return boundary_coordinates
	

def closest_exit(boundary_coordinates,point):
	disp = np.array(boundary_coordinates[0]) -point
	min = np.dot(disp,disp)
	min_x = np.array(boundary_coordinates[0])
	for pos in boundary_coordinates[1:]:
		disp = np.array(pos) - point
		dist = np.dot(disp,disp)
		if dist < min:
			min = dist
			min_x = pos
	return min_x
	

def delete(swig):
	swig.__swig_destroy__(swig)

def feval(function,points):
	point = np.array(points)
	values = np.zeros(len(points))
	function.eval(values,points)
	return values

values2 = np.zeros(2)
def feval2(function,points):
	point = np.array(points)
	function.eval(values2,points)
	return values2


def feval_p(function,point):
	return feval(function,points)[0]

def cell_index(mesh,apoint):
	point = Point(*apoint)
	cells = ArrayUInt()
	mesh.intersection(point,cells)
	result = cells[0]
	delete(cells)
	return result

def vert_index(mesh,point):
	indices = mesh.cells()[cell_index(mesh,point)]
	coord = mesh.coordinates()
	#return indices[0]
	disp = point-coord[indices[0]]
	min = np.dot(disp,disp)
	min_id = indices[0]
	for id in indices[1:]:
		disp = point-coord[id]
		dist = np.dot(disp,disp)
		if dist < min:
			min = dist
			min_id=id
	return min_id

def set_cell(mesh,function,point,value):
	index = vert_index(mesh,point)
	func_array = function.vector().array()
	func_array[index] = float(value)
	function.vector().set(func_array)

def set_cellid(mesh,function,id,value):
	func_array = function.vector().array()
	func_array[id] = float(value)
	function.vector().set(func_array)

def get_cell(mesh,function,point):
	id = vert_index(mesh,point)
	array = function.vector().array()
	return function.vector().array()[id]

def get_vec(mesh,function,point):
	#print id,function.vector().array().size
	return feval(function,point)

def get_vec2(mesh,function,point):
	#print id,function.vector().array().size
	return feval2(function,point)


def get_vecs(mesh,function,points):
	#print id,function.vector().array().size
	return feval(function,points)


def alter_cell(mesh,function,point,delta):
	value = get_cell(mesh,function,point)
	set_cell(mesh,function,point,value+delta)

def alter_cellid(mesh,function,id,delta):
	value = function.vector().array()[id]
	set_cellid(mesh,function,id,value+delta)

def out_of_bounds(mesh,pos):
	inter =  ArrayUInt()
	mesh.intersection(Point(*pos),inter)
	if(inter.size() == 0):
		delete(inter)
		return True
	delete(inter)
	return False
