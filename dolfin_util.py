import numpy as np
from dolfin import *
import traceback
from random import randint as rndi

def delete(swig):
	swig.__swig_destroy__(swig)

def feval(function,points):
	point = np.array(points)
	values = np.zeros(len(points))
	function.eval(values,points)
	return values

def feval_p(function,point):
	return feval(function,points)[0]

def cell_index(mesh,point):
	point = Point(*point)
	cells = ArrayUInt()
	mesh.intersection(point,cells)
	result = cells[0]
	delete(cells)
	return result

def vert_index(mesh,point):
	indices = mesh.cells()[cell_index(mesh,point)]
	return indices[0]
	#return indices[rndi(0,len(indices)-1)]

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
