from dolfin import *
import os, inspect
from os import path
from numpy import zeros, array, float

class Probe: 

    def __init__(self, x, V, probe_id): #x in a numpy array which stores a coordinate

        #1. Find cell that contains probe
        self.probe_id = probe_id
        self.x=x
        self.dimension=V.element().space_dimension()
        x_point = Point(*x) #Probe position
        #2. Figure out which cell we are in
        tree=V.mesh().bounding_box_tree()
        ident = tree.compute_first_entity_collision(x_point) #element id
        if ident<V.mesh().num_cells(): #Make sure to only do this on procs that have the cell
            print(ident, MPI.rank(MPI.comm_world))
            self.ident=ident
            dolfin_cell = Cell(V.mesh(), ident)
            self.dolfin_cell = dolfin_cell	
            coordinate_dofs = dolfin_cell.get_vertex_coordinates()
            self.coordinate_dofs=coordinate_dofs
            cell_orientation = dolfin_cell.orientation()
            basis = V.element().evaluate_basis_all(x, coordinate_dofs, cell_orientation) 
            self.basis = basis
            #print(basis)	
        else:
            self.dolfin_cell = None
            print('Cell = ', ident, 'Not on processor:', MPI.rank(MPI.comm_world))
	        	
    def __call__(self, x_, t, newfolder, V):
        #When I call this, I want the probe to be evaluated
        
        if self.dolfin_cell is not None:
            x_eval = self.evaluate_func(x_, V)  
            #2. Call print_snapshot_to_file to print the probe value and its ID number to its own file
            self.print_snapshot_to_file(x_eval, t, newfolder)
        
    def print_snapshot_to_file(self, x_eval, t, newfolder): 
        fsnap = open(path.join(newfolder, '{}_probe_snapshots.txt'.format(self.probe_id)), 'a')
        print(format(t, '.4f'), x_eval, file=fsnap, sep=',')  
        fsnap.close()  
            
    def evaluate_func(self, x_, V):
        #1. Evaluate the function at the point in the element
        #u(vertex1)*basisFunc[1]+u(vertex2)+ etc. and make sure we only return the one value that comes from the relevant
        #print(Vector(x_['u0'])[4])
        vert_indices = V.dofmap().cell_dofs(self.ident)
        #print(vert_indices)
        #print(x_['u0'].get_local(vert_indices))
        u0_vertices = x_['u0'].get_local(vert_indices)
        u1_vertices = x_['u1'].get_local(vert_indices)
        u2_vertices = x_['u2'].get_local(vert_indices)
        #print(u0_vertices)
        x_eval1 = 0
        x_eval2 = 0
        x_eval3 = 0
        for i in range(self.dimension):
            x_eval1 = u0_vertices[i]*self.basis[i]+x_eval1
            x_eval2 = u1_vertices[i]*self.basis[i]+x_eval2
            x_eval3 = u2_vertices[i]*self.basis[i]+x_eval3
        x_eval = array([x_eval1, x_eval2, x_eval3])
        return x_eval 
