import numpy as np
# cystal structure parameter MoS2
a = 3.161			# A
thickness = 3.01227 # A distance between S-S
c = 12.295			# A
z_spacing = c/2

alpha = 0. #degree
k = 1-alpha/360.

def upper_spiral(range_x, layer_number, element, ABC_site, XMX_position,z_shift):
				# range_x: number of unit cells along x direction
				# Layer_number: the nth layer of spirals
				# Mo, W, S, Se, etc.
				# ABC_site: in hexagonal unit cell A = 0 B = 1 C = -1
				# XMX_position M = 0 X = 1 and -1
				# z_shift: shift in z position
	# this function convert a hexagonal cell to orthorombic cell 
	
	total_atoms_temp = 0
	for i in range(-range_x,range_x+1):
		for j in range(1,range_x+1):
			x = a*(i+0.5*ABC_site)
			y = a*(j*np.sqrt(3)+ABC_site*np.sqrt(3.)/6.)
			
			if ((y+np.sqrt(3)*x)<(a*(range_x-0.05)*np.sqrt(3))):
				if  ((y-np.sqrt(3)*x)<(a*(range_x-0.05)*np.sqrt(3))):
					z = z_spacing*(np.arctan2(x,y)/np.pi-layer_number)+thickness*XMX_position/2.+z_shift
					r = np.sqrt(x**2+y**2)+np.sqrt(1-k**2)*thickness*XMX_position/2.
					t = np.arctan2(y,x)+np.pi*layer_number

					r_cone_projection = k*r
					t_cone_projection = t/k
					x_cone_projection = r_cone_projection*np.cos(t_cone_projection)
					y_cone_projection = r_cone_projection*np.sin(t_cone_projection)
					z_cone_projection = z - np.sqrt(x**2+y**2)*np.sqrt(1-k**2) +k*thickness*XMX_position/2.

					print element,"    ", x_cone_projection,"    ", y_cone_projection,"    ", z_cone_projection
					total_atoms_temp = total_atoms_temp+1
					
			x = a*(i+0.5+0.5*ABC_site)
			y = a*((j+0.5)*np.sqrt(3)+ABC_site*np.sqrt(3)/6)
			
			if ((y+np.sqrt(3)*x)<(a*(range_x-0.05)*np.sqrt(3))):
				if  ((y-np.sqrt(3)*x)<(a*(range_x-0.05)*np.sqrt(3))):
					z = z_spacing*(np.arctan2(x,y)/np.pi-layer_number)+thickness*XMX_position/2.+z_shift		
					r = np.sqrt(x**2+y**2)+np.sqrt(1-k**2)*thickness*XMX_position/2.
					t = np.arctan2(y,x)+np.pi*layer_number

					r_cone_projection = k*r
					t_cone_projection = t/k
					x_cone_projection = r_cone_projection*np.cos(t_cone_projection)
					y_cone_projection = r_cone_projection*np.sin(t_cone_projection)
					z_cone_projection = z - np.sqrt(x**2+y**2)*np.sqrt(1-k**2) +k*thickness*XMX_position/2.

					# print element,"    ", x,"    ", y,"    ", z
					print element,"    ", x_cone_projection,"    ", y_cone_projection,"    ", z_cone_projection
					total_atoms_temp = total_atoms_temp+1
	return total_atoms_temp

def lower_spiral(range_x, range_y, range_y_d, layer_number, element, ABC_site, XMX_position,z_shift):
				# range_x: number of unit cells along x direction
				# Layer_number: the nth layer of spirals
				# Mo, W, S, Se, etc.
				# ABC_site: in hexagonal unit cell A=0 B = 1 C = -1
				# XMX_position M = 0 X = 1 and -1
				# z_shift: shift in z position
	# this function convert a hexagonal cell to orthorombic cell

	total_atoms_temp = 0
	for i in range(-(range_x+range_y),(range_x+range_y+range_y_d+1)):
		for j in range(-range_y,1):
			x = a*(i+0.5*ABC_site)
			y = a*(j*np.sqrt(3)+ABC_site*np.sqrt(3)/6)
			if ((y+np.sqrt(3)*x)<(a*(range_x+range_y_d-0.05)*np.sqrt(3))):
				if  ((y-np.sqrt(3)*x)<(a*(range_x-0.05)*np.sqrt(3))):
					z = z_spacing*(np.arctan2(x,y)/np.pi-layer_number)+thickness*XMX_position/2.+z_shift
					if (x>=0):
						z = z_spacing*(np.arctan2(x,y)/np.pi-layer_number-1-1)+thickness*XMX_position/2.+z_shift
					r = np.sqrt(x**2+y**2)+np.sqrt(1-k**2)*thickness*XMX_position/2.
					if y <0:
						t = np.arctan2(y,x)+np.pi*layer_number+2*np.pi
					else:
						if x > 0:
							t = np.arctan2(y,x)+np.pi*layer_number+alpha*np.pi/180.
						else:
							t = np.arctan2(y,x)+np.pi*layer_number
	
					r_cone_projection = k*r
					t_cone_projection = t/k
					x_cone_projection = r_cone_projection*np.cos(t_cone_projection)
					y_cone_projection = r_cone_projection*np.sin(t_cone_projection)
					z_cone_projection = z - np.sqrt(x**2+y**2)*np.sqrt(1-k**2) +k*thickness*XMX_position/2.

					# print element,"    ", x,"    ", y,"    ", z
					print element,"    ", x_cone_projection,"    ", y_cone_projection,"    ", z_cone_projection
					total_atoms_temp = total_atoms_temp+1

			x = a*(i+0.5+0.5*ABC_site)
			y = a*((j+0.5)*np.sqrt(3)+ABC_site*np.sqrt(3)/6)
			if ((y+np.sqrt(3)*x)<(a*(range_x+range_y_d-0.05)*np.sqrt(3))):
				if  ((y-np.sqrt(3)*x)<(a*(range_x-0.05)*np.sqrt(3))):
					z = z_spacing*(np.arctan2(x,y)/np.pi-layer_number)+thickness*XMX_position/2.+z_shift
					if (x>=0):
						z = z_spacing*(np.arctan2(x,y)/np.pi-layer_number-1-1)+thickness*XMX_position/2.+z_shift
					r = np.sqrt(x**2+y**2)+np.sqrt(1-k**2)*thickness*XMX_position/2.
					if y <0:
						t = np.arctan2(y,x)+np.pi*layer_number+2*np.pi
					else:
						if x > 0:
							t = np.arctan2(y,x)+np.pi*layer_number+alpha*np.pi/180.
						else:
							t = np.arctan2(y,x)+np.pi*layer_number

					r_cone_projection = k*r
					t_cone_projection = t/k
					x_cone_projection = r_cone_projection*np.cos(t_cone_projection)
					y_cone_projection = r_cone_projection*np.sin(t_cone_projection)
					z_cone_projection = z - np.sqrt(x**2+y**2)*np.sqrt(1-k**2) +k*thickness*XMX_position/2.

					# print element,"    ", x,"    ", y,"    ", z
					print element,"    ", x_cone_projection,"    ", y_cone_projection,"    ", z_cone_projection
					total_atoms_temp = total_atoms_temp+1
	return total_atoms_temp

initial_range_x = 10
spacing_A = 10
spacing_B = spacing_A/2
z_shift = 0

total_atoms = 0
total_atoms_W = 0
total_atoms_S = 0 

for kk in range(0,5):
	# upper_spiral(range_x, layer_number, element, ABC_site, XMX_position,z_shift)
	W_number_u = upper_spiral(initial_range_x, 2*kk, "W", 0, 0, z_shift) 	#"W"
	S_number_u1 = upper_spiral(initial_range_x, 2*kk, "S", -1, -1, z_shift) 	#"S"
	S_number_u2 = upper_spiral(initial_range_x, 2*kk, "S", -1, 1, z_shift) 	#"S"

	range_y_temp = spacing_B*(kk+1)
	# lower_spiral(range_x, range_y, range_y_d, layer_number, element, ABC_site, XMX_position,z_shift)
	W_number_l = lower_spiral(initial_range_x, range_y_temp, spacing_A, 2*kk, "W", 0, 0, z_shift)	#W
	S_number_l1 = lower_spiral(initial_range_x, range_y_temp, spacing_A, 2*kk, "S", -1, -1, z_shift)	#S1
	S_number_l2 = lower_spiral(initial_range_x, range_y_temp, spacing_A, 2*kk, "S", -1, 1, z_shift)	#S2
	
	initial_range_x = initial_range_x+spacing_A

	total_atoms_W = total_atoms_W + W_number_u +W_number_l
	total_atoms_S = total_atoms_S + S_number_u1+S_number_u2+S_number_l1+S_number_l2

total_atoms = total_atoms_W+total_atoms_S
print total_atoms
