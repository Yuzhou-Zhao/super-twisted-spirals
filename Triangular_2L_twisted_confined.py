import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
matplotlib.rcParams['figure.figsize'] = [10, 10]
matplotlib.rcParams['font.size'] = 20

Twist_angles = np.array([-15,-10,10,15])
layer_number = 30

### second layer is on top of the first one
### L2_angle_offset needs to be between [0,120)
L2_angle_offset = 59.9

def x3_is_min_x1_x2_with_offet(	x3,
								x1,
								x2,
								idx_x1_initial,
								idx_x1_end,
								offset
								):
	x3[idx_x1_initial:idx_x1_end+1] = np.minimum(x1[idx_x1_initial:idx_x1_end+1],x2[idx_x1_initial+offset:idx_x1_end+offset+1])	


for alpha in Twist_angles:
	ax = plt.subplot(111, projection='polar')
	# ax.grid(False)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.spines['polar'].set_visible(False)
	ax.set_rlim(0,2.2*layer_number)

	k = 1.-alpha/360.

	beta_initial = 0.5*np.pi/k
	beta_end = (2*layer_number*np.pi+4*np.pi/3.)/k
	### theta_idx_spacing determines the smallest angle difference.
	theta_idx_spacing = 3600
	### idx_offset_alpha is used to compensate the missing angle between a full period to 2pi when compare r
	idx_offset_alpha = int(theta_idx_spacing*alpha*1.0/360.)
	print idx_offset_alpha
	theta_spacing = 2.*np.pi/theta_idx_spacing
	L2_theta_offset_k = (L2_angle_offset)*np.pi/(k*180.)
	L2_idx_offset_k = int(round(L2_theta_offset_k/theta_spacing)) 


	### theta is the angle of the bottom layer
	### as well as the original angle of other layer before rotation
	theta = np.arange(beta_initial, beta_end, theta_spacing)
	### the same idx from theta was used in theta2 to realize the rotation of second layer
	theta2 = np.arange(beta_initial+L2_theta_offset_k, beta_end+L2_theta_offset_k, theta_spacing)
	### r is used to store boundary confined r
	r = np.full(theta.shape, 0.)
	r2 = np.full(theta.shape, 0.)
	### r_temp is used to store boundary unconfined r
	r_temp = np.full(theta.shape, 0.)
	r2_temp = np.full(theta.shape, 0.)


	### Set up the bottom layer
	n = layer_number+1
	two_n_pi = 2*n*np.pi
	beta0 = (-np.pi+two_n_pi)/k
	beta1 = (-np.pi/3.+two_n_pi)/k
	idx_0 = (np.abs(theta - beta0)).argmin()
	idx_1 = (np.abs(theta - beta1)).argmin()
	idx_end = idx_1 
	r_temp[idx_0:idx_1+1] = n*k/np.cos(k*theta[idx_0:idx_1+1]+2*np.pi/3.+two_n_pi)
	
	beta2_0 = (-np.pi+two_n_pi)/k
	beta2_1 = (-np.pi/3.+two_n_pi)/k
	### Note theta is used here instead of theta2
	idx2_0 = (np.abs(theta - beta2_0)).argmin()
	idx2_1 = (np.abs(theta - beta2_1)).argmin()
	idx2_end = idx_1 
	### Note theta is used here instead of theta2
	r2_temp[idx2_0:idx2_1+1] = (n-0.5)*k/np.cos(k*theta[idx2_0:idx2_1+1]+2*np.pi/3.+two_n_pi)


	### in bilayer: r1_temp is confinced by r2, r2_temp is confinced by r1
	### in trilayer: r1_temp by r2, r2_temp by r3, r3_temp by r1
	### ...
	r[idx_0:idx_1+1] = r_temp[idx_0:idx_1+1]
	r2[idx_1-L2_idx_offset_k:idx_1+1] = r2_temp[idx_1-L2_idx_offset_k:idx_1+1]

	x3_is_min_x1_x2_with_offet(	x3 = r2,
								x1 = r2_temp,
								x2 = r,
								idx_x1_initial = idx_0,
								idx_x1_end = idx_1-L2_idx_offset_k,
								offset = L2_idx_offset_k
								)
	n = layer_number
	two_n_pi = 2*n*np.pi
	beta0 = (-np.pi+two_n_pi)/k
	beta1 = (-np.pi/3.+two_n_pi)/k
	beta2 = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n))+two_n_pi)/k
	beta3 = (np.pi+two_n_pi)/k
	idx_0 = (np.abs(theta - beta0)).argmin()
	idx_1 = (np.abs(theta - beta1)).argmin()
	idx_2 = (np.abs(theta - beta2)).argmin()
	idx_3 = (np.abs(theta - beta3)).argmin()
	r_temp[idx_2:idx_3+1] = (n+1)*k/np.cos(k*theta[idx_2:idx_3+1]-2*np.pi/3.+two_n_pi)
	r_temp[idx_1:idx_2+1] = n*k/np.cos(k*theta[idx_1:idx_2+1]+two_n_pi)
	r_temp[idx_0:idx_1+1] = n*k/np.cos(k*theta[idx_0:idx_1+1]+2*np.pi/3.+two_n_pi)
	

	beta2_0 = (-np.pi+two_n_pi)/k
	beta2_1 = (-np.pi/3.+two_n_pi)/k
	beta2_2 = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n-1.5))+two_n_pi)/k
	beta2_3 = (np.pi+two_n_pi)/k
	### Note theta is used here instead of theta2
	idx2_0 = (np.abs(theta - beta2_0)).argmin()
	idx2_1 = (np.abs(theta - beta2_1)).argmin()
	idx2_2 = (np.abs(theta - beta2_2)).argmin()
	idx2_3 = (np.abs(theta - beta2_3)).argmin()
	### Note theta is used here instead of theta2
	r2_temp[idx2_2:idx2_3+1] = (n+1-0.5)*k/np.cos(k*theta[idx2_2:idx2_3+1]-2*np.pi/3.+two_n_pi)
	r2_temp[idx2_1:idx2_2+1] = (n-0.5)*k/np.cos(k*theta[idx2_1:idx2_2+1]+two_n_pi)
	r2_temp[idx2_0:idx2_1+1] = (n-0.5)*k/np.cos(k*theta[idx2_0:idx2_1+1]+2*np.pi/3.+two_n_pi)


	n_t = layer_number-1
	two_n_t_pi = 2*n_t*np.pi
	beta2_t = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n))+two_n_t_pi)/k
	beta3_t = (np.pi+two_n_t_pi)/k
	idx_2_t = (np.abs(theta - beta2_t)).argmin()
	idx_3_t = (np.abs(theta - beta3_t)).argmin()
	r_temp[idx_2_t:idx_3_t+1] = (n_t+1)*k/np.cos(k*theta[idx_2_t:idx_3_t+1]-2*np.pi/3.+two_n_t_pi)
	
	beta2_2_t = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n-1.5))+two_n_t_pi)/k
	beta2_3_t = (np.pi+two_n_t_pi)/k
	### Note theta is used here instead of theta2
	idx2_2_t = (np.abs(theta - beta2_2_t)).argmin()
	idx2_3_t = (np.abs(theta - beta2_3_t)).argmin()
	### Note theta is used here instead of theta2
	r2_temp[idx2_2_t:idx2_3_t+1] = (n_t+1-0.5)*k/np.cos(k*theta[idx2_2_t:idx2_3_t+1]-2*np.pi/3.+two_n_t_pi)



	r[idx_end+1-theta_idx_spacing:idx_3+1] = r_temp[idx_end+1-theta_idx_spacing:idx_3+1]
	x3_is_min_x1_x2_with_offet(	x3 = r,
								x1 = r_temp,
								x2 = r,
								idx_x1_initial = idx_end-theta_idx_spacing-L2_idx_offset_k,
								idx_x1_end = idx_end-theta_idx_spacing,
								offset = theta_idx_spacing
								)
	
	x3_is_min_x1_x2_with_offet(	x3 = r2,
								x1 = r2_temp,
								x2 = r,
								idx_x1_initial = idx_1,
								idx_x1_end = idx_3,
								offset = L2_idx_offset_k
								)
	x3_is_min_x1_x2_with_offet(	x3 = r,
								x1 = r_temp,
								x2 = r2,
								idx_x1_initial = idx_0-10,
								idx_x1_end = idx_end-theta_idx_spacing-L2_idx_offset_k,
								offset = theta_idx_spacing-L2_idx_offset_k
								)
	x3_is_min_x1_x2_with_offet(	x3 = r2,
								x1 = r2_temp,
								x2 = r,
								idx_x1_initial = idx_0-10,
								idx_x1_end = idx_1,
								offset = L2_idx_offset_k
								)
	
	### Other layers
	for n in range(layer_number-1,0,-1):
		two_n_pi = 2*n*np.pi
		beta0 = (-np.pi+two_n_pi)/k
		beta1 = (-np.pi/3.+two_n_pi)/k
		beta2 = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n))+two_n_pi)/k
		beta3 = (np.pi+two_n_pi)/k
		idx_0 = (np.abs(theta - beta0)).argmin()
		idx_1 = (np.abs(theta - beta1)).argmin()
		idx_2 = (np.abs(theta - beta2)).argmin()
		idx_3 = (np.abs(theta - beta3)).argmin()
		r_temp[idx_2:idx_3+1] = (n+1)*k/np.cos(k*theta[idx_2:idx_3+1]-2*np.pi/3.+two_n_pi)
		r_temp[idx_1:idx_2+1] = n*k/np.cos(k*theta[idx_1:idx_2+1]+two_n_pi)
		r_temp[idx_0:idx_1+1] = n*k/np.cos(k*theta[idx_0:idx_1+1]+2*np.pi/3.+two_n_pi)
		
		beta2_0 = (-np.pi+two_n_pi)/k
		beta2_1 = (-np.pi/3.+two_n_pi)/k
		beta2_2 = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n-1.5))+two_n_pi)/k
		beta2_3 = (np.pi+two_n_pi)/k
		### Note theta is used here instead of theta2
		idx2_0 = (np.abs(theta - beta2_0)).argmin()
		idx2_1 = (np.abs(theta - beta2_1)).argmin()
		idx2_2 = (np.abs(theta - beta2_2)).argmin()
		idx2_3 = (np.abs(theta - beta2_3)).argmin()
		### Note theta is used here instead of theta2
		r2_temp[idx2_2:idx2_3+1] = (n+1-0.5)*k/np.cos(k*theta[idx2_2:idx2_3+1]-2*np.pi/3.+two_n_pi)
		r2_temp[idx2_1:idx2_2+1] = (n-0.5)*k/np.cos(k*theta[idx2_1:idx2_2+1]+two_n_pi)
		r2_temp[idx2_0:idx2_1+1] = (n-0.5)*k/np.cos(k*theta[idx2_0:idx2_1+1]+2*np.pi/3.+two_n_pi)


		n_t = n-1
		two_n_t_pi = 2*n_t*np.pi
		beta2_t = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n))+two_n_t_pi)/k
		beta3_t = (np.pi+two_n_t_pi)/k
		idx_2_t = (np.abs(theta - beta2_t)).argmin()
		idx_3_t = (np.abs(theta - beta3_t)).argmin()
		r_temp[idx_2_t:idx_3_t+1] = (n_t+1)*k/np.cos(k*theta[idx_2_t:idx_3_t+1]-2*np.pi/3.+two_n_t_pi)
		
		beta2_2_t = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n-1.5))+two_n_t_pi)/k
		beta2_3_t = (np.pi+two_n_t_pi)/k
		### Note theta is used here instead of theta2
		idx2_2_t = (np.abs(theta - beta2_2_t)).argmin()
		idx2_3_t = (np.abs(theta - beta2_3_t)).argmin()
		### Note theta is used here instead of theta2
		r2_temp[idx2_2_t:idx2_3_t+1] = (n_t+1-0.5)*k/np.cos(k*theta[idx2_2_t:idx2_3_t+1]-2*np.pi/3.+two_n_t_pi)
		

		x3_is_min_x1_x2_with_offet(	x3 = r,
								x1 = r_temp,
								x2 = r2,
								idx_x1_initial = idx_0+L2_idx_offset_k+idx_offset_alpha,
								idx_x1_end = idx_3,
								offset = theta_idx_spacing-L2_idx_offset_k
								)

		x3_is_min_x1_x2_with_offet(	x3 = r2,
								x1 = r2_temp,
								x2 = r,
								idx_x1_initial = idx_1,
								idx_x1_end = idx_3,
								offset = L2_idx_offset_k
								)
		### for alpha>0, +idx_offset_alpha fills up a missing part 
		x3_is_min_x1_x2_with_offet(	x3 = r,
								x1 = r_temp,
								x2 = r2,
								idx_x1_initial = idx_0-10,
								idx_x1_end = idx_0+L2_idx_offset_k+idx_offset_alpha,
								offset = theta_idx_spacing-L2_idx_offset_k
								)
		x3_is_min_x1_x2_with_offet(	x3 = r2,
								x1 = r2_temp,
								x2 = r,
								idx_x1_initial = idx_0-10,
								idx_x1_end = idx_1,
								offset = L2_idx_offset_k
								)



	r_temp[0:idx_0+1] = k/np.cos(k*theta[0:idx_0+1]-2*np.pi/3.)
	r2_temp[0:idx2_0+1] = 0.5*k/np.cos(k*theta[0:idx_0+1]-2*np.pi/3.)
	
	x3_is_min_x1_x2_with_offet(	x3 = r,
								x1 = r_temp,
								x2 = r2,
								idx_x1_initial = 0,
								idx_x1_end = idx_0,
								offset = theta_idx_spacing-L2_idx_offset_k
								)
	x3_is_min_x1_x2_with_offet(	x3 = r2,
								x1 = r2_temp,
								x2 = r,
								idx_x1_initial = 0,
								idx_x1_end = idx_0,
								offset = L2_idx_offset_k
								)

	theta_rotate = theta-beta_end -np.pi/2.
	theta2_rotate = theta2-beta_end -np.pi/2.

	theta_rotate = np.insert(theta_rotate, 0, theta_rotate[0])
	r = np.insert(r, 0, 0.00)
	theta2_rotate = np.insert(theta2_rotate, 0, theta2_rotate[0])
	r2 = np.insert(r2, 0, 0.00)

	ax.plot(theta_rotate, r,color = "black")
	ax.plot(theta2_rotate, r2,color = "red")

	# plt.fill_between(theta_rotate, r,  alpha=0.2)


	ax.set_title("Twist Angles = "+str(alpha), va='bottom')
	plt.savefig("Tri_BC2_60offset_" +str(alpha) + ".png", format='png', dpi=150, transparent=True)
	plt.clf()
	plt.close()
