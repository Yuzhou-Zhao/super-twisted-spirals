import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = [10, 10]
matplotlib.rcParams['font.size'] = 20

Twist_angles = np.array([-10,0,10])
layer_number = 5

for alpha in Twist_angles:
	ax = plt.subplot(111, projection='polar')
	ax.grid(False)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.spines['polar'].set_visible(False)
	ax.set_rlim(0,3*layer_number)

	k = 1.-alpha/360.

	beta_initial = 0.5*np.pi/k
	beta_end = (2*layer_number*np.pi+4*np.pi/3.)/k
	theta_idx_spacing = 3600
	idx_offset_alpha = int(theta_idx_spacing*alpha*1.1/360.)
	# print idx_offset_alpha
	theta_spacing = 2.*np.pi/theta_idx_spacing
	

	theta = np.arange(beta_initial, beta_end, theta_spacing)
	r = np.full(theta.shape, 0.)
	r_temp = np.full(theta.shape, 0.)


	### Bottom layer
	n = layer_number+1
	two_n_pi = 2*n*np.pi
	beta0 = (-np.pi+two_n_pi)/k
	beta1 = (-np.pi/3.+two_n_pi)/k

	idx_0 = (np.abs(theta - beta0)).argmin()
	idx_1 = (np.abs(theta - beta1)).argmin()
	idx_end = idx_1 

	r[idx_0:idx_1+1] = n*k/np.cos(k*theta[idx_0:idx_1+1]+2*np.pi/3.+two_n_pi)
	
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

	### Take care near idx0 in case of Alpha <0
	r[idx_0:idx_end+1-theta_idx_spacing] = np.minimum(r[idx_0+theta_idx_spacing:idx_end+1], r_temp[idx_0:idx_end+1-theta_idx_spacing])

	r[idx_end+1-theta_idx_spacing:idx_3+1] = r_temp[idx_end+1-theta_idx_spacing:idx_3+1]
	r[idx_0:idx_3+1-theta_idx_spacing] = np.minimum(r[idx_0+theta_idx_spacing:idx_3+1], r_temp[idx_0:idx_3+1-theta_idx_spacing])

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

		### For Alpha > 0 case, fill up the part for boundary check 
		### For Alpha < 0 case, this is redundant  	
		r[idx_3-idx_offset_alpha:idx_3+1] = r_temp[idx_3-idx_offset_alpha:idx_3+1]
		r[idx_3-idx_offset_alpha:idx_3+1] = np.minimum(r[idx_3-idx_offset_alpha+theta_idx_spacing:idx_3+1+theta_idx_spacing], r_temp[idx_3-idx_offset_alpha:idx_3+1])

		idx_0_offset = idx_0+theta_idx_spacing
		idx_3_offset = idx_3+1+theta_idx_spacing
		r[idx_0:idx_3+1] = np.minimum(r[idx_0_offset:idx_3_offset], r_temp[idx_0:idx_3+1])


	r[0:idx_0+1] = k/np.cos(k*theta[0:idx_0+1]-2*np.pi/3.)

	theta_rotate = theta-beta_end -np.pi/2.
	theta_rotate = np.insert(theta_rotate, 0, theta_rotate[0])
	r = np.insert(r, 0, 0.00)
	ax.plot(theta_rotate, r,color = "black")
	plt.fill_between(theta_rotate, r,  alpha=0.2)


	ax.set_title("Twist Angles = "+str(alpha), va='bottom')
	plt.savefig("Tri_BC2_5L_" +str(alpha) + ".png", format='png', dpi=150, transparent=True)
	plt.clf()
	plt.close()
