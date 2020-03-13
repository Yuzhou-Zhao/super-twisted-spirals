import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = [10, 10]
matplotlib.rcParams['font.size'] = 20

Twist_angles = np.array([-15,-10,-5,-3,-1,0,1,3,5,10,15])
layer_number = 10

for alpha in Twist_angles:
	ax = plt.subplot(111, projection='polar')
	ax.grid(False)
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.spines['polar'].set_visible(False)
	ax.set_rlim(0,2*layer_number)

	k = 1.-alpha/360.
	for n in range(1,layer_number):
		two_n_pi = 2*n*np.pi
		beta0 = (-np.pi+two_n_pi)/k
		beta1 = (-np.pi/3.+two_n_pi)/k
		beta2 = (np.arctan(np.sqrt(3)+2*np.sqrt(3)/(3*n))+two_n_pi)/k
		beta3 = (np.pi+two_n_pi)/k

		theta1 = np.arange(beta0, beta1, 0.01)
		r1 = n*k/np.cos(k*theta1+2*np.pi/3.+two_n_pi)

		theta2 = np.arange(beta1, beta2, 0.01)
		r2 = n*k/np.cos(k*theta2+two_n_pi)

		theta3 = np.arange(beta2, beta3, 0.01)
		r3 = (n+1)*k/np.cos(k*theta3-2*np.pi/3.+two_n_pi)

		ax.plot(theta1, r1,color = "black")
		ax.plot(theta2, r2,color = "black")
		ax.plot(theta3, r3,color = "black")

	ax.set_title("Twist Angles = "+str(alpha), va='bottom')
	plt.savefig("Twist_angles_10_"+str(alpha) + ".png", format='png', dpi=150, transparent=True)
	plt.clf()
	plt.close()
