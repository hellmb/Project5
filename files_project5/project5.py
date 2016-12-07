from pylab import *
from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes

# load text-files
exp_euler10 = loadtxt('explicit_euler10.txt')
exp_euler100 = loadtxt('explicit_euler100.txt')
exp_euler1000 = loadtxt('explicit_euler1000.txt')
exp_euler10000 = loadtxt('explicit_euler10000.txt')

imp_euler10 = loadtxt('implicit_euler10.txt')
imp_euler100 = loadtxt('implicit_euler100.txt')
imp_euler1000 = loadtxt('implicit_euler1000.txt')
imp_euler10000 = loadtxt('implicit_euler10000.txt')

imp_CN10 = loadtxt('implicit_CN10.txt')
imp_CN100 = loadtxt('implicit_CN100.txt')
imp_CN1000 = loadtxt('implicit_CN1000.txt')
imp_CN10000 = loadtxt('implicit_CN10000.txt')


n = 100
time = linspace(0, 1, n)
x = linspace(0, 1, n)
nn = linspace(1, 100, 100)

# analytical solution
analytical_solution = x + (2.0/pi) * sum( (-1.0**n/n) * sin(n * pi * x) * exp( -(n * pi)**2 * time) )

explicit_euler = False
if explicit_euler:

	# plot all curves in same figure
	plot(time, exp_euler10, label=r'Timestep $= 10$')
	plot(time, exp_euler100, label=r'Timestep $= 100$')
	plot(time, exp_euler1000, label=r'Timestep $= 1000$')
	plot(time, exp_euler10000, label=r'Timestep $= 10000$')
	title(r'Explicit Euler scheme over different timesteps', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u(x, t)$', fontsize=18)
	legend(loc='best', fontsize=15)
	show()

implicit_euler = False
if implicit_euler:

	# plot all curves in same figure
	plot(time, imp_euler10, label=r'Timestep $= 10$')
	plot(time, imp_euler100, label=r'Timestep $= 100$')
	plot(time, imp_euler1000, label=r'Timestep $= 1000$')
	plot(time, imp_euler10000, label=r'Timestep $= 10000$')
	title(r'Implicit Euler scheme over different timesteps', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u(x,t)$', fontsize=18)
	legend(loc='best', fontsize=15)
	show()

implicit_CN = False
if implicit_CN:

	# plot all curves in same figure
	plot(time, imp_CN10, label=r'Timestep $= 10$')
	plot(time, imp_CN100, label=r'Timestep $= 100$')
	plot(time, imp_CN1000, label=r'Timestep $= 1000$')
	plot(time, imp_CN10000, label=r'Timestep $= 10000$')
	title(r'Implicit Crank-Nicolson scheme over different timesteps', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u(x, t)$', fontsize=18)
	legend(loc='best', fontsize=15)
	show()

plot_analytical_solution = True
if plot_analytical_solution:

	fig = figure(figsize=(25,8), facecolor='white')
	ax = fig.add_subplot(122)

	subplot(121)
	plot(time, analytical_solution, 'm', label=r'Analytical solution')
	plot(time, exp_euler1000, 'r', label=r'Explicit Euler')
	plot(time, imp_euler1000, 'b', label=r'Implicit Euler')
	plot(time, imp_CN1000, 'g', label=r'Implicit Crank-Nicolson')

	title(r'Compare analytical solution to all three schemes at $t_1 = 1000$', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u(x, t)$', fontsize=18)
	legend(loc='best', fontsize=15)

	subplot(122)
	plot(time, analytical_solution, 'm', label=r'Analytical solution')
	plot(time, exp_euler10000, 'r', label=r'Explicit Euler')
	plot(time, imp_euler10000, 'b', label=r'Implicit Euler')
	plot(time, imp_CN10000, 'g', label=r'Implicit Crank-Nicolson')

	title(r'Compare analytical solution to all three schemes at $t_2 = 10000$', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u(x, t)$', fontsize=18)
	legend(loc='best', fontsize=15)


	# subplot inside figure
	inset_axes = inset_axes(ax, width='30%', height=1.0, loc=5)
	plot(time, analytical_solution)
	plot(time, exp_euler10000)
	plot(time, imp_euler10000)
	plot(time, imp_CN10000)
	ylim([0.5,0.51])
	xlim([0.5,0.51])

	xticks([0.5, 0.509],)
	yticks([0.5, 0.509])
	#tight_layout()

	show()

	# find values to put in table
	value_eeuler1 = exp_euler10000[10]
	value_ieuler1 = imp_euler10000[10]
	value_cn1 = imp_CN10000[10]
	value_as1 = analytical_solution[10]

	value_eeuler2 = exp_euler10000[50]
	value_ieuler2 = imp_euler10000[50]
	value_cn2 = imp_CN10000[50]
	value_as2 = analytical_solution[50]

	value_eeuler3 = exp_euler10000[90]
	value_ieuler3 = imp_euler10000[90]
	value_cn3 = imp_CN10000[90]
	value_as3 = analytical_solution[90]

	print 'Explicit euler1: %g\nImplicit euler1: %g\nImplicit CN1: %g\nAnalytical solution1: %g' % (value_eeuler1, value_ieuler1, value_cn1, value_as1)

	'''
	Explicit euler1: 0.0997185
	Implicit euler1: 0.102649
	Implicit CN1: 0.10052
	Analytical solution1: 0.10101
	'''

	print 'Explicit euler2: %g\nImplicit euler2: %g\nImplicit CN2: %g\nAnalytical solution2: %g' % (value_eeuler2, value_ieuler2, value_cn2, value_as2)

	'''
	Explicit euler value: 0.500912
	Implicit euler value: 0.502422
	Implicit CN value: 0.501323
	Analytical solution value: 0.505051
	'''

	print 'Explicit euler3: %g\nImplicit euler3: %g\nImplicit CN3: %g\nAnalytical solution3: %g' % (value_eeuler3, value_ieuler3, value_cn3, value_as3)

	'''
	Explicit euler3: 0.907925
	Implicit euler3: 0.908191
	Implicit CN3: 0.907997
	Analytical solution3: 0.909091
	'''








