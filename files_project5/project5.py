from pylab import *
from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes

time = linspace(0, 1, 100)
x = linspace(0, 1, 100)
n = 100.0

explicit_euler = False
if explicit_euler:

	# load text-files
	exp_euler10 = loadtxt('explicit_euler10.txt')
	exp_euler100 = loadtxt('explicit_euler100.txt')
	exp_euler1000 = loadtxt('explicit_euler1000.txt')
	exp_euler10000 = loadtxt('explicit_euler10000.txt')

	# plot all curves in same figure
	plot(time, exp_euler10, label=r'$u(x, t_{10})$')
	plot(time, exp_euler100, label=r'$u(x, t_{100})$')
	plot(time, exp_euler1000, label=r'$u(x, t_{1000})$')
	plot(time, exp_euler10000, label=r'$u(x, t_{10000})$')
	title(r'Explicit Euler scheme over different timesteps', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u_{i, j+1}$', fontsize=18)
	legend(loc='best', fontsize=15)
	show()

implicit_euler = False
if implicit_euler:

	# load text-files
	imp_euler10 = loadtxt('implicit_euler10.txt')
	imp_euler100 = loadtxt('implicit_euler100.txt')
	imp_euler1000 = loadtxt('implicit_euler1000.txt')
	imp_euler10000 = loadtxt('implicit_euler10000.txt')

	# plot all curves in same figure
	plot(time, imp_euler10, label=r'$u(x, t_{10})$')
	plot(time, imp_euler100, label=r'$u(x, t_{100})$')
	plot(time, imp_euler1000, label=r'$u(x, t_{1000})$')
	plot(time, imp_euler10000, label=r'$u(x, t_{10000})$')
	title(r'Implicit Euler scheme over different timesteps', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u_{i, j}$', fontsize=18)
	legend(loc='best', fontsize=15)
	show()

implicit_CN = False
if implicit_CN:

	# load text-files
	imp_CN10 = loadtxt('implicit_CN10.txt')
	imp_CN100 = loadtxt('implicit_CN100.txt')
	imp_CN1000 = loadtxt('implicit_CN1000.txt')
	imp_CN10000 = loadtxt('implicit_CN10000.txt')

	# plot all curves in same figure
	plot(time, imp_CN10, label=r'$u(x, t_{10})$')
	plot(time, imp_CN100, label=r'$u(x, t_{100})$')
	plot(time, imp_CN1000, label=r'$u(x, t_{1000})$')
	plot(time, imp_CN10000, label=r'$u(x, t_{10000})$')
	title(r'Implicit Crank-Nicolson scheme over different timesteps', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u_{i, j}$', fontsize=18)
	legend(loc='best', fontsize=15)
	show()

# analytical solution
analytical_solution = x + (2.0/pi) * sum( (-1.0**n/n) * sin(n * pi * x) * exp( -(n * pi)**2 * time) )

plot_analytical_solution = True
if plot_analytical_solution:

	exp_euler10000 = loadtxt('explicit_euler10000.txt')
	imp_euler10000 = loadtxt('implicit_euler10000.txt')
	imp_CN10000 = loadtxt('implicit_CN10000.txt')

	fig = figure(figsize=(25,8), facecolor='white')
	ax = fig.add_subplot(121)
	plot(time, analytical_solution, label=r'Analytical solution')
	plot(time, exp_euler10000, label=r'Explicit Euler')
	plot(time, imp_euler10000, label=r'Implicit Euler')
	plot(time, imp_CN10000, label=r'Implicit Crank-Nicolson')
	title(r'Compare analytical solution to all three schemes', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u_{i, j}$', fontsize=18)
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
	value1 = exp_euler10000[50]
	value2 = imp_euler10000[50]
	value3 = imp_CN10000[50]
	value4 = analytical_solution[50]
	print 'Explicit euler value: %g\nImplicit euler value: %g\nImplicit CN value: %g\nAnalytical solution value: %g' % (value1, value2, value3, value4)

	'''
	Explicit euler value: 0.500912
	Implicit euler value: 0.502422
	Implicit CN value: 0.501323
	Analytical solution value: 0.505051
	'''










