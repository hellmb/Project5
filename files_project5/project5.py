from pylab import *
from numpy import zeros
from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.mplot3d import Axes3D

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

exp_euler10_unst = loadtxt('explicit_euler_unstable10.txt')
exp_euler100_unst = loadtxt('explicit_euler_unstable100.txt')

solution_10 = loadtxt('solution10.txt')
solution_1000 = loadtxt('solution1000.txt')
solution_10000 = loadtxt('solution10000.txt')

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


explicit_euler_unstable = False
if explicit_euler_unstable:

	plot(time, exp_euler10_unst, label=r'Explicit Euler')
	plot(time, analytical_solution, label=r'Analytical solution')
	title(r'Explicit Euler for unstable $\alpha$ at $10$ timesteps against analytical solution', fontsize=20)
	xlabel(r'$Time$', fontsize=18)
	ylabel(r'$u(x, t)$', fontsize=18)
	legend(loc='best', fontsize=15)
	show()

	plot(time, exp_euler100_unst, label=r'Explicit Euler')
	plot(time, analytical_solution, label=r'Analytical solution')
	title(r'Explicit Euler for unstable $\alpha$ at $100$ timesteps against analytical solution', fontsize=20)
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

plot_analytical_solution = False
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
	plot(time, analytical_solution, 'm')
	plot(time, exp_euler10000, 'r')
	plot(time, imp_euler10000, 'b')
	plot(time, imp_CN10000, 'g')
	ylim([0.5,0.51])
	xlim([0.5,0.51])

	xticks([0.5, 0.51],)
	yticks([0.5, 0.51])
	#tight_layout()

	show()

	# find values to put in table
	value_eeuler1 = exp_euler10000[10]
	value_ieuler1 = imp_euler10000[10]
	value_cn1 = imp_CN10000[10]
	value_as1 = analytical_solution[10]

	value_eeuler2 = exp_euler10000[25]
	value_ieuler2 = imp_euler10000[25]
	value_cn2 = imp_CN10000[25]
	value_as2 = analytical_solution[25]

	value_eeuler3 = exp_euler10000[50]
	value_ieuler3 = imp_euler10000[50]
	value_cn3 = imp_CN10000[50]
	value_as3 = analytical_solution[50]

	value_eeuler4 = exp_euler10000[75]
	value_ieuler4 = imp_euler10000[75]
	value_cn4 = imp_CN10000[75]
	value_as4 = analytical_solution[75]

	value_eeuler5 = exp_euler10000[90]
	value_ieuler5 = imp_euler10000[90]
	value_cn5 = imp_CN10000[90]
	value_as5 = analytical_solution[90]

	print 'Explicit euler1: %g\nImplicit euler1: %g\nImplicit CN1: %g\nAnalytical solution1: %g' % (value_eeuler1, value_ieuler1, value_cn1, value_as1)

	print 'Explicit euler2: %g\nImplicit euler2: %g\nImplicit CN2: %g\nAnalytical solution2: %g' % (value_eeuler2, value_ieuler2, value_cn2, value_as2)

	print 'Explicit euler3: %g\nImplicit euler3: %g\nImplicit CN3: %g\nAnalytical solution3: %g' % (value_eeuler3, value_ieuler3, value_cn3, value_as3)

	print 'Explicit euler4: %g\nImplicit euler4: %g\nImplicit CN4: %g\nAnalytical solution4: %g' % (value_eeuler4, value_ieuler4, value_cn4, value_as4)

	print 'Explicit euler5: %g\nImplicit euler5: %g\nImplicit CN5: %g\nAnalytical solution5: %g' % (value_eeuler5, value_ieuler5, value_cn5, value_as5)


plot_initial = False
if plot_initial:

	dim = solution_10.shape[0]

	x = linspace(0,1, dim)
	y = linspace(0,1, dim)

	X, Y = meshgrid(x,y)
	Z = sin( 2 * pi * X ) * sin( 2 * pi * Y )

	fig = figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X,Y,Z, cmap=cm.rainbow, rstride=5, cstride=5)
	ax.set_title(r'Initial condition $ u(x,t,0) = sin( 2 \pi x ) sin( 2 \pi y ) $', fontsize=20)
	ax.set_xlabel('x', fontsize=18)
	ax.set_ylabel('y', fontsize=18)
	ax.set_zlabel('z', fontsize=18)
	show()

plot_solution5 = False
if plot_solution5:

	dim = solution_10.shape[0]

	x = linspace(0,1, dim)
	y = linspace(0,1, dim)

	X, Y = meshgrid(x,y)

	fig = figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X,Y,solution_10, cmap=cm.rainbow, rstride=5, cstride=5)
	ax.set_title('Numerical result of the diffusion equation for $t = 10$', fontsize=20)
	ax.set_xlabel('x', fontsize=18)
	ax.set_ylabel('y', fontsize=18)
	ax.set_zlabel('z', fontsize=18)
	show()


plot_solution1000 = False
if plot_solution1000:

	dim = solution_1000.shape[0]

	x = linspace(0,1, dim)
	y = linspace(0,1, dim)

	X, Y = meshgrid(x,y)

	fig = figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X,Y,solution_1000, cmap=cm.rainbow, rstride=5, cstride=5)
	ax.set_title(r'Numerical result of the diffusion equation for $t = 1000$', fontsize=20)
	ax.set_xlabel('x', fontsize=18)
	ax.set_ylabel('y', fontsize=18)
	ax.set_zlabel('z', fontsize=18)
	show()

plot_solution10000 = False
if plot_solution10000:

	dim = solution_10000.shape[0]

	x = linspace(0,1, dim)
	y = linspace(0,1, dim)

	X, Y = meshgrid(x,y)

	fig = figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X,Y,solution_10000, cmap=cm.rainbow, rstride=5, cstride=5)
	ax.set_title(r'Numerical result of the diffusion equation for $t = 10000$', fontsize=20)
	ax.set_xlabel('x', fontsize=18)
	ax.set_ylabel('y', fontsize=18)
	ax.set_zlabel('z', fontsize=18)
	show()

plot_analytical = False
if plot_analytical:

	dim = solution_1000.shape[0]

	x = linspace(0,1, dim)
	y = linspace(0,1, dim)
	t = linspace(0,1, dim)

	X, Y = meshgrid(x,y)

	Z = sin(2 * pi * X) * sin(2 * pi * Y) * exp(-4 * pi**2 * 1000)

	fig = figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X,Y,Z, cmap=cm.rainbow, rstride=5, cstride=5)
	ax.set_title('Analytical solution to the two-dimensional diffusion equation', fontsize=20)
	ax.set_xlabel('x', fontsize=18)
	ax.set_ylabel('y', fontsize=18)
	ax.set_zlabel('z', fontsize=18)
	show()


















