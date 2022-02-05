using Plots: plot, savefig
using Ipopt
using JuMP

"""
	Histogram(xs::Vector{Float64}, ys::Vector{Float64},N,dx)

Histogram structure storing ``N`` x values ``xs`` with uniform bin width ``dx`` and ``N`` y values ``ys``.
"""
struct Histogram
	xs::Vector{Float64}
	ys::Vector{Float64}
	N::Int64
	dx::Float64
end

"""
	velocity_histogram(sys::ParticleSystem; v_max = 0, N = 10)

Building the histrogram of 2D velocities (norms) with ``v_max`` the maximum velocity in the histogram and ``N`` bins.
"""
function velocity_histogram(sys::ParticleSystem; v_max = 0.0, N = 100)::Histogram
	if v_max == 0.0 # if v_max = 0, find the maximum velocity of the particles
		for k in 1:length(sys.particles) 
			v = norm(sys.particles[k].v)
			if v > v_max
				v_max = v
			end
		end
	end

	# Find the heights of the histogram bins
	dv = v_max/N # velocity increment between the bins
	vs = 0.:dv:v_max
	ns = zeros(length(vs))
	for k in 1:length(sys.particles)
		v = norm(sys.particles[k].v)
		n = Int64(round(v/dv))
		if 1 <= n <= length(ns)
			ns[n] += 1.0/(dv*length(sys.particles))
		end
	end

	return Histogram(vs,ns,100,dv)
end



"""
	kB

Boltzmann constant (in the SI units)
"""
const kB = 1.380649e-23

"""
	entropy(fMB::Histogram)::Float64

Calculate Boltzmann entropy of a 2D Maxwell-Boltzmann distribution approximated by an ``fMB`` histogram.
"""
function entropy_2D_MB(fMB::Histogram)::Float64
	@assert(fMB.xs[1] == 0) # Assuming that the histogram starts at zero velocity

	S = 0.0

	# Approximating the reduced entropy near v=0, where a numerical singularity could appear
	fMBder = (fMB.ys[2]-fMB.ys[1])/fMB.dx
	if fMBder > 0
		S = - fMB.ys[1] * (log(fMBder)*fMB.dx - fMBder*(fMB.dx^3)/6)
	end

	# Approximating the rest of entropy
	for k in 2:length(fMB.xs)
		if fMB.xs[k] != 0
			if fMB.ys[k] > 0
				S += -fMB.ys[k] * log(fMB.ys[k]/fMB.xs[k]) * fMB.dx
			end
		end
	end

	return S
end

""" 
	plot_velocity_distr(sys::ParticleSystem, name::String)

Plots the distribution of velocity magnitudes among the particles of ``sys`` and saves the resulting pdf,
together with a fit of the Maxwell-Boltzmann distribution, to file ``name``. ``m`` is the mass of each particle.
Returns the temperature.
"""
function plot_velocity_distr(sys::ParticleSystem, m::Float64, name::String; v_max = 0.0)::Float64
	distr = velocity_histogram(sys, v_max=v_max, N = 100)
	
	# fitting the histogram to a 2D Maxwell-Boltzmann distribution
	model = Model(Ipopt.Optimizer)
	@variable(model, beta)
	@NLobjective(
			model,
			Min,
			sum((distr.ys[i] - m*beta*distr.xs[i]*exp(-0.5*m*beta*distr.xs[i]^2))^2 for i in 1:length(distr.xs)),
		) 
	optimize!(model)
	beta = value(beta)

	# Plotting both the actual histogram and the fitted Maxwell-Boltzmann distribution
	ns_boltz = zeros(length(distr.xs))
	for i in 1:length(ns_boltz)
		ns_boltz[i] = m*beta*distr.xs[i]*exp(-0.5*m*beta*distr.xs[i]^2)
	end
	@show(beta)
	T = 1/(beta*kB)
	@show(T)

	p = plot(distr.xs, [distr.ys ns_boltz], label = ["data" "Maxwell-Boltzmann, T="*string(T)])
	savefig(p, name)
	return T
end

#"""
#	plot_velocity_distr(path::String, name:: String)
#
#Plot velocity distribution of particles loaded from a VTK file.
#"""
#function plot_velocity_distr(path::String, name::String)
#	@error "Not working, sorry"
#	domain = Rectangle(-box_width, -box_width, 2*box_width, 3*box_height) 
#	sys = ParticleSystem(Particle, domain, h)
#	read_vtk!(sys, path, x -> Particle(x = x, type = 0.0))
#	plot_velocity_distr(sys, name)  #not yet working
#end
