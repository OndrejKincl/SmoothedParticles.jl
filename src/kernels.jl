#all kernel functions are normalized such that their support radius is h

@fastmath function pos(x::Float64)::Float64
	return (x > 0.0 ? x : 0.0)
end

"""
    spline23(h::Float64, r::Float64)::Float64

Returns ``w(r)``, the value of a 2d cubic spline ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function spline23(h::Float64, r::Float64)::Float64
	x = r/h
	if (x < 0.5)
		# 40/7pi = 1.8189136353359467
		return 1.8189136353359467*(1.0 - 6.0*x^2 + 6.0*x^3)/h^2
	elseif (x < 1.0)
		# 80/7pi = 3.6378272706718935
		return 3.6378272706718935*(1.0 - x)^3/h^2
	end
	return 0.0
end

"""
    Dspline23(h::Float64, r::Float64)::Float64

Returns ``\\frac{\\text{d}w}{\\text{d}r}(r)``, the first derivative of a 2d cubic spline ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function Dspline23(h::Float64, r::Float64)::Float64
	x = r/h
	if (x < 0.5)
		# 240/7pi = 10.91348181201568
		return -10.91348181201568*(2.0*x - 3.0*x^2)/h^3
	elseif (x < 1.0)
		return -10.91348181201568*(1.0 - x)^2/h^3
	end
	return 0.0
end

"""
    rDspline23(h::Float64, r::Float64)::Float64

Returns ``\\frac{1}{r}\\,\\frac{\\text{d}w}{\\text{d}r}(r)``, the reduced first derivative of a 2d cubic spline ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function rDspline23(h::Float64, r::Float64)::Float64
	x = r/h
	if (x < 0.5)
		# 240/7pi = 10.91348181201568
		return -10.91348181201568*(2.0 - 3.0*x)/h^4
	elseif (x < 1.0)
		return -10.91348181201568*(1.0 - x)^2/(x*h^4)
	end
	return 0.0
end

"""
    spline24(h::Float64, r::Float64)::Float64

Returns ``w(r)``, the value of a 2d quartic spline ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function spline24(h::Float64, r::Float64)::Float64
	x = r/h
	return 6.222175110452539*(pos(1.0 - x)^4 - 5*pos(0.6 - x)^4 + 10*pos(0.2 - x)^4)/h^2
end

"""
    Dspline24(h::Float64, r::Float64)::Float64

Returns ``\\frac{\\text{d}w}{\\text{d}r}(r)``, the first derivative of a 2d quartic spline ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function Dspline24(h::Float64, r::Float64)::Float64
	x = r/h
	return -24.888700441810155*(pos(1.0 - x)^3 - 5*pos(0.6 - x)^3 + 10*pos(0.2 - x)^3)/h^3
end

"""
    rDspline24(h::Float64, r::Float64)::Float64

Returns ``\\frac{1}{r}\\,\\frac{\\text{d}w}{\\text{d}r}(r)``, the reduced first derivative of a 2d quartic spline ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function rDspline24(h::Float64, r::Float64)::Float64
	x = r/h
	if x > 0.2
		return -24.888700441810155*(pos(1.0 - x)^3 - 5*pos(0.6 - x)^3)/(x*h^4)
	end
	return -24.888700441810155*(1.2 - 6.0*x^2)/h^4
end

"""
    wendland2(h::Float64, r::Float64)::Float64

Returns ``w(r)``, the value of a 2d quintic Wendland kernel ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function wendland2(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	# 7/pi ≐ 2.228169203286535
	return 2.228169203286535*((1.0 - x)^4)*(1.0 + 4.0*x)/h^2
end

"""
    Dwendland2(h::Float64, r::Float64)::Float64

Returns ``\\frac{\\text{d}w}{\\text{d}r}(r)``, the first derivative of a 2d quintic Wendland kernel ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function Dwendland2(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	# 140/pi ≐ 44.563384065730695
	return -44.563384065730695*x*((1.0 - x)^3)/h^3
end

"""
    rDwendland2(h::Float64, r::Float64)::Float64

Returns ``\\frac{1}{r}\\,\\frac{\\text{d}w}{\\text{d}r}(r)``, the reduced first derivative of a 2d quintic Wendland kernel ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function rDwendland2(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	# 140/pi ≐ 44.563384065730695
	return -44.563384065730695*((1.0 - x)^3)/h^4
end

"""
    wendland3(h::Float64, r::Float64)::Float64

Returns ``w(r)``, the value of a 3d quintic Wendland kernel ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function wendland3(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	# 21/2pi ≐ 3.3422538049298023
	return 3.3422538049298023*((1.0 - x)^4)*(1.0 + 4.0*x)/h^3
end

"""
    Dwendland3(h::Float64, r::Float64)::Float64

Returns ``\\frac{\\text{d}w}{\\text{d}r}(r)``, the first derivative of a 3d quintic Wendland kernel ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function Dwendland3(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	# 210/pi ≐ 66.84507609859604
	return -66.84507609859604*x*((1.0 - x)^3)/h^4
end

"""
    rDwendland3(h::Float64, r::Float64)::Float64

Returns ``\\frac{1}{r}\\,\\frac{\\text{d}w}{\\text{d}r}(r)``, the reduced first derivative of a 3d quintic Wendland kernel ``w`` with support radius `h`.
Integrates to unity.

"""
@fastmath function rDwendland3(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	# 210/pi ≐ 66.84507609859604
	return -66.84507609859604*((1.0 - x)^3)/h^5
end

@fastmath function DDwendland3(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	# 630/pi ≐ 200.53522829578813
	return 200.53522829578813*((1.0 - x)^2)/h^5
end