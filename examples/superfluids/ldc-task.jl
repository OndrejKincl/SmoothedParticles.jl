include("ldc_fast.jl")
#ldc.main(336)

for n in 1:6
	N = round(Int64, 50*1.1^n)
	ldc.main(N)
end
