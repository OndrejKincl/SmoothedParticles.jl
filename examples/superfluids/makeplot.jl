
module makeplot
using Plots
using CSV
using DataFrames
using LaTeXStrings
using Interpolations
using LinearAlgebra
using Plots.PlotMeasures

const fontsize = 20

function lp_average(a::Array{Float64}, p::Int = 2)::Float64
    out = 0.
    dx = 1/length(a)
    for x in a
        out += dx*abs(x)^p
    end
    return out^(1/p)
end

function waves(name = "waves_nonlinear")
    linf2_error = Float64[]
    Ns = Int64[]
    p1 = plot(
        xlabel = L"\hat{t}", ylabel = L"\epsilon_1",
        legend = :none,
        linewidth = 3,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    bottom_margin = 10mm,
    )
    for N in 1:1000
        path = "results/"*name*"/error_data"*string(N)*".csv"
        if ispath(path)
            csvfile = CSV.read(path, DataFrame)
            plot!(p1, csvfile.t, csvfile.error, label = latexstring("N = "*string(N)))
            push!(linf2_error, maximum(csvfile.error))
            push!(Ns, N)
        end
    end
    savefig(p1, "results/"*name*"/error_data.pdf")

    X = -log10.(Ns)
    Y = log10.(linf2_error)
    p2 = plot(
        X, Y,
        xlabel = L"\log_{10} (1/N)", ylabel = L"\log_{10} (\epsilon_1)",
        legend = :none,
        linewidth = 3,
        markershape = :diamond,
        ms = 5,
        aspect_ratio = 1,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    bottom_margin = 10mm,
    )
    #linear regression
    (a,b) = [X ones(length(X))]\Y
    println("l^inf,2 error regression: ",a)
    plot!(p2, X, a*X .+ b, line = :dash)
    savefig(p2, "results/waves/convergence.pdf")

    Ns = Int64[]
    energy_error = Float64[]
    for N in 1:1000
        path = "results/"*name*"/energy_data"*string(N)*".csv"
        if ispath(path)
            csvfile = CSV.read(path, DataFrame)
            push!(Ns, N)
            push!(energy_error, maximum(abs.(csvfile.total .- csvfile.total[1])))
        end
    end
    csvpath = "results/"*name*"/energy_data"*string(maximum(Ns))*".csv"
    csvfile = CSV.read(csvpath, DataFrame)
    p3 = plot(
        csvfile.t, 
        [csvfile.kinetic csvfile.bulk csvfile.heat csvfile.total],
        legend = :outertopright,
        label = ["kinetic" "compression" "heat" "total"],
        xlabel = L"\hat{t}",
        ylabel = L"\hat{E}",
        linewidth = 3,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    bottom_margin = 10mm,
    )
    savefig(p3, "results/"*name*"/energy_balance.pdf")

    p4 = plot(
        csvfile.t, csvfile.total .- csvfile.total[1],
        legend = :outertopright,
        label = :none,
        xlabel = L"\hat{t}",
        ylabel = L"\hat{E}(t)",
        linewidth = 3,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    )
    savefig(p4, "results/"*name*"/energy_error_over_time.pdf")
    X = -log10.(Ns)
    Y = log10.(energy_error)
    p5 = plot(
        X, Y,
        xlabel = L"\log_{10} (1/N)", ylabel = L"\log_{10}(\epsilon_2)",
        legend = :none,
        linewidth = 3,
        markershape = :diamond,
        ms = 5,
        aspect_ratio = 1,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    )
    #linear regression
    (a,b) = [X ones(length(X))]\Y
    println("energy error regression: ",a)
    plot!(p2, X, a*X .+ b, line = :dash)
    savefig(p5, "results/"*name*"/energy_error.pdf")
end


function ldc()
    for Re in [100, 400, 1000, 3200, 5000, 7500, 10000]
        path = string("results/ldc/Re",Re)
        if !ispath(path)
            continue
        end
        ref_x2vy = CSV.read("ref/ldc-x2vy.csv", DataFrame)
        ref_y2vx = CSV.read("ref/ldc-y2vx.csv", DataFrame)
        propertyname = Symbol("Re",Re)
        ref_vy = getproperty(ref_x2vy, propertyname)
        ref_vx = getproperty(ref_y2vx, propertyname)
        ref_x = ref_x2vy.x
        ref_y = ref_y2vx.y
        @show Re
        #remove 6th point for Re400
        
        if Re == 400
            deleteat!(ref_x,6)
            deleteat!(ref_vy,6)
        end
        
        es = Float64[]
        Ns = Int64[]
        for N in 1:1000
            if ispath(path*"/fluxes"*string(N)*".csv")
                @show N
                df = CSV.read(path*"/fluxes"*string(N)*".csv", DataFrame)
                push!(Ns, N)
                fun1 = interpolate(df.s, df.v1, LinearMonotonicInterpolation())
                fun2 = interpolate(df.s, df.v2, LinearMonotonicInterpolation())
                e1 = maximum(abs.(fun2.(ref_x) - ref_vy)[2:end-1])
                e2 = maximum(abs.(fun1.(ref_y) - ref_vx)[2:end-1])
                push!(es, max(e1,e2))
            end
        end
        if isempty(Ns)
            continue
        end
        if length(Ns) > 1
            X = -log10.(Ns) 
            Y = log10.(es)
            p1 = plot(
                X, Y,
                xlabel = L"\log_{10} (1/N)", ylabel = L"\log_{10} (\epsilon)",
                legend = :none,
                linewidth = 4,
                markershape = :diamond,
                ms = 4,
                aspect_ratio = 1,
                tickfontsize = fontsize,
		        labelfontsize = fontsize,
		        legendfontsize = fontsize,
                bottom_margin = 10mm,
            )
            #linear regression
            (a,b) = [X ones(length(X))]\Y
            println("x2vy error regression: ",a)
            #plot!(p1, X, a*X .+ b, line = :dash, linewidth = 4)
            savefig(p1, path*"/convergence.pdf")
        end

        N = maximum(Ns)
        df = CSV.read(path*"/fluxes"*string(N)*".csv", DataFrame)
        p3 = plot(
            df.s, df.v2,
            xlabel = L"x",
            ylabel = L"v_y",
            label = "SPH",
            linewidth = 4,
            legend = :bottomleft,
            color = :orange,
            tickfontsize = fontsize,
            labelfontsize = fontsize,
            legendfontsize = fontsize,
            aspect_ratio = 1,
            bottom_margin = 10mm,
        )
        scatter!(p3, ref_x, ref_vy, label = "REF", markershape = :diamond, ms = 5, color = :blue)
        savefig(p3, path*"/ldc-x2vy.pdf")    
        p4 = plot(
            df.v1, df.s,
            xlabel = L"v_x",
            ylabel = L"y",
            label = "SPH",
            linewidth = 4,
            legend = :bottomright,
            color = :orange,
            tickfontsize = fontsize,
            labelfontsize = fontsize,
            legendfontsize = fontsize,
            aspect_ratio = 1,
            bottom_margin = 10mm,
        )
        scatter!(p4, ref_vx, ref_y,  label = "REF", markershape = :diamond, ms = 5, color = :blue)
        savefig(p4, path*"/ldc-y2vx.pdf")
    end
end

function irreversible()
    csvfile = CSV.read("results/irreversible/energy_data.csv", DataFrame)
    p = plot(
        csvfile.t, 
        [csvfile.kinetic csvfile.bulk csvfile.heat csvfile.total],
        legend = :outertopright,
        label = ["kinetic" "compression" "heat" "total"],
        xlabel = L"\hat{t}",
        ylabel = L"\hat{E}",
        linewidth = 3,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    )
    savefig(p, "results/irreversible/energy_balance.pdf")
    p = plot(
        csvfile.t, csvfile.total .- csvfile.total[1],
        legend = :outertopright,
        label = :none,
        xlabel = L"\hat{t}",
        ylabel = L"\hat{E}(t)",
        linewidth = 3,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,    )
    savefig(p, "results/irreversible/energy_error.pdf")
end

function moving_average(a::Array, R::Integer)::Array
    b = zeros(length(a))
    for i in eachindex(a)
        n = 0
        for j in max(i - R, 1) : min(i + R, length(a))
            b[i] += a[j]
            n += 1
        end
        b[i] = b[i]/n
    end
    return b
end

function fountain(name = "fountain100", w_dot::Float64 = 100.0, R = 10)
    csvfile = CSV.read("results/"*name*"/fountain.csv", DataFrame)
    T = 1.65
    S = 335
    r_cap = 1.55e-3/2
    rho0 = 145.2352
    v_theory = w_dot/(2.0*T*S*rho0*r_cap)
    p1 = plot(
        csvfile.t, 
        moving_average(csvfile.T_cell .+ T, R),
        legend = :none,
        xlabel = L"t \; [\mathrm{s}]",
        ylabel = L"T \; [\mathrm{K}]",
        linewidth = 3,
        tickfontsize = fontsize,
        labelfontsize = fontsize,
        legendfontsize = fontsize,
        bottom_margin = 5mm,
    )
    plot!(0.231*ones(10), range(minimum(csvfile.T_cell) .+ T , maximum(csvfile.T_cell) .+ T, 10),
        style = :dash,
        linewidth = 3,
        tickfontsize = fontsize,
        labelfontsize = fontsize,
        legendfontsize = fontsize,
        bottom_margin = 5mm,
        )
    savefig(p1, "results/"*name*"/cell_temperature.pdf")
    p2 = plot(
        csvfile.t, 
        [moving_average(csvfile.v_jet, R) v_theory*ones(length(csvfile.v_jet))],
        label = ["SPH" "THEORY"],
        legend = :bottomright,
        style = [:solid :dash],
        xlabel = L"t \; [\mathrm{s}]",
        ylabel = L"v \; [\mathrm{m/s}]",
        linewidth = 3,
        ylims = (0.0, 1.0),
        tickfontsize = fontsize,
        labelfontsize = fontsize,
        legendfontsize = fontsize,
        bottom_margin = 5mm,
    )
    savefig(p2, "results/"*name*"/jet_speed.pdf")
    K = div(length(csvfile.v_jet), 2)
    v_computed = sum(csvfile.v_jet[end-K+1:end])/K
    @show v_computed
    @show v_theory
    println("efficiency = ", 100*v_computed/v_theory, "%")
    println("deviation = ", 100*(1.0 - v_computed/v_theory), "%")
end

function fountain_3d(name = "fountain_3d", w_dot::Float64 = 0.05, K = 10)
    csvfile = CSV.read("results/"*name*"/fountain.csv", DataFrame)
    t = range(0.0, 0.4, length(csvfile.v_jet)) 
    T = 1.65 # + sum(csvfile.T_cell[end-K+1:end])/K
    S = 335.0 # + 11.4*(T - 1.65)
    r_cap = 1.55e-3/2
    rho0 = 145.2352
    v_theory = w_dot/(pi*T*S*rho0*r_cap*r_cap)
    p1 = plot(
        t,
        csvfile.T_cell .+ T,
        legend = :none,
        xlabel = L"t \; [\mathrm{s}]",
        ylabel = L"T \; [\mathrm{K}]",
        linewidth = 3,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    )
    savefig(p1, "results/"*name*"/cell_temperature.pdf")
    p2 = plot(
        t, 
        [moving_average(csvfile.v_jet, 5) v_theory*ones(length(csvfile.v_jet))],
        label = ["SPH" "THEORY"],
        legend = :bottomright,
        style = [:solid :dash],
        xlabel = L"t \; [\mathrm{s}]",
        ylabel = L"v \; [\mathrm{m/s}]",
        linewidth = 3,
        ylims = (0.0, 1.0),
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    )
    savefig(p2, "results/"*name*"/jet_speed.pdf")
    v_computed = sum(csvfile.v_jet[end-K+1:end])/K
    @show v_computed
    @show v_theory
    println("efficiency = ", 100*v_computed/v_theory, "%")
end


function wave_midpoint()
    csvfile = CSV.read("results/waves_ctest_fine/midpoint_data_t5_A1.csv", DataFrame)
    p1 = plot(
        1000.0*csvfile.t, 
        csvfile.s_exact,
        label = "exact",
        xlabel = L"t \; [\mathrm{ms}]",
        ylabel = L"s - s_0 \; [\mathrm{J/(K*kg)}]",
        linewidth = 3,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    )
    scatter!(p1,
        1000.0*csvfile.t[1:2:end],
        csvfile.s_computed[1:2:end],
        label = "SPH",
        markershape = :diamond,
        legend = :bottomleft,
        ms = 5,
	tickfontsize = fontsize,
	labelfontsize = fontsize,
	legendfontsize = fontsize,
    )
    savefig(p1, "results/waves/midpoint_s.pdf")

    p2 = plot(
        1000.0*csvfile.t, 
        zeros(length(csvfile.p)),
        label = "exact",
        xlabel = L"t \; [\mathrm{ms}]",
        ylabel = L"p \; [Pa]",
        linewidth = 3,
        tickfontsize = fontsize,
        labelfontsize = fontsize,
        legendfontsize = fontsize,
    )
    scatter!(p2,
        1000.0*csvfile.t[1:2:end],
        csvfile.p[1:2:end],
        label = "SPH",
        markershape = :diamond,
        legend = :bottomleft,
        ms = 5,
        tickfontsize = fontsize,
        labelfontsize = fontsize,
        legendfontsize = fontsize,
    )
    savefig(p2, "results/waves/midpoint_p.pdf")
end

function wave_ctest()
    As = Float64[]
    errs = Float64[]
    for k in 1:10
        csvfile = CSV.read("results/waves_ctest/error_data_t5_A"*string(k)*".csv", DataFrame)
        push!(As, k/100)
        push!(errs, maximum(abs.(csvfile.error)))
    end
    errs_coarse = Float64[]
    for k in 1:10
        csvfile = CSV.read("results/waves_ctest_coarse/error_data_t5_A"*string(k)*".csv", DataFrame)
        push!(errs_coarse, maximum(abs.(csvfile.error)))
    end
    errs_fine = Float64[]
    for k in 1:10
        csvfile = CSV.read("results/waves_ctest_fine/error_data_t5_A"*string(k)*".csv", DataFrame)
        push!(errs_fine, maximum(abs.(csvfile.error)))
    end
    p1 = plot(
        As,
        [errs_coarse errs errs_fine],
        label = [L"N = 100" L"N = 200" L"N = 300"],
        xlabel = L"A",
        ylabel = L"\epsilon_1",
        ms = 5,
        markershape = :diamond,
        legend = :outertopright,
        lw = 4,
        xticks = [0.02, 0.06, 0.10],
        yticks = [0.01, 0.015, 0.02],
        ylims = (0.005, 0.025),
        tickfontsize = 18,
        labelfontsize = 20,
        legendfontsize = 20,
        bottom_margin = 10mm,
        top_margin = 5mm
    )
        #linear regression
    #(a,b) = [As ones(length(As))]\errs
    #println("x2vy error regression: ",a)
    #plot!(p1, As, a*As .+ b, line = :dash, lw = 3)
    savefig(p1, "results/waves_ctest/ctest.pdf")
end

function wave_etest()
    ts = Float64[]
    errs = Float64[]
    for k in 1:10
        csvfile = CSV.read("results/waves_etest/energy_data_t"*string(k)*"_A10.csv", DataFrame)
        push!(ts, k/100)
        push!(errs, abs(maximum(csvfile.total) - 1.0))
    end
    errs_linear = Float64[]
    for k in 1:10
        csvfile = CSV.read("results/waves_etest_linear/energy_data_t"*string(k)*"_A1.csv", DataFrame)
        #push!(ts, k/100)
        push!(errs_linear, abs(maximum(csvfile.total) - 1.0))
    end
    X = log10.(ts)
    Y = log10.(errs)
    p1 = plot(legend = :none)#, aspect_ratio = 1.0)
        #linear regression
    (a,b) = [X ones(length(X))]\Y
    println("x2vy error regression: ",a)
    plot!(p1, 
        log10.(ts),
        [log10.(errs) log10.(errs_linear)],
        xlabel = L"\log_{10} \left( \frac{1}{M} \right)",
        ylabel = L"\log_{10} ( \epsilon_2)",
        label = [L"A = 0.1" L"A = 0.01"],
        linewidth = 4,
        markershape = :diamond,
        ms = 4,
        yticks = collect(-7.0 : 0.6 : -4.0),
        xticks = [-1.8, -1.5, -1.2],
        legend = :outertopright,
        tickfontsize = 18,
        labelfontsize = 20,
        legendfontsize = 20,
        bottom_margin = 10mm,
        top_margin = 5mm
    )
    #plot!(p1, X,a*X .+ b, line = :dash)
    savefig(p1, "results/waves_etest/etest.pdf")
end


function fountain_viscosity_test(w_dot::Float64 = 80.0, K = 20)
    T = 1.652 #sum(csvfile.T_cell[end-K+1:end])/K
    S = 335.0 + 11.4*(T - 1.65)
    r_cap = 1.55e-3/2
    rho0 = 145.2352
    v_theory = w_dot/(2.0*T*S*rho0*r_cap)
    Res = [145,452,1447,4824]
    etas = Float64[]
    for Re in Res
        path = "results/viscosity_test/Re"*string(Re)
        csvfile = CSV.read(path*"/fountain.csv", DataFrame)
        v_jet = sum(csvfile.v_jet[end-K+1:end])/K
        push!(etas, v_jet/v_theory)
    end
    p1 = plot(
        Res,
        etas,
        xlabel = L"\mathrm{Re}",
        ylabel = L"\hat{v}",
        linewidth = 3,
        markershape = :diamond,
        ms = 5,
        xaxis = :log,
        legend = :none,
        xticks = [10^2.5, 10^3, 10^3.5, 10^4, 10^4.5],
        ylims = (0.8, 1.0),
        style = :dash,
        tickfontsize = 18,
        labelfontsize = fontsize,
        legendfontsize = fontsize,
        bottom_margin = 5mm
    )
    savefig(p1, "results/viscosity_test/graph.pdf")
end

end
