using Plots, LinearAlgebra, Polynomials
gr()

const tau = 2pi
const t_last_plot = 5tau
const t_last_error = tau

function f(t::Float64, y::Array{Float64})
    [y[2], -y[1]]
end

function euler(t_ary, v_c_ary::Matrix{Float64}, dt::Float64)
    for i = 2:length(t_ary)
        v_c_ary[i, :] = next(t_ary, v_c_ary, i, dt)
    end
end

function next(t_ary, v_c_ary::Matrix{Float64}, i::Int, dt::Float64)
    return v_c_ary[i - 1, :] + dt * f(0.0, v_c_ary[i - 1, :])
end

# main
function main()
    p_list = 3:18
    dt_list = [tau / 2^p for p in p_list]
    log2E_r = Float64[]

    @time for dt in dt_list
        t_ary = range(0, stop = t_last_error, step = dt)
        v_c_ary = Matrix{Float64}(undef, length(t_ary), 2)

        # make v_a_ary
        v_a_ary = Matrix{Float64}(undef, length(t_ary), 2)
        for i = 1:length(t_ary)
            v_a_ary[i,:] = [sin(t_ary[i]), cos(t_ary[i])]
        end

        # initial
        v_c_ary[1, :] = [0.0, 1.0]

        euler(t_ary, v_c_ary, dt)

        e_r = Array{Float64}(undef, length(t_ary))
        for i = 1:length(t_ary)
            e_r[i] = norm(v_c_ary[i,:] - v_a_ary[i,:])
        end
        push!(log2E_r, log2(maximum(e_r)))
    end

    p = scatter(p_list, log2E_r, legend = :none, xlabel = "p", ylabel = "log2(E_r)",
    title = "The maximum value of the error", xticks = [p for p in p_list if p % 3 == 0])
    png(p, "img/3-2-2.png")
    a, b = coeffs(polyfit(p_list, log2E_r, 1))
    println(b)
    # scatter(p_list, log2E_r, legend=:none, xlabel="p", ylabel="log2(E_r)",
    # title="The maximum value of the error",
    # xticks=[p for p in p_list if p % 3 == 0])
end

main()
