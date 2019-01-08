using Plots, LinearAlgebra, Polynomials
gr()

const tau = 2pi
const t_last_plot = 5tau
const t_last_error = tau

function f(t::Float64, y::Array{Float64,1})
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
    dt = tau / 64
    t_ary = range(0, stop = t_last_plot, step = dt)
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

    p1 = plot(t_ary, e_r, legend = :none, xlabel = "time", ylabel = "error",
    title = "Time evolution of the error")
    p2 = scatter(v_c_ary[:, 1], v_c_ary[:, 2], xlabel = "v_x", ylabel = "v_y",
    title = "Trajectory of v_c", aspect_ratio = 1, size = (400, 450), label = "Euler")
    plot!(v_a_ary[:, 1], v_a_ary[:, 2], label = "analytical", legend = :topright)

    # plot(p1, p2)

    png(p1, "img/3-2-1-a.png")
    png(p2, "img/3-2-1-b.png")
end

@time main()
