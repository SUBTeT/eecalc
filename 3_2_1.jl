using Plots, LinearAlgebra, Polynomials
gr()

const tau = 2pi
const t_last_plot = 5tau
const t_last_error = tau

function f(t, y)
    [y[2], -y[1]]
end

function euler(t_ary, v_c_ary, dt)
    for i = 2:length(t_ary)
        v_c_ary[i, :] = next(t_ary, v_c_ary, i, dt)
    end
end

function next(t_ary, v_c_ary, i, dt)
    return v_c_ary[i-1, :] + dt * f(undef, v_c_ary[i-1, :])
end

# main
function main()
    dt = tau/64
    t_ary = range(0, stop=t_last_plot, step=dt)
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

    p1 = plot(t_ary, e_r, legend=:none, xlabel="time", ylabel="error",
    title="Time evolution of the error", aspect_ratio=10)
    p2 = scatter(v_c_ary[:, 1], v_c_ary[:, 2], legend=:none, xlabel="v_x", ylabel="v_y",
    title="Trajectory of v_c", aspect_ratio=1)

    p = plot(p1, p2)
    png(p, "img/3-2-1.png")
end

main()
