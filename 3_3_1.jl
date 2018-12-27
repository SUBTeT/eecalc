using Plots, LinearAlgebra, Polynomials
gr()

const sigma = 10
const b = 8/3
const r = 28
const dt = 0.001


function f(t::Float64, y::Array{Float64})
    [sigma*(y[2] - y[1]), r*y[1]-y[2]-y[1]y[3], y[1]y[2]-b*y[3]]
end

function runge_kutta(t_ary, v_c_ary, dt::Float64)
    for i = 2:length(t_ary)
        v_c_ary[i, :] = next(t_ary, v_c_ary, i, dt)
    end
end

function next(t_ary, v_c_ary, i::Int, dt::Float64)
    k1 = f(0.0, v_c_ary[i-1, :])
    k2 = f(0.0, v_c_ary[i-1, :] + dt / 2 * k1)
    k3 = f(0.0, v_c_ary[i-1, :] + dt / 2 * k2)
    k4 = f(0.0, v_c_ary[i-1, :] + dt * k3)
    return v_c_ary[i-1, :] + dt * (k1 + 2*k2 + 2*k3 + k4)/6
end

function main()
    t_last = 10

    t_ary = range(0, stop=t_last, step=dt)
    y_ary = Matrix{Float64}(undef, length(t_ary), 3)

    # initial
    y_ary[1, :] = [1.0, 0.0, 20.0]

    runge_kutta(t_ary, y_ary, dt)

    p1 = plot(t_ary, y_ary[:, 1], legend=:none, xlabel="time", ylabel="y_1",
    title="Time evolution of y_1")
    p2 = plot(y_ary[:, 1], y_ary[:, 2], legend=:none, xlabel="y_1", ylabel="y_2",
    title="Trajectory of (y_1, y_2)")

    png(p1, "img/3-3-1-a.png")
    png(p2, "img/3-3-1-b.png")
    # plot(p1, p2)
    #=
    @time anim = @animate for i=1:50:length(y_ary[:, 1])
        plot(y_ary[1:i, 1], y_ary[1:i, 2], y_ary[1:i, 3], legend=:none, title="Lorenz attractor",
        xlabel="x", ylabel="y", zlabel="z",
        xlims=(-20, 20), ylims=(-25, 25), zlims=(0, 50))
    end
    gif(anim, "gif/lorentz.gif", fps = 25)
    =#
end

main()
