using Plots, LinearAlgebra, Polynomials
gr()
# if you want to make 3d plot gif, then use pyplot()

sigma = 10
b = 8 / 3
r = 28
dt = 0.001

function f(t::Float64, y::Array{Float64})
    [sigma * (y[2] - y[1]), r * y[1] - y[2] - y[1]y[3], y[1]y[2] - b * y[3]]
end

function runge_kutta(t_ary, v_c_ary, dt::Float64)
    for i = 2:length(t_ary)
        v_c_ary[i, :] = next(t_ary, v_c_ary, i, dt)
    end
end

function next(t_ary, v_c_ary, i::Int, dt::Float64)
    k1 = f(0.0, v_c_ary[i - 1, :])
    k2 = f(0.0, v_c_ary[i - 1, :] + dt / 2 * k1)
    k3 = f(0.0, v_c_ary[i - 1, :] + dt / 2 * k2)
    k4 = f(0.0, v_c_ary[i - 1, :] + dt * k3)
    return v_c_ary[i - 1, :] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
end

function savegif(y_ary_of_list::Array{Float64,3})
    anim = @animate for i = 1:80:length(y_ary_of_list[1, :, 1])
        plot(y_ary_of_list[1, 1:i, 1], y_ary_of_list[1, 1:i, 2], y_ary_of_list[1, 1:i, 3],
        legend = :none, title = "difference by the initial value",
        xlabel = "x", ylabel = "y", zlabel = "z",
        xlims = (-25, 25), ylims = (-30, 30), zlims = (0, 55),
        c = :red)
        plot!(y_ary_of_list[2, 1:i, 1], y_ary_of_list[2, 1:i, 2], y_ary_of_list[2, 1:i, 3], c = :green)
        plot!(y_ary_of_list[3, 1:i, 1], y_ary_of_list[3, 1:i, 2], y_ary_of_list[3, 1:i, 3], c = :blue)
    end

    gif(anim, "gif/dif_by_init.gif", fps = 15)
end

function main()
    t_last = 10
    q_list = [0, 10^(-2), 10^(-3)]

    t_ary = range(0, stop = t_last, step = dt)
    y_ary_of_list = Array{Float64,3}(undef, length(q_list), length(t_ary), 3)

    for i = 1:length(q_list)
        y_ary = @view y_ary_of_list[i, :, :]
        # initial
        y_ary[1, :] = [1.0 + q_list[i], 0.0, 20.0]

        runge_kutta(t_ary, y_ary, dt)
    end

    error_list = Array{Float64,2}(undef, 2, length(t_ary))
    for i = 1:2
        for j = 1:length(t_ary)
            error_list[i, j] = norm(y_ary_of_list[1, j, :] - y_ary_of_list[1 + i, j, :])
        end
    end

    p1 = plot(t_ary, error_list[1, :], legend = :none, xlabel = "time", ylabel = "error",
    title = "Time evolution of error", label = "q=10^(-2)") # 0, 10^(-2)
    plot!(t_ary, error_list[2, :], legend = :topleft, xlabel = "time", ylabel = "error",
    title = "Time evolution of error", label = "q=10^(-3)") # 0, 10^(-3)
    png(p1, "img/3-3-2.png")
    # plot(p1)

    # savegif(y_ary_of_list)
end

@time main()
