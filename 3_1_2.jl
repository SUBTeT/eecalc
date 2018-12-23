#=
3_1_2:
- Julia version: 
- Author: Tetsushi
- Date: 2018-12-20
=#
@time using Reactive, Interact, Plots

clibrary(:colorcet)
gr()
ENV["PLOTS_TEST"] = "true"

N = 100
dx = 1.0 / N
dt = 0.001
c = 5
T_D = 0.1
mu = c^2 * dt^2 / dx^2

limit = T_D / dt * 4  # period

x = range(0, stop = 1 - dx, step = dx)
t = range(0, stop = limit, step = dt)

function calc_next(u)
    u_k_next = zeros(Float64, N)

    last = length(u[:, 1])
    # fixed end
    u_k_next[N] = -u[last - 1][N] + mu * u[last][N - 1] + 2 * (1 - mu) * u[last][N]
    # free end
    # u_k_next[N - 1] = -u[-2][N - 1] + mu * u[-1][N - 2] + (2 - mu) * u[-1][N - 1]

    if last * dt < T_D
        u_k_next[1] = 1.
    end

    for i = 2:N - 1
        u_k_next[i] = -u[last - 1][i] + mu * u[last][i - 1] + 2 * (1 - mu) * u[last][i] + mu * u[last][i + 1]
    end

    return u_k_next
end

function make_u(u)
    for i = 1:limit
        u_k_next = calc_next(u)
        push!(u, u_k_next)
    end
    return u
end

# initial
u = []
a = zeros(Float64, N)
a[1] = 1.0
push!(u, a)
push!(u, a)

@time u = make_u(u)

@time anim = @animate for i in 1:5:length(u)
    plot(u[i, :])
    ylims!(-1.3, 1.3)
end

gif(anim, "gif/3-1-2.gif", fps = 10)
