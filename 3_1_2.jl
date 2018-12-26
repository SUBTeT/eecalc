@time using Reactive, Interact, Plots
gr()

N = 100
dx = 1.0 / N
dt = 0.001
c = 5
T_D = 0.1
mu = c^2 * dt^2 / dx^2

limit = trunc(Int, T_D / dt * 4)  # period

x = range(0, stop = 1 - dx, step = dx)
t = range(0, stop = limit, step = dt)

function calc_next(u, i)
    u_k_next = zeros(Float64, N)

    # fixed end
    u_k_next[N] = -u[i-2, N] + mu * u[i-1, N - 1] + 2 * (1 - mu) * u[i-1, N]
    # free end
    # u_k_next[N - 1] = -u[-2][N - 1] + mu * u[-1][N - 2] + (2 - mu) * u[-1][N - 1]

    if (i-1) * dt < T_D
        u_k_next[1] = 1.0
    end

    for j = 2:N - 1
        u_k_next[j] = -u[i-2, j] + mu * u[i-1, j - 1] + 2 * (1 - mu) * u[i-1, j] + mu * u[i-1, j + 1]
    end

    return u_k_next
end

function make_u(u)
    for i = 3:limit
        u_k_next = calc_next(u, i)
        u[i,:] = u_k_next
    end
    return u
end

function main()
    u = Array{Float64, 2}(undef, limit, N)
    for i = 1:2
        u[i, 1] = 1.0
        for j in 2:N
            u[i, j] = 0.0
        end
    end
    u = make_u(u)

    @time anim = @animate for i in 1:2:limit
        plot(u[i, :], legend=:none, xlabel="x", ylabel="u",
        title="Time evolution of u")
        ylims!(-1.3, 1.3)
    end
    gif(anim, "gif/3-1-2.gif", fps = 30)


end

main()
