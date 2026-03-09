using Plots

l = 2
nx = 101
dx = l/(nx-1)
nt = 100
c = 0.2
u_amp = 1
total_time = 0

x = LinRange(0,l,nx)
centre = l/2  #centre of domain
width = 0.2   #spread of wave
u = zeros(nx) + u_amp .* exp.(-(x .- centre).^2 / (2 * width^2))

function convective(un,i,dx)

F_e = un[i+1]^2
F_w = un[i-1]^2

return (F_e - F_w)*(-1/dx)
end

function BC(u)
u[end] = u[end-1]
u[1] = u[end]
end

un = zeros(nx)

anim = @animate for n in 1:nt
un .= u
u_max = maximum(abs.(u))
dt = (c*dx)/u_max
@inbounds @simd for i in 2:nx-1
C = convective(un,i,dx)
u[i] = 0.5*(un[i+1]+un[i-1]) + dt*C
end
global total_time += dt
BC(u)
plot(x,u, title = "Time = $n", xlabel = "X", ylabel = "velocity", ylims = (0,2))
println("max u:", maximum(abs.(u)),"  at --> n:  ", n)
end
println("total time:", total_time)
gif(anim,"1D_convection_Lax_Friedrich.gif",fps = 30)
