using Plots

l = 2
nx = 101
dx = l/(nx-1)
nt = 300
c = 0.4
total_time = 0

x = LinRange(0,l,nx)
centre = l/2  #centre of domain
width = 0.2   #spread of wave
amplitude = 1.0  #amplitude added to u = 1
u = zeros(nx) + amplitude .* exp.(-(x .- centre).^2 / (2 * width^2))

function convection(un,i,dx)
u_e_avg = 0.5*(un[i+1]+un[i])
u_w_avg = 0.5*(un[i]+un[i-1])

if u_e_avg >= 0 
F_e = un[i]
else
F_e = un[i+1]
end

if u_w_avg >= 0 
F_w = un[i-1]
else
F_w = un[i]
end

conv = -(F_e - F_w)
return conv
end

function BC(un)
u[end] = u[end-1]
u[1] = u[end]
end

anim = @animate for n in 1:nt
un = copy(u)
u_max = maximum(abs.(u))
dt = (c*dx)/u_max

@inbounds @simd for i in 2:nx-1
C = convection(un,i,dx)
u[i] = un[i] + c*C
end
global total_time += dt
BC(un)

plot(x,u, title = "Time = $n", xlabel = "X", ylabel = "velocity", ylims = (0,2))
println("max u:", maximum(abs.(u)),"  at --> n:  ", n)
end
println("total time:", total_time)
gif(anim,"1D_convection_upwind_linear.gif",fps = 30)
