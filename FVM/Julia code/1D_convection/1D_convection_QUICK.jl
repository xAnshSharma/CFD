using Plots

l = 2
nx = 201
dx = l/(nx-1)
nt = 200
c = 0.1
u_amp = 1

x = LinRange(0, l, nx)
#Gaussian Pulse
centre = l/2  #centre of domain
width = 0.2   #spread of wave
u = zeros(nx) + u_amp .* exp.(-(x .- centre).^2 / (2 * width^2))

function convective(un,dx,i)

u_e = 0.5*(un[i+1] + un[i])
u_w = 0.5*(un[i] + un[i-1])

if u_e >=0
F_e = (1/8)*(6*(un[i]*un[i]) + 3*(un[i+1]*un[i+1]) - (un[i-1]*un[i-1]))
else
F_e = (1/8)*(6*(un[i+1]*un[i+1]) + 3*(un[i]*un[i]) - (un[i+2]*un[i+2]))
end

if u_w >=0
F_w = (1/8)*(6*(un[i-1]*un[i-1]) + 3*(un[i]*un[i]) - (un[i-2]*un[i-2]))
else
F_w = (1/8)*(6*(un[i]*un[i]) + 3*(un[i-1]*un[i-1]) - (un[i+1]*un[i+1]))
end

conv_u = (-1/dx)*(F_e - F_w)
return conv_u
end

function BC(u)

u[end-1] = u[end-2]
u[end] = u[end-1]
u[1] = u[end]
u[2] = u[1]

end

un = zeros(nx)
anim = @animate for n in 1:nt
un .= u
u_max = maximum(abs.(u))
dt = ((c*dx)/u_max)
@inbounds @simd for i in 3:nx-2
Cu = convective(un,dx,i)
u[i] = un[i] + dt*Cu
end
BC(u)
plot(x,u, title = "Time = $n", xlabel = "X", ylabel = "velocity", ylims = (0,2))
println("max u:", maximum(abs.(u)),"  at --> n:  ", n)
end
gif(anim, "1D_convection_QUICK.gif", fps = 30)
