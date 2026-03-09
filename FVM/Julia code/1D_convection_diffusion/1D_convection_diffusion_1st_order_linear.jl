using Plots   #to include contour plots and animations

l = 2   #dimensions of domain
nx = 101   #number of cell nodes in domain
c = 0.4   #courant number
Fo = 0.5    #fourier number
u_amp = 1   #disturbance amplitude
dx = l/(nx-1)    #distance between cell nodes
nt = 200    #total number of time steps
nu = 1e-3   #kinematic viscosity
total_time = 0	#total simulated time
vel = 1		#pulse velocity

x = LinRange(0,l,nx)   #defining axes
u = zeros(nx)    #initialising velocity fields

xc = l/2
sigma = 0.1

@inbounds @simd for i in 1:nx
    pulse = exp(-((x[i]-xc)^2 )/(2*sigma^2))

    u[i] = u_amp*pulse
end

#function to calculate advective term in x
function convective(un,dx,i)

u_e = 0.5*(un[i+1] + un[i])
u_w = 0.5*(un[i] + un[i-1])

if u_e >=0
F_e = un[i]
else
F_e = un[i+1]
end

if u_w >=0
F_w = un[i-1]
else
F_w = un[i]
end

conv_u = (-1/dx)*(F_e - F_w)
return conv_u
end

#function to calculate diffusion term
function diffusion_u(un,dx,i,nu)
diff_u = (un[i+1] + un[i-1] -2*un[i])*(nu/(dx^2))
return diff_u
end

#Boundary conditions
function BC(u)

u[end] = u[end-1]   #u at x = Lx
u[1] = u[end]   #u at x = 0


end

#Solution
anim = @animate for n in 1:nt
un = copy(u)

u_max = maximum(abs.(u))
dt = min(((c*dx)/u_max),((Fo*(dx^2))/nu))

@inbounds @simd for i in 2:nx-1
Cu = convective(un,dx,i)
Du = diffusion_u(un,dx,i,nu)

u[i] = un[i] + (vel*dt)*(Cu + Du)

end
BC(u)
global total_time += dt

#resultant velocityhat 
plot(x,u, title = "Time step = $n", xlabel = "X", ylabel = "velocity", ylims = (0,2))
u_max = maximum(abs.(u))
println("max u:", u_max,"  at --> n:  ", n)

end
#generating animation
gif(anim, "1D_convection_diffusion_1st_order_linear.gif", fps = 30)
