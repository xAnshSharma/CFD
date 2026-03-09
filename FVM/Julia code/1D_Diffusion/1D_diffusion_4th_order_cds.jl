
using Plots   #to include contour plots and animations

l = 2   #dimensions of domain
nx = 101   #number of cell nodes in domain
Fo = 0.4    #fourier number
u_amp = 1   #disturbance amplitude
dx = l/(nx-1)    #distance between cell nodes
nt = 200    #total number of time steps
nu = 1e-1   #kinematic viscosity
total_time = 0	#total simulated time

x = LinRange(0,l,nx)   #defining axes
u = zeros(nx)    #initialising velocity fields

xc = l/2
sigma = 0.1

@inbounds @simd for i in 1:nx
    pulse = exp(-((x[i]-xc)^2)/(2*sigma^2))

    u[i] = u_amp*pulse
end

#function to calculate diffusion term in x (4th order)
function diffusion_u(un,dx,i,nu)
diff_u = (-un[i+2] + 16*un[i+1] -30*un[i] + 16*un[i-1] -un[i-2])*(nu/(dx^2))*(1/12)
return diff_u
end


#Boundary conditions
function BC(u)

u[1:2] .= 0   #u at x = 0
u[end-1:end] .= 0   #u at x = Lx

end

#Solution
anim = @animate for n in 1:nt
un = copy(u)

dt = ((Fo*(dx^2))/nu)


@inbounds @simd for i in 3:nx-2
Du = diffusion_u(un,dx,i,nu)

u[i] = un[i] + dt*Du

end
BC(u)
global total_time += dt

#resultant velocityhat 
plot(x,u, title = "Time step = $n", xlabel = "X", ylabel = "velocity", ylims = (0,2))
u_max = maximum(abs.(u))
println("max u:", u_max,"  at --> n:  ", n)

end
println("total simulated time:", total_time)
#generating animation
gif(anim, "1D_diffusion_4th_order_cds.gif", fps = 30)
