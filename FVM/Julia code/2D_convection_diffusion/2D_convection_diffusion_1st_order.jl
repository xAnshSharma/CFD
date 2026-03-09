# ======================================================================
# This solver solves 2D viscid burger's equation using FVM
# It uses 1st and 2nd order schemes to discretize various terms in the burger's equation
# Temporal term : 1st order Explicit Euler scheme
# Advective term: 1st Order UPWIND scheme
# DIffusion term: 2nd order Central Difference Scheme
# ======================================================================

using Plots   #to include contour plots and animations

lx,ly = 2,2   #dimensions of domain
nx,ny = 401,401   #number of cell nodes in domain
c = 0.1   #courant number
Fo = 0.1    #fourier number
u_amp,v_amp = 1,1   #disturbance magnitude
dx,dy = lx/(nx-1), ly/(ny-1)    #distance between cell nodes
nt = 1000    #total number of time steps
nu = 1e-3   #kinematic viscosity
total_time = 0  #total simulated time 

x,y = LinRange(0,lx,nx),LinRange(0,ly,ny)   #defining axes
u,v,u_res = zeros(nx,ny),zeros(nx,ny),zeros(nx,ny)    #initialising velocity fields

#Creating Gaussian Pulse disturbance
xc, yc = lx/2, ly/2  #defining centre 
sigma = 0.1  #spread of gaussian pulse

@inbounds for j in 1:ny  
@simd for i in 1:nx
    pulse = exp(-((x[i]-xc)^2 + (y[j]-yc)^2)/(2*sigma^2))

    u[i,j] = u_amp*pulse
    v[i,j] = v_amp*pulse
end
end

#function to calculate advective term in x
function convective(un,vn,dx,dy,i,j)

#computing face velocities using simple average
u_e = 0.5*(un[i+1,j] + un[i,j]) 
u_w = 0.5*(un[i,j] + un[i-1,j])
v_n = 0.5*(vn[i,j+1] + vn[i,j])
v_s = 0.5*(vn[i,j] + vn[i,j-1])

#Defining Flux at each face based on face velocity direction (UPWIND)
if u_e >=0
F_e = un[i,j]*un[i,j]
G_e = un[i,j]*vn[i,j]
else
F_e = un[i+1,j]*un[i+1,j]
G_e = un[i+1,j]*vn[i+1,j]
end

if u_w >=0
F_w = un[i-1,j]*un[i-1,j]
G_w = un[i-1,j]*vn[i-1,j]
else
F_w = un[i,j]*un[i,j]
G_w = un[i,j]*vn[i,j]
end

if v_n >=0
G_n = un[i,j]*vn[i,j]
H_n = vn[i,j]*vn[i,j]
else
G_n = un[i,j+1]*vn[i,j+1]
H_n = vn[i,j+1]*vn[i,j+1]
end

if v_s >=0
G_s = un[i,j-1]*vn[i,j-1]
H_s = vn[i,j-1]*vn[i,j-1]
else
G_s = un[i,j]*vn[i,j]
H_s = vn[i,j]*vn[i,j]
end

conv_u = -((1/dx)*(F_e - F_w)) - ((1/dy)*(G_n - G_s))
conv_v = -((1/dx)*(G_e - G_w)) - ((1/dy)*(H_n - H_s))
return conv_u,conv_v
end

#function to calculate diffusion term in x
function diffusion_u(un,vn,dx,dy,i,j,nu)
diff_u = (un[i+1,j] + un[i-1,j] -2*un[i,j])*(nu/(dx^2)) + (un[i,j+1] + un[i,j-1] -2*un[i,j])*(nu/(dy^2))
return diff_u
end

#function to calculate diffusion term in y
function diffusion_v(un,vn,dx,dy,i,j,nu)
diff_v = (vn[i+1,j] + vn[i-1,j] -2*vn[i,j])*(nu/(dx^2)) + (vn[i,j+1] + vn[i,j-1] -2*vn[i,j])*(nu/(dy^2))
return diff_v
end

#Boundary conditions
function BC(u,v)

u[1,:] .= 0   #u at x = 0
u[end,:] .= 0   #u at x = Lx
u[:,1] .= 0   #u at y = 0
u[:,end] .= 0   #u at y = Ly

v[1,:] .= 0   #v at x = 0
v[end,:] .= 0   #v at x = Lx
v[:,1] .= 0   #v at y = 0
v[end,:] .= 0   #v at y = Ly

end

#Solution
anim = @animate for n in 1:nt
un = copy(u)
vn = copy(v)

u_max = maximum(abs.(u))
v_max = maximum(abs.(v))

#Dynamic time stepping based on CFL and DIffusion stability to for stable solver
dt_conv = min(((c*dx)/u_max),((c*dy)/v_max))  #dt based on CFL
dt_diff = min(((Fo*(dx^2))/nu),((Fo*(dy^2))/nu))  #dt based on Diffusion stability criterion
dt = min(dt_conv,dt_diff)  #final dt selection based on both stability criterion

#spatial loop
@inbounds for j in 2:ny-1
@simd for i in 2:nx-1
Cu,Cv = convective(un,vn,dx,dy,i,j)
Du = diffusion_u(un,vn,dx,dy,i,j,nu)
Dv = diffusion_v(un,vn,dx,dy,i,j,nu)

u[i,j] = un[i,j] + dt*(Cu + Du)
v[i,j] = vn[i,j] + dt*(Cv + Dv)

end
end
BC(u,v)
global total_time += dt #to compute the total simulated time

#resultant velocity 
u_res .= sqrt.(u.^2 .+ v.^2)
heatmap(y,x,u_res, title = "Time step = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
u_max = maximum(abs.(u))
println("max u:", u_max,"  at --> n:  ", n)  #shows variation of u_max throughout simulation to check for blowups
end
println("total time:", total_time)
#generating animation
gif(anim, "2D_convection_diffusion.gif", fps = 30)
