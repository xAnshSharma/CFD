using Plots    #to generate plots and gif

lx,ly = 2,2	#Defining extent of domain
nx,ny = 401,401		#number of nodes in each direction in domain
c = 0.1		#CFL
u_amp,v_amp = 1,1	#disturbance amplitude
dx,dy = lx/(nx-1), ly/(ny-1)	
nt = 1000	#total number of time steps
nu = 1e-3	#kinematic viscosity
Fo = 0.1	#Fourier number
total_time = 0    #total simulated time

x,y = LinRange(0, lx, nx),LinRange(0, ly, ny)	#Defining axes
u,v,u_res = zeros(nx,ny),zeros(nx,ny),zeros(nx,ny)    #initialising velocity fields
#Gaussian pulse disturbance	
xc, yc = lx/2, ly/2	#defining centre
sigma = 0.1	#spread of gaussian pulse

#creating 2D gaussian pulse at centre of domain
@inbounds for j in 1:ny
@simd for i in 1:nx
    pulse = exp(-((x[i]-xc)^2 + (y[j]-yc)^2)/(2*sigma^2))

    u[i,j] = u_amp*pulse
    v[i,j] = v_amp*pulse
end
end

#function to calculate advective term in x
function smart_face_x(phi,i,j,sgn)
#Defining normalized variable according to flow direction
if sgn > 0
U = phi[i-1,j]		#upwind node
C = phi[i,j]		#central node (node under consideration)
D = phi[i+1,j]		#downwind node
else
U = phi[i+2,j]
C = phi[i+1,j]
D = phi[i,j]
end

denom = D - U
if abs(denom) < 1e-12	#to avoid divison by 0
return C
end
phi_c_t = (C - U)/denom	#normalized variable 

phi_c_t = clamp(phi_c_t,0.0,1.0)	#for robustness
#defining piecewise function of SMART
if phi_c_t <= (1/6)	
return U + (3*phi_c_t*denom)
elseif phi_c_t <= (5/6)
phi_f = (3/8) + (0.75*phi_c_t)
return U + phi_f*denom
else
return D
end 
end

#similarly normalized varible in y
function smart_face_y(phi,i,j,sgn)
if sgn > 0
U = phi[i,j-1]
C = phi[i,j]
D = phi[i,j+1]
else
U = phi[i,j+2]
C = phi[i,j+1]
D = phi[i,j]
end

denom = D - U
if abs(denom) < 1e-12
return C
end
phi_c_t = (C - U)/denom

phi_c_t = clamp(phi_c_t,0.0,1.0)
if phi_c_t <= (1/6)
return U + (3*phi_c_t*denom)
elseif phi_c_t <= (5/6)
phi_f = (3/8) + (0.75*phi_c_t)
return U + phi_f*denom
else
return D
end 
end
function convective(un,vn,dx,dy,i,j)

#averaging face velocities to decide flow direction at face
u_e_avg = 0.5*(un[i+1,j] + un[i,j])
u_w_avg = 0.5*(un[i,j] + un[i-1,j])
v_n_avg = 0.5*(vn[i,j+1] + vn[i,j])
v_s_avg = 0.5*(vn[i,j] + vn[i,j-1])

if u_e_avg >=0
sgn_e = 1
else
sgn_e = -1
end

if u_w_avg >=0
sgn_w = 1
else
sgn_w = -1
end

if v_n_avg >=0
sgn_n = 1
else
sgn_n = -1
end

if v_s_avg >=0
sgn_s = 1
else
sgn_s = -1
end

#calculating face velocity using SMART function
u_e = smart_face_x(un,i,j,sgn_e)
u_w = smart_face_x(un,i-1,j,sgn_w)
v_n = smart_face_y(vn,i,j,sgn_n)
v_s = smart_face_y(vn,i,j-1,sgn_s)

v_e = smart_face_x(vn,i,j,sgn_e)
v_w = smart_face_x(vn,i-1,j,sgn_w)
u_n = smart_face_y(un,i,j,sgn_n)
u_s = smart_face_y(un,i,j-1,sgn_s)

#Calculating fluxes for momentum eqn
F_e = u_e^2
F_w = u_w^2
G_n = u_n*v_n
G_s = u_s*v_s
G_e = v_e*u_e
G_w = v_w*u_w
H_n = v_n^2
H_s = v_s^2

conv_u = -((1/dx)*(F_e - F_w)) - ((1/dy)*(G_n - G_s))
conv_v = -((1/dx)*(G_e - G_w)) - ((1/dy)*(H_n - H_s))
return conv_u,conv_v
end

#function to calculate diffusion term in x (4th order)
function diffusion_u(un,vn,dx,dy,i,j,nu)
diff_u = ((-un[i+2,j] + 16*un[i+1,j] -30*un[i,j] + 16*un[i-1,j] -un[i-2,j])*(nu/(dx^2)) + (-un[i,j+2] + 16*un[i,j+1] -30*un[i,j] + 16*un[i,j-1] -un[i,j-2])*(nu/(dy^2)))*(1/12)
return diff_u
end

#function to calculate diffusion term in y (4th order)
function diffusion_v(un,vn,dx,dy,i,j,nu)
diff_v = ((-vn[i+2,j] + 16*vn[i+1,j] -30*vn[i,j] + 16*vn[i-1,j] -vn[i-2,j])*(nu/(dx^2)) + (-vn[i,j+2] + 16*vn[i,j+1] -30*vn[i,j] + 16*vn[i,j-1] -vn[i,j-2])*(nu/(dy^2)))*(1/12)
return diff_v
end

#Computing RHS of momentum eqn in x and similarly in y
function residual_u(un,vn,i,j,dx,dy,nu)
conv,ign = convective(un,vn,dx,dy,i,j)
diff = diffusion_u(un,vn,dx,dy,i,j,nu)
return conv + diff
end

function residual_v(un,vn,i,j,dx,dy,nu)
ign,conv = convective(un,vn,dx,dy,i,j)
diff = diffusion_v(un,vn,dx,dy,i,j,nu)
return conv + diff
end

#Defining Boundary conditions
function BC(u,v)

u[1:3,:] .= 0	#u near x = 0
u[end-2:end,:] .= 0	#u near x = Lx

u[:,1:3] .= 0	#u near y = 0
u[:,end-2:end] .= 0	#u near y = Ly

v[1:3,:] .= 0	#v near x = 0
v[end-2:end,:] .= 0	#u near v = Lx

v[:,1:3] .= 0	#v near y = 0
v[:,end-2:end] .= 0	#v near y = Ly

end

#initialising residual matrices
Ru = zeros(nx,ny)
Rv = zeros(nx,ny)
Ru_1 = zeros(nx,ny)
Rv_1 = zeros(nx,ny)
Ru_2 = zeros(nx,ny)
Rv_2 = zeros(nx,ny)
#initialising velocity vectors
u_1 = copy(u)
v_1 = copy(v)
u_2 = copy(u)
v_2 = copy(v)

#Main simulation
anim = @animate for n in 1:nt
un = copy(u)
vn = copy(v)

#Define maximum velocity in each direction
u_max = maximum(abs.(u))
v_max = maximum(abs.(v))

dt_conv = min(((c*dx)/u_max),((c*dy)/v_max))	#compute dt using CFL
dt_diff = min(((Fo*(dx^2))/nu),((Fo*(dy^2))/nu))	#compute dt using diffusion stability
dt = min(dt_conv,dt_diff)	#compute dt using both CFL and diffusion stability

#RK3 stage one
@inbounds for j in 4:ny-3
@simd for i in 4:nx-3
Ru[i,j] = residual_u(un,vn,i,j,dx,dy,nu)
Rv[i,j] = residual_v(un,vn,i,j,dx,dy,nu)

u_1[i,j] = un[i,j] + dt*Ru[i,j]
v_1[i,j] = vn[i,j] + dt*Rv[i,j]

end
end
BC(u_1,v_1)

#RK3 stage two
@inbounds for j in 4:ny-3
@simd for i in 4:nx-3
Ru_1[i,j] = residual_u(u_1,v_1,i,j,dx,dy,nu)
Rv_1[i,j] = residual_v(u_1,v_1,i,j,dx,dy,nu)

u_2[i,j] = (3/4)*un[i,j] + 0.25*(u_1[i,j] + (dt*Ru_1[i,j]))
v_2[i,j] = (3/4)*vn[i,j] + 0.25*(v_1[i,j] + (dt*Rv_1[i,j]))

end
end

BC(u_2,v_2)

#RK3 stage three
@inbounds for j in 4:ny-3
@simd for i in 4:nx-3
Ru_2[i,j] = residual_u(u_2,v_2,i,j,dx,dy,nu)
Rv_2[i,j] = residual_v(u_2,v_2,i,j,dx,dy,nu)

u[i,j] = (1/3)*un[i,j] + (2/3)*(u_2[i,j] + (dt*Ru_2[i,j]))
v[i,j] = (1/3)*vn[i,j] + (2/3)*(v_2[i,j] + (dt*Rv_2[i,j]))

end
end
BC(u,v)

global total_time += dt	#compute total simulated time
u_res .= sqrt.(u.^2 .+ v.^2)	#resultant velocity
heatmap(y,x,u_res, title = "Time step = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
u_max = maximum(abs.(u))
println("max u:", u_max,"  at --> n:  ", n)    #shows variation of u_max throughout simulation to check for blowups

end
println("total time:", total_time)
#Animation
gif(anim, "2D_convection_diffusion_SMART_RK3.gif", fps = 30)
