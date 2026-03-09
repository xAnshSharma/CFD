# ======================================================================
# This solver solves 2D viscid burger's equation using FVM
# It uses 2nd order schemes to discretize various terms in the burger's equation
# Temporal term : Two-stage Runge-Kutta method
# Advective term: 2nd order MUSCL with minmod flux limiter
# DIffusion term: 2nd order Central Difference Scheme
# ======================================================================

using Plots	#importing plots

lx,ly = 2,2	#defining extent of domain
nx,ny = 401,401		#defining nodes in each direction in domain
dx,dy = lx/(nx-1),ly/(ny-1)	#computing dx,dy
c = 0.1		#courant number
Fo = 0.1	#fourier number
nu = 1e-3	#kinematic viscosity
u_amp,v_amp = 1,1	#disturbance amplitude
nt = 1000	#total number of time steps
total_time = 0    #total simulated time 

x,y = LinRange(0,lx,nx),LinRange(0,ly,ny)   #defining axes
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

#Defining function for minmod flux limiter
function minmod(a,b)
if a*b <= 0.0	#extrema
return 0.0	#upwind
else 
return sign(a)*min(abs(a),abs(b))	#limited slope
end
end

#reconstructing face velocity using MUSCL with minmod flux limiter
function velocity_components(un,vn,i,j)
#East Face
#x - component
slope_x_ij = minmod(((un[i,j] - un[i-1,j])/dx),((un[i+1,j] - un[i,j])/dx))	#calculating slope inside cell i,j
u_e_L = un[i,j] + 0.5*dx*slope_x_ij	#computing velocity from left (i,j) of face e

slope_x_ip1j = minmod(((un[i+1,j] - un[i,j])/dx),((un[i+2,j] - un[i+1,j])/dx))
u_e_R = un[i+1,j] - 0.5*dx*slope_x_ip1j

#y - component
slope_x_ij = minmod(((vn[i,j] - vn[i-1,j])/dx),((vn[i+1,j] - vn[i,j])/dx))
v_e_L = vn[i,j] + 0.5*dx*slope_x_ij

slope_x_ip1j = minmod(((vn[i+1,j] - vn[i,j])/dx),((vn[i+2,j] - vn[i+1,j])/dx))
v_e_R = vn[i+1,j] - 0.5*dx*slope_x_ip1j

#Applying upwind at east face
if (u_e_L + u_e_R)/2 > 0
u_e = u_e_L
v_e = v_e_L
else
u_e = u_e_R
v_e = v_e_R
end

#West face
#x component
slope_x_im1j = minmod(((un[i-1,j] - un[i-2,j])/dx),((un[i,j] - un[i-1,j])/dx))
u_w_L = un[i-1,j] + 0.5*dx*slope_x_im1j

slope_x_ij = minmod(((un[i,j] - un[i-1,j])/dx),((un[i+1,j] - un[i,j])/dx))
u_w_R = un[i,j] - 0.5*dx*slope_x_ij

#y - component

slope_x_im1j = minmod(((vn[i-1,j] - vn[i-2,j])/dx),((vn[i,j] - vn[i-1,j])/dx))
v_w_L = vn[i-1,j] + 0.5*dx*slope_x_im1j

slope_x_ij = minmod(((vn[i,j] - vn[i-1,j])/dx),((vn[i+1,j] - vn[i,j])/dx))
v_w_R = vn[i,j] - 0.5*dx*slope_x_ij

#Applying upwind at west face
if (u_w_L + u_w_R)/2 > 0
v_w = v_w_L
u_w = u_w_L
else
v_w = v_w_R
u_w = u_w_R
end

#North Face
#y - component
slope_y_ij = minmod(((vn[i,j] - vn[i,j-1])/dy),((vn[i,j+1] - vn[i,j])/dy))
v_n_B = vn[i,j] + 0.5*dy*slope_y_ij

slope_y_ijp1 = minmod(((vn[i,j+1] - vn[i,j])/dy),((vn[i,j+2] - vn[i,j+1])/dy))
v_n_U = vn[i,j+1] - 0.5*dy*slope_y_ijp1

#x component
slope_y_ij = minmod(((un[i,j] - un[i,j-1])/dy),((un[i,j+1] - un[i,j])/dy))
u_n_B = un[i,j] + 0.5*dy*slope_y_ij

slope_y_ijp1 = minmod(((un[i,j+1] - un[i,j])/dy),((un[i,j+2] - un[i,j+1])/dy))
u_n_U = un[i,j+1] - 0.5*dy*slope_y_ijp1

#Applying upwind at north face
if (v_n_B + v_n_U)/2 > 0
v_n = v_n_B
u_n = u_n_B
else
v_n = v_n_U
u_n = u_n_U
end

#South face
#x component
slope_y_ijm1 = minmod(((un[i,j-1] - un[i,j-2])/dy),((un[i,j] - un[i,j-1])/dy))
u_s_B = un[i,j-1] + 0.5*dy*slope_y_ijm1

slope_y_ij = minmod(((un[i,j] - un[i,j-1])/dy),((un[i,j+1] - un[i,j])/dy))
u_s_U = un[i,j] - 0.5*dy*slope_y_ij

#y component
slope_y_ijm1 = minmod(((vn[i,j-1] - vn[i,j-2])/dy),((vn[i,j] - vn[i,j-1])/dy))
v_s_B = vn[i,j-1] + 0.5*dy*slope_y_ijm1

slope_y_ij = minmod(((vn[i,j] - vn[i,j-1])/dy),((vn[i,j+1] - vn[i,j])/dy))
v_s_U = vn[i,j] - 0.5*dy*slope_y_ij

#Applying upwind at south face
if (v_s_B + v_s_U)/2 > 0
v_s = v_s_B
u_s = u_s_B
else
v_s = v_s_U
u_s = u_s_U
end
return u_e,u_w,v_e,v_w,u_n,u_s,v_n,v_s
end

#computing flux to later use in momentum eqns in x direction
function conv_u(un,vn,dx,dy,i,j)
u_e,u_w,v_e,v_w,u_n,u_s,v_n,v_s = velocity_components(un,vn,i,j)

Fe = u_e^2
Fw = u_w^2
Gn = u_n*v_n
Gs = u_s*v_s

conv_u = (Fe - Fw)/dx + (Gn - Gs)/dy	#return convective flux term in x
return conv_u
end

#computing flux to later use in momentum eqns in y direction
function conv_v(un,vn,dx,dy,i,j)
u_e,u_w,v_e,v_w,u_n,u_s,v_n,v_s = velocity_components(un,vn,i,j)

Ge = u_e*v_e
Gw = u_w*v_w
Hn = v_n^2
Hs = v_s^2

conv_v = (Ge - Gw)/dx + (Hn - Hs)/dy	#return convective flux in y
return conv_v
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

#Defining boundary conditions
function BC(u,v)

u[1:2,:] .= 0	#u at x = 0,1
u[end-1,:] .= u[end-2,:]	#neumen condn for o/p
u[end,:] .= u[end-1,:]

u[:,1:2] .= 0	#u at y = 0,1
u[:,end] .= u[:,end-1]	#neumen condn for o/p
u[:,end-1] .= u[:,end-2]

v[1:2,:] .= 0	#v at x = 0,1
v[end-1,:] .= v[end-2,:]	#neumen condn for o/p	
v[end,:] .= v[end-1,:]

v[:,1:2] .= 0	#v at y = 0,1
v[:,end-1] .= v[:,end-2]	#neumen cond for o/p
v[:,end] .= v[:,end-1]

end

#Computing RHS of momentum eqn in x and similarly in y
function residual_u(un,vn,dx,dy,i,j)
conv = conv_u(un,vn,dx,dy,i,j)
diff = diffusion_u(un,vn,dx,dy,i,j,nu)
return -conv  + diff
end

function residual_v(un,vn,dx,dy,i,j)
conv = conv_v(un,vn,dx,dy,i,j)
diff = diffusion_v(un,vn,dx,dy,i,j,nu)
return -conv  + diff
end

#initialising residual matrices
Ru = zeros(nx,ny)
Rv = zeros(nx,ny)

#Main simulation
anim = @animate for n in 1:nt
#Define maximum velocity in each direction
u_max = maximum(abs.(u))
v_max = maximum(abs.(v))

#computing dt inside time loop for robust stability
dt_conv = min(((c*dx)/u_max),((c*dy)/v_max))		#compute dt using CFL
dt_diff = min(((Fo*(dx^2))/nu),((Fo*(dy^2))/nu))	#compute dt using Diffusive stability
dt = min(dt_conv,dt_diff)	#select dt using both CFL and diffusive stability
un = copy(u)
vn = copy(v)
u_star = copy(u)
v_star = copy(v)
@inbounds for j in 3:ny-2
@simd for i in 3:nx-2
Ru[i,j] = residual_u(un,vn,dx,dy,i,j)
Rv[i,j] = residual_v(un,vn,dx,dy,i,j)


#Predictor step for RK2
u_star[i,j] = un[i,j] + dt*Ru[i,j]
v_star[i,j] = vn[i,j] + dt*Rv[i,j]
end 
end

BC(u_star,v_star)

@inbounds for j in 3:ny-2
@simd for i in 3:nx-2
Ru_star = residual_u(u_star,v_star,dx,dy,i,j)
Rv_star = residual_v(u_star,v_star,dx,dy,i,j)

#Corrector step for RK2
u[i,j] = 0.5*(un[i,j] + u_star[i,j] + dt*Ru_star)
v[i,j] = 0.5*(vn[i,j] + v_star[i,j] + dt*Rv_star)

end
end
global total_time += dt    #to compute the total simulated time
BC(u,v)
u_res .= sqrt.(u.^2 .+ v.^2)	#Compute resultant velocity
heatmap(y,x,u_res, title = "Time step = $n",xlabel = "x position", ylabel = "y position", legend = false)
u_max = maximum(abs.(u))
println("max u:", u_max,"  at --> n:  ", n)    #shows variation of u_max throughout simulation to check for blowups

end
#Animation
gif(anim,"2D_Burger's_MUSCL_RK2.gif",fps = 30)
