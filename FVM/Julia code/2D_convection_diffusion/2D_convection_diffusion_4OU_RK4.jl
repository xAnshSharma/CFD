using Plots	#to generate plots and gifs

lx,ly = 2,2	#domain length in each direction
nx,ny = 401,401	#number of nodes in each direction
c = 0.1	#CFL
u_amp,v_amp = 1,1	#disturbance amplitude
dx,dy = lx/(nx-1), ly/(ny-1)	
nt = 1000	#total number of time steps
nu = 1e-3	#kinematic viscosity
Fo = 0.1	#Fourier number
total_time = 0	#total simulated time

x,y = LinRange(0, lx, nx),LinRange(0, ly, ny)	#defining axes
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

#face velocity function using 4th order UPWIND
function face4_x(phi,i,j,vel)
if vel >= 0	#positive velocity
phi_f = (-phi[i-2,j] + 7*phi[i-1,j] + 7*phi[i,j] - phi[i+1,j])/12
else	#negative velocity
phi_f = (-phi[i+2,j] + 7*phi[i+1,j] + 7*phi[i,j] - phi[i-1,j])/12
end
return phi_f
end

function face4_y(phi,i,j,vel)	
if vel >= 0
phi_f = (-phi[i,j-2] + 7*phi[i,j-1] + 7*phi[i,j] - phi[i,j+1])/12
else
phi_f = (-phi[i,j+2] + 7*phi[i,j+1] + 7*phi[i,j] - phi[i,j-1])/12
end
return phi_f
end

#Defining function to compute the convective term
function convective(dx,dy,i,j,un,vn)

#averaging face velocities for direction
u_e_avg = 0.5*(un[i+1,j] + un[i,j])
u_w_avg = 0.5*(un[i,j] + un[i-1,j])
v_n_avg = 0.5*(vn[i,j+1] + vn[i,j])
v_s_avg = 0.5*(vn[i,j] + vn[i,j-1])

#computing face velocity using 4th order UPWIND
u_e = face4_x(un,i,j,u_e_avg)
u_w = face4_x(un,i-1,j,u_w_avg)
v_n = face4_y(vn,i,j,v_n_avg)
v_s = face4_y(vn,i,j-1,v_s_avg)

v_e = face4_x(vn,i,j,u_e_avg)
v_w = face4_x(vn,i-1,j,u_w_avg)
u_n = face4_y(un,i,j,v_n_avg)
u_s = face4_y(un,i,j-1,v_s_avg)

#computing convective fluxes
F_e = u_e_avg*u_e
F_w = u_w_avg*u_w

G_n = v_n_avg*u_n
G_s = v_s_avg*u_s

G_e = u_e_avg*v_e
G_w = u_w_avg*v_w

H_n = v_n_avg*v_n	
H_s = v_s_avg*v_s

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

#calculating RHS of mom eqn in x and y
function residual_u(un,vn,i,j,dx,dy,nu)
conv,ign = convective(dx,dy,i,j,un,vn)
diff = diffusion_u(un,vn,dx,dy,i,j,nu)
return conv + diff
end

function residual_v(un,vn,i,j,dx,dy,nu)
ign,conv = convective(dx,dy,i,j,un,vn)
diff = diffusion_v(un,vn,dx,dy,i,j,nu)
return conv + diff
end

#defining boundary conditions
function BC(u,v)

u[1:3,:] .= 0	#u near x = 0
u[end-2:end,:] .= 0	#u near x = L

u[:,1:3] .= 0	#u near y = 0
u[:,end-2:end] .= 0	#u near y = L

v[1:3,:] .= 0	#v near x = 0
v[end-2:end,:] .= 0	#v near x = L

v[:,1:3] .= 0	#v near y = 0
v[:,end-2:end] .= 0	#v near y = L

end

#initialising residual matrices
Ru = zeros(nx,ny)
Rv = zeros(nx,ny)
Ru_1 = zeros(nx,ny)
Rv_1 = zeros(nx,ny)
Ru_2 = zeros(nx,ny)
Rv_2 = zeros(nx,ny)
Ru_3 = zeros(nx,ny)
Rv_3 = zeros(nx,ny)

#intialising velocity vectors
u_1 = copy(u)
v_1 = copy(v)
u_2 = copy(u)
v_2 = copy(v)
u_3 = copy(u)
v_3 = copy(v)

file = open("simulation_data_3.csv", "w")
println(file, "time,KE,u_max")

#Main simulation
anim = @animate for n in 1:nt
un = copy(u)
vn = copy(v)
#computing maximum velocities
u_max = maximum(abs.(u))
v_max = maximum(abs.(v))

dt_conv = min(((c*dx)/u_max),((c*dy)/v_max))	#computing dt using CFL
dt_diff = min(((Fo*(dx^2))/nu),((Fo*(dy^2))/nu))	#computing dt using diffusive stability
dt = min(dt_conv,dt_diff)	#computing dt using both CFL and diffusive stability

#RK4 stage one
@inbounds for j in 4:ny-3
@simd for i in 4:nx-3
Ru[i,j] = residual_u(un,vn,i,j,dx,dy,nu)
Rv[i,j] = residual_v(un,vn,i,j,dx,dy,nu)

u_1[i,j] = un[i,j] + 0.5*dt*Ru[i,j]
v_1[i,j] = vn[i,j] + 0.5*dt*Rv[i,j]

end
end
BC(u_1,v_1)

#RK4 stage two
@inbounds for j in 4:ny-3
@simd for i in 4:nx-3
Ru_1[i,j] = residual_u(u_1,v_1,i,j,dx,dy,nu)
Rv_1[i,j] = residual_v(u_1,v_1,i,j,dx,dy,nu)

u_2[i,j] = un[i,j] + 0.5*dt*Ru_1[i,j]
v_2[i,j] = vn[i,j] + 0.5*dt*Rv_1[i,j]

end
end

BC(u_2,v_2)

#RK4 stage three
@inbounds for j in 4:ny-3
@simd for i in 4:nx-3
Ru_2[i,j] = residual_u(u_2,v_2,i,j,dx,dy,nu)
Rv_2[i,j] = residual_v(u_2,v_2,i,j,dx,dy,nu)

u_3[i,j] = un[i,j] + dt*Ru_2[i,j]
v_3[i,j] = vn[i,j] + dt*Rv_2[i,j]

end
end
BC(u_3,v_3)

#RK4 stage four
@inbounds for j in 4:ny-3
@simd for i in 4:nx-3
Ru_3[i,j] = residual_u(u_3,v_3,i,j,dx,dy,nu)
Rv_3[i,j] = residual_v(u_3,v_3,i,j,dx,dy,nu)

u[i,j] = un[i,j] + (dt/6)*(Ru[i,j] + 2*Ru_1[i,j] + 2*Ru_2[i,j] + Ru_3[i,j])
v[i,j] = vn[i,j] + (dt/6)*(Rv[i,j] + 2*Rv_1[i,j] + 2*Rv_2[i,j] + Rv_3[i,j])

end
end

global total_time += dt	#computing total simulated time
BC(u,v)
u_res .= sqrt.(u.^2 .+ v.^2)	#resultant velocity
heatmap(y,x,u_res, title = "Time = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
KE = 0.5 * sum(u.^2 + v.^2) * dx * dy
u_max = maximum(abs.(u))
println("KE:", KE,"  at --> n:  ", n)
println("max u:", u_max,"  at --> n:  ", n)	#shows variation of u_max throughout simulation to check for blowups
println(file, "$total_time,$KE,$u_max")

if n == 300
p = heatmap(y,x,u_res, title = "Time step = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
savefig(p,"300_Time_4.png")
elseif n == 600
p = heatmap(y,x,u_res, title = "Time step = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
savefig(p,"600_Time_4.png")
elseif n == 900
p = heatmap(y,x,u_res, title = "Time step = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
savefig(p,"900_Time_4.png")
else
end

end
println("total time:", total_time)
close(file)
#Animation
gif(anim, "2D_convection_diffusion_4O'U_RK4.gif", fps = 30)
