using Plots

lx,ly = 2,2
nx,ny = 201,201
c = 0.1
u_amp,v_amp = 1,1
dx,dy = lx/(nx-1), ly/(ny-1)
nt = 200
vel = 1
total_time = 0

x,y = LinRange(0, lx, nx),LinRange(0, ly, ny)
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

function convective(un,vn,dx,dy,i,j)

u_e = 0.5*(un[i+1,j] + un[i,j])
u_w = 0.5*(un[i,j] + un[i-1,j])
v_n = 0.5*(vn[i,j+1] + vn[i,j])
v_s = 0.5*(vn[i,j] + vn[i,j-1])

if u_e >=0
F_e = (1/8)*(6*(un[i,j]) + 3*(un[i+1,j]) - un[i-1,j])
G_e = (1/8)*(6*(vn[i,j]) + 3*(vn[i+1,j]) - vn[i-1,j])
else
F_e = (1/8)*(6*(un[i+1,j]) + 3*(un[i,j]) - un[i+2,j])
G_e = (1/8)*(6*(vn[i+1,j]) + 3*(vn[i,j]) - vn[i+2,j])
end

if u_w >=0
F_w = (1/8)*(6*(un[i-1,j]) + 3*(un[i,j]) - un[i-2,j])
G_w = (1/8)*(6*(vn[i-1,j]) + 3*(vn[i,j]) - vn[i-2,j])
else
F_w = (1/8)*(6*(un[i,j]) + 3*(un[i-1,j]) - un[i+1,j])
G_w = (1/8)*(6*(vn[i,j]) + 3*(vn[i-1,j]) - vn[i+1,j])
end

if v_n >=0
G_n = (1/8)*(6*(vn[i,j]) + 3*(vn[i,j+1]) - vn[i,j-1])
H_n = (1/8)*(6*(vn[i,j]) + 3*(vn[i,j+1]) - vn[i,j-1])
else
G_n = (1/8)*(6*(vn[i,j+1]) + 3*(vn[i,j]) - vn[i,j+2])
H_n = (1/8)*(6*(vn[i,j+1]) + 3*(vn[i,j]) - vn[i,j+2])
end

if v_s >=0
G_s = (1/8)*(6*(vn[i,j-1]) + 3*(vn[i,j]) - vn[i,j-2])
H_s = (1/8)*(6*(vn[i,j-1]) + 3*(vn[i,j]) - vn[i,j-2])
else
G_s = (1/8)*(6*(vn[i,j]) + 3*(vn[i,j-1]) - vn[i,j+1])
H_s = (1/8)*(6*(vn[i,j]) + 3*(vn[i,j-1]) - vn[i,j+1])
end

conv_u = -((1/dx)*(F_e - F_w)) - ((1/dy)*(G_n - G_s))
conv_v = -((1/dx)*(G_e - G_w)) - ((1/dy)*(H_n - H_s))
return conv_u,conv_v
end

function BC(u,v)

u[end-1,:] .= u[end-2,:]
u[end,:] .= u[end-1,:]
u[1,:] .= u[end,:]
u[2,:] .= u[1,:]

u[:,end-1] .= u[:,end-2]
u[:,end] .= u[:,end-1]
u[:,1] .= u[:,end]
u[:,2] .= u[:,1]

v[end-1,:] .= v[end-2,:]
v[end,:] .= v[end-1,:]
v[1,:] .= v[end,:]
v[2,:] .= v[1,:]

v[:,end-1] .= v[:,end-2]
v[:,end] .= v[:,end-1]
v[:,1] .= v[:,end]
v[:,2] .= v[:,1]

end

anim = @animate for n in 1:nt
un = copy(u)
vn = copy(v)
u_max = maximum(abs.(u))
v_max = maximum(abs.(v))
dt = min(((c*dx)/u_max),((c*dy)/v_max))
@inbounds for j in 3:ny-2
@simd for i in 3:nx-2
Cu,Cv = convective(un,vn,dx,dy,i,j)
u[i,j] = un[i,j] + vel*dt*Cu
v[i,j] = vn[i,j] + vel*dt*Cv

end
end
global total_time += dt
BC(u,v)
u_res = sqrt.((u.^2) + (v.^2))
heatmap(y,x,u_res, title = "Time = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
println("max u:", maximum(abs.(u)),"  at --> n:  ", n)
end
println("total time:", total_time)
gif(anim, "2D_convection_QUICK_linear.gif", fps = 30)
