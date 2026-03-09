using Plots

lx,ly = 2,2
nx,ny = 401,401
c = 0.1
u_amp,v_amp = 1,1
dx,dy = lx/(nx-1), ly/(ny-1)
nt = 200
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

anim = @animate for n in 1:nt
un = copy(u)
vn = copy(v)

u_max = maximum(abs.(u))
v_max = maximum(abs.(v))
dt = min(((c*dx)/u_max),((c*dy)/v_max))

@inbounds for j in 2:ny-1
@simd for i in 2:nx-1
Cu,Cv = convective(un,vn,dx,dy,i,j)

u[i,j] = un[i,j] + dt*Cu
v[i,j] = vn[i,j] + dt*Cv

end
end
BC(u,v)

#resultant velocity
u_res = sqrt.((u.^2) + (v.^2))
heatmap(y,x,u_res, title = "Time step = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
println("max u:", maximum(abs.(u)),"  at --> n:  ", n)
end
println("total time:", total_time)
gif(anim, "2D_convection_upwind.gif", fps = 30)
