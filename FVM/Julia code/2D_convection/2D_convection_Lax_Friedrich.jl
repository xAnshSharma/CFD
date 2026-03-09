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

function convective(un,vn,i,j,dx)

F_e = un[i+1,j]^2
F_w = un[i-1,j]^2

G_n = vn[i,j+1]*un[i,j+1]
G_s = vn[i,j-1]*un[i,j-1]

G_e = un[i+1,j]*vn[i+1,j]
G_w = un[i-1,j]*vn[i-1,j]

H_n = vn[i.j+1]^2
H_s = vn[i,j-1]^2

conv_u = -((1/dx)*(F_e - F_w)) - ((1/dy)*(G_n - G_s))
conv_v = -((1/dx)*(G_e - G_w)) - ((1/dy)*(H_n - H_s))
return conv_u,conv_v
end

function BC(u)
u[1,:] .= 0
u[end,:] .= 0
u[:,1] .= 0
u[end,:] .= 0

v[1,:] .= 0
v[end,:] .= 0
v[:,1] .= 0
v[end,:] .= 0

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

u[i,j] = 0.25*(un[i+1,j] + un[i-1,j] + un[i,j+1] + un[i,j-1]) + dt*Cu
v[i,j] = 0.25*(vn[i+1,j] + vn[i-1,j] + vn[i,j+1] + vn[i,j-1]) + dt*Cv

end
end
global total_time += dt
BC(u,v)

u_res = sqrt.((u.^2) + (v.^2))
heatmap(y,x,u_res, title = "Time = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
println("max u:", maximum(abs.(u)),"  at --> n:  ", n)
end
println("total time:", total_time)
gif(anim, "2D_convection_Lax_Friedrich.gif", fps = 30)
