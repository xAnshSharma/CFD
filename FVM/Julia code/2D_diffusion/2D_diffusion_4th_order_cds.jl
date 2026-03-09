using Plots   #to include contour plots and animations

lx,ly = 2,2   #dimensions of domain
nx,ny = 401,401   #number of cell nodes in domain
Fo = 0.1    #fourier number
u_amp,v_amp = 1,1   #disturbance amplitude
dx,dy = lx/(nx-1), ly/(ny-1)    #distance between cell nodes
nt = 200    #total number of time steps
nu = 1e-3   #kinematic viscosity
total_time = 0

x,y = LinRange(0,lx,nx),LinRange(0,ly,ny)   #defining axes
u,v,u_res = zeros(nx,ny),zeros(nx,ny),zeros(nx,ny)    #initialising velocity fields

xc, yc = lx/2, ly/2
sigma = 0.1

@inbounds for j in 1:ny
@simd for i in 1:nx
    pulse = exp(-((x[i]-xc)^2 + (y[j]-yc)^2)/(2*sigma^2))

    u[i,j] = u_amp*pulse
    v[i,j] = v_amp*pulse
end
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

#Boundary conditions
function BC(u,v)

u[1:2,:] .= 0   #u at x = 0
u[end-1:end,:] .= 0   #u at x = Lx
u[:,1:2] .= 0   #u at y = 0
u[:,end-1:end] .= 0   #u at y = Ly

v[1:2,:] .= 0   #v at x = 0
v[end-1:end,:] .= 0   #v at x = Lx
v[:,1:2] .= 0   #v at y = 0
v[:,end-1:end] .= 0   #v at y = Ly

end

#Solution
anim = @animate for n in 1:nt
un = copy(u)
vn = copy(v)

dt = min(((Fo*(dx^2))/nu),((Fo*(dy^2))/nu))

@inbounds for j in 3:ny-2
@simd for i in 3:nx-2
Du = diffusion_u(un,vn,dx,dy,i,j,nu)
Dv = diffusion_v(un,vn,dx,dy,i,j,nu)

u[i,j] = un[i,j] + dt*Du
v[i,j] = vn[i,j] + dt*Dv

end
end
BC(u,v)
global total_time += dt

#resultant velocityhat 
u_res .= sqrt.(u.^2 .+ v.^2)
heatmap(y,x,u_res, title = "Time step = $n", xlabel = "X", ylabel = "Y", legend = false, xlims = (0,lx), ylims = (0,ly))
u_max = maximum(abs.(u))
println("max u:", u_max,"  at --> n:  ", n)

end

#generating animation
println("total simulated time:", total_time)
gif(anim, "2D_diffusion_4th_order_cds.gif", fps = 30)
