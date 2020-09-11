using Plots;
using Printf;
gr();
anim = Animation();
f(x) = 1/2*cos(8pi*(x-1/2));

function forward(u_old::Array, i::Int64)
	return -(u[ir2[i]] - 2u[ir[i]] + 9u[i] - 10u[il[i]] + 2u[il2[i]])C / 6;
end

function backward(u_old::Array, i::Int64)
	return -(-2u[ir2[i]] + 10u[ir[i]] - 9u[i] + 2[il[i]] - u[il2[i]])C / 6;
end


function rk4()
	##First Substep
	phi1 = zeros(Float64, nx);
	@inbounds for i = 1:nx
		phi1[i] = u[i] + dt * forward(u, i);
	end

	##Second Substep
	phi2 = zeros(Float64, nx);
	@inbounds for i = 1:nx
		phi2[i] = u[i] + dt * forward((u + phi1) ./ 2, i);
	end

	##Third Substep
	phi3 = zeros(Float64, nx);
	@inbounds for i = 1:nx
		phi3[i] = u[i] + dt * forward((u + phi2) ./ 2, i);
	end

	##Forth Substep
	phi4 = zeros(Float64, nx);
	@inbounds for i = 1:nx
		phi4[i] = u[i] +  dt * forward(phi3, i);
	end

	@inbounds for i = 1:nx
		global u[i] = (phi1[i] + 2phi2[i] + 2phi3[i] + phi4[i]) / 6.0;
	end
end

function main()
#=********Numerical Set Up*********=#
	global c = 1.0;
	global nx = 101;
	global lx = 100;
	global x = range(0.0, stop=lx, length=nx);
	global dx = lx /(nx-1);
	global u  = zeros(Float64, nx);
	@inbounds for i = 1:nx
		if 40.0 < x[i] < 60.0; 
			u[i] = f(x[i]);
		end
	end
	global t = 0.0;
	global dt = 5e-1;
	global tlims = 100;
	global C = dt * c / dx;
	global ir = zeros(Int64, nx);
	global ir2 = zeros(Int64, nx);
	global il = zeros(Int64, nx);
	global il2 = zeros(Int64, nx);
	@inbounds for i = 1:nx
		ir[i] = i + 1;
		ir2[i] = i + 2;
		il[i] = i - 1;
		il2[i] = i - 2;
	end
	ir[end] = 1;
	ir2[end - 1] = 1;
	ir2[end] = 2;
	il[begin] = nx;
	il2[begin + 1] = nx;
	il2[begin] = nx - 1;
#=******Finish Numerical Set Up*************=#

	while t < tlims
		rk4();
		plt = plot(x, u, xlims=(0, lx), xticks=0:lx/5:lx, ylims=(-0.2, 1.2), yticks=-0.2:0.2:1.2,
		xlabel="x", ylabel="u", color="blue", title="KK scheme by using RK4", legend=false);
		frame(anim, plt);
		t += dt;
		@printf("t = %.4f\n", t);
	end
	gif(anim, "KK_rk4.gif", fps=5);
end

main();
