local class = require 'ext.class'

local function fill2d(dst,src)
	for i=1,#dst do
		for j=1,#dst[i] do
			dst[i][j] = src[i][j]
		end
	end
end

local ForwardEuler = class()
ForwardEuler.name = 'F.E.'
function ForwardEuler:integrate(qs, dt, dq_dts)
	return qs + dt * dq_dts()
end

local RungeKutta4 = class()
RungeKutta4.name = 'RK4'
function RungeKutta4:integrate(qs, dt, dq_dts)
	local orig = qs * 1	-- copy?
	local k1 = dt * dq_dts()
	fill2d(qs, orig + k1 * .5)
	local k2 = dt * dq_dts()
	fill2d(qs, orig + .5 * k2)
	local k3 = dt * dq_dts()
	fill2d(qs, orig + k3)
	local k4 = dt * dq_dts()
	return orig + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
end

local BackwardEuler = class()
BackwardEuler.name = 'B.E.'
function BackwardEuler:integrate(qs, dt, dq_dts)
	return require 'solver.gmres'{
		A = function(u)
			-- TODO this is only sampling dq_dts() from the u(t) timestep
			-- so you can either store it once outside the solve loop
			-- or you can try to recalculate it each time based on u instead of qs
			-- (but that might mean rearranging dq_dts to accept input parameters of its state)
			return u - dt * dq_dts()
		end,
		x = qs:clone(),
		b = qs:clone(),
		clone = qs.clone,
		dot = qs.dot,
	}
end

return {
	ForwardEuler = ForwardEuler,
	RungeKutta4 = RungeKutta4,
	BackwardEuler = BackwardEuler,
}
