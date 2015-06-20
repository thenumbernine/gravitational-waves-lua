local class = require 'ext.class'

local ForwardEuler = class()

function ForwardEuler:integrate(qs, dt, dq_dts)
	return qs + dt * dq_dts()
end

local RungeKutta4 = class()

function RungeKutta4:integrate(qs, dt, dq_dts)
	local k1 = dq_dts(qs)
	local k2 = dq_dts(qs + .5 * k1)
	local k3 = dq_dts(qs + .5 * k2)
	local k4 = dq_dts(qs + k3)
	return qs + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
end

return {
	ForwardEuler = ForwardEuler,
	RungeKutta4 = RungeKutta4,
}

