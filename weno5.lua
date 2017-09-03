-- based on code at https://www.mathworks.com/matlabcentral/fileexchange/44639-weighted-essentially-non-oscillatory--weno--scheme
local class = require 'ext.class'
local SolverFV = require 'solverfv'

local WENO = class(SolverFV)

function WENO:init(args)
	self.equation = assert(args.equation or self.equation)

	self.name = self.equation.name .. ' WENO'

	WENO.super.init(self, args)
end

function WENO:calcFluxes(dt)
--[[
	-- Lax-Friedrichs Flux Splitting
	local a=max(abs(dflux(w)))
	local v=0.5*(flux(w)+a*w)
	local u=circshift(0.5*(flux(w)-a*w),[0,-1])

	-- Right Flux
	-- Choose the positive fluxes, 'v', to compute the left cell boundary flux:
	-- $u_{i+1/2}^{-}$
	local vmm = circshift(v,2)
	local vm  = circshift(v,1)
	local vp  = circshift(v,-1)
	local vpp = circshift(v,-2)

	-- Polynomials
	p0n = (2*vmm - 7*vm + 11*v)/6
	p1n = ( -vm  + 5*v  + 2*vp)/6
	p2n = (2*v   + 5*vp - vpp )/6

	-- Smooth Indicators (Beta factors)
	B0n = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2
	B1n = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2
	B2n = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2

	-- Constants
	d0n = 1/10
	d1n = 6/10
	d2n = 3/10
	epsilon = 1e-6

	-- Alpha weights
	alpha0n = d0n./(epsilon + B0n).^2
	alpha1n = d1n./(epsilon + B1n).^2
	alpha2n = d2n./(epsilon + B2n).^2
	alphasumn = alpha0n + alpha1n + alpha2n

	-- ENO stencils weigths
	w0n = alpha0n./alphasumn
	w1n = alpha1n./alphasumn
	w2n = alpha2n./alphasumn

	-- Numerical Flux at cell boundary, $u_{i+1/2}^{-}$
	hn = w0n.*p0n + w1n.*p1n + w2n.*p2n

	---- Left Flux
	-- Choose the negative fluxes, 'u', to compute the left cell boundary flux:
	-- $u_{i-1/2}^{+}$
	umm = circshift(u,[0 2])
	um  = circshift(u,[0 1])
	up  = circshift(u,[0 -1])
	upp = circshift(u,[0 -2])

	-- Polynomials
	p0p = ( -umm + 5*um + 2*u  )/6
	p1p = ( 2*um + 5*u  - up   )/6
	p2p = (11*u  - 7*up + 2*upp)/6

	-- Smooth Indicators (Beta factors)
	B0p = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2
	B1p = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2
	B2p = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2

	-- Constants
	d0p = 3/10
	d1p = 6/10
	d2p = 1/10
	epsilon = 1e-6

	-- Alpha weights
	alpha0p = d0p./(epsilon + B0p).^2
	alpha1p = d1p./(epsilon + B1p).^2
	alpha2p = d2p./(epsilon + B2p).^2
	alphasump = alpha0p + alpha1p + alpha2p

	-- ENO stencils weigths
	w0p = alpha0p./alphasump
	w1p = alpha1p./alphasump
	w2p = alpha2p./alphasump

	-- Numerical Flux at cell boundary, $u_{i-1/2}^{+}$
	hp = w0p.*p0p + w1p.*p1p + w2p.*p2p

	---- Compute finite volume residual term, df/dx.
	res = (hp-circshift(hp,[0,1])+hn-circshift(hn,[0,1]))/dx - S(w)
--]]
end

return WENO
