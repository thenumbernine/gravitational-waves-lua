--[[
[rho],t =  [0											1						0			] [rho],x
[ m ],t = -[(gamma-3)/2*m^2/rho^2						(3-gamma)*m/rho			gamma-1		] [ m ],x
[ E ],t =  [(-gamma*m*E/rho^2 + (gamma-1)*m^3/rho^3)		H+(1-gamma)*m^2/rho^2	gamma*m/rho	] [ E ],x

[rho],t = -m,x
[ m ],t = (3-gamma)/2*m^2/rho^2 * rho,x							+ (gamma-3)*m/rho * m,x			+ (1-gamma) * E,x
[ E ],t = (gamma*m*E/rho^2 + (1-gamma)*m^3/rho^3) * rho,x		+ ((gamma-1)*m^2/rho^2-H) * m,x	+ (-gamma*m/rho) * E,x
	for H = E/rho + (gamma-1)*rho*(E/rho - .5*m^2/rho^2)

spatial derivative:
let q_j = 1/N sum{k=0,N-1 of qhat_k exp(i 2 pi j k / N)}
then dq/dj = 1/N sum{k=0,N-1 of i 2 pi k / N qhat_k exp(i 2 pi j k / N)}

x = (j-1)/(N-1) * (xmax-xmin) + xmin
dx/dj = (xmax - xmin) / (N-1)
dj/dx = (N-1)/(xmax - xmin)
dq/dx = dq/dj * dj/dx
dq/dx = dq/dj * (N-1)/(xmax-xmin)

time derivative:
let q(t)_j = 1/N sum{k=0,N-1 of qhat(t)_k exp(i 2 pi j k / N)}
then q'(t)_j = 1/N sum{k=0,N-1 of qhat'(t)_k exp(i 2 pi j k / N)}
for qhat(t)_j = sum{k=0,N-1 of q(t)_k exp(-i 2 pi j k / N)}
and qhat'(t)_j = sum{k=0,N-1 of q'(t)_k exp(-i 2 pi j k / N)}

q(t+dt)_j = q(t)_j + dt * q'(t)_j
		 = q(t)_j + dt * 1/N sum{k=0,N-1 of qhat'(t)_k exp(i 2 pi j k / N)}			<- based on spatial derivatives of dft of q

rho(t+dt)_j = rho_j + dt * rho_j,t
		= rho_j + dt * (-rhoFlux_j,x)
		= rho_j + dt * 1/N sum{k=0,N-1 of (i 2 pi k / N) (-rhoFluxHat_k) exp(i 2 pi j k / N)}	<- based on spatial derivatives of dft of q

--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Solver = require 'solver'

local EulerDFT = class(Solver)

EulerDFT.fixed_dt = 1/100

local function real(x)
	if type(x) == 'table' then return x[1] end
	return x
end

local function imag(x)
	if type(x) == 'table' then return x[2] end
	return 0
end

local function dft(dir, xs)
	local hats = {}
	local N = #xs
	for i=1,N do
		local re = 0
		local im = 0
		for k=1,N do
			local theta = dir * 2 * math.pi * (i-1) * (k-1) / N
			local costheta = math.cos(theta)
			local sintheta = math.sin(theta)
			local xre = real(xs[k])
			local xim = imag(xs[k])
			re = re + xre * costheta - xim * sintheta
			im = im + xre * sintheta + xim * costheta
		end
		hats[i] = {re,im}
	end
	return hats
end

-- x-derivative of forward transform
local function dft_x(xs)
	local hats = {}
	local N = #xs
	for i=1,N do
		local re = 0
		local im = 0
		for k=1,N do
			local theta = 2 * math.pi * (i-1) * (k-1) / N
			local expre = math.cos(theta)
			local expim = math.sin(theta)
			local xre = real(xs[k])
			local xim = imag(xs[k])
			-- apply derivative
			local s = 2 * math.pi * (k-1) / N
			-- why does there need to be a -1 here?
			s = s * -1
			-- multiply by 'i'...
			xre, xim = -s * xim, s * xre
			re = re + xre * expre - xim * expim
			im = im + xre * expim + xim * expre
		end
		hats[i] = {re,im}
	end
	return hats
end

function EulerDFT:iterate()
	self:applyBoundary()
	
	local dt = self.fixed_dt
	
	local N = self.gridsize
	local xmin = self.domain.xmin
	local xmax = self.domain.xmax

-- [=[ dft on flux then calc the dft's x-deriv
	local rhoFluxes = {}
	local jFluxes = {}
	local EFluxes = {}
	for i=1,self.gridsize do
		rhoFluxes[i], jFluxes[i], EFluxes[i] = unpack(self.equation:calcFluxForState(self.qs[i]))
	end
	
	local rhoFluxHats = dft(-1, rhoFluxes)		-- -1 = inverse dft
	local jFluxHats = dft(-1, jFluxes)		-- -1 = inverse dft
	local EFluxHats = dft(-1, EFluxes)		-- -1 = inverse dft
	
	local dx_rhoFluxes = dft_x(rhoFluxHats)
	local dx_jFluxes = dft_x(jFluxHats)
	local dx_EFluxes = dft_x(EFluxHats)
	for i=1,N do
		for j=1,2 do
			dx_rhoFluxes[i][j] = dx_rhoFluxes[i][j] * (N-1) / (xmax - xmin) / N^2
			dx_jFluxes[i][j] = dx_jFluxes[i][j] * (N-1) / (xmax - xmin) / N^2
			dx_EFluxes[i][j] = dx_EFluxes[i][j] * (N-1) / (xmax - xmin) / N^2
		end
	end

	for i=1,N do
		self.qs[i][1] = self.qs[i][1] - dt * real(dx_rhoFluxes[i])
		self.qs[i][2] = self.qs[i][2] - dt * real(dx_jFluxes[i])
		self.qs[i][3] = self.qs[i][3] - dt * real(dx_EFluxes[i])
	end
--]=]

--[=[ test
	local rhos = {}
	local js = {}
	local Es = {}
	for i=1,N do
		rhos[i], js[i], Es[i] = unpack(self.qs[i])
	end

	local rhoHats = dft(-1, rhos)
	local jHats = dft(-1, js)
	local EHats = dft(-1, Es)

	local rhos = dft(1, rhoHats)
	local js = dft(1, jHats)
	local Es = dft(1, EHats)

	-- replace channel 2 with derivative of channel 1 - for verification 
	local drhos = dft_x(rhoHats)
	for i=1,N do
		-- where did I lose a minus sign?
		drhos[i][1] = drhos[i][1] * (N-1) / (xmax - xmin) / N
		drhos[i][2] = drhos[i][2] * (N-1) / (xmax - xmin) / N
	end
	js = drhos
	
	for i=1,N do
		self.qs[i] = {real(rhos[i])/N, real(js[i])/N, real(Es[i])/N}
	end
--]=]

	self.t = self.t + dt
end

return EulerDFT
