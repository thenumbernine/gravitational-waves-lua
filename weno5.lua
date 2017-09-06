-- based on code at https://www.mathworks.com/matlabcentral/fileexchange/44639-weighted-essentially-non-oscillatory--weno--scheme
local class = require 'ext.class'
local matrix = require 'matrix'

-- technically a Godunov solver and not a Roe solver
-- and that's where I should put the eigenbasis stuff
local SolverFV = require 'solverfv'

local WENO5 = class(SolverFV)

function WENO5:init(args)
	self.equation = assert(args.equation or self.equation)

	self.name = self.equation.name .. ' WENO5'

	WENO5.super.init(self, args)
end

function WENO5:calcFluxes(dt)

	local numStates = self.numStates
	
	-- max eigenvalue at cell center
	local A = self:newState()
	for i=1,self.gridsize do
		local lambdaMin, lambdaMax = self.equation:calcCellMinMaxEigenvalues(self, i)
		A[i] = math.max(
			math.abs(lambdaMin), 
			math.abs(lambdaMax))
	end

	-- cell centered prims
	local ps = self:newState()
	for i=1,self.gridsize do
		ps[i] = matrix{self.equation:calcPrimFromCons(
			table.unpack(self.qs[i])
		)}
	end

	-- cell-centered fluxes
	local F = self:newState()
	for i=1,self.gridsize do
		F[i] = matrix{self.equation:calcFluxForState(self.qs[i])}
	end

	for i=3,self.gridsize-3 do

		-- interface prims based on cell averaging
		local ip = .5 * (ps[i] + ps[i+1])

		-- interface cons based on interface prims
		local iq = matrix{self.equation:calcConsFromPrim(ip:unpack())}

		-- interface eigensystem based on interface cons
		local ilambda = matrix.zeros{numStates}
		local ievr = matrix.zeros{numStates, numStates}
		local ievl = matrix.zeros{numStates, numStates}
		self.equation:calcEigenBasisFromCons(ilambda, ievr, ievl, nil, iq)

		-- max wavespeed in local stencil
		--[[
		local ml = math.max(A:unpack(i-2, i+3))
		--]]
		-- [[
		local ml = -math.huge
		for j=1,6 do
			local m = i + j - 3
			ml = math.max(ml, A[m])
		end
		--]]

		-- fp[curve sample index][state variable]
		local fp = matrix()
		local fm = matrix()
		for j=1,6 do
			local m = i + j - 3
			
			local Fp = matrix()
			local Fm = matrix()
			for q=1,numStates do
				Fp[q] = .5 * (F[m][q] + ml * self.qs[m][q])
				Fm[q] = .5 * (F[m][q] - ml * self.qs[m][q])
			end
		
			-- matrix mul -- apply left eigenvalues
			--[[
			fp[j] = ievl * Fp
			fm[j] = ievl * Fm
			--]]
			-- [[
			fp[j] = matrix()
			fm[j] = matrix()
			for q=1,numStates do
				fp[j][q] = 0
				fm[j][q] = 0
				for r=1,numStates do
					fp[j][q] = fp[j][q] + ievl[q][r] * Fp[r]
					fm[j][q] = fm[j][q] + ievl[q][r] * Fm[r]
				end
			end
			--]]
		end
	
		-- reconstruct ...
		local function sqr(x) return x*x end
		local function weno5(v, c, d)
			local vs = {
				c[1][1]*v[3] + c[1][2]*v[4] + c[1][3]*v[5],
				c[2][1]*v[2] + c[2][2]*v[3] + c[2][3]*v[4],
				c[3][1]*v[1] + c[3][2]*v[2] + c[3][3]*v[3],
			}
			
			local B = { -- smoothness indicators
				(13/12)*sqr(v[ 3] - 2*v[ 4] + v[ 5]) + (1/4)*sqr(3*v[ 3] - 4*v[ 4] +   v[ 5]),
				(13/12)*sqr(v[ 2] - 2*v[ 3] + v[ 4]) + (1/4)*sqr(  v[ 2]           -   v[ 4]),
				(13/12)*sqr(v[ 1] - 2*v[ 2] + v[ 3]) + (1/4)*sqr(  v[ 1] - 4*v[ 2] + 3*v[ 3]),
			}

			local w = matrix()
			--[[
			if (IS_mode == ImprovedBorges08) {
				eps = eps_prime = 1e-14; -- Borges uses 1e-40, but has Matlab
				const double tau5 = fabs(B[0] - B[2]);

				-- Calculate my weights with new smoothness indicators accoding to Borges
				w[0] = d[0] * (1.0 + (tau5 / (B[0] + eps)));
				w[1] = d[1] * (1.0 + (tau5 / (B[1] + eps)));
				w[2] = d[2] * (1.0 + (tau5 / (B[2] + eps)));
			}
			else if (IS_mode == ImprovedShenZha10) {
				eps = 1e-6;
				eps_prime = 1e-10;
				const double A = shenzha10_A; -- [0 (less aggressive) -> ~100 (more aggressive)]
				const double minB = min3(B), maxB = max3(B);
				const double R0 = minB / (maxB + eps_prime);
				B[0] = R0*A*minB + B[0];
				B[1] = R0*A*minB + B[1];
				B[2] = R0*A*minB + B[2];
				w[0] = d[0] / sqr(eps_prime + B[0]);
				w[1] = d[1] / sqr(eps_prime + B[1]);
				w[2] = d[2] / sqr(eps_prime + B[2]);
			}
			else -- Use OriginalJiangShu96
			]] do	
				 -- recommended value by Jiang and Shu
				local eps = 1e-6
				local eps_prime = 1e-6
				w[1] = d[1] / sqr(eps_prime + B[1])
				w[2] = d[2] / sqr(eps_prime + B[2])
				w[3] = d[3] / sqr(eps_prime + B[3])
			end

			local wtot = w[1] + w[2] + w[3]
			return (w[1]*vs[1] + w[2]*vs[2] + w[3]*vs[3])/wtot
	
		end

		local cls = { {11/6, -7/6,  1/3 }, { 1/3,  5/6, -1/6 }, {-1/6,  5/6,  1/3 } }
		local crs = { { 1/3,  5/6, -1/6 }, {-1/6,  5/6,  1/3 }, { 1/3, -7/6, 11/6 } }
		local dls = { 0.1, 0.6, 0.3 }
		local drs = { 0.3, 0.6, 0.1 }
	
		-- fpT[state variable][curve sample index]
		local fpT = fp:T()
		local fmT = fm:T()
		local f = matrix()
		for q=1,numStates do
			f[q] = weno5(fpT[q], cls, dls)
					+ weno5(fmT[q], crs, drs)
		end

		--[[
		self.fluxes[i] = ievr * f	
		--]]
		for j=1,numStates do
			self.fluxes[i][j] = 0
			for k=1,numStates do
				self.fluxes[i][j] = self.fluxes[i][j] + ievr[j][k] * f[k]
			end
		end
	end

--[[
	-- TODO characteristic variables
	local w = self.qs
	
	-- flux(w) returns the flux for the state
	local F = flux(w)

	-- ... what kind of flux function?
	-- dflux(w) = dF/dx
	local dF = dflux(w)

	-- Lax-Friedrichs Flux Splitting
	local a=dF:map(mat.abs):sup()
	local v=0.5*(F+a*w)
	local u = circshift(0.5*(F-a*w),-1)

	-- Right Flux
	-- Choose the positive fluxes, 'v', to compute the left cell boundary flux:
	-- $u_{i+1/2}^{-}$
	local vmm = circshift(v,2)
	local vm  = circshift(v,1)
	local vp  = circshift(v,-1)
	local vpp = circshift(v,-2)

	-- Polynomials
	local p0n = (2*vmm - 7*vm + 11*v)/6
	local p1n = ( -vm  + 5*v  + 2*vp)/6
	local p2n = (2*v   + 5*vp - vpp )/6

	-- Smooth Indicators (Beta factors)
	local B0n = 13/12*(vmm-2*vm+v  )^2 + 1/4*(vmm-4*vm+3*v)^2
	local B1n = 13/12*(vm -2*v +vp )^2 + 1/4*(vm-vp)^2
	local B2n = 13/12*(v  -2*vp+vpp)^2 + 1/4*(3*v-4*vp+vpp)^2

	-- Constants
	local d0n = 1/10
	local d1n = 6/10
	local d2n = 3/10
	local epsilon = 1e-6

	-- Alpha weights
	local alpha0n = d0n/((epsilon + B0n)^2)
	local alpha1n = d1n/((epsilon + B1n)^2)
	local alpha2n = d2n/((epsilon + B2n)^2)
	local alphasumn = alpha0n + alpha1n + alpha2n

	-- ENO stencils weigths
	local w0n = alpha0n / alphasumn
	local w1n = alpha1n / alphasumn
	local w2n = alpha2n / alphasumn

	-- Numerical Flux at cell boundary, $u_{i+1/2}^{-}$
	local hn = w0n*p0n + w1n*p1n + w2n*p2n

	---- Left Flux
	-- Choose the negative fluxes, 'u', to compute the left cell boundary flux:
	-- $u_{i-1/2}^{+}$
	local umm = circshift(u, 2)
	local um  = circshift(u, 1)
	local up  = circshift(u, -1)
	local upp = circshift(u, -2)

	-- Polynomials
	local p0p = ( -umm + 5*um + 2*u  )/6
	local p1p = ( 2*um + 5*u  - up   )/6
	local p2p = (11*u  - 7*up + 2*upp)/6

	-- Smooth Indicators (Beta factors)
	local B0p = 13/12*(umm-2*um+u  )^2 + 1/4*(umm-4*um+3*u)^2
	local B1p = 13/12*(um -2*u +up )^2 + 1/4*(um-up)^2
	local B2p = 13/12*(u  -2*up+upp)^2 + 1/4*(3*u -4*up+upp)^2

	-- Constants
	local d0p = 3/10
	local d1p = 6/10
	local d2p = 1/10
	local epsilon = 1e-6

	-- Alpha weights
	local alpha0p = d0p/((epsilon + B0p)^2)
	local alpha1p = d1p/((epsilon + B1p)^2)
	local alpha2p = d2p/((epsilon + B2p)^2)
	local alphasump = alpha0p + alpha1p + alpha2p

	-- ENO stencils weigths
	local w0p = alpha0p / alphasump
	local w1p = alpha1p / alphasump
	local w2p = alpha2p / alphasump

	-- Numerical Flux at cell boundary, $u_{i-1/2}^{+}$
	local hp = w0p*p0p + w1p*p1p + w2p*p2p

	---- Compute finite volume residual term, df/dx.
	return (hp-circshift(hp,1)+hn-circshift(hn,1))/dx - S(w)
--]]
end

return WENO5
