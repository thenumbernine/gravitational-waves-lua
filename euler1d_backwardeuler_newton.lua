--[[

m = rho u
E = rho (e + 1/2 u^2)
E/rho = e + 1/2 u^2
e = E/rho - 1/2 u^2
p = (gamma - 1) rho e
p = (gamma - 1) rho (E/rho - 1/2 u^2)
p = (gamma - 1) (E - 1/2 rho u^2)
p = (gamma - 1) (E - 1/2 m^2 / rho)

diff_t rho = diff_x (-m)
diff_t m = diff_x ((1 - gamma) (E - 1/2 m^2 / rho) - m^2 / rho)
diff_t E = diff_x ((1/2 (gamma - 1) m^2 / rho - gamma E) m/rho)

derivative approximation:

diff_t rho[x] = -(m[x+1] - m[x-1])/(2*dx)
diff_t m[x] =   ( ( (1-gamma) E[x+1] + (gamma-3)/2 m^2[x+1]/rho[x+1])
				- ( (1-gamma) E[x-1] + (gamma-3)/2 m^2[x-1]/rho[x-1]) )/(2*dx)
diff_t E[x] = ( (1/2(gamma-1) m[x+1]^3/rho[x+1]^2 - gamma E[x+1] m[x+1]/rho[x+1])
			  - (1/2(gamma-1) m[x-1]^3/rho[x-1]^2 - gamma E[x-1] m[x-1]/rho[x-1]) )/(2*dx)

first derivatives:
	
	d/d(m[x+1]) diff_t rho[x] = -1/(2*dx)
	d/d(m[x-1]) diff_t rho[x] = 1/(2*dx)
	
	d/d(rho[x+1]) diff_t m[x] = -(gamma-3)/2 m[x+1]^2/rho[x+1]^2/(2*dx)
	d/d(rho[x-1]) diff_t m[x] = (gamma-3)/2 m[x-1]^2/rho[x-1]^2/(2*dx)
	
	d/d(m[x+1]) diff_t m[x] = (gamma-3) m[x+1]/rho[x+1]/(2*dx)
	d/d(m[x-1]) diff_t m[x] = -(gamma-3) m[x-1]/rho[x-1]/(2*dx)
	
	d/d(E[x+1]) diff_t m[x] = (1-gamma)/(2*dx)
	d/d(E[x-1]) diff_t m[x] = -(1-gamma)/(2*dx)
	
	d/d(rho[x+1]) diff_t E[x] = ((1-gamma) m[x+1]^3/rho[x+1]^3 + gamma E[x+1] m[x+1]/rho[x+1]^2)/(2*dx)
	d/d(rho[x-1]) diff_t E[x] = -((1-gamma) m[x-1]^3/rho[x-1]^3 + gamma E[x-1] m[x-1]/rho[x-1]^2)/(2*dx)
	
	d/d(m[x+1]) diff_t E[x] = -(3/2(1-gamma) m[x+1]^2/rho[x+1]^2 + gamma E[x+1]/rho[x+1])/(2*dx)
	d/d(m[x-1]) diff_t E[x] = (3/2(1-gamma) m[x-1]^2/rho[x-1]^2 + gamma E[x-1]/rho[x-1])/(2*dx)
		
	d/d(E[x+1]) diff_t E[x] = -(gamma m[x+1]/rho[x+1])/(2*dx)
	d/d(E[x-1]) diff_t E[x] = (gamma m[x-1]/rho[x-1])/(2*dx)

second derivatives:

	d/d(rho[x+1]) d/d(rho[x+1]) diff_t m[x] = (gamma-3) m[x+1]^2/rho[x+1]^3/(2*dx)
	d/d(rho[x-1]) d/d(rho[x-1]) diff_t m[x] = -(gamma-3) m[x-1]^2/rho[x-1]^3/(2*dx)
	
	d/d(m[x+1]) d/d(rho[x+1]) diff_t m[x] = -(gamma-3) m[x+1]/rho[x+1]^2/(2*dx)
	d/d(m[x-1]) d/d(rho[x-1]) diff_t m[x] = (gamma-3) m[x-1]/rho[x-1]^2/(2*dx)
	
	d/d(m[x+1]) d/d(m[x+1]) diff_t m[x] = (gamma-3)/rho[x+1]/(2*dx)
	d/d(m[x-1]) d/d(m[x-1]) diff_t m[x] = -(gamma-3)/rho[x-1]/(2*dx)
	
	d/d(rho[x+1]) d/d(rho[x+1]) diff_t E[x] = -(2 gamma E[x+1] m[x+1]/rho[x+1]^3 + 3(1-gamma) m[x+1]^3/rho[x+1]^4)/(2*dx)
	d/d(rho[x-1]) d/d(rho[x-1]) diff_t E[x] = (2 gamma E[x-1] m[x-1]/rho[x-1]^3 + 3(1-gamma) m[x-1]^3/rho[x-1]^4)/(2*dx)

	d/d(m[x+1]) d/d(rho[x+1]) diff_t E[x] = (gamma E[x+1]/rho[x+1]^2 + 3(1-gamma) m[x+1]^2/rho[x+1]^3)/(2*dx)
	d/d(m[x-1]) d/d(rho[x-1]) diff_t E[x] = -(gamma E[x-1]/rho[x-1]^2 + 3(1-gamma) m[x-1]^2/rho[x-1]^3)/(2*dx)
	
	d/d(E[x+1]) d/d(rho[x+1]) diff_t E[x] = gamma m[x+1]/rho[x+1]^2/(2*dx)
	d/d(E[x-1]) d/d(rho[x-1]) diff_t E[x] = -gamma m[x-1]/rho[x-1]^2/(2*dx)

	d/d(m[x+1]) d/d(m[x+1]) diff_t E[x] = 3(gamma-1) m[x+1]/rho[x+1]^2/(2*dx)
	d/d(m[x-1]) d/d(m[x-1]) diff_t E[x] = -3(gamma-1) m[x-1]/rho[x-1]^2/(2*dx)
	
	d/d(E[x+1]) d/d(m[x+1]) diff_t E[x] = -gamma/rho[x+1]/(2*dx)
	d/d(E[x-1]) d/d(m[x-1]) diff_t E[x] = gamma/rho[x-1]/(2*dx)

--]]

local class = require 'ext.class'
local table = require 'ext.table'

local Solver = require 'solver'

local EulerBackwardEulerNewton = class(Solver)

-- implicit needs a default fixed dt
-- I suppose I could use wavespeeds...
EulerBackwardEulerNewton.fixed_dt = 1/512

function EulerBackwardEulerNewton:init(args)
	args = table(args)
	args.equation = require 'euler1d'(args)
	EulerBackwardEulerNewton.super.init(self, args)
end

function EulerBackwardEulerNewton:iterate()
	local dt = self.fixed_dt

	local gamma = self.equation.gamma
	local q = self.qs		-- q(t)
	local nq = q:clone()	-- q(t+dt)

	-- hack for now ...
	local dx = self.ixs[2] - self.ixs[1]

	function dot(a,b)
		local y = 0
		for i=1,self.gridsize do
			for j=1,self.numStates do
				y = y + a[i][j] * b[i][j]
			end
		end
		return y
	end

	function norm(a) return dot(a,a) end

	-- newton iteration
	local newtonMaxIter = 50
	for newtonIter=1,newtonMaxIter do
		-- conjugate residual to solve the inverse of the hessian times the gradient
		-- solving for x with A*x = b
		-- aka solving x = A^-1 * b
		-- A = hessian matrix = d/dq_i d/dq_j e
		-- b = d/dq_i e
		
		-- f_i = partial_t q_i(t+dt)
		local f = self:newState() 
		for i=1,self.gridsize do
			-- diff_t rho[x] = -(m[x+1] - m[x-1])/(2*dx)
			-- diff_t m[x] =   ( ( (1-gamma) E[x+1] + (gamma-3)/2 m^2[x+1]/rho[x+1])
			-- 				   - ( (1-gamma) E[x-1] + (gamma-3)/2 m^2[x-1]/rho[x-1]) )/(2*dx)
			-- diff_t E[x] = ( (1/2(gamma-1) m[x+1]^3/rho[x+1]^2 - gamma E[x+1] m[x+1]/rho[x+1])
			-- 				 - (1/2(gamma-1) m[x-1]^3/rho[x-1]^2 - gamma E[x-1] m[x-1]/rho[x-1]) )/(2*dx)
			if i>1 then
				f[i][1] = f[i][1] + nq[i-1][2]/(2*dx)
				f[i][2] = f[i][2] - ((1-gamma) * nq[i-1][3] + (gamma-3)/2 * nq[i-1][2]^2/nq[i-1][1])/(2*dx)
				f[i][3] = f[i][3] - (.5*(gamma-1) * nq[i-1][2]^3/nq[i-1][1]^2 - gamma * nq[i-1][3]*nq[i-1][2]/nq[i-1][1])/(2*dx)
			end
			if i<self.gridsize then
				f[i][1] = f[i][1] - nq[i+1][2]/(2*dx)
				f[i][2] = f[i][2] + ((1-gamma) * nq[i+1][3] + (gamma-3)/2 * nq[i+1][2]^2/nq[i+1][1])/(2*dx)
				f[i][3] = f[i][3] + (.5*(gamma-1) * nq[i+1][2]^3/nq[i+1][1]^2 - gamma * nq[i+1][3]*nq[i+1][2]/nq[i+1][1])/(2*dx)
			end
		end
		
		local g = nq - dt * f - q		-- f returns the time deriv vector
		
		-- df_dq_ij = partial_q_i f_j = partial_q_i (partial_t q_j)
		-- returns y_i = df_dq_ij x_j = partial_q_i partial_t q_j x_j
		local function df_dq(x)
			local y = self:newState()
			for i=1,self.gridsize do
				if i>1 then
					-- d/d(m[x+1]) diff_t rho[x] = -1/(2*dx)
					y[i][2] = y[i][2] + x[i-1][1] * -1
					
					-- d/d(rho[x+1]) diff_t m[x] = -(gamma-3)/2 (m[x+1]^2/rho[x+1])^2/(2*dx)
					y[i][1] = y[i][1] + x[i-1][2] * -(gamma-3)/2 * nq[i][2]^2/nq[i][1]^2
					-- d/d(m[x+1]) diff_t m[x] = (gamma-3) m[x+1]/rho[x+1]/(2*dx)
					y[i][2] = y[i][2] + x[i-1][2] * (gamma-3) * nq[i][2]/nq[i][1]
					-- d/d(E[x+1]) diff_t m[x] = (1-gamma)/(2*dx)
					y[i][3] = y[i][3] + x[i-1][2] * (1-gamma)
					
					-- d/d(rho[x+1]) diff_t E[x] = ((1-gamma) m[x+1]^3/rho[x+1]^3 + gamma E[x+1] m[x+1]/rho[x+1]^2)/(2*dx)
					y[i][1] = y[i][1] + x[i-1][3] * ((1-gamma) * nq[i][2]^3/nq[i][1]^3 + gamma * nq[i][3] * nq[i][2]/nq[i][1]^2)
					-- d/d(m[x+1]) diff_t E[x] = -(3/2(1-gamma) m[x+1]^2/rho[x+1]^2 + gamma E[x+1]/rho[x+1])/(2*dx)
					y[i][2] = y[i][2] + x[i-1][3] * -(3/2*(1-gamma) * nq[i][2]^2/nq[i][1]^2 + gamma * nq[i][3]/nq[i][1])
					-- d/d(E[x+1]) diff_t E[x] = -(gamma m[x+1]/rho[x+1])/(2*dx)
					y[i][3] = y[i][3] + x[i-1][3] * -gamma * nq[i][2]/nq[i][1]
				end
				if i<self.gridsize then
					-- d/d(m[x-1]) diff_t rho[x] = 1/(2*dx)
					y[i][2] = y[i][2] - x[i+1][1] * -1
					
					-- d/d(rho[x-1]) diff_t m[x] = (gamma-3)/2 (m[x-1]^2/rho[x-1])^2/(2*dx)
					y[i][1] = y[i][1] - x[i+1][2] * -(gamma-3)/2 * nq[i][2]^2/nq[i][1]^2
					-- d/d(m[x-1]) diff_t m[x] = -(gamma-3) m[x-1]/rho[x-1]/(2*dx)
					y[i][2] = y[i][2] - x[i+1][2] * (gamma-3) * nq[i][2]/nq[i][1]
					-- d/d(E[x-1]) diff_t m[x] = -(1-gamma)/(2*dx)
					y[i][3] = y[i][3] -	x[i+1][2] * (1-gamma)
					
					-- d/d(rho[x-1]) diff_t E[x] = -((1-gamma) m[x-1]^3/rho[x-1]^3 + gamma E[x-1] m[x-1]/rho[x-1]^2)/(2*dx)
					y[i][1] = y[i][1] - x[i+1][3] * ((1-gamma) * nq[i][2]^3/nq[i][1]^3 + gamma * nq[i][3] * nq[i][2]/nq[i][1]^2)
					-- d/d(m[x-1]) diff_t E[x] = -(3/2(gamma-1) m[x-1]^2/rho[x-1]^2 - gamma E[x-1]/rho[x-1])/(2*dx)
					y[i][2] = y[i][2] - x[i+1][3] * -(3/2*(1-gamma) * nq[i][2]^2/nq[i][1]^2 + gamma * nq[i][3]/nq[i][1])
					-- d/d(E[x-1]) diff_t E[x] = -(gamma m[x-1]/rho[x-1])/(2*dx)
					y[i][3] = y[i][3] - x[i+1][3] * -gamma * nq[i][2]/nq[i][1]
				end
				y[i][1] = y[i][1] / (2 * dx)
				y[i][2] = y[i][2] / (2 * dx)
				y[i][3] = y[i][3] / (2 * dx)
			end
			return y
		end
		
		-- y_i = partial_q_j partial_t q_i x_j
		local function df_dq_t(x)
			local y = self:newState()
			for i=1,self.gridsize do
				if i<self.gridsize then
					-- d/d(m[x+1]) diff_t rho[x] = -1/(2*dx)
					y[i][1] = y[i][1] + x[i+1][2] * -1

					y[i][2] = y[i][2]
					-- d/d(rho[x+1]) diff_t m[x] = -(gamma-3)/2 m[x+1]^2/rho[x+1]^2/(2*dx)
						+ x[i+1][1] * -(gamma-3)/2 * nq[i+1][2]^2/nq[i+1][1]^2
					-- d/d(m[x+1]) diff_t m[x] = (gamma-3) m[x+1]/rho[x+1]/(2*dx)
						+ x[i+1][2] * (gamma-3) * nq[i+1][2]/nq[i+1][1]
					-- d/d(E[x+1]) diff_t m[x] = (1-gamma)/(2*dx)
						+ x[i+1][3] * (1-gamma)

					y[i][3] = y[i][3]
					-- d/d(rho[x+1]) diff_t E[x] = ((1-gamma) m[x+1]^3/rho[x+1]^3 + gamma E[x+1] m[x+1]/rho[x+1]^2)/(2*dx)
						+ x[i+1][1] * ((1-gamma) * nq[i+1][2]^3/nq[i+1][1]^3 + gamma * nq[i+1][3] * nq[i+1][2] / nq[i+1][1]^2)
					-- d/d(m[x+1]) diff_t E[x] = -(3/2(1-gamma) m[x+1]^2/rho[x+1]^2 + gamma E[x+1]/rho[x+1])/(2*dx)
						+ x[i+1][2] * -(3/2*(1-gamma) * nq[i+1][2]^2/nq[i+1][1]^2 + gamma * nq[i+1][3]/nq[i+1][1])
					-- d/d(E[x+1]) diff_t E[x] = -(gamma m[x+1]/rho[x+1])/(2*dx)
						+ x[i+1][3] * -(gamma * nq[i+1][2]/nq[i+1][1])
				end
				if i>1 then
					-- d/d(m[x-1]) diff_t rho[x] = 1/(2*dx)
					y[i][1] = y[i][1] - x[i-1][2] * -1
					
					y[i][2] = y[i][2]
					-- d/d(rho[x-1]) diff_t m[x] = (gamma-3)/2 m[x-1]^2/rho[x-1]^2/(2*dx)
						- x[i-1][1] * -(gamma-3)/2 * nq[i-1][2]^2/nq[i-1][1]^2
					-- d/d(m[x-1]) diff_t m[x] = -(gamma-3) m[x-1]/rho[x-1]/(2*dx)
						- x[i-1][2] * (gamma-3) * nq[i-1][2]/nq[i-1][1]
					-- d/d(E[x-1]) diff_t m[x] = -(1-gamma)/(2*dx)
						- x[i-1][3] * (1-gamma)
					
					y[i][3] = y[i][3]
					-- d/d(rho[x-1]) diff_t E[x] = -((1-gamma) m[x-1]^3/rho[x-1]^3 + gamma E[x-1] m[x-1]/rho[x-1]^2)/(2*dx)
						- x[i-1][1] * ((1-gamma) * nq[i-1][2]^3/nq[i-1][1]^3 + gamma * nq[i-1][3] * nq[i-1][2] / nq[i-1][1]^2)
					-- d/d(m[x-1]) diff_t E[x] = (3/2(1-gamma) m[x-1]^2/rho[x-1]^2 + gamma E[x-1]/rho[x-1])/(2*dx)
						+ x[i-1][2] * -(3/2*(1-gamma) * nq[i-1][2]^2/nq[i-1][1]^2 + gamma * nq[i-1][3]/nq[i-1][1])
					-- d/d(E[x-1]) diff_t E[x] = (gamma m[x-1]/rho[x-1])/(2*dx)
						+ x[i-1][3] * -(gamma * nq[i-1][2]/nq[i-1][1])
				end
			end
			return y
		end
	
		-- A_ij = d/dq_i d/dq_j e
		-- A_ij = delta_ij - dt * d/dq_j f_i - dt * sum_k (d/dq_j g_k * d/dq_i f_k + g_k * d/dq_i d/dq_j f_k )
		-- A_ij = delta_ij - dt * d/dq_j f_i - dt * sum_k (d/dq_j (q_k(t+dt) - dt * f_k - q_k(t)) * d/dq_i f_k + g_k * d/dq_i d/dq_j f_k )
		-- A_ij = delta_ij - dt * (d/dq_j f_i + d/dq_i f_j) + dt * sum_k (dt * d/dq_j f_k * d/dq_i f_k) - dt * sum_k (g_k * d/dq_i d/dq_j f_k)
		-- y_i = A_ij x_j
		function A(x)
			local y = self:newState()
			local df_dq_x = df_dq(x)
			local df_dq_t_x = df_dq_t(x)
			-- delta_ij
			y = y + x
			-- - dt * (d/dq_j f_i + d/dq_i f_j)
			y = y - dt * (df_dq_x + df_dq_t_x)
			-- + dt * sum_k (dt * d/dq_i f_k * d/dq_j f_k)
			y = y + dt * df_dq(df_dq_t_x)
			--  - dt * sum_k (g_k * d/dq_i d/dq_j f_k)
			for i=1,self.gridsize do
				if i>1 then
					-- d/d(rho[x+1]) d/d(rho[x+1]) diff_t m[x] = (gamma-3) m[x+1]^2/rho[x+1]^3/(2*dx)
					y[i][1] = y[i][1] - dt * x[i][1] * g[i-1][2] * (gamma-3)*nq[i][2]^2/nq[i][1]^3/(2*dx)
					-- d/d(m[x+1]) d/d(rho[x+1]) diff_t m[x] = -(gamma-3) m[x+1]/rho[x+1]^2/(2*dx)
					y[i][1] = y[i][1] - dt * x[i][2] * g[i-1][2] * -(gamma-3)*nq[i][2]/nq[i][1]^2/(2*dx)
					y[i][2] = y[i][2] - dt * x[i][1] * g[i-1][2] * -(gamma-3)*nq[i][2]/nq[i][1]^2/(2*dx)
					-- d/d(m[x+1]) d/d(m[x+1]) diff_t m[x] = (gamma-3)/rho[x+1]/(2*dx)
					y[i][2] = y[i][2] - dt * x[i][2] * g[i-1][2] * (gamma-3)/nq[i][1]/(2*dx)
					
					-- d/d(rho[x+1]) d/d(rho[x+1]) diff_t E[x] = -(2 gamma E[x+1] m[x+1]/rho[x+1]^3 + 3(1-gamma) m[x+1]^3/rho[x+1]^4)/(2*dx)
					y[i][1] = y[i][1] - dt * x[i][1] * g[i-1][3] * -(2*gamma * nq[i][3]*nq[i][2]/nq[i][1]^3 + 3*(1-gamma)*nq[i][2]^3/nq[i][1]^4)/(2*dx)
					-- d/d(m[x+1]) d/d(rho[x+1]) diff_t E[x] = (gamma E[x+1]/rho[x+1]^2 + 3(1-gamma) m[x+1]^2/rho[x+1]^3)/(2*dx)
					y[i][2] = y[i][2] - dt * x[i][1] * g[i-1][3] * (gamma*nq[i][3]/nq[i][1]^2 + 3*(1-gamma)*nq[i][2]^2/nq[i][1]^3)/(2*dx)
					y[i][1] = y[i][1] - dt * x[i][2] * g[i-1][3] * (gamma*nq[i][3]/nq[i][1]^2 + 3*(1-gamma)*nq[i][2]^2/nq[i][1]^3)/(2*dx)
					-- d/d(E[x+1]) d/d(rho[x+1]) diff_t E[x] = gamma m[x+1]/rho[x+1]^2/(2*dx)
					y[i][3] = y[i][3] - dt * x[i][1] * g[i-1][3] * (gamma*nq[i][2]/nq[i][1]^2)/(2*dx)
					y[i][1] = y[i][1] - dt * x[i][3] * g[i-1][3] * (gamma*nq[i][2]/nq[i][1]^2)/(2*dx)
					-- d/d(m[x+1]) d/d(m[x+1]) diff_t E[x] = 3(gamma-1) m[x+1]/rho[x+1]^2/(2*dx)
					y[i][2] = y[i][2] - dt * x[i][2] * g[i-1][3] * 3*(gamma-1)*nq[i][2]/nq[i][1]^2/(2*dx)
					-- d/d(E[x+1]) d/d(m[x+1]) diff_t E[x] = -gamma/rho[x+1]/(2*dx)
					y[i][3] = y[i][3] - dt * x[i][2] * g[i-1][3] * -gamma/nq[i][1]/(2*dx)
					y[i][2] = y[i][2] - dt * x[i][3] * g[i-1][3] * -gamma/nq[i][1]/(2*dx)
				end
				if i<self.gridsize then
					-- d/d(rho[x-1]) d/d(rho[x-1]) diff_t m[x] = -(gamma-3) m[x-1]^2/rho[x-1]^3/(2*dx)
					y[i][1] = y[i][1] + dt * x[i][1] * g[i+1][2] * (gamma-3)*nq[i][2]^2/nq[i][1]^3/(2*dx)
					-- d/d(m[x-1]) d/d(rho[x-1]) diff_t m[x] = (gamma-3) m[x-1]/rho[x-1]^2/(2*dx)
					y[i][1] = y[i][1] + dt * x[i][2] * g[i+1][2] * -(gamma-3)*nq[i][2]/nq[i][1]^2/(2*dx)
					y[i][2] = y[i][2] + dt * x[i][1] * g[i+1][2] * -(gamma-3)*nq[i][2]/nq[i][1]^2/(2*dx)
					-- d/d(m[x-1]) d/d(m[x-1]) diff_t m[x] = -(gamma-3)/rho[x-1]/(2*dx)
					y[i][2] = y[i][2] + dt * x[i][2] * g[i+1][2] * (gamma-3)/nq[i][1]/(2*dx)
	
					-- d/d(rho[x-1]) d/d(rho[x-1]) diff_t E[x] = (2 gamma E[x-1] m[x-1]/rho[x-1]^3 + 3(1-gamma) m[x-1]^3/rho[x-1]^4)/(2*dx)
					y[i][1] = y[i][1] + dt * x[i][1] * g[i+1][3] * -(2*gamma * nq[i][3]*nq[i][2]/nq[i][1]^3 + 3*(1-gamma)*nq[i][2]^3/nq[i][1]^4)/(2*dx)
					-- d/d(m[x-1]) d/d(rho[x-1]) diff_t E[x] = -(gamma E[x-1]/rho[x-1]^2 + 3(1-gamma) m[x-1]^2/rho[x-1]^3)/(2*dx)
					y[i][2] = y[i][2] + dt * x[i][1] * g[i+1][3] * (gamma*nq[i][3]/nq[i][1]^2 + 3*(1-gamma)*nq[i][2]^2/nq[i][1]^3)/(2*dx)
					y[i][1] = y[i][1] + dt * x[i][2] * g[i+1][3] * (gamma*nq[i][3]/nq[i][1]^2 + 3*(1-gamma)*nq[i][2]^2/nq[i][1]^3)/(2*dx)
					-- d/d(E[x-1]) d/d(rho[x-1]) diff_t E[x] = -gamma m[x-1]/rho[x-1]^2/(2*dx)
					y[i][3] = y[i][3] + dt * x[i][1] * g[i+1][3] * (gamma*nq[i][2]/nq[i][1]^2)/(2*dx)
					y[i][1] = y[i][1] + dt * x[i][3] * g[i+1][3] * (gamma*nq[i][2]/nq[i][1]^2)/(2*dx)
					-- d/d(m[x-1]) d/d(m[x-1]) diff_t E[x] = -3(gamma-1) m[x-1]/rho[x-1]^2/(2*dx)
					y[i][2] = y[i][2] + dt * x[i][2] * g[i+1][3] * 3*(gamma-1)*nq[i][2]/nq[i][1]^2/(2*dx)
					-- d/d(E[x-1]) d/d(m[x-1]) diff_t E[x] = gamma/rho[x-1]/(2*dx)
					y[i][3] = y[i][3] + dt * x[i][2] * g[i+1][3] * -gamma/nq[i][1]/(2*dx)
					y[i][2] = y[i][2] + dt * x[i][3] * g[i+1][3] * -gamma/nq[i][1]/(2*dx)
				end
			end
			return y
		end
		
		local b = g - dt * df_dq(g)	-- b = de/dq_i, df_dq returns a matrix, the state deriv of the time deriv vector
-- [=[ gradient descent
		local step = b
--]=]
--[=[ scale by inverse hessian using conjugate residual ... which is diverging ...
		-- TODO better initial guess for H^-1 * de/dq	
		--local hq = nq:clone()
		local hq = b:clone()
		local r = b - A(hq)
		local Ar = A(r)
		local rAr = dot(r, Ar)
		local p = r:clone()
		local Ap = A(p)
		local cgMaxIter = 50
		for cgIter=1,cgMaxIter do
			local alpha = dot(r, Ar) / norm(Ap)
			hq = hq + alpha * p
			local nr = r - alpha * Ap
			local Anr = A(nr)
			local nrAr = dot(nr, Anr)
			local beta = nrAr / rAr
			rAr = nrAr
			Ar = Anr
			p = r + beta * p
			Ap = Ar + beta * Ap
		end
		-- hq = inverse hessian times gradient
		local step = hq
--]=]
		print('newton iter step norm',norm(step))

		-- TODO stop on min hq 
		nq = nq - step
	end

	self.qs = nq
end

return EulerBackwardEulerNewton

