#!/usr/bin/env luajit

--DEBUG_PPM=true

require 'ext'

-- stack operations
function first(...) return select(1, ...) end
function last(...) return select(select('#', ...), ...) end
function firstAndLast(...) return first(...), last(...) end
function fill(t, ...)
	for i=1,select('#', ...) do
		t[i] = select(i, ...)
	end
end

-- useful
function assertfinite(x, msg)
	assert(math.isfinite(x), msg or 'is not finite!')
end

local limiter = require 'limiter' 
local boundaryMethods = require 'boundary'
local integrators = require 'integrators'

-- here's some globals I have to get rid of

-- meta __index operation -- equivalent of operator[]() in C++
function index(k,v) return k[v] end

local symmath = require 'symmath'

-- solvers:
local HLL = require 'hll'
local Roe = require 'roe'

local MUSCLBehavior = require 'muscl'
local RoeMUSCL = MUSCLBehavior(Roe)
local HLLMUSCL = MUSCLBehavior(HLL)

local PLMBehavior = require 'plm'
local RoePLM = PLMBehavior(Roe)
local HLLPLM = PLMBehavior(HLL)

local RoeImplicitLinearized = require 'roe_implicit_linearized'
-- equations:
local Maxwell = require 'maxwell'
local Euler1D = require 'euler1d'
local Euler3D = require 'euler3d'
local MHD = require 'mhd'
local EMHD = require 'emhd'
local ADM1D3Var = require'adm1d3var'
local ADM1D3to5Var = require'adm1d3to5var'
local ADM1D5Var = require'adm1d5var'
local BSSNOK1D = require 'bssnok1d'
local ADM2DSpherical = require'adm2dspherical'
local ADM3D = require 'adm3d'
local Z41D = require 'z4-1d'
local Z41Dv2 = require 'z4-1d-v2'

-- setup
local sims = table()

-- [[	1D Gaussian curve perturbation / shows coordinate shock waves in 1 direction
do
	local x = symmath.var'x'
	local alpha = symmath.var'alpha'
	local xc = 150
	local H = 5
	local sigma = 10
	local h = H * symmath.exp(-(x - xc)^2 / sigma^2)
	local gamma_xx = 1 - h:diff(x)^2
	local K_xx = -h:diff(x,x) / gamma_xx^.5
	local kappa = 1
	local args = {
		gridsize = 200,
		domain = {xmin=0, xmax=300},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = limiter.donorCell,
		linearSolver = require 'linearsolvers'.gmres,
		linearSolverEpsilon = 1e-10,
		linearSolverMaxIter = 100,
	}
	local equationArgs = {
		-- the symbolic math driving it:
		x = x,	-- variable
		alpha = 1,
		--[=[
		alpha = 1/2 * (
			(1 + symmath.sqrt(1 + kappa))
			* symmath.sqrt((1-h:diff(x))/(1+h:diff(x)))
			- 
			kappa / (1 + symmath.sqrt(1 + kappa))
			* symmath.sqrt((1+h:diff(x))/(1-h:diff(x)))
		),
		--]=]
		gamma_xx = gamma_xx,	-- gamma_xx
		-- a_x = d/dx alpha
		-- d_xxx = 1/2 d/dx gamma_xx
		K_xx = K_xx,	-- K_xx
		-- Bona-Masso slicing conditions:
		f_param = alpha,
		--f = .49,
		--f = .5,
		f = 1,
		--f = 1.5,
		--f = 1.69,
		--f = 1 + kappa/alpha^2,
	}

	-- [=[ compare different equations/formalisms 
	-- these two match:
	sims:insert(Roe(table(args, {equation = ADM1D3Var(equationArgs)})))		-- \_ these two are identical
	--sims:insert(Roe(table(args, {equation = ADM1D3to5Var(equationArgs)})))	-- /
	-- these two match, but differ from the first two:
	sims:insert(Roe(table(args, {equation = ADM1D5Var(equationArgs)})))		--> this one, for 1st iter, calcs a_x half what it should (or the others calculate it double what it should be ...)
	--sims:insert(Roe(table(args, {equation = ADM3D(equationArgs)})))
	-- this one is similar to the last two, but off by just a bit (and has an asymmetric evolution of alpha)
	--sims:insert(Roe(table(args, {equation = BSSNOK1D(equationArgs)})))
	--sims:insert(Roe(table(args, {equation = Z41D(equationArgs)})))
	--sims:insert(Roe(table(args, {equation = Z41Dv2(equationArgs)})))
	
	-- ... and plm (was working before when I was using the Athena paper implementation, but I broke it when trying to use something more simple):
	--sims:insert(RoePLM(table(args, {equation=ADM1D3Var(equationArgs), fluxLimiter=limiter.donorCell})))
	--sims:insert(RoePLM(table(args, {equation=ADM1D3to5Var(equationArgs), fluxLimiter=limiter.donorCell})))
	--sims:insert(RoePLM(table(args, {equation=ADM1D5Var(equationArgs), fluxLimiter=limiter.donorCell})))
	--sims:insert(RoePLM(table(args, {equation=ADM3D(equationArgs), fluxLimiter=limiter.donorCell})))
	--sims:insert(RoePLM(table(args, {equation=BSSNOK1D(equationArgs), fluxLimiter=limiter.donorCell})))
	--sims:insert(RoePLM(table(args, {equation=Z41D(equationArgs), fluxLimiter=limiter.donorCell})))
	--sims:insert(RoePLM(table(args, {equation=Z41Dv2(equationArgs), fluxLimiter=limiter.donorCell})))
	
	-- and here's the start of my looking into implicit solvers.
	--sims:insert(RoeImplicitLinearized(table(args, {equation = ADM1D3Var(equationArgs)})))
	--sims:insert(RoeImplicitLinearized(table(args, {equation = ADM1D3to5Var(equationArgs)})))
	--sims:insert(RoeImplicitLinearized(table(args, {equation = ADM1D5Var(equationArgs)})))
	--sims:insert(RoeImplicitLinearized(table(args, {equation = ADM3D(equationArgs)})))
	--sims:insert(require'bssnok1d_backwardeuler_linear'(table(args, equationArgs)))
	--sims:insert(require'bssnok1d_original_backwardeuler_linear'(table(args, equationArgs)))
	--sims:insert(require'bssnok1d_backwardeuler_newton'(args))
	--sims:insert(RoeImplicitLinearized(table(args, {equation = Z41D(equationArgs)})))
	--sims:insert(RoeImplicitLinearized(table(args, {equation = Z41Dv2(equationArgs)})))
	--]=]

	--[=[
	for _,sim in ipairs(sims) do
		sim.stopAtTimes = {100}
		sim.fixed_dt = 0.125
	end
	--]=]
end
--]]

--[[	1D profile of 2D spherical Gaussian curve perturbation / coordinate shock wave
do
	-- r - eta(rs) = M ln(((rs + eta(r)) / (rs - eta(rs)))
	-- eta(rs) = sqrt(rs^2 - 2 M rs)
	local rc = 300
	local r = symmath.var'r'
	local alpha = symmath.var'alpha'
	local h = 5 * symmath.exp(-((r - rc) / 10)^2)
	sims:insert(Roe{
		gridsize = 200,
		domain = {xmin=100, xmax=500},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = limiter.donorCell,
		equation = ADM2DSpherical{
			-- the symbolic math driving it:
			r = r,
			h = h,
			gamma_rr = 1 - h:diff(r)^2,
			gamma_hh = r^2,
			alpha = 1,
			-- Bona-Masso slicing conditions:
			f_param = alpha,
			f = 1,
			--f = 1.69,
			--f = .49,
			--f = 1 + 1/alpha^2,
		},
	})
end
--]]

--[[	-- 1D profile of 3D Gaussian curve perturbation / coordinate shock wave in 3 directions
do
	local x = symmath.var'x'
	local y = symmath.var'y'
	local z = symmath.var'z'
	local xc = 150
	local yc = 0
	local zc = 0
	local alpha = symmath.var'alpha'
	local sigma = 10
	local h = 5 * symmath.exp(-((x - xc)^2 + (y - yc)^2 + (z - zc)^2) / sigma^2)
	-- simplifying the math:
	-- dh = h',i'	<- convert expression to rank-1 covariant tensor
	--  this means overriding __newindex for handling tensor derivatives for *all* Expressions
	--  it also means specifying a default coordinate basis
	local hx = h:diff(x)
	local hy = h:diff(y)
	local hz = h:diff(z)
	-- sqrt_det_g = symmath.sqrt(1 - dh'^k' * dh'_k')
	local sqrt_det_g = symmath.sqrt(1 - hx*hx - hy*hy - hz*hz)
	-- d2h = dh',i'	<- convert rank-1 covariant to rank-1 covariant
	local hxx = hx:diff(x)
	local hxy = hx:diff(y)
	local hxz = hx:diff(z)
	local hyy = hy:diff(y)
	local hyz = hy:diff(z)
	local hzz = hz:diff(z)
	sims:insert(Roe{
		gridsize = 100,
		domain = {xmin=0, xmax=300},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = limiter.donorCell,
		equation = ADM3D{
			-- the symbolic math driving it:
			x = x,
			y = y,
			z = z,
			alpha = 1,
			-- g = delta'_ij' - h',i' * h',j',
			gamma_xx = 1 - hx * hx,
			gamma_xy = -hx * hy,
			gamma_xz = -hx * hz,
			gamma_yy = 1 - hy * hy,
			gamma_yz = -hy * hz,
			gamma_zz = 1 - hz * hz,
			-- K = -h',ij' / sqrt_det_g,
			K_xx = -hxx / sqrt_det_g,
			K_xy = -hxy / sqrt_det_g,
			K_xz = -hxz / sqrt_det_g,
			K_yy = -hyy / sqrt_det_g,
			K_yz = -hyz / sqrt_det_g,
			K_zz = -hzz / sqrt_det_g,
			-- Bona-Masso slicing conditions:
			f_param = alpha,	
			--f = 1,
			--f = 1.69,
			--f = .49,
			--f = 1/3,
			f = 1 + 1/alpha^2,
		},
	})
end
--]]

--[[	-- Alcubierre warp drive
do
	local x = symmath.var'x'
	local y = symmath.var'y'
	local z = symmath.var'z'
	local alpha = symmath.var'alpha'
	local xs = 1	-- bubble location?
	local vs = symmath.Constant(1)	-- bubble speed = d/dt xs
	local rs = symmath.sqrt((x - xs)^2 + y^2 + z^2)
	local R = 1	-- parameter
	local sigma = 1 -- parameter
	local f = (symmath.tanh(sigma * (rs + R)) - symmath.tanh(sigma * (rs - R))) / (2 * symmath.tanh(sigma * R))
	sims:insert(Roe{
		gridsize = 100,
		domain = {xmin=-10, xmax=10},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = limiter.donorCell,
		equation = ADM3D{
			-- the symbolic math driving it:
			x = x,
			y = y,
			z = z,
			alpha = 1,
			-- hmm... beta is important ... so I need to incorporate lapse into my 3D ADM
			-- beta^x = -vs * f(rs(t))
			-- g = delta'_ij'
			gamma_xx = 1,
			gamma_xy = 0,
			gamma_xz = 0,
			gamma_yy = 1,
			gamma_yz = 0,
			gamma_zz = 1,
			-- K_ij = -alpha Gamma^t_ij, which I have precomputed for the Alcubierre warp bubble
			K_xx = -(vs * f:diff(x) + vs:diff(x) * f),
			K_xy = -(vs * f:diff(y) + vs:diff(y) * f) / 2,
			K_xz = -(vs * f:diff(z) + vs:diff(z) * f) / 2,
			K_yy = 0,
			K_yz = 0,
			K_zz = 0,
			-- Bona-Masso slicing conditions:
			f_param = alpha,	
			--f = 1,
			--f = 1.69,
			--f = .49,
			--f = 1/3,
			f = 1 + 1/alpha^2,
		},
	})
end
--]]


--[[	shockwave test via Roe (or Brio-Wu for the MHD simulation)
do
	local args = {
		equation = Euler1D(),
		--stopAtTimes = {.1},
		gridsize = 256,
		domain = {xmin=-1, xmax=1},
		boundaryMethod = boundaryMethods.freeFlow,
		--boundaryMethod = boundaryMethods.freeFlow,
		--linearSolver = require 'linearsolvers'.jacobi,
		--linearSolver = require 'linearsolvers'.conjgrad,
		--linearSolver = require 'linearsolvers'.conjres,
		linearSolver = require 'linearsolvers'.gmres,
		--linearSolver = require 'linearsolvers'.bicgstab,	-- working on this ...
		--fluxLimiter = limiter.donorCell,
		fluxLimiter = limiter.superbee,
		integrator = integrators.ForwardEuler,
		--integrator = integrators.RungeKutta4,
		--[=[ TODO broken
		scheme = require 'euler1d_muscl'{
			baseScheme = require 'euler1d_burgers'(),
			slopeLimiter = limiter.Fromm,
		}
		--]=]
	}
	
	-- [=[ compare schemes
	--sims:insert(require 'euler1d_burgers'(args))
	--sims:insert(require 'euler1d_godunov'(table(args, {godunovMethod='exact', sampleMethod='alt'})))
	--sims:insert(require 'euler1d_godunov'(table(args, {godunovMethod='exact'})))
	--sims:insert(require 'euler1d_godunov'(table(args, {godunovMethod='pvrs'})))
	--sims:insert(require 'euler1d_godunov'(table(args, {godunovMethod='twoshock'})))
	--sims:insert(require 'euler1d_godunov'(table(args, {godunovMethod='adaptive'})))
	--sims:insert(HLL(args))
	--sims:insert(Roe(args))
	--sims:insert(RoePLM(table(args, {fluxLimiter=limiter.donorCell})))
	--sims:insert(HLLPLM(args))
	--sims:insert(Roe(table(args, {equation = require 'euler1d_quasilinear'()})))
	--sims:insert(require 'euler1d_selfsimilar'(table(args, {gridsize=50, domain={xmin=-5, xmax=5}})))
	--sims:insert(Roe(table(args, {usePPM=true})))
	--sims:insert(RoeImplicitLinearized(args))
	--sims:insert(RoeImplicitLinearized(table(args, {fixed_dt = .005})))
	--sims:insert(require 'euler1d_backwardeuler_newton'(args))
	--sims:insert(require 'euler1d_backwardeuler_linear'(args))
	--sims:insert(require 'euler1d_dft'(args))
	
	-- mhd:
	--sims:insert(Roe(table(args, {equation=MHD()})))
	--sims:insert(HLL(table(args, {equation=MHD()})))
	--sims:insert(RoePLM(table(args, {equation=MHD(), fluxLimiter=limiter.donorCell})))
	--sims:insert(RoeImplicitLinearized(table(args, {equation=MHD()})))

	-- srhd Marti & Muller 2003 problem #1
	--sims:insert(require 'srhd1d_roe'(table(args, {stopAtTimes={.4249}, gridsize=400, domain={xmin=0, xmax=1}, equation=require 'srhd1d'()})))
	-- srhd Marti & Muller 2003 problem #2
	--sims:insert(require 'srhd1d_roe'(table(args, {stopAtTimes={.43}, gridsize=2000, domain={xmin=0, xmax=1}, equation=require 'srhd1d'()})))
	-- srhd Marti & Muller 2003 blast wave interaction
	--sims:insert(require 'srhd1d_roe'(table(args, {stopAtTimes={.1, .26, .426}, gridsize=4000, domain={xmin=0, xmax=1}, equation=require 'srhd1d'()})))
	-- srhd versus euler
	--sims:insert(Roe(table(args, {gridsize=2000})))
	
	-- problem #1 with PLM ... not working
	--sims:insert(PLMBehavior(require 'srhd1d_roe')(table(args, {stopAtTimes={.4249}, gridsize=400, domain={xmin=0, xmax=1}, equation=require 'srhd1d'()})))
	
	-- emhd
	--sims:insert(Roe(table(args, {equation=EMHD()})))
	--]=]

	--[=[ compare various flux limiters
	sims:insert(Roe(table(args, {fluxLimiter = limiter.donorCell})))
	sims:insert(Roe(table(args, {fluxLimiter = limiter.LaxWendroff})))
	sims:insert(Roe(table(args, {fluxLimiter = limiter.MC})))
	sims:insert(Roe(table(args, {fluxLimiter = limiter.superbee})))
	--]=]

	--[=[ compare flux vs plm slope limiter
	--sims:insert(require 'sod_exact'(table(args, {gridsize=2000})))
	sims:insert(Roe(table(args, {fluxLimiter=limiter.superbee})))
	sims:insert(Roe(table(args, {fluxLimiter=limiter.donorCell}))) 	-- eliminate the flux limiter, so only the PLM slope limiter is applied
	sims:insert(RoePLM(table(args, {fluxLimiter=limiter.donorCell}))) 	-- eliminate the flux limiter, so only the PLM slope limiter is applied
	--sims:insert(RoePLM(table(args, {fluxLimiter=limiter.superbee}))) 	-- this is applying the limiter twice: the flux and the slope
	--sims:insert(RoeMUSCL(table(args, {fluxLimiter=limiter.superbee}))) 	-- this is applying the limiter twice: the flux and the slope
	--sims:insert(RoeMUSCL(table(args, {fluxLimiter=limiter.donorCell}))) 	-- eliminate the flux limiter, so only the MUSCL slope limiter is applied
	--]=]

	--[=[ compare schemes
	sims:insert(require 'euler1d_burgers'(args))
	sims:insert(Roe(args))
	sims:insert(HLL(args))
	--sims:insert(require 'euler1d_muscl'(table(args, {baseScheme = Roe(args)})))	-- TODO I broke this
	--sims:insert(require 'euler1d_muscl'(table(args, {baseScheme = HLL()})))	-- TODO I broke this
	--sims:insert(Roe(table(args, {scheme = schemes.HLLC()})))	-- TODO 
	--]=]

	--[=[ compare integrators
	args.stopAtTimes = {.4}
	sims:insert(Roe(table(args, {integrator=integrators.ForwardEuler()})))
	sims:insert(Roe(table(args, {integrator=integrators.RungeKutta4()})))
	--]=]
end
--]]

--[[
sims:insert(Roe{
	equation = Maxwell(),
	gridsize = 200,
	domain = {xmin=-1, xmax=1},
	boundaryMethod = boundaryMethods.freeFlow,
	fluxLimiter = limiter.donorCell,
})
--]]


for _,sim in ipairs(sims) do
	do
		local r, g, b = math.random(), math.random(), math.random()
		--local l = math.sqrt(r^2 + g^2 + b^2)
		local l = math.max(r,g,b)
		sim.color = {r / l, g / l, b / l}
		--print(sim.name, unpack(sim.color))
	end
	sim:reset()
end
if #sims >= 1 then sims[1].color = {.2,1,1} end
if #sims >= 2 then sims[2].color = {1,.4,1} end

--[[ text
local printState = function()
	for _,sim in ipairs(sims) do
		for infoIndex,info in ipairs(sim.graphInfos) do
			io.write(info.name)
			for i=1,sim.gridsize do
				local y = info.getter(i)
				io.write('\t', y)
			end
			print()
		end
	end
end
printState()
for iter=1,1 do
	for _,sim in ipairs(sims) do
		sim:iterate()
	end
	printState()
end
os.exit()
--]]

-- [=[ graphics
local ffi = require 'ffi'
local gl = require 'ffi.OpenGL'
local sdl = require 'ffi.sdl'
local GLTex2D = require 'gl.tex2d'
local Font = require 'gui.font'

-- [[ with ImGui
local ig = require 'ffi.imgui'
local ImGuiApp = require 'imguiapp'
local TestApp = class(ImGuiApp)
--]]
--[[ disable ImGui
local ig = setmetatable({
}, {
	__index = function() 
		return function() end
	end,
})
local GLApp = require 'glapp'
local TestApp = class(GLApp)
--]]


TestApp.width = 800
TestApp.height = 600
TestApp.showFPS = false

function TestApp:initGL(...)
	if TestApp.super.initGL then TestApp.super.initGL(self, ...) end
	self.doIteration = false
	-- [[ need to get image loading working
	local fonttex = GLTex2D{
		filename = 'font.png',
		minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
		magFilter = gl.GL_LINEAR,
	}
	if not pcall(function()
		gl.glGenerateMipmap(gl.GL_TEXTURE_2D) 
	end) then
		gl.glTexParameteri(fonttex.target, gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST)
		gl.glTexParameteri(fonttex.target, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR) 
	end
	self.font = Font{tex=fonttex}
	--]]
end

function TestApp:event(event)
	if event.type == sdl.SDL_KEYDOWN then
		local callback = self.keyDownCallbacks[event.key.keysym.sym]
		if callback then
			callback(self)
		end
	end
end

TestApp.keyDownCallbacks = {
	[sdl.SDLK_r] = function(self)
		for _,sim in ipairs(sims) do
			sim:reset()
		end
	end,
	[sdl.SDLK_e] = function(self)
		self.reportError = not self.reportError	
	end,
	[sdl.SDLK_SPACE] = function(self)
		self.doIteration = not self.doIteration
	end,
	[sdl.SDLK_u] = function(self)
		self.doIteration = 'once'
	end,
}

local graphNamesEnabled = table()
graphNamesEnabled:insert{
	name = 'all',
	ptr = ffi.new('bool[1]', true),
}

local solverGens 
do
	local numRel1DArgs	-- used for NR simulations
	do
		local x = symmath.var'x'
		local alpha = symmath.var'alpha'
		local xc = 150
		local H = 5
		local sigma = 10
		local h = H * symmath.exp(-(x - xc)^2 / sigma^2)
		local gamma_xx = 1 - h:diff(x)^2
		local K_xx = -h:diff(x,x) / gamma_xx^.5
		local kappa = 1
		local args = {
			gridsize = 200,
			domain = {xmin=0, xmax=300},
			boundaryMethod = boundaryMethods.freeFlow,
			fluxLimiter = limiter.donorCell,
			linearSolver = require 'linearsolvers'.gmres,
			linearSolverEpsilon = 1e-10,
			linearSolverMaxIter = 100,
		}
		numRel1DArgs = {
			-- the symbolic math driving it:
			x = x,	-- variable
			alpha = 1,
			--[=[
			alpha = 1/2 * (
				(1 + symmath.sqrt(1 + kappa))
				* symmath.sqrt((1-h:diff(x))/(1+h:diff(x)))
				- 
				kappa / (1 + symmath.sqrt(1 + kappa))
				* symmath.sqrt((1+h:diff(x))/(1-h:diff(x)))
			),
			--]=]
			gamma_xx = gamma_xx,	-- gamma_xx
			-- a_x = d/dx alpha
			-- d_xxx = 1/2 d/dx gamma_xx
			K_xx = K_xx,	-- K_xx
			-- Bona-Masso slicing conditions:
			f_param = alpha,
			--f = 1,
			--f = 1.69,
			--f = .49,
			f = 1 + kappa/alpha^2,
		}
	end

	local numRel2DSphericalArgs
	do
		-- r - eta(rs) = M ln(((rs + eta(r)) / (rs - eta(rs)))
		-- eta(rs) = sqrt(rs^2 - 2 M rs)
		local rc = 300
		local r = symmath.var'r'
		local alpha = symmath.var'alpha'
		local h = 5 * symmath.exp(-((r - rc) / 10)^2)	
		numRel2DSphericalArgs = {
			-- the symbolic math driving it:
			r = r,
			h = h,
			gamma_rr = 1 - h:diff(r)^2,
			gamma_hh = r^2,
			alpha = 1,
			-- Bona-Masso slicing conditions:
			f_param = alpha,
			f = 1,
			--f = 1.69,
			--f = .49,
			--f = 1 + 1/alpha^2,
		}
	end

	local numRel3DArgs
	do
		local x = symmath.var'x'
		local y = symmath.var'y'
		local z = symmath.var'z'
		local xc = 150
		local yc = 0
		local zc = 0
		local alpha = symmath.var'alpha'
		local sigma = 10
		local h = 5 * symmath.exp(-((x - xc)^2 + (y - yc)^2 + (z - zc)^2) / sigma^2)
		-- simplifying the math:
		-- dh = h',i'	<- convert expression to rank-1 covariant tensor
		--  this means overriding __newindex for handling tensor derivatives for *all* Expressions
		--  it also means specifying a default coordinate basis
		local hx = h:diff(x)
		local hy = h:diff(y)
		local hz = h:diff(z)
		-- sqrt_det_g = symmath.sqrt(1 - dh'^k' * dh'_k')
		local sqrt_det_g = symmath.sqrt(1 - hx*hx - hy*hy - hz*hz)
		-- d2h = dh',i'	<- convert rank-1 covariant to rank-1 covariant
		local hxx = hx:diff(x)
		local hxy = hx:diff(y)
		local hxz = hx:diff(z)
		local hyy = hy:diff(y)
		local hyz = hy:diff(z)
		local hzz = hz:diff(z)
		numRel3DArgs = {
			-- the symbolic math driving it:
			x = x,
			y = y,
			z = z,
			alpha = 1,
			-- g = delta'_ij' - h',i' * h',j',
			gamma_xx = 1 - hx * hx,
			gamma_xy = -hx * hy,
			gamma_xz = -hx * hz,
			gamma_yy = 1 - hy * hy,
			gamma_yz = -hy * hz,
			gamma_zz = 1 - hz * hz,
			-- K = -h',ij' / sqrt_det_g,
			K_xx = -hxx / sqrt_det_g,
			K_xy = -hxy / sqrt_det_g,
			K_xz = -hxz / sqrt_det_g,
			K_yy = -hyy / sqrt_det_g,
			K_yz = -hyz / sqrt_det_g,
			K_zz = -hzz / sqrt_det_g,
			-- Bona-Masso slicing conditions:
			f_param = alpha,	
			--f = 1,
			--f = 1.69,
			--f = .49,
			--f = 1/3,
			f = 1 + 1/alpha^2,
		}
	end

	solverGens = table{
		-- accepts varying equations
			-- TODO PLMBehavior applied to any of those
			-- Euler1D:
		{name='Euler1D HLL', gen=function(args) return HLL(table(args, {equation=Euler1D()})) end},
		{name='Euler1D Roe', gen=function(args) return Roe(table(args, {equation=Euler1D()})) end},
		{name='Euler1D Roe Implicit Linearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=Euler1D()})) end},
			-- MHD
		{name='MHD HLL', gen=function(args) return HLL(table(args, {equation=MHD()})) end},
		{name='MHD Roe', gen=function(args) return Roe(table(args, {equation=MHD()})) end},
		{name='MHD Roe Implicit Linearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=MHD()})) end},
			-- Maxwell
		{name='Maxwell HLL', gen=function(args) return HLL(table(args, {equation=Maxwell()})) end},
		{name='Maxwell Roe', gen=function(args) return Roe(table(args, {equation=Maxwell()})) end},
		{name='Maxwell Roe Implicit Linearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=Maxwell()})) end},
		-- equations that only work with Roe:
		{name='SRHD Roe', gen=function(args) return require 'srhd1d_roe'(table(args, {equation=require'srhd1d'()})) end},
		{name='ADM 1D 3-var Roe', gen=function(args) return Roe(table(args, {equation=ADM1D3Var(numRel1DArgs)})) end},
		{name='ADM 1D 3-to-5-var Roe', gen=function(args) return Roe(table(args, {equation=ADM1D3to5Var(numRel1DArgs)})) end},
		{name='ADM 1D 5-var Roe', gen=function(args) return Roe(table(args, {equation=ADM1D5Var(numRel1DArgs)})) end},
		{name='BSSNOK 1D Roe', gen=function(args) return Roe(table(args, {equation=BSSNOK1D(numRel1DArgs)})) end},
		{name='ADM 2D Spherical Roe', gen=function(args) return Roe(table(args, {equation=ADM2DSpherical(numRel2DSphericalArgs)})) end},
		{name='ADM 3D Roe', gen=function(args) return Roe(table(args, {equation=ADM3D(numRel3DArgs)})) end},
			-- ... and their implicit linearized versions ...
		{name='ADM 1D 3-var RoeImplicitLinearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=ADM1D3Var(numRel1DArgs)})) end},
		{name='ADM 1D 3-to-5-var RoeImplicitLinearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=ADM1D3to5Var(numRel1DArgs)})) end},
		{name='ADM 1D 5-var RoeImplicitLinearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=ADM1D5Var(numRel1DArgs)})) end},
		{name='BSSNOK 1D RoeImplicitLinearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=BSSNOK1D(numRel1DArgs)})) end},
		{name='ADM 2D Spherical RoeImplicitLinearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=ADM2DSpherical(numRel2DSphericalArgs)})) end},
		{name='ADM 3D RoeImplicitLinearized', gen=function(args) return RoeImplicitLinearized(table(args, {equation=ADM3D(numRel3DArgs)})) end},
		-- stuff I've never finished:
		{name='BSSNOK 1D Backward Euler Linear', gen=function(args) return require 'bssnok1d_backwardeuler_linear'(table(args, numRel1DArgs)) end},
		{name='BSSNOK 1D Original Backward Euler Linear', gen=function(args) return require 'bssnok1d_original_backwardeuler_linear'(table(args, numRel1DArgs)) end},
		{name='BSSNOK 1D Backward Euler Newton', gen=function(args) return require 'bssnok1d_backwardeuler_newton'(table(args, numRel1DArgs)) end},
		-- only implemented for Euler1D
		{name='Euler1D Burgers', gen=require 'euler1d_burgers'},
		{name='Euler1D Self-Similar', gen=require 'euler1d_selfsimilar'},
		{name='Euler1D Godunov Exact', gen=function(args) return require 'euler1d_godunov'(table(args, {godunovMethod='exact'})) end},
		{name='Euler1D Godunov Exact Alt Sampling', gen=function(args) return require 'euler1d_godunov'(table(args, {godunovMethod='exact', sampleMethod='alt'})) end},
		{name='Euler1D Godunov PVRS', gen=function(args) return require 'euler1d_godunov'(table(args, {godunovMethod='pvrs'})) end},
		{name='Euler1D Godunov Two-Shock', gen=function(args) return require 'euler1d_godunov'(table(args, {godunovMethod='twoshock'})) end},
		{name='Euler1D Godunov Adaptive', gen=function(args) return require 'euler1d_godunov'(table(args, {godunovMethod='adaptive'})) end},
		{name='Euler1D Backwards Euler Newton', gen=require 'euler1d_backwardeuler_newton'},
		{name='Euler1D Backwards Euler Linear', gen=require 'euler1d_backwardeuler_linear'},
		{name='Euler1D DFT', gen=require 'euler1d_dft'},
	}
end
local solverGenIndex = ffi.new('int[1]', 0)

function TestApp:updateGUI()
	ig.igText('simulations:')
	local toRemove
	for i,sim in ipairs(sims) do
		ig.igPushIdStr(i..' '..sim.name)
		sim.visiblePtr = sim.visiblePtr or ffi.new('bool[1]', true)
		ig.igPushIdStr('visible')
		ig.igCheckbox('', sim.visiblePtr)
		ig.igPopId()
		ig.igSameLine()
		if ig.igButton('X') then
			toRemove = toRemove or table()
			toRemove:insert(1,i)	-- insert last-to-first
		end
		ig.igSameLine()
		ig.igText(sim.name)
		ig.igPopId()
	end
	if toRemove then
		for _,i in ipairs(toRemove) do
			sims:remove(i)
		end
	end

	ig.igCombo('new sim type', solverGenIndex, solverGens:map(function(solverGen) return solverGen.name end))
	if ig.igButton('Add New...') then
		print('adding new sim '..solverGenIndex[0]..' named '..solverGens[solverGenIndex[0]+1].name)
		local sim = solverGens[solverGenIndex[0]+1].gen{
			gridsize = 100,
			domain = {xmin=-1, xmax=1},
			boundaryMethod = boundaryMethods.freeFlow,
			linearSolver = require 'linearsolvers'.gmres,
			--linearSolver = require 'linearsolvers'.conjres,
			fluxLimiter = limiter.superbee,
			integrator = integrators.ForwardEuler,
		}
		local r, g, b = math.random(), math.random(), math.random()
		--local l = math.sqrt(r^2 + g^2 + b^2)
		local l = math.max(r,g,b)
		sim.color = {r / l, g / l, b / l}
		sim:reset()
		sims:insert(sim)
	end

	for _, graphNameEnabled in ipairs(graphNamesEnabled) do
		graphNameEnabled.found = false
	end
	graphNamesEnabled[1].found = true
	for _,sim in ipairs(sims) do
		for _,graphInfo in ipairs(sim.equation.graphInfos) do
			local _, graphNameEnabled = graphNamesEnabled:find(nil, function(graphName)
				return graphName.name == graphInfo.name
			end)
			if graphNameEnabled then
				graphNameEnabled.found = true
			else
				graphNamesEnabled:insert{
					name = graphInfo.name,
					ptr = ffi.new('bool[1]', true),
					found = true,
				}
			end
		end
	end
	for i=#graphNamesEnabled,1,-1 do
		if not graphNamesEnabled[i].found then
			graphNamesEnabled:remove(i)
		end
	end
	--[[ if you do sort, keep 'all' at the top
	graphNamesEnabled = graphNamesEnabled:sort(function(a,b)
		if #a.name ~= #b.name then return #a.name < #b.name end
		return a.name < b.name
	end)
	--]]

	local allBefore = graphNamesEnabled[1].ptr[0]
	ig.igText('variables:')
	for i,graphNameEnabled in ipairs(graphNamesEnabled) do
		ig.igPushIdStr(tostring(i))
		ig.igCheckbox(graphNameEnabled.name, graphNameEnabled.ptr)
		ig.igPopId()
	end
	local allAfter = graphNamesEnabled[1].ptr[0]
	if allAfter ~= allBefore then
		for i=2,#graphNamesEnabled do
			graphNamesEnabled[i].ptr[0] = allAfter
		end
	end

	if self.doIteration then
		if ig.igButton('pause') then
			self.doIteration = nil
		end
	else
		if ig.igButton('start') then
			self.doIteration = true
		end
	end
	if ig.igButton('step') then
		self.doIteration = 'once'
	end
	if ig.igButton('reset') then
		for _,sim in ipairs(sims) do
			sim:reset()
		end
	end
end

function TestApp:update(...)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	
	for _,sim in ipairs(sims) do
		if sim.stopAtTimes then
			local t = sim.t
			if self.oldt then
				for _,t2 in ipairs(sim.stopAtTimes) do
					if t >= t2 and self.oldt < t2 then
						self.doIteration = false
					end
				end
			end
			self.oldt = t
		end
	end

	if self.doIteration then
		-- [[ iterate the furthest back
		local oldestSim = sims:inf(function(a,b)
			return a.t < b.t
		end)
		if oldestSim then oldestSim:iterate() end
		--]]
		--[[ iterate all
		for _,sim in ipairs(sims) do
			sim:iterate()
		end
		--]]
	end

	if self.doIteration == 'once' then
		self.doIteration = false
	end

	local w, h = self:size()

	local numEnabled = 0
	for i=2,#graphNamesEnabled do
		if graphNamesEnabled[i].ptr[0] then
			numEnabled = numEnabled + 1
		end
	end
	local graphsWide = math.ceil(math.sqrt(numEnabled))
	local graphsHigh = math.ceil(numEnabled/graphsWide)
	local graphCol = 0
	local graphRow = 0

	-- just use the first sim's infos
	--for name,info in pairs(sims[1].equation.graphInfoForNames) do
	if #sims > 0 then
		for j=2,#graphNamesEnabled do
			local graphNameEnabled = graphNamesEnabled[j]
			if graphNameEnabled.ptr[0] then		
				local name = graphNameEnabled.name

				local xmin, xmax, ymin, ymax
				for _,sim in ipairs(sims) do
					if sim.visiblePtr and sim.visiblePtr[0] then
						sim.ys = {}
						local simymin, simymax
						for i=3,sim.gridsize-2 do
							local siminfo = sim.equation.graphInfoForNames[name]
							if siminfo then
								
								local y = siminfo.getter(sim,i)
								if not y then 
									--error("failed to get for getter "..name)
								else
									sim.ys[i] = y
									if y == y and math.abs(y) < math.huge then
										if not simymin or y < simymin then simymin = y end
										if not simymax or y > simymax then simymax = y end
									end
								end
							end
						end
						if self.reportError then
							print(name, 'min', simymin, 'max', simymax)
						end

						if not simymin or not simymax or simymin ~= simymin or simymax ~= simymax then
							--simymin = -1
							--simymax = 1
						--elseif math.abs(simymin) == math.huge or math.abs(simymax) == math.huge then
						else
							local base = 10	-- round to nearest base-10
							local scale = 10 -- ...with increments of 10
							simymin, simymax = 1.1 * simymin - .1 * simymax, 1.1 * simymax - .1 * simymin	
							local newymin = (simymin<0 and -1 or 1)*(math.abs(simymin)==math.huge and 1e+100 or base^math.log(math.abs(simymin),base))
							local newymax = (simymax<0 and -1 or 1)*(math.abs(simymax)==math.huge and 1e+100 or base^math.log(math.abs(simymax),base))
							simymin, simymax = newymin, newymax
							do
								local minDeltaY = 1e-5
								local deltaY = simymax - simymin
								if deltaY < minDeltaY then
									simymax = simymax + .5 * minDeltaY
									simymin = simymin - .5 * minDeltaY
								end
							end
						end

						local simxmin, simxmax = sim.domain.xmin, sim.domain.xmax
						simxmin, simxmax = 1.1 * simxmin - .1 * simxmax, 1.1 * simxmax - .1 * simxmin

						xmin = xmin or simxmin
						xmax = xmax or simxmax
						ymin = ymin or simymin
						ymax = ymax or simymax
							
						if xmin and simxmin then xmin = math.min(xmin, simxmin) end
						if xmax and simxmax then xmax = math.max(xmax, simxmax) end
						if ymin and simymin then ymin = math.min(ymin, simymin) end
						if ymax and simymax then ymax = math.max(ymax, simymax) end
					end
				end
				
				if not xmin or not xmax or xmin ~= xmin or xmax ~= xmax then
					xmin = -1
					xmax = 1
				end
				if not ymin or not ymax or ymin ~= ymin or ymax ~= ymax then
					ymin = -1
					ymax = 1
				end

				-- display
				-- TODO viewports per variable and maybe ticks too
				gl.glViewport(
					graphCol / graphsWide * w,
					(1 - (graphRow + 1) / graphsHigh) * h,
					w / graphsWide,
					h / graphsHigh)
				gl.glMatrixMode(gl.GL_PROJECTION)
				gl.glLoadIdentity()
				gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)
				gl.glMatrixMode(gl.GL_MODELVIEW)
				gl.glLoadIdentity()

				gl.glColor3f(.1, .1, .1)
				local xrange = xmax - xmin
				local xstep = 10^math.floor(math.log(xrange, 10) - .5)
				local xticmin = math.floor(xmin/xstep)
				local xticmax = math.ceil(xmax/xstep)
				gl.glBegin(gl.GL_LINES)
				for x=xticmin,xticmax do
					gl.glVertex2f(x*xstep,ymin)
					gl.glVertex2f(x*xstep,ymax)
				end
				gl.glEnd()
				local yrange = ymax - ymin
				local ystep = 10^math.floor(math.log(yrange, 10) - .5)
				local yticmin = math.floor(ymin/ystep)
				local yticmax = math.ceil(ymax/ystep)
				gl.glBegin(gl.GL_LINES)
				for y=yticmin,yticmax do
					gl.glVertex2f(xmin,y*ystep)
					gl.glVertex2f(xmax,y*ystep)
				end
				gl.glEnd()
					
				gl.glColor3f(.5, .5, .5)
				gl.glBegin(gl.GL_LINES)
				gl.glVertex2f(xmin, 0)
				gl.glVertex2f(xmax, 0)
				gl.glVertex2f(0, ymin)
				gl.glVertex2f(0, ymax)
				gl.glEnd()
			
				-- should I show ghost cells? for some derived values it causes errors...
				for _,sim in ipairs(sims) do
					if sim.visiblePtr and sim.visiblePtr[0] then
						gl.glColor3f(unpack(sim.color))
						gl.glPointSize(2)
						if #sim.ys > 0 then
							for _,mode in ipairs{
								gl.GL_LINE_STRIP,
								DEBUG_PPM and gl.GL_POINTS
							} do
								gl.glBegin(mode)
								for i=3,sim.gridsize-2 do
									gl.glVertex2f(sim.xs[i], sim.ys[i])
								end
								gl.glEnd()
							end
						end
			-- [[ special PPM hack
			if DEBUG_PPM then
						local channel = 2
						local ppmCount = 0
						local ppmYs = table()
						gl.glColor3f(1,1,0)
						gl.glBegin(gl.GL_LINE_STRIP)
						for n=0,#sim.xs*10 do
							local x = (xmax - xmin) / (#sim.xs*10) * n + xmin
							-- getter ... abstracts the index of the state variable ...
							local y = sim:getPPM(x,channel)
							if y then
								ppmYs:insert(y)
								ppmCount = ppmCount + 1
								gl.glVertex2f(x,y)
							end
						end
						gl.glEnd()
						--print(unpack(ppmYs,1,10))
						--print(ppmCount) 
						if sim.ppm_iqs then
							gl.glColor3f(0,1,0)
							gl.glBegin(gl.GL_LINES)
							for i=3,sim.gridsize-2 do
								gl.glVertex2f(sim.ixs[i], sim.ppm_qLs[i][channel])
								gl.glVertex2f(sim.ixs[i+1], sim.ppm_qRs[i][channel])
							end
							gl.glEnd()
							gl.glColor3f(1,0,1)
							gl.glBegin(gl.GL_POINTS)
							for i=3,sim.gridsize-1 do
								gl.glVertex2f(sim.ixs[i],sim.ppm_iqs[i][channel])
							end
							gl.glEnd()
						end
			end
			--]]
						gl.glPointSize(1)
						
						if self.font then
							local fontSizeX = (xmax - xmin) * .05
							local fontSizeY = (ymax - ymin) * .05
							local ystep = ystep * 2
							for y=math.floor(ymin/ystep)*ystep,math.ceil(ymax/ystep)*ystep,ystep do
								self.font:draw{
									pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
									text=tostring(y),
									color = {1,1,1,1},
									fontSize={fontSizeX, -fontSizeY},
									multiLine=false,
								}
							end
							self.font:draw{
								pos={xmin, ymax},
								text=name,
								color = {1,1,1,1},
								fontSize={fontSizeX, -fontSizeY},
								multiLine=false,
							}
						end
					end
				end

				gl.glViewport(0,0,w,h)
				gl.glMatrixMode(gl.GL_PROJECTION)
				gl.glLoadIdentity()
				gl.glOrtho(0, w/h, 0, 1, -1, 1)
				gl.glMatrixMode(gl.GL_MODELVIEW)
				gl.glLoadIdentity()

				if self.font then
					local simNames = sims:map(function(sim)
						return {
							text = ('(%.3f) '):format(sim.t)..sim.name,
							color = sim.color,
						}
					end)
					local fontSizeX = .02
					local fontSizeY = .02
					local maxlen = simNames:map(function(simName)
						return self.font:draw{
							text = simName.text,
							fontSize = {fontSizeX, -fontSizeY},
							dontRender = true,
							multiLine = false,
						}
					end):inf()
					for i,simName in ipairs(simNames) do
						self.font:draw{
							pos = {w/h-maxlen,fontSizeY*(i+1)},
							text = simName.text,
							color = {simName.color[1], simName.color[2], simName.color[3],1},
							fontSize = {fontSizeX, -fontSizeY},
							multiLine = false,
						}
					end
				end
			
				graphCol = graphCol + 1
				if graphCol == graphsWide then
					graphCol = 0
					graphRow = graphRow + 1
				end
			end
		end
		if self.reportError then
			self.reportError = false
		end
		gl.glViewport(0,0,w,h)
	end	

	if TestApp.super.update then
		return TestApp.super.update(self, ...)
	else
		self:updateGUI()
	end
end

TestApp():run()
--]=]

