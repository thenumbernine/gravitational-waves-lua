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


local fluxLimiters = require 'limiter' 
local boundaryMethods = require 'boundary'
local integrators = require 'integrators'

-- here's some globals I have to get rid of

-- meta __index operation -- equivalent of operator[]() in C++
function index(k,v) return k[v] end

-- makes math easier!
for k,v in pairs(math) do _G[k] = v end

local symmath = require 'symmath'

-- solvers:
local HLL = require 'hll'
local Roe = require 'roe'
local PLMBehavior = require 'plm'
local RoeImplicitLinearized = require 'roe_implicit_linearized'
-- equations:
local Maxwell = require 'maxwell'
local Euler1D = require 'euler1d'
local MHD = require 'mhd'
local ADM1D3Var = require'adm1d3var'
local ADM1D3to5Var = require'adm1d3to5var'
local ADM1D5Var = require'adm1d5var'
local BSSNOK1D = require 'bssnok1d'
local ADM2DSpherical = require'adm2dspherical'
local ADM3D = require 'adm3d'

-- setup
local sims = table()

--[[	1D Gaussian curve perturbation / shows coordinate shock waves in 1 direction
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
		fluxLimiter = fluxLimiters.donorCell,
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
		-- A_x = d/dx alpha
		-- D_xxx = 1/2 d/dx gamma_xx
		K_xx = K_xx,	-- K_xx
		-- Bona-Masso slicing conditions:
		f_param = alpha,
		--f = 1,
		--f = 1.69,
		--f = .49,
		f = 1 + kappa/alpha^2,
	}

	-- [=[ compare different equations/formalisms 
	-- these two match:
	sims:insert(Roe(table(args, {equation = ADM1D3Var(equationArgs)})))		-- \_ these two are identical
	sims:insert(Roe(table(args, {equation = ADM1D3to5Var(equationArgs)})))	-- /
	
	-- these two match, but differ from the first two:
	--sims:insert(Roe(table(args, {equation = ADM1D5Var(equationArgs)})))		--> this one, for 1st iter, calcs A_x half what it should
	--sims:insert(Roe(table(args, {equation = ADM3D(equationArgs)})))
	
	-- this is similar to the last two, but off by just a bit (and has an asymmetric evolution of alpha)
	--sims:insert(Roe(table(args, {equation = BSSNOK1D(equationArgs)})))
	
	-- and here's the start of my looking into implicit solvers.  they're broke at the moment:
	--sims:insert(RoeImplicitLinearized(table(args, {equation = ADM1D3to5Var(equationArgs)})))
	--sims:insert(RoeImplicitLinearized(table(args, {equation = ADM1D5Var(equationArgs)})))
	--sims:insert(RoeImplicitLinearized(table(args, {equation = ADM3D(equationArgs)})))
	--sims:insert(require'bssnok1d_backwardeuler_linear'(table(args, equationArgs)))
	--sims:insert(require'bssnok1d_original_backwardeuler_linear'(table(args, equationArgs)))
	--sims:insert(require'bssnok1d_backwardeuler_newton'(args))
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
		fluxLimiter = fluxLimiters.donorCell,
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
		fluxLimiter = fluxLimiters.donorCell,
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
		fluxLimiter = fluxLimiters.donorCell,
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


-- [[	shockwave test via Roe (or Brio-Wu for the MHD simulation)
do
	local args = {
		equation = Euler1D(),
		--stopAtTimes = {.1},
		gridsize = 100,
		domain = {xmin=-1, xmax=1},
		boundaryMethod = boundaryMethods.freeFlow,
		--boundaryMethod = boundaryMethods.freeFlow,
		--linearSolver = require 'linearsolvers'.jacobi,
		--linearSolver = require 'linearsolvers'.conjgrad,
		linearSolver = require 'linearsolvers'.conjres,
		--linearSolver = require 'linearsolvers'.bicgstab,	-- working on this ...
		--fluxLimiter = fluxLimiters.donorCell,
		fluxLimiter = fluxLimiters.superbee,
		integrator = integrators.ForwardEuler,
		--integrator = integrators.RungeKutta4,
		--[=[ TODO broken
		scheme = require 'euler1d_muscl'{
			baseScheme = require 'euler1d_burgers'(),
			slopeLimiter = fluxLimiters.Fromm,
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
	--sims:insert(require 'euler1d_selfsimilar'(table(args, {gridsize=50, domain={xmin=-5, xmax=5}})))
	--sims:insert(Roe(table(args, {usePPM=true})))
	--sims:insert(RoeImplicitLinearized(table(args, {fixed_dt = .01})))
	--sims:insert(require 'euler1d_backwardeuler_newton'(args))
	--sims:insert(require 'euler1d_backwardeuler_linear'(args))
	--sims:insert(require 'euler1d_dft'(args))
	sims:insert(Roe(table(args, {equation = MHD()})))

	-- srhd Marti & Muller 2003 problem #1
	--sims:insert(require 'srhd1d_roe'(table(args, {stopAtTimes={.4249], gridsize=400, domain={xmin=0, xmax=1}, equation=require 'srhd1d'()})))
	-- srhd Marti & Muller 2003 problem #2
	--sims:insert(require 'srhd1d_roe'(table(args, {stopAtTimes={.43}, gridsize=2000, domain={xmin=0, xmax=1}, equation=require 'srhd1d'()})))
	-- srhd Marti & Muller 2003 blast wave interaction
	--sims:insert(require 'srhd1d_roe'(table(args, {stopAtTimes={.1, .26, .426}, gridsize=4000, domain={xmin=0, xmax=1}, equation=require 'srhd1d'()})))
	-- srhd versus euler
	--sims:insert(Roe(table(args, {gridsize=2000})))
	--]=]

	--[=[ compare flux limiters
	sims:insert(Roe(table(args, {fluxLimiter = fluxLimiters.donorCell})))
	sims:insert(Roe(table(args, {fluxLimiter = fluxLimiters.LaxWendroff})))
	sims:insert(Roe(table(args, {fluxLimiter = fluxLimiters.MC})))
	sims:insert(Roe(table(args, {fluxLimiter = fluxLimiters.superbee})))
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

	--[=[ compare constant vs piecewise linear
	args.stopAtTimes = {.4}
	sims:insert(Roe(args))
	local Roe_PLM = class(PLMBehavior(Roe))
	sims:insert(Roe_PLM(args)) 
	--]=]
end
--]]

--[[
sims:insert(Roe{
	equation = Maxwell(),
	gridsize = 200,
	domain = {xmin=-1, xmax=1},
	boundaryMethod = boundaryMethods.freeFlow,
	fluxLimiter = fluxLimiters.donorCell,
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
local GLApp = require 'glapp'
local gl = require 'ffi.OpenGL'
local GLTex2D = require 'gl.tex2d'
local sdl = require 'ffi.sdl'
local Font = require 'gui.font'

local TestApp = class(GLApp)

TestApp.width = 800
TestApp.height = 600
TestApp.showFPS = false

function TestApp:initGL()
	self.doIteration = false
	-- [[ need to get image loading working
	local fonttex = GLTex2D{
		filename = 'font.png',
		minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
		magFilter = gl.GL_LINEAR,
	}
	gl.glGenerateMipmap(gl.GL_TEXTURE_2D)
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
		sims:inf(function(a,b)
			return a.t < b.t
		end):iterate()
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
	
	-- just use the first sim's infos
	for name,info in pairs(sims[1].equation.graphInfoForNames) do
			
		local xmin, xmax, ymin, ymax
		for _,sim in ipairs(sims) do
			sim.ys = {}
			local simymin, simymax
			for i=3,sim.gridsize-2 do
				local siminfo = sim.equation.graphInfoForNames[name]
				if siminfo then
					local y = siminfo.getter(sim,i)
					if not y then 
						--error("failed to get for getter "..info.name)
					else
						sim.ys[i] = y
						if y == y and abs(y) < huge then
							if not simymin or y < simymin then simymin = y end
							if not simymax or y > simymax then simymax = y end
						end
					end
				end
			end
			if self.reportError then
				print(info.name, 'min', simymin, 'max', simymax)
			end

			if not simymin or not simymax or simymin ~= simymin or simymax ~= simymax then
				--simymin = -1
				--simymax = 1
			--elseif abs(simymin) == huge or abs(simymax) == huge then
			else
				local base = 10	-- round to nearest base-10
				local scale = 10 -- ...with increments of 10
				simymin, simymax = 1.1 * simymin - .1 * simymax, 1.1 * simymax - .1 * simymin	
				local newymin = (simymin<0 and -1 or 1)*(abs(simymin)==huge and 1e+100 or base^log(abs(simymin),base))
				local newymax = (simymax<0 and -1 or 1)*(abs(simymax)==huge and 1e+100 or base^log(abs(simymax),base))
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
		
		if not ymin or not ymax or ymin ~= ymin or ymax ~= ymax then
			ymin = -1
			ymax = 1
		end

		-- display
		-- TODO viewports per variable and maybe ticks too
		gl.glViewport(
			info.viewport[1] * w,
			(1 - info.viewport[2] - info.viewport[4]) * h,
			info.viewport[3] * w,
			info.viewport[4] * h)
		gl.glMatrixMode(gl.GL_PROJECTION)
		gl.glLoadIdentity()
		gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)
		gl.glMatrixMode(gl.GL_MODELVIEW)
		gl.glLoadIdentity()

		gl.glColor3f(.1, .1, .1)
		local xrange = xmax - xmin
		local xstep = 10^floor(log(xrange, 10) - .5)
		local xticmin = floor(xmin/xstep)
		local xticmax = ceil(xmax/xstep)
		gl.glBegin(gl.GL_LINES)
		for x=xticmin,xticmax do
			gl.glVertex2f(x*xstep,ymin)
			gl.glVertex2f(x*xstep,ymax)
		end
		gl.glEnd()
		local yrange = ymax - ymin
		local ystep = 10^floor(log(yrange, 10) - .5)
		local yticmin = floor(ymin/ystep)
		local yticmax = ceil(ymax/ystep)
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
				for y=floor(ymin/ystep)*ystep,ceil(ymax/ystep)*ystep,ystep do
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
					text=info.name,
					color = {1,1,1,1},
					fontSize={fontSizeX, -fontSizeY},
					multiLine=false,
				}
			end
		
		end

		gl.glViewport(0,0,w,h)
		gl.glMatrixMode(gl.GL_PROJECTION)
		gl.glLoadIdentity()
		gl.glOrtho(0, w/h, 0, 1, -1, 1)
		gl.glMatrixMode(gl.GL_MODELVIEW)
		gl.glLoadIdentity()

		if self.font then
			local strings = sims:map(function(sim)
				return {
					text = ('(%.3f) '):format(sim.t)..sim.name,
					color = sim.color,
				}
			end)
			local fontSizeX = .02
			local fontSizeY = .02
			local maxlen = strings:map(function(string)
				return self.font:draw{
					text = string.text,
					fontSize = {fontSizeX, -fontSizeY},
					dontRender = true,
					multiLine = false,
				}
			end):inf()
			for i,string in ipairs(strings) do
				self.font:draw{
					pos = {w/h-maxlen,fontSizeY*(i+1)},
					text = string.text,
					color = {string.color[1],string.color[2],string.color[3],1},
					fontSize = {fontSizeX, -fontSizeY},
					multiLine = false,
				}
			end
		end
	end
	if self.reportError then
		self.reportError = false
	end
	gl.glViewport(0,0,w,h)
end
TestApp():run()
--]=]

