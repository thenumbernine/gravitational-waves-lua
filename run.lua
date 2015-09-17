#!/usr/bin/env luajit

local class = require 'ext.class'
local table = require 'ext.table'

local fluxLimiters = require 'limiter' 
local boundaryMethods = require 'boundary'

-- here's some globals I have to get rid of

-- meta __index operation -- equivalent of operator[]() in C++
function index(k,v) return k[v] end

-- makes math easier!
for k,v in pairs(math) do _G[k] = v end


local symmath = require 'symmath'

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
	local g_xx = 1 - h:diff(x)^2
	local K_xx = -h:diff(x,x) / g_xx^.5
	local kappa = 1
	local args = {
		gridsize = 200,
		domain = {xmin=0, xmax=300},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = fluxLimiters.donorCell,
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
		g_xx = g_xx,	-- g_xx
		-- A_x = d/dx alpha
		-- D_xxx = 1/2 d/dx g_xx
		K_xx = K_xx,	-- K_xx
		-- Bona-Masso slicing conditions:
		f_param = alpha,
		--f = 1,
		--f = 1.69,
		--f = .49,
		f = 1 + kappa/alpha^2,
	}

	-- [=[ compare different equations/formalisms 
	--sims:insert(require'adm1d3var_roe'(args))		-- \_ these two are identical
	sims:insert(require'adm1d3to5var_roe'(args))	-- /
	--sims:insert(require'adm1d5var_roe'(args))		--> this one, for 1st iter, calcs A_x half what it should
	--sims:insert(require'bssnok1d_roe'(args))
	--sims:insert(require'adm3d_roe'(args))
	--sims:insert(require'bssnok1d_backwardeuler_newton'(args))
	--]=]

	--[=[
	for _,sim in ipairs(sims) do
		sim.stopAtTime = 100
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
	sims:insert(require'adm2dspherical_roe'{
		gridsize = 200,
		domain = {xmin=100, xmax=500},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = fluxLimiters.donorCell,
		-- the symbolic math driving it:
		r = r,
		h = h,
		g_rr = 1 - h:diff(r)^2,
		g_hh = r^2,
		alpha = 1,
		-- Bona-Masso slicing conditions:
		f_param = alpha,
		f = 1,
		--f = 1.69,
		--f = .49,
		--f = 1 + 1/alpha^2,
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
	sims:insert(require'adm3d'{
		gridsize = 100,
		domain = {xmin=0, xmax=300},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = fluxLimiters.donorCell,
		-- the symbolic math driving it:
		x = x,
		y = y,
		z = z,
		alpha = 1,
		-- g = delta'_ij' - h',i' * h',j',
		g_xx = 1 - hx * hx,
		g_xy = -hx * hy,
		g_xz = -hx * hz,
		g_yy = 1 - hy * hy,
		g_yz = -hy * hz,
		g_zz = 1 - hz * hz,
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
	sims:insert(require'adm3d'{
		gridsize = 100,
		domain = {xmin=-10, xmax=10},
		boundaryMethod = boundaryMethods.freeFlow,
		fluxLimiter = fluxLimiters.donorCell,
		-- the symbolic math driving it:
		x = x,
		y = y,
		z = z,
		alpha = 1,
		-- hmm... beta is important ... so I need to incorporate lapse into my 3D ADM
		-- beta^x = -vs * f(rs(t))
		-- g = delta'_ij'
		g_xx = 1,
		g_xy = 0,
		g_xz = 0,
		g_yy = 1,
		g_yz = 0,
		g_zz = 1,
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
	})
end
--]]


-- [[	shockwave test via Roe (or Brio-Wu for the MHD simulation)
do
	local solverClass = require 'euler1d_roe'
	--local solverClass = require 'euler1d_hll'
	--local solverClass = require 'euler1d_burgers'
	
	local args = {
		gridsize = 200,
		domain = {xmin=-1, xmax=1},
		boundaryMethod = boundaryMethods.mirror,
		--boundaryMethod = boundaryMethods.freeFlow,
		--linearSolver = require 'linearsolvers'.conjres,		-- actually works
		linearSolver = require 'linearsolvers'.conjgrad,	-- not so well
		--linearSolver = require 'linearsolvers'.jacobi,	-- nope
		--fluxLimiter = fluxLimiters.donorCell,
		fluxLimiter = fluxLimiters.superbee,
		--[=[
		scheme = require 'euler1d_muscl'{
			baseScheme = require 'euler1d_burgers'(),
			slopeLimiter = fluxLimiters.Fromm,
		}
		--]=]
	}
	
	-- [=[ compare schemes
	--sims:insert(require 'euler1d_burgers'(args))
	--sims:insert(require 'euler1d_hll'(args))
	--sims:insert(require 'euler1d_roe'(args))
	sims:insert(require 'euler1d_roe_backwardeuler_linear'(args))
	--sims:insert(require 'euler1d_backwardeuler_newton'(args))
	--sims:insert(require 'euler1d_backwardeuler_linear'(args))
	--sims:insert(require 'euler1d_dft'(args))
	--sims:insert(require 'mhd_roe'(args))
	--]=]

	--[=[ compare flux limiters
	sims:insert(solverClass(table(args, {fluxLimiter = fluxLimiters.donorCell})))
	sims:insert(solverClass(table(args, {fluxLimiter = fluxLimiters.LaxWendroff})))
	sims:insert(solverClass(table(args, {fluxLimiter = fluxLimiters.MC})))
	sims:insert(solverClass(table(args, {fluxLimiter = fluxLimiters.superbee})))
	--]=]

	--[=[ compare schemes
	--sims:insert(solverClass(table(args, {scheme = schemes.EulerBurgers()})))
	sims:insert(solverClass(table(args, {scheme = schemes.Roe()})))
	--sims:insert(solverClass(table(args, {scheme = schemes.HLL()})))
	sims:insert(solverClass(table(args, {scheme = schemes.EulerMUSCL{baseScheme = schemes.Roe()}})))
	--sims:insert(solverClass(table(args, {scheme = schemes.EulerMUSCL{baseScheme = schemes.HLL()}})))
	--sims:insert(solverClass(table(args, {scheme = schemes.HLLC()})))	-- TODO 
	--]=]

	--[=[
	for _,sim in ipairs(sims) do
		sim.fixed_dt = 1/512
	end
	--]=]
end
--]]

--[[
sims:insert(require'maxwell'{
	gridsize = 200,
	domain = {xmin=-1, xmax=1},
	boundaryMethod = boundaryMethods.freeFlow,
	fluxLimiter = fluxLimiters.donorCell,
})
--]]


for _,sim in ipairs(sims) do
	do
		local r, g, b = math.random(), math.random(), math.random()
		local l = math.sqrt(r^2 + g^2 + b^2)
		sim.color = {r / l, g / l, b / l}
		print(sim.name, unpack(sim.color))
	end
	sim:reset()
end

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
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
	}
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
		if sim.stopAtTime then
			local t = sim.t
			if self.oldt then
				if t >= sim.stopAtTime and self.oldt < sim.stopAtTime then
					self.doIteration = false
				end
			end
			self.oldt = t
		end
	end

	for _,sim in ipairs(sims) do
		if self.doIteration then
			sim:iterate()
		end
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
			for i=1,sim.gridsize do
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

		gl.glColor3f(.25, .25, .25)
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

		for _,sim in ipairs(sims) do
			gl.glColor3f(unpack(sim.color))
			gl.glPointSize(2)
			if #sim.ys > 0 then
				for _,mode in ipairs{
					--gl.GL_POINTS,
					gl.GL_LINE_STRIP
				} do
					gl.glBegin(mode)
					for i=1,sim.gridsize do
						gl.glVertex2f(sim.xs[i], sim.ys[i])
					end
					gl.glEnd()
				end
			end
			gl.glPointSize(1)
			
			if self.font then
				local fontSizeX = (xmax - xmin) * .05
				local fontSizeY = (ymax - ymin) * .1
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
	end
	if self.reportError then
		self.reportError = false
	end
	gl.glViewport(0,0,w,h)
end
TestApp():run()
--]=]

