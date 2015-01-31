#!/usr/bin/env luajit

require 'ext'
local slopeLimiters = require 'relativity.limiter' 
local boundaryMethods = require 'relativity.boundary'

-- here's some globals I have to get rid of

-- meta __index operation -- equivalent of operator[]() in C++
function index(k,v) return k[v] end

-- makes math easier!
for k,v in pairs(math) do _G[k] = v end


local ADM1D3VarSim = require 'relativity.adm1d3var'
local ADM1D5VarSim = require 'relativity.adm1d5var'
local ADM2D = require 'relativity.adm2d'
local EulerSim = require 'relativity.euler'
local symmath = require 'symmath'

-- setup
-- [[
--local sim = ADM1D3VarSim{
local sim
do
	local sigma = 10
	local xc = 150
	local H = 5

	-- dependent vars
	local x = symmath.var'x'
	local alpha = symmath.var'alpha'
	
	local h = H * symmath.exp(-(x - xc)^2 / sigma^2)
	sim = ADM1D5VarSim{
		gridsize = 1200,
		domain = {xmin=0, xmax=300},
		boundaryMethod = boundaryMethods.freeFlow,	-- still reflecting despite freeflow ...
		slopeLimiter = slopeLimiters.donorCell,
		-- the symbolic math driving it:
		x = x,
		h = h,
		g = 1 - h:diff(x)^2,
		alpha = 1,
		-- Bona-Masso slicing conditions:
		f_var = alpha,	
		f = 1,
		--f = 1.69,
		--f = .49,
		--f = (1 + 1/alpha^2),
	}
end
sim.fixed_dt = nil	--.125 
sim.stopAtTime = 70
--]]

--[[
local sim = EulerSim{
	gridsize = 200,
	domain = {xmin=-1, xmax=1},
	boundaryMethod = boundaryMethods.mirror,
	slopeLimiter = slopeLimiters.superbee,
}
--]]

sim:reset()


--[[ text
for iter=1,100 do
	sim:iterate()
	for j=1,sim.numStates do
		io.write(({'alpha','g','A','D','K'})[j])
		for i=1,sim.gridsize do
			io.write('\t', qs[i][j])
		end
		print()
	end
end
--]]

-- [=[ graphics
local GLApp = require 'glapp'
local gl = require 'ffi.OpenGL'
local sdl = require 'ffi.sdl'


local TestApp = class(GLApp)

TestApp.width = 640
TestApp.height = 640
TestApp.showFPS = false

function TestApp:init()
	TestApp.super.init(self)
	self.doIteration = false
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
		sim:reset()
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

	if sim.stopAtTime then
		local t = sim.t
		if self.oldt then
			if t >= sim.stopAtTime and self.oldt < sim.stopAtTime then
				self.doIteration = false
			end
		end
		self.oldt = t
	end

	if self.doIteration then
		if self.doIteration == 'once' then
			self.doIteration = false
		end
		sim:iterate()
	end
	
	local xs = sim.xs
	local w, h = self:size()
	for infoIndex,info in ipairs(sim.graphInfos) do
		local ys = {}
		local ymin, ymax
		for i=1,sim.gridsize do
			local y = info.getter(i)
			ys[i] = y
			if y == y and abs(y) < huge then
				if not ymin or y < ymin then ymin = y end
				if not ymax or y > ymax then ymax = y end
			end
		end
		if self.reportError then
			print(info.name, 'min', ymin, 'max', ymax)
		end
	
		if info.range then
			ymin, ymax = unpack(info.range)
		else
			if not ymin or not ymax or ymin ~= ymin or ymax ~= ymax then
				ymin = -1
				ymax = 1
			--elseif abs(ymin) == huge or abs(ymax) == huge then
			else
				local base = 10	-- round to nearest base-10
				local scale = 10 -- ...with increments of 10
				ymin = (ymin<0 and -1 or 1)*(abs(ymin)==huge and 1e+100 or base^((log(abs(1.1*ymin-.1*ymax),base)*scale)/scale))
				ymax = (ymax<0 and -1 or 1)*(abs(ymax)==huge and 1e+100 or base^((log(abs(1.1*ymax-.1*ymin),base)*scale)/scale))
			end
		end
		
		local xmin, xmax = sim.domain.xmin, sim.domain.xmax
		xmin, xmax = 1.1 * xmin - .1 * xmax, 1.1 * xmax - .1 * xmin	
		
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

-- [[
		gl.glColor3f(.25, .25, .25)
		local xrange = xmax - xmin
		local xstep = 10^floor(log(xrange, 10) - .5)
		gl.glBegin(gl.GL_LINES)
		for x=floor(xmin/xstep)*xstep,ceil(xmax/xstep)*xstep,xstep do
			gl.glVertex2f(x,ymin)
			gl.glVertex2f(x,ymax)
		end
		gl.glEnd()
		local yrange = ymax - ymin
		local ystep = 10^floor(log(yrange, 10) - .5)
		gl.glBegin(gl.GL_LINES)
		for y=floor(ymin/ystep)*ystep,ceil(ymax/ystep)*ystep,ystep do
			gl.glVertex2f(xmin,y)
			gl.glVertex2f(xmax,y)
		end
		gl.glEnd()
		gl.glColor3f(.5, .5, .5)
		gl.glBegin(gl.GL_LINES)
		gl.glVertex2f(xmin, 0)
		gl.glVertex2f(xmax, 0)
		gl.glVertex2f(0, ymin)
		gl.glVertex2f(0, ymax)
		gl.glEnd()
--]]

		gl.glColor3f(unpack(info.color))
		gl.glPointSize(2)
		for _,mode in ipairs{gl.GL_POINTS, gl.GL_LINE_STRIP} do
			gl.glBegin(mode)
			for i=1,sim.gridsize do
				gl.glVertex2f(xs[i], ys[i])
			end
			gl.glEnd()
		end
		gl.glPointSize(1)
	end
	if self.reportError then
		self.reportError = false
	end
	gl.glViewport(0,0,w,h)
end
TestApp():run()
--]=]

