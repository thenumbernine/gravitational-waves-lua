#!/usr/bin/env luajit

require 'ext'
local slopeLimiters = require 'limiter' 
local boundaryMethods = require 'boundary'

-- here's some globals I have to get rid of

-- meta __index operation -- equivalent of operator[]() in C++
function index(k,v) return k[v] end

-- makes math easier!
for k,v in pairs(math) do _G[k] = v end


local symmath = require 'symmath'

-- setup

-- [[
local sim
do
	local xc = 150
	local x = symmath.var'x'
	local alpha = symmath.var'alpha'
	local h = 5 * symmath.exp(-((x - xc) / 10)^2)
	--sim = require'adm1d3var'{
	sim = require'adm1d5var'{
		gridsize = 2000,
		domain = {xmin=0, xmax=300},
		boundaryMethod = boundaryMethods.freeFlow,
		slopeLimiter = slopeLimiters.donorCell,
		-- the symbolic math driving it:
		x = x,
		h = h,
		g = 1 - h:diff(x)^2,
		alpha = 1,
		-- Bona-Masso slicing conditions:
		f_param = alpha,	
		--f = 1,
		--f = 1.69,
		--f = .49,
		f = 1 + 1/alpha^2,
	}
end
--]]

--[[
local sim
do
	-- r - eta(rs) = M ln(((rs + eta(r)) / (rs - eta(rs)))
	-- eta(rs) = sqrt(rs^2 - 2 M rs)
	local rc = 300
	local r = symmath.var'r'
	local alpha = symmath.var'alpha'
	local h = 5 * symmath.exp(-((r - rc) / 10)^2)
	sim = require'adm2dspherical'{
		gridsize = 200,
		domain = {xmin=100, xmax=500},
		boundaryMethod = boundaryMethods.freeFlow,
		slopeLimiter = slopeLimiters.donorCell,
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
	}
end
--]]

--[[
local sim = require'euler1d'{
	gridsize = 200,
	domain = {xmin=-1, xmax=1},
	boundaryMethod = boundaryMethods.mirror,
	slopeLimiter = slopeLimiters.superbee,
}
--]]

--[[
local sim = require'maxwell'{
	gridsize = 200,
	domain = {xmin=-1, xmax=1},
	boundaryMethod = boundaryMethods.freeFlow,
	slopeLimiter = slopeLimiters.donorCell,
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

		if not ymin or not ymax or ymin ~= ymin or ymax ~= ymax then
			ymin = -1
			ymax = 1
		--elseif abs(ymin) == huge or abs(ymax) == huge then
		else
			local base = 10	-- round to nearest base-10
			local scale = 10 -- ...with increments of 10
			ymin, ymax = 1.1 * ymin - .1 * ymax, 1.1 * ymax - .1 * ymin	
			local newymin = (ymin<0 and -1 or 1)*(abs(ymin)==huge and 1e+100 or base^log(abs(ymin),base))
			local newymax = (ymax<0 and -1 or 1)*(abs(ymax)==huge and 1e+100 or base^log(abs(ymax),base))
			ymin, ymax = newymin, newymax
		end
		do
			local minDeltaY = 1e-5
			local deltaY = ymax - ymin
			if deltaY < minDeltaY then
				ymax = ymax + .5 * minDeltaY
				ymin = ymin - .5 * minDeltaY
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

		gl.glColor3f(unpack(info.color))
		gl.glPointSize(2)
		for _,mode in ipairs{
			--gl.GL_POINTS,
			gl.GL_LINE_STRIP
		} do
			gl.glBegin(mode)
			for i=1,sim.gridsize do
				gl.glVertex2f(xs[i], ys[i])
			end
			gl.glEnd()
		end
		gl.glPointSize(1)
		
		if self.font then
			local fontSizeX = (xmax - xmin) * .05
			local fontSizeY = (ymax - ymin) * .1
			local ystep = ystep * 10
			for y=floor(ymin/ystep)*ystep,ceil(ymax/ystep)*ystep,ystep do
				self.font:draw{
					pos={xmin * .9 + xmax * .1, y + fontSizeY * .5},
					text=tostring(y),
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

