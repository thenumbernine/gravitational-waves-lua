local class = require 'ext.class'

local State = class()

function State:init(h, w)
	for i=1,h do
		self[i] = {}
		for j=1,w do
			self[i][j] = 0
		end
	end
end

function State.__add(a,b)
	local c = State(#a, #a[1])
	for i=1,#a do
		for j=1,#a[1] do
			c[i][j] = a[i][j] + b[i][j]
		end
	end
	return c
end

function State.__mul(a,b)
	local src = State.is(a) and a or b
	local c = State(#src, #src[1])
	for i=1,#src do
		for j=1,#src[1] do
			local aij = type(a) == 'number' and a or a[i][j]
			local bij = type(b) == 'number' and b or b[i][j]
			c[i][j] = aij * bij
		end
	end
	return c
end

function State.__div(a,b)
	assert(type(b) == 'number')
	return a * (1/b)
end

function State.__unm(a)
	return -1 * a
end

function State.__sub(a,b)
	return a + -b
end

function State:clone()
	local h = #self
	local w = #self[1]
	local new = State(h, w)
	for i=1,h do
		for j=1,w do
			new[i][j] = self[i][j]
		end
	end
	return new
end

function State.dot(a,b)
	local s = 0
	for i=1,#a do
		for j=1,#a[1] do
			s = s + a[i][j] * b[i][j]
		end
	end
	return s
end

function State.perElementMultiply(a,b)
	local c = State(#a, #a[1])
	for i=1,#a do
		for j=1,#a[i] do
			c[i][j] = a[i][j] * b[i][j]
		end
	end
	return c
end

function State.perElementDivide(a,b)
	local c = State(#a, #a[1])
	for i=1,#a do
		for j=1,#a[i] do
			c[i][j] = a[i][j] / b[i][j]
		end
	end
	return c
end

return State
