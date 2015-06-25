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
	local function is(x) return type(x) == 'table' and x.isa and x:isa(State) end
	local src = is(a) and a or b
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

return State

