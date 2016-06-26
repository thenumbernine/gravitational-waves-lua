local limiters = {
	{name='donorCell', func=function(r) return 0 end},
	{name='LaxWendroff', func=function(r) return 1 end},
	{name='MC', func=function(r) return math.max(0, math.min(2, .5 * (1 + r), 2 * r)) end},
	{name='superbee', func=function(r) return math.max(0, math.min(1, 2*r), math.min(2,r)) end},
	{name='BeamWarming', func=function(r) return r end},
	{name='Fromm', func=function(r) return .5 * (1 + r) end},
	{name='vanLeer', func=function(r) return (r + math.abs(r)) / (1 + math.abs(r)) end},
}

for _,limiter in ipairs(limiters) do
	limiters[limiter.name] = limiter
end

return limiters
