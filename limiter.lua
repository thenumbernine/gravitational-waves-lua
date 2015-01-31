return {
	donorCell = function(r) return 0 end,
	LaxWendroff = function(r) return 1 end,
	MC = function(r) return math.max(0, math.min(2, .5 * (1 + r), 2 * r)) end,
	superbee = function(r) return math.max(0, math.min(1, 2*r), math.min(2,r)) end,
	BeamWarming = function(r) return r end,
	Fromm = function(r) return .5 * (1 + r) end,
	vanLeer = function(r) return (r + math.abs(r)) / (1 + math.abs(r)) end,
}
