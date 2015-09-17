return {
	mirror = function(qs)
		local gridsize = #qs
		local numStates = #qs[1]
		for i=1,numStates do
			qs[1][i] = qs[4][i]
			qs[2][i] = qs[3][i]
			qs[gridsize-1][i] = qs[gridsize-2][i]
			qs[gridsize][i] = qs[gridsize-3][i]
		end
		-- ... and negative the wavespeed variable ... which varies per-equation 
		qs[1][2] = -qs[4][2]
		qs[2][2] = -qs[3][2]
		qs[gridsize-1][2] = -qs[gridsize-2][2]
		qs[gridsize][2] = -qs[gridsize-3][2]
		-- special case for MHD ...
		if #qs[1] == 8 then
			qs[1][5] = -qs[1][5]
			qs[2][5] = -qs[2][5]
			qs[gridsize-1][5] = -qs[gridsize-2][5]
			qs[gridsize][5] = -qs[gridsize-3][5]
		end
	end,

	freeFlow = function(qs)
		local gridsize = #qs
		local numStates = #qs[1]
		for i=1,numStates do
			qs[1][i] = qs[3][i]
			qs[2][i] = qs[3][i]
			qs[gridsize][i] = qs[gridsize-2][i]
			qs[gridsize-1][i] = qs[gridsize-2][i]
		end
	end,
}

