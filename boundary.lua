return {
	mirror = function(self)
		for i=1,self.numStates do
			self.qs[1][i] = self.qs[4][i]
			self.qs[2][i] = self.qs[3][i]
			self.qs[self.gridsize-1][i] = self.qs[self.gridsize-2][i]
			self.qs[self.gridsize][i] = self.qs[self.gridsize-3][i]
		end
		-- ... and negative the wavespeed variable ... which varies per-equation 
		self.qs[1][2] = -self.qs[4][2]
		self.qs[2][2] = -self.qs[3][2]
		self.qs[self.gridsize-1][2] = -self.qs[self.gridsize-2][2]
		self.qs[self.gridsize][2] = -self.qs[self.gridsize-3][2]
		-- special case for MHD ...
		if #self.qs[1] == 8 then
			self.qs[1][5] = -self.qs[1][5]
			self.qs[2][5] = -self.qs[2][5]
			self.qs[self.gridsize-1][5] = -self.qs[self.gridsize-2][5]
			self.qs[self.gridsize][5] = -self.qs[self.gridsize-3][5]
		end
	end,

	freeFlow = function(self)
		for i=1,self.numStates do
			self.qs[1][i] = self.qs[3][i]
			self.qs[2][i] = self.qs[3][i]
			self.qs[self.gridsize][i] = self.qs[self.gridsize-2][i]
			self.qs[self.gridsize-1][i] = self.qs[self.gridsize-2][i]
		end
	end,
}

