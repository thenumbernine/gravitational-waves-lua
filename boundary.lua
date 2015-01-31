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

