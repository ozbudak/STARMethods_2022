"""
Define a slice
Copyright (C) 2016 Ahmet Ay, Dong Mai, Soo Bin Kwon

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy

class Slice:
	def __init__(self, top, bottom, top_left_xpos, top_right_xpos, bottom_left_xpos, bottom_right_xpos, last_slice):
		self.top = top
		self.bottom = bottom
		self.top_left_xpos = top_left_xpos
		self.top_right_xpos = top_right_xpos
		self.bottom_left_xpos = bottom_left_xpos
		self.bottom_right_xpos = bottom_right_xpos
		self.slopeL = float("inf") if top_left_xpos == bottom_left_xpos else (top - bottom)/(top_left_xpos - bottom_left_xpos)
		self.slopeR = float("inf") if top_right_xpos == bottom_right_xpos else (top - bottom)/(top_right_xpos - bottom_right_xpos)
		self.last_slice = last_slice # whether this slice is the last (rightmost) slice in the region
		self.cells = []
		self.her1_levels = []
		self.her7_levels = []
		self.keep_cells = []
		self.her1_bgNlevels = [] 	# background normalized her1 level
		self.her7_bgNlevels = []	# background normalized her7 level
		self.num_cells = 0
		self.num_keep_cells = 0
		self.valid = True # whether this slice is usable (i.e. has 3 or more cells)

	def identify_cells(self, cell_list):
		for cell in cell_list:
			left_boundary = self.bottom_left_xpos + (cell.ypos - self.bottom)/self.slopeL
			right_boundary = self.bottom_right_xpos + (cell.ypos - self.bottom)/self.slopeR
			if (left_boundary <= cell.xpos if self.last_slice else left_boundary < cell.xpos) and cell.xpos <= right_boundary and self.bottom <= cell.ypos and cell.ypos <= self.top:
				# Remain information from all cell for heatmap plotting visualization purpose
				self.cells.append(cell)
				self.her1_levels.append(cell.her1)
				self.her7_levels.append(cell.her7)
				
				# Write background subtracted her value to slices.xls
				self.keep_cells.append(cell)
				self.her1_bgNlevels.append(cell.her1_bgN)
				self.her7_bgNlevels.append(cell.her7_bgN)
		self.num_cells = len(self.cells)
		self.num_keep_cells = len(self.keep_cells)
		
		# Ignore slices with fewer than 3 cells
		if self.num_cells <= 2:
			self.valid = False
			return False

		return True
		
	def slice_variance_her1(self):
		return numpy.var(self.her1_levels) 

	def slice_variance_her7(self):
		return numpy.var(self.her7_levels) 

	def slice_mean_her1(self):
		return numpy.mean(self.her1_levels)

	def slice_mean_her7(self):
		return numpy.mean(self.her7_levels)

	def slice_variance_her1_bgN(self):
		return numpy.var(self.her1_bgNlevels) 

	def slice_variance_her7_bgN(self):
		return numpy.var(self.her7_bgNlevels) 

	def slice_mean_her1_bgN(self):
		return numpy.mean(self.her1_bgNlevels)

	def slice_mean_her7_bgN(self):
		return numpy.mean(self.her7_bgNlevels)
