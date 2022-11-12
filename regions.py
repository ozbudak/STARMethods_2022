"""
Modify by leeyy, on 10 Feb 2018
	- adding in new function create_slices_dynamics() to alter the slope gradually based on given function
	- run new function in __init__()
	
Define a region/half
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
import math
import itertools
from slices import Slice

class Region:
	def __init__(self, num_sec, cell_lists, name, angle, delta_angle):
		self.secs = []
		self.name = name
		self.num_sec = num_sec
		self.cell_list = []
		for i in range(num_sec):
			self.cell_list += cell_lists[i]		
		self.slice_width = 8
		self.radian = angle/180 * math.pi
		self.delta_radian = delta_angle/180 * math.pi
		self.slope = math.tan(self.radian)
		self.xs_ys_calc()
		self.single_cell_boundaries()
		self.left_corners_calc()
		self.right_corners_calc()
		#self.create_slices()	# Old function to create fix angle slices
		self.create_dynamic_slice()
		

	def xs_ys_calc(self):	# to read x and y position for every cell
		self.xs = []
		self.ys = []
		for cell in self.cell_list:
			self.xs.append(cell.xpos)
			self.ys.append(cell.ypos)

	def single_cell_boundaries(self):	# find bottom, top, leftmost and rightmost cell position 
		self.bottom = self.cell_list[self.ys.index(min(self.ys))] # bottom cell
		self.top = self.cell_list[self.ys.index(max(self.ys))] # top cell
		self.left = self.cell_list[self.xs.index(min(self.xs))] # leftmost cell
		self.right = self.cell_list[self.xs.index(max(self.xs))] # rightmost cell

	def left_corners_calc(self):	# Define left border with fix angle
		bottom_candidates = []
		top_candidates = []
		for cell in self.cell_list:
			bottom_candidates.append(cell.xpos + (self.bottom.ypos - cell.ypos)/self.slope)
			top_candidates.append(cell.xpos + (self.top.ypos - cell.ypos)/self.slope)
		self.bottom_left_xpos = min(bottom_candidates)
		self.top_left_xpos = min(top_candidates)

						
	def right_corners_calc(self):	# Calculate right border for angle change algorithm
		bottom_candidates = []
		top_candidates = []
		radian = self.radian
		delta_radian = self.delta_radian
		baseline_left_xpos = max(self.bottom_left_xpos, self.top_left_xpos)
		for cell in self.cell_list:
			approx_slope = math.tan(radian + (cell.xpos - baseline_left_xpos) * delta_radian)
			bottom_candidates.append(cell.xpos + (self.bottom.ypos - cell.ypos)/approx_slope)
			top_candidates.append(cell.xpos + (self.top.ypos - cell.ypos)/approx_slope)
		self.bottom_right_xpos = max(bottom_candidates)
		self.top_right_xpos = max(top_candidates)	
		
	def create_slices(self):		# Create slices with fix angle - not in used for angle change version
		self.slices = []
		width = self.bottom_right_xpos - self.bottom_left_xpos			
		self.num_slices = int(width/self.slice_width) if width%self.slice_width == 0 else (int(width/self.slice_width) + 1)
		# first slice
		bottom_left_xpos = self.bottom_left_xpos
		bottom_right_xpos = bottom_left_xpos + self.slice_width	
		top_left_xpos = self.top_left_xpos
		top_right_xpos = top_left_xpos + self.slice_width
		for j in range(self.num_slices):
			if j == (self.num_slices-1): # last slice
				self.slices.append(Slice(self.top.ypos, self.bottom.ypos, top_left_xpos, self.top_right_xpos, bottom_left_xpos, self.bottom_right_xpos, True))			
			else: # not the last slice
				self.slices.append(Slice(self.top.ypos, self.bottom.ypos, top_left_xpos, top_right_xpos, bottom_left_xpos, bottom_right_xpos, False))
			bottom_left_xpos += self.slice_width
			bottom_right_xpos += self.slice_width
			top_left_xpos += self.slice_width
			top_right_xpos += self.slice_width
	
	def create_dynamic_slice(self):	# Create slice with increasing/decreasing angle based on delta_angle		
		self.slices = []
		self.height = self.top.ypos - self.bottom.ypos
		delta_radian = self.delta_radian * self.slice_width		# angle change for every slice 
		cur_radian = self.radian
		
		if self.name == 'L':				
			width = self.bottom_right_xpos - self.bottom_left_xpos	# For left PSM, use bottom x-width to define border		
			self.num_slices = int(width/self.slice_width) if width%self.slice_width == 0 else (int(width/self.slice_width) + 1)

			# Define first slice in left PSM					
			bottom_left_xpos = self.bottom_left_xpos
			bottom_right_xpos = bottom_left_xpos + self.slice_width	# taking constant step at bottom right
			top_left_xpos = self.top_left_xpos
			top_right_xpos = top_left_xpos + self.slice_width + abs(self.height/math.tan(cur_radian+delta_radian) - self.height/math.tan(cur_radian))	# modify each step based on angle change at top right
			
			for j in range(self.num_slices):	
				if j == (self.num_slices-1): # last slice
					self.slices.append(Slice(self.top.ypos, self.bottom.ypos, top_left_xpos, top_right_xpos, bottom_left_xpos, bottom_right_xpos, True))			
				else: # not the last slice
					self.slices.append(Slice(self.top.ypos, self.bottom.ypos, top_left_xpos, top_right_xpos, bottom_left_xpos, bottom_right_xpos, False))				

				cur_radian += delta_radian
				bottom_left_xpos += self.slice_width
				bottom_right_xpos += self.slice_width
				top_left_xpos = top_right_xpos	# replace next top left xpos with current top right xpos
				top_right_xpos += (self.slice_width + abs(self.height/math.tan(cur_radian+delta_radian) - self.height/math.tan(cur_radian)))
				
				
		elif self.name == 'R':								
			width = self.top_right_xpos - self.top_left_xpos	# For right PSM, use top x-width to define border
			self.num_slices = int(width/self.slice_width) if width%self.slice_width == 0 else (int(width/self.slice_width) + 1)			

			# Define first slice in right PSM
			top_left_xpos = self.top_left_xpos
			top_right_xpos = top_left_xpos + self.slice_width	# taking constant step at top right
			bottom_left_xpos = self.bottom_left_xpos
			bottom_right_xpos = bottom_left_xpos + self.slice_width	+ abs(self.height/math.tan(cur_radian+delta_radian) - self.height/math.tan(cur_radian))	# modify each step based on angle change at bottom right

			for j in range(self.num_slices):
				if j == (self.num_slices-1): # last slice
					self.slices.append(Slice(self.top.ypos, self.bottom.ypos, top_left_xpos, top_right_xpos, bottom_left_xpos, bottom_right_xpos, True))			
				else: # not the last slice
					self.slices.append(Slice(self.top.ypos, self.bottom.ypos, top_left_xpos, top_right_xpos, bottom_left_xpos, bottom_right_xpos, False))

				cur_radian += delta_radian 
				top_left_xpos += self.slice_width
				top_right_xpos += self.slice_width
				bottom_left_xpos = bottom_right_xpos	# replace next bottom left xpos with current bottom right xpos
				bottom_right_xpos += (self.slice_width + abs(self.height/math.tan(cur_radian+delta_radian) - self.height/math.tan(cur_radian)))
						