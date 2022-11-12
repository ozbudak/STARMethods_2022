"""
Shared functionality for Python scripts
Copyright (C) 2017 Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi, Dong Mai, Soo Bin Kwon, Ha Vu

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
import os

def ensureDir(directory): # ensure the given directory exists
	if directory[-1] == '/':
		directory = directory[:-1]
	if not os.path.exists(directory):
		os.makedirs(directory)
	return directory

def isInt(s): # check if a given string is an integer
    try: 
        int(s)
        return True
    except ValueError:
        return False
        
def isFloat(s): # check if a given string is a float
    try: 
       	float(s)
        return True
    except ValueError:
        return False
        
# try to open the file specified by the given filename
def openFile(filename, mode):
	try:
		return open(filename, mode)
	except IOError:
		print("Couldn't open '" + filename + "'!")
		exit(1)

# try to cast the given value to an integer
def toInt(val):
	try:
		val = int(val)
	except ValueError:
		print("'" + val + "' is not an integer!")
		exit(2)
	return val

# try to cast the given value to a float
def toFlo(val):
	try:
		val = float(val)
	except ValueError:
		print("'" + val + "' is not a number!")
		exit(2)
	return val

