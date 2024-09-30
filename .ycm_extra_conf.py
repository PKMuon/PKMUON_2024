'''
.ycm_extra_conf.py - YouCompleteMe customizing options.

Copyright (C) 2024  @lyazj <lyazj@github.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
import os

# Split path of current file.
filename = os.path.abspath(__file__)
basename = os.path.basename(filename)
dirname = os.path.dirname(filename)
assert dirname != '/'

# Find next configuration script in higher level directories.
dirname_org = dirname
dirname = os.path.dirname(dirname)
while dirname != '/':
    if os.path.exists(os.path.join(dirname, basename)): break
    dirname = os.path.dirname(dirname)
filename = os.path.join(dirname, basename)

# Execute parent configuration script.
exec(open(filename).read())

# Add include paths for this project.
flags.append('-I' + os.path.join(dirname_org, 'include'))
flags.append('-I' + os.path.join(dirname_org, 'cry', 'src'))
