#!/usr/bin/env python
from __future__ import print_function
from sys import stdout, argv
import re
import bisect
from math import radians, cos, sqrt
#from random import random #!! (for random colors in debug)


# Path simplification thresholds:
eps = 0.02 # deviation (original distance units)
ang = cos(radians(180 - 45)) # kink angle (degrees)

# Color of electrodes:
meshcolor = '0.549 0.431 0.314'
# Color translations:
recolor = {
  meshcolor          : '0 0 0', # electrodes -> black
  '0.737 0.000 0.000': '0.8 0 0', # red
  '0.000 0.737 0.000': '0 0.6 0', # green
  '0.000 0.000 0.737': '0 0 1',   # blue
  '0.988 0.784 0.784': '1 0.5 0.5', # light red
  '0.627 0.988 0.627': '0.7 1 0.7', # light green
  '0.784 0.784 0.992': '0.6 0.6 1', # light blue
}

# Output line width [pt]:
linewidth = 0.5
# Time marker radius [pt]:
markersize = 1


# Regular expressions for color switch, line segment and other:
recol = re.compile(r'([\.\d]+\s+[\.\d]+\s+[\.\d]+)\s+col')
relin = re.compile(r'\s*([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)\s+drw1')
rebox = re.compile(r'\s*([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\.\d]+)\s+boxfill')
reend = re.compile(r'%(-|=|%Trailer)')

# (for regexps in if ... elif ...)
class Var(object):
  '''
  Helper class for assignments within expressions.
  var = Var() creates a storage variable.
  var(expr) returns the expr value and saves it for future use.
            (Note: if expr == None, the value is not saved (see below).
            If needed, use var.set(expr) instead.)
  var() returns the saved value.
  Any attribute attr can be accessed directly as var.attr.
  '''
  def __init__(self, val = None):
    self.__val = val
  def __getattr__(self, attr):
    return getattr(self.__val, attr)
  def __call__(self, arg = None):
    if arg is not None:
      self.__val = arg
    return self.__val
  def set(self, arg):
    self.__val = arg
    return self.__val

# Global array for line segments read from input (or constructed):
lines = []

# Convert electrodes mesh to boundaries and "ideal grids".
# (lines -> lines, linesg)
def process_mesh():
  global lines
  global linesg # SIMION "grids" (non-solid electrodes)

  print('  Constructing mesh...')
  # Recover the underlying mesh from lines:
  # make sorted lists of all x and y coordinates:
  xs, ys = [], []
  def add_coord(arr, c):
    i = bisect.bisect_left(arr, c)
    if i == len(arr) or arr[i] != c:
      arr.insert(i, c)        
  for line in lines:
    add_coord(xs, line[0])
    add_coord(ys, line[1])
    add_coord(xs, line[2])
    add_coord(ys, line[3])
  # add guarding edges:
  add_coord(xs, xs[ 0] - 1)
  add_coord(xs, xs[-1] + 1)
  add_coord(ys, ys[ 0] - 1)
  add_coord(ys, ys[-1] + 1)
  # final dimensions:
  Nx, Ny = len(xs), len(ys)
  print('   Mesh size: {} X {}.'.format(Nx, Ny))
  # nodes have bit flags according to present segments (0x0 = no segments):
  nodes = [[0x0 for _ in range(Ny)] for _ in range(Nx)]
      
  # !! Mesh printing (for debug).
  def print_mesh():
    for iy in reversed(range(Ny)):
      for ix in range(Nx):
        stdout.write('{:X}'.format(nodes[ix][iy]))
      stdout.write('\n')

  # Fill the mesh (nodes) with line segments:
  print('   Filling mesh...')
  for line in lines:
    x0, y0, x1, y1 = line
    # horizontal:
    if y0 == y1:
      iy = bisect.bisect_left(ys, y0)
      ix0, ix1 = bisect.bisect_left(xs, x0), bisect.bisect_left(xs, x1)
      if ix0 > ix1:
        ix0, ix1 = ix1, ix0
      for ix in range(ix0, ix1): # all segments
        nodes[ix][iy] |= 0x1 # to right
    # vertical:
    else:
      ix = bisect.bisect_left(xs, x0)
      iy0, iy1 = bisect.bisect_left(ys, y0), bisect.bisect_left(ys, y1)
      if iy0 > iy1:
        iy0, iy1 = iy1, iy0
      for iy in range(iy0, iy1): # all segments
        nodes[ix][iy] |= 0x2 # to top
#  print_mesh()

  # Close open cells (truncated lines) in mesh data:
  print('   Closing open cells...')
  # Sides:
  closed_l = ''
  closed_r = ''
  closed_b = ''
  closed_t = ''
  for ix in range(1, Nx - 1):
    # bottom:
    iy = 1
    if (nodes[ix][iy] == 0x2 and
        nodes[ix + 1][iy] == 0x2 and
        nodes[ix][iy + 1] & 0x1):
      nodes[ix][iy] |= 0x1
      closed_b = 'B'
    # top:
    iy = Ny - 2
    if (nodes[ix][iy] == 0x0 and
        nodes[ix][iy - 1] == 0x3 and
        nodes[ix + 1][iy - 1] & 0x2):
      nodes[ix][iy] = 0x1
      closed_t = 'T'
  for iy in range(1, Ny - 1):
    # left:
    ix = 1
    if (nodes[ix][iy] == 0x1 and
        nodes[ix][iy + 1] == 0x1 and
        nodes[ix + 1][iy] & 0x2):
      nodes[ix][iy] |= 0x2
      closed_l = 'L'
    # right:
    ix = Nx - 2
    if (nodes[ix][iy] == 0x0 and
        nodes[ix - 1][iy] == 0x3 and
        nodes[ix - 1][iy + 1] & 0x1):
      nodes[ix][iy] = 0x2
      closed_r = 'R'
  # Corners:
  # left, bottom:
  if (closed_l and closed_b and
      nodes[1][2] & 0x1 and
      nodes[2][1] & 0x2):
    nodes[1][1] = 0x3
  # right, bottom:
  if (closed_r and closed_b and
      nodes[Nx - 3][1] & 0x2 and
      nodes[Nx - 3][2] & 0x1):
    nodes[Nx - 3][1] |= 0x1
    nodes[Nx - 2][1] = 0x2
  # left, top:
  if (closed_l and closed_t and
      nodes[1][Ny - 3] & 0x1 and
      nodes[2][Ny - 3] & 0x2):
    nodes[1][Ny - 3] |= 0x2
    nodes[1][Ny - 2] = 0x1
  # right, top:
  if (closed_r and closed_t and
      nodes[Nx - 3][Ny - 3] == 0x3):
    nodes[Nx - 3][Ny - 2] = 0x1
    nodes[Nx - 2][Ny - 3] = 0x2
  print('    closed: {}.'.format(closed_l + closed_r + closed_b + closed_t or 'none'))
#  print_mesh()

  # Collect boundary lines for solid and separately non-solid:
  print('  Searching boundaries...')
  lines = []  # solid
  linesg = [] # non-solid (SIMION "grids")
  for ix in range(1, Nx - 1):
    for iy in range(1, Ny - 1):
      f  = nodes[ix][iy]
      # horizontal (to right):
      if f & 0x1:
        b = nodes[ix][iy - 1] & nodes[ix + 1][iy - 1] & 0x2 # cell below
        a = f & nodes[ix + 1][iy] & 0x2 # cell above
        if b != a: # boundary
          lines.append((xs[ix], ys[iy], xs[ix + 1], ys[iy]))
        elif not b: # (= not a) not inner
          linesg.append((xs[ix], ys[iy], xs[ix + 1], ys[iy]))
      # vertical (to top):
      if f & 0x2:
        l = nodes[ix - 1][iy] & nodes[ix - 1][iy + 1] & 0x1 # cell left
        r = f & nodes[ix][iy + 1] & 0x1 # cell right
        if l != r: # boundary
          lines.append((xs[ix], ys[iy], xs[ix], ys[iy + 1]))
        elif not l: # (= not r) not inner
          linesg.append((xs[ix], ys[iy], xs[ix], ys[iy + 1]))
  print('   {} boundary segments.'.format(len(lines)))
  print('   {} "ideal grid" segments.'.format(len(linesg)))


# Collect separate line segments into continuous paths and output them.
def process_paths(do_sort = True, fill_closed = False):
  global lines

  print('  Collecting paths...')
  # (preliminary construction of paths from segments works much faster than
  # joining huge ammount of short paths)
  if do_sort:
    print('   (sort)')
    lines.sort() # speeds up and improves the processing.
                 # (?? a throughout sorting and binary search should help
                 #  !! but can create "wrong" loops for intersecting paths
                 #  ?? find the "most straight" continuation?)
  paths = []
  for l in lines:
    p0, p1 = l[0:2], l[2:4]
    joined = False
    # try to find existing path (exactly one) to add to:
    for path in paths:
      if   p0 == path[0]:  # beginning at path beginning
        path.insert(0, p1) # insert end there
        joined = True; break
      elif p0 == path[-1]: # beginning at path end
        path.append(p1)    # insert end there
        joined = True; break
      elif p1 == path[0]:  # end at path beginning
        path.insert(0, p0) # insert beginning there
        joined = True; break
      elif p1 == path[-1]: # end at path end
        path.append(p0)    # insert beginning there
        joined = True; break
    # otherwise -- create new path:
    if not joined:
      paths.append([p0, p1])
  print('   {} paths.'.format(len(paths)))
  lines = [] # not needed any more; must be emptied for the next group

  # Join all paths that can be joined:
  # Try to find two paths and join them.
  def joinpaths():
    N = len(paths)
    for i in range(N):
      for j in range(i + 1, N):
        # end of i == start of j:
        if paths[i][-1] == paths[j][0]:
          paths[i] += paths[j][1:]
          del paths[j]
          #print('   Joined {} and {}.'.format(i, j))
        # start of i == end of j:
        elif paths[i][0] == paths[j][-1]:
          paths[j] += paths[i][1:]
          del paths[i]
          #print('   Joined {} and {}.'.format(j, i))
        # start of i == start of j:
        elif paths[i][0] == paths[j][0]:
          paths[j].reverse()
          paths[j] += paths[i][1:]
          del paths[i]
          #print('   Joined reversed {} and {}.'.format(j, i))
        # end of i == end of j:
        elif paths[i][-1] == paths[j][-1]:
          paths[j].reverse()
          paths[i] += paths[j][1:]
          del paths[j]
          #print('   Joined and {} reversed {}.'.format(i, j))
        else: # not joined -- try next
          continue
        return True # joined -- return
    return False # nothing found for joining
  # try to fing all possibilities:
  stdout.write('  Joining paths') # start progress indicator
  while joinpaths():
    stdout.write('.') # progress indicator
    stdout.flush()
  stdout.write('\n') # end progress indicator
  print('   {} paths.'.format(len(paths)))

  # Remove intermediate points with small deviations from the line:
  stdout.write('  Simplifying paths') # (start progress indicator)
  for path in paths:
    stdout.write('.') # (progress indicator)
    stdout.flush()
    # Check maximal deviation (False if distance > eps or sharp turn).
    def checkpath(b, e):
      pb, pe = path[b], path[e]
      dX, dY = pe[0] - pb[0], pe[1] - pb[1]
      dR = sqrt(dX**2 + dY**2)
      for i in range(b + 1, e):
        p = path[i]
        dx, dy = p[0] - pb[0], p[1] - pb[1]
        # Check that angle at p is not less obtuse than threshold.
        def kink():
          a, b = path[i - 1], path[i + 1]
          ax, ay = a[0] - p[0], a[1] - p[1]
          bx, by = b[0] - p[0], b[1] - p[1]
          return (ax * bx + ay * by) > \
                 ang * sqrt((ax * ax + ay * ay) * (bx * bx + by * by))
        if abs(dx * dY - dX * dy) > eps * dR or kink():
          return False
      return True
    # Remove longest chunks (greedy):
    # loop over possible beginnings:
    i = 0
    while i < len(path) - 2:
      # loop over possible ends:
      # (search with progressive steps)
      j = i + 2
      dj = 1
      l = len(path)
      while j < l:
        if checkpath(i, j):
          if j == l - 1: # fine upto the path end
            del path[i + 1:-1]
            i = l # set i > new len(path) to stop the beginnings loop
            break
          # increase step, advance possible end:
          dj *= 2
          if j + dj >= l: # too far
            dj = (l - 1) - j
          j += dj
        else: # large deviation
          if dj == 1: # the step was minimal ==> j - 1 was the last good point
            del path[i + 1:j - 1]
            i += 1 # set the next possible beginning to the found end
            break
          else: # step was large
            j = j - dj + 1 # return to the previous good point + 1
            dj = 1         # reset step size to minimum
      """ # direct scan (leads to total time ~ O(N^3)!)
      for j in range(i + 2, len(path)):
        if not checkpath(i, j):
          # previous end was fine ==> remove the chunk before it
          del path[i + 1:j - 1]
          i += 1 # set the next possible beginning to the "previous end"
          break # return to the beginnings loop
      else: # last segment was fine entirely ==> leave the last point only
        del path[i + 1:-1]
      #"""
  stdout.write('\n') # (end progress indicator)

  # Output of each path (in differential form):
  print('  Writting paths...')
  # If requested, fill closed paths (before drawing open ones):
  if fill_closed:
    found_closed = False
    i = 0
    while i < len(paths):
      path = paths[i]
      if path[0] == path[-1]: # closed
        found_closed = True
#        fout.write('{:g} {:g} {:g} col\n'.format(random(), random(), random()))
        p0 = path[0]
        fout.write('{:g} {:g} m\n'.format(*path[0]))
        for p in path[1:-1]: # without the last point (== first)
          fout.write('{:g} {:g} r\n'.format(p[0] - p0[0], p[1] - p0[1]))
          p0 = p
        # (path will be closed automatically)
        del paths[i] # delete it (only open should remain for drawing)
      else:
        i += 1
    # even-odd fill all closed paths together (for nested regions):
    if found_closed:
      fout.write('f\n')
  # Draw all open paths separately:
  for path in paths:
#    fout.write('{:g} {:g} {:g} col\n'.format(random(), random(), random()))
    p0 = path[0]
    fout.write('{:g} {:g} m\n'.format(*path[0]))
    for p in path[1:]:
      fout.write('{:g} {:g} r\n'.format(p[0] - p0[0], p[1] - p0[1]))
      p0 = p
    fout.write('s\n')

colors = [] # list of colors in current group
colorlines = [] # line segments in current group, by color
colorboxes = [] # boxes in current group, by color
color, colorn = '', -1 # current color and its number in colors[]

# Check command line:
if (len(argv) != 2):
  print('Clean-up SIMION output (EPS print-out).')
  print('Usage: cleanSI.py in.eps')
  print('  in.eps  EPS file created by SIMION')
  print('  (output is written to out.eps)')
  print('Look inside the program source for adjustable parameters.')
  exit()
# Open input and output:
try:
  fin = open(argv[1], 'rU')
except:
  print('Cannot open input file "{}"!'.format(argv[1]))
  exit()
try:
  fout = open('out.eps', 'w')
except:
  print('Cannot open output file "out.eps"!')
  exit()
# Processing:
mesh = False # mesh data mode (switched by color)
for textline in fin: # go over all input lines...
  m = Var() # (variable to store MatchObject)
  # line segments:
  if m(relin.match(textline)):
    line = [float(x) for x in m.groups()]
    # add only if not degenerate:
    if line[0:2] != line[2:4]:
      colorlines[colorn].append(line)
  # boxes (time markers):
  elif m(rebox.match(textline)):
    colorboxes[colorn].append([float(x) for x in m.groups()])
  # colors:
  elif m(recol.match(textline)):
    color = m.group(1)
    # select color number:
    try:
      colorn = colors.index(color)
    except: # not found ==> add
      colorn = len(colors)
      colors.append(color)
      colorlines.append([])
      colorboxes.append([])
  # end of line groups ("%%Trailer" or manually inserted additional
  # "%-" or "%=" (latter means do not reorder segments)):
  elif m(reend.match(textline)):
    stdout.write('Group ended with {}'.format(textline))
    do_sort = m.group(1) != '='
    # go over each found color:
    for colorn in range(len(colors)):
      # change color:
      color = colors[colorn]
      print(' Color "{}":'.format(color))
      if color in recolor:
        color = recolor[color]
        print(' Changed to "{}".'.format(color))
      fout.write('{} col\n'.format(color))
      # take its lines:
      lines = colorlines[colorn]
      print(' {} line segments.'.format(len(lines)))
      if len(lines):
        # process and output:
        # if electrodes color -- as mesh:
        if colors[colorn] == meshcolor:
          # convert solid to outlines and separate non-solid:
          process_mesh()
          process_paths(do_sort, True) # filling closed paths
          # process non-solid:
          lines = linesg
          process_paths(do_sort)
        # otherwise -- as regular curves:
        else:
          process_paths(do_sort)
      else: # no lines
        print('  Skipped.')
      # and boxes:
      boxes = colorboxes[colorn]
      if len(boxes):
        print(' {} boxes.'.format(len(boxes)))
        for b in boxes:
          # put marker at box center:
          xc = (b[0] + b[2]) / 2
          yc = (b[1] + b[3]) / 2
          fout.write('{} {} o\n'.format(xc, yc))
    # clear data for the next group:
    colors = []
    colorlines = []
    # output the original input line:
    fout.write(textline)
  # end of EPS prolog:
  elif textline == '%%EndProlog\n':
    # append new definitions:
    fout.write('/m {moveto} def\n/r {rlineto} def\n')
    fout.write('/s {stroke} def\n/f {eofill} def\n')
    fout.write('/o {{{} 0 360 arc fill}} def\n'.format(markersize))
    fout.write('/wid {pop} def\n') # disable SIMION's line width
    # output the "EndProlog" itself:
    fout.write(textline)
  # line stroke settings:
  elif textline == '1 setlinewidth 2 setlinecap 0 setgray\n':
    # replace with round ends and joints:
    fout.write('{} setlinewidth 1 setlinecap 1 setlinejoin 0 setgray\n'.
               format(linewidth))  
  # all other text lines:
  else:
    # just copy to output:
    fout.write(textline)
# Close input and output:
fin.close()
fout.close()
