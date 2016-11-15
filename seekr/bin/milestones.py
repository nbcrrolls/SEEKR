#!/usr/bin/python

'''
milestones.py

Contains functions to be called by seekr.py that create anchor grids and their associated milestones.


'''

import numpy as np
from numpy import arange
from math import sqrt, sin, cos, ceil, log10, atan2, pi
import cmath
import unittest

import xml.etree.cElementTree as ET
from xml.dom.minidom import Document
from xml.dom import minidom

from bd import dict2xml # here is the problem I think

import positions_orient
from copy import deepcopy



verbose = True

'''class Anchor():
  'represents a milestoning anchor'
  def __init__(self,coord):
    self.coord = coord'''


'''
quaternions_tesseract = [[0.5, 0.5, 0.5, 0.5],
                        [-0.5, 0.5, 0.5, 0.5],
                        [ 0.5,-0.5, 0.5, 0.5],
                        [ 0.5, 0.5,-0.5, 0.5],
                        [ 0.5, 0.5, 0.5,-0.5],
                        [-0.5,-0.5, 0.5, 0.5],
                        [ 0.5,-0.5,-0.5, 0.5],
                        [-0.5, 0.5,-0.5, 0.5]]  # quaternions for the cells of a tesseract; useful for uniform SO(3) rotation
#tesseract_neighbors = {0:(1,2,3,4), 1:(0, 5, 6, 7), 2:(0, 5, 6, 7)}

#quaternions_hyperoctahedron =

quaternions_simplex = normalize([[1/sqrt(10.0), 1/sqrt(6.0), 1/sqrt(3.0), 1.0],
                       [1/sqrt(10.0), 1/sqrt(6.0), 1/sqrt(3.0),-1.0],
                       [1/sqrt(10.0), 1/sqrt(6.0),-2/sqrt(3.0), 0.0],
                       [1/sqrt(10.0), -sqrt(3.0/2.0), 0.0, 0.0],
                       [-2*sqrt(2.0/5.0), 0.0, 0.0, 0.0]]) # quaternions for a pentachoron: a 4-simplex (tetrahedron in 4D)

simplex_neighbors = {0:[1,2,3,4], 1:[0,2,3,4], 2:[0,1,3,4], 3:[0,1,2,4], 4:[0,1,2,3]}
'''

class Milestone():
  '''represents a surface in phase space that is monitored for crossings'''
  def __init__(self, anchor, dimensions, index, siteid, center_type, absolute='False', shape='sphere', md=True, bd=False):
    self.shape = shape # can be 'sphere', 'plane', 'ellipse'?
    self.fullname = ""
    self.directory = ""
    self.anchor = anchor # the location of this milestone's anchor
    self.rotation = [1,0,0,0]
    self.neighbors = []
    self.index = index
    self.siteid = siteid
    self.center_type = center_type
    self.dimensions = dimensions
    self.absolute = absolute
    self.md = md # whether md simulations are run fro this milestone
    self.bd = bd # whether bd simulations are run from this milestone
    self.bd_adjacent = None
    self.end = False

def in_ellipsoid(a, b, c, x, y, z, i, j, k):
  '''determine whether point i,j,k are inside the ellipse described by a,b,c and x,y,z
'''
  for argument in ['a','b','c','x','y','z','i','j','k']: # convert all variables to floats
    command = "%s = float(%s)" % (argument, argument)
    exec command
  #print sqrt( ((i-x)**2 / a**2) + ((j-y)**2 / b**2) + ((k-z)**2 / c**2) )
  if sqrt( ((i-x)**2 / a**2) + ((j-y)**2 / b**2) + ((k-z)**2 / c**2) ) <= 1.0:
    return True # then the point i, j, k falls within the ellipsoid
  else:
    return False

def quadratic_formula(a,b,c):
  '''calculates the roots of the equation:
       a*x**2 + b*x + c = 0
     using the quadratic formula
  '''
  results = []
  square_root = b*b - 4.0 * a * c
  for i in [-1.0, 1.0]:
    if square_root >= 0.0: # then the roots will be real
      root = (-b + (i * sqrt(square_root))) / (2.0 * a)
    else: # then the roots will be complex
      root = (-b + (i * cmath.sqrt(square_root))) / (2.0 * a)
    results.append(root)
  return results


def find_rad_of_skewed_vectors(A, B, r):
  '''Finds x such that norm(A + xB) = r
  input:
    A: the vector pointing to the first radial milestone
    B: the vector pointing from A to the milestone of interest
    r: the radius we are trying to obtain
  output:
    x: the length to scale a normalized B in order to make norm(A + B) == r
  '''
  assert np.linalg.norm(A) <= r, "The milestone radius must be larger than the minimum milestone radius"
  if np.linalg.norm(B) > 0.0: B_prime = B / np.linalg.norm(B)
  a = np.dot(B_prime,B_prime) # 'a' for the quadratic formula
  b = 2.0 * np.dot(A,B_prime) # 'b' for the quadratic formula
  c = np.dot(A,A) - r**2 # 'c' for the quadratic formula
  roots = quadratic_formula(a,b,c)
  x = min(map(abs,roots)) # get the root with the smallest absolute value from the quadratic formula
  return x

def generate_cubic_positions(a, b, c, x, y, z, increment, ellipsoid=False, neighbors_in_ellipsoid_only=False):
  '''given a box with dimensions a,b,c centered at x,y,z, will generate a list of all positions within an
ellipse with the same dimensions. will also return a list of neighbors of every corresponding grid point in the ellipse'''
  neighbor_dict = {} # a dictionary with an index for each grid point, keeping track of all neighbors
  num_pos = 0 # the total number of positions within the ellipse
  milestones = [] # a list of of anchors within the ellipse
  # first loop through every point of our cube
  for i in arange(x-a, x+a+increment, increment):
    for j in arange(y-b, y+b+increment, increment):
      for k in arange(z-c, z+c+increment, increment):
        if ellipsoid: # then we have to make sure its within an ellipsoid
          if not in_ellipsoid(a, b, c, x, y, z, i, j, k): continue # skip if not in ellipsoid
        #anchorlist.append((i,j,k))
        #neighbors[(i,j,k)] = [] # create a space for the neighbors of this point
        #num_pos += 1
        anchor = None #Anchor((i,j,k)) # NOTE: this entire function needs to be fixed or deleted
        milestone = Milestone(shape='plane',anchor=anchor, index=num_pos)
        milestones.append(milestone)
        neighbor_dict[(i,j,k)] = num_pos # create a space for the neighbors of this point
        num_pos += 1

  # now run through and find adjacent neighbors
  for counter in range(num_pos):
    for i1 in [-increment,0,increment]:
      for j1 in [-increment,0,increment]:
        for k1 in [-increment,0,increment]:
          if i1==j1==k1==0: continue
          #[i,j,k] = anchorlist[counter]
          if neighbors_in_ellipsoid_only and not in_ellipsoid(a, b, c, x, y, z, i+i1, j+j1, k+k1): # # if we only care about neighbors within our ellipse  then check if this neighbor is within the ellipse
            continue
          #neighbors[(i,j,k)].append((i+i1,j+j1,k+k1))

          milestone[counter].neighbors.append(neighbor_dict[(i+i1,j+j1,k+k1)])

  return milestones

def generate_elliptic_positions(a,b,c,x,y,z, increment, neighbors_in_ellipsoid_only=False):
  '''returns an elliptic grid and neighbors'''
  return generate_cubic_positions(a,b,c,x,y,z, increment, ellipsoid=True, neighbors_in_ellipsoid_only=neighbors_in_ellipsoid_only)

def generate_spherical_positions(r,x,y,z,increment,neighbors_in_ellipsoid_only=False):
  '''returns a spherical grid and neighbors'''
  return generate_elliptic_positions(r,r,r,x,y,z, increment, neighbors_in_ellipsoid_only=neighbors_in_ellipsoid_only)

def generate_planar_positions(origx,origy,origz,normx, normy, normz, lowest, highest, increment, siteid, site_index, k_off, absolute='False'):
  '''returns a series of planar anchors and their neighbors'''
  vector = np.array([normx,normy,normz])
  vector = vector/np.linalg.norm(vector) # normalize the vector
  origin = np.array([origx,origy,origz])
  milestones = []
  planes = np.arange(lowest, highest+increment, increment)
  print "planes: ", planes
  planes_in_site = len(planes)
  for i in range(planes_in_site): # create a position along each increment
    plane = planes[i]
    center = origin + vector * plane
    dimensions = {'centerx':center[0], 'centery':center[1], 'centerz':center[2], 'distance':plane, 'normal':str(tuple(vector)).replace(' ','').replace(',',' ').replace(')','').replace('(','')}
    anchor = tuple(center)
    milestone = Milestone(anchor=anchor, dimensions=dimensions, index=i, siteid=siteid, absolute=absolute, center_type='coord', shape='plane')
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(planes_in_site*int(site_index)), i, siteid, center[0], center[1], center[2])
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    else:
      milestone.end = True

    if i < len(planes) - 1:
      milestone.neighbors.append(i+1)
    else:
      milestone.end = True

    milestones.append(milestone)

  # generate one rotational milestone
  anchor = [1,0,0,0]
  dimensions = {}
  index = len(planes)
  milestone = Milestone(anchor=anchor, dimensions=dimensions, index=index, siteid=siteid, center_type="irrelevant", shape='rotational')
  milestone.fullname = "rotation_%d" % index
  milestones.append(milestone)
  print "number of nonrotational milestones:", index
  return milestones

def generate_concentric_spheres(r,x,y,z,vx,vy,vz,increment, siteid, site_index, k_off, startvx = 0.0, startvy = 0.0, startvz = 0.0, absolute='False', r_low=1):
  '''returns grid/neighbors on concentric spheres along a vector''' # NOTE: not including the origin itself
  vector = np.array([vx,vy,vz])
  startvector = np.array([startvx, startvy, startvz]) # the direction in which we go up to the first radius
  milestones = []
  if np.linalg.norm(vector) > 0.0: vector = vector/np.linalg.norm(vector) # normalize the vector
  if np.linalg.norm(startvector) > 0.0: startvector = startvector/np.linalg.norm(startvector) # normalize the startvector
  vector = vector/np.linalg.norm(vector) # normalize the vector
  origin = np.array([x,y,z])
  sphere_radii = np.arange(r_low,r+increment,increment)
  lowest_radius = sphere_radii[0] # the lowest radius
  spheres_in_site = len(sphere_radii)
  for i in range(spheres_in_site): # create a position along each increment
    radius = sphere_radii[i]
    dimensions = {'centerx':x, 'centery':y, 'centerz':z, 'radius':radius}
    diff_rad = find_rad_of_skewed_vectors(startvector*lowest_radius, vector, radius) # radius - lowest_radius # the difference between the current radius and the lowest radius
    anchor =  tuple(origin + startvector*lowest_radius + vector*diff_rad) # choose the location on the next surface.
    milestone = Milestone(anchor=anchor, dimensions=dimensions, index=i, siteid=siteid, absolute=absolute, center_type='coord')
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(spheres_in_site*int(site_index)), i, siteid, anchor[0], anchor[1], anchor[2])
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    else:
      if not k_off: milestone.end = True

    if i < len(sphere_radii) - 1:
      milestone.neighbors.append(i+1)
    else: # then this is the outermost one, so set BD to true
      milestone.bd = True
      milestone.bd_adjacent = milestones[-1] # make the adjacent milestone the previously created md milestone
      milestone.end = True
      milestone.md = False
    #neighbors = [i-1, i+1] # the milestones to either side

    milestones.append(milestone)

  # generate one rotational milestone
  anchor = [1,0,0,0]
  dimensions = {}
  index = len(sphere_radii)
  milestone = Milestone(anchor=anchor, dimensions=dimensions, index=index, siteid=siteid, center_type="irrelevant", shape='rotational')
  milestone.fullname = "rotation_%d" % index
  milestones.append(milestone)
  return milestones

def generate_concentric_spheres_atom(r, atomid, x, y, z, vx, vy, vz, increment, siteid, site_index, k_off, startvx = 0.0, startvy = 0.0, startvz = 0.0, absolute='False', r_low=1, radius_list=None):
  vector = np.array([vx,vy,vz])
  vector = vector/np.linalg.norm(vector) # normalize the vector
  startvector = np.array([startvx, startvy, startvz]) # the direction in which we go up to the first radius
  milestones = []
  if np.linalg.norm(vector) > 0.0: vector = vector/np.linalg.norm(vector) # normalize the vector
  if np.linalg.norm(startvector) > 0.0: startvector = startvector/np.linalg.norm(startvector) # normalize the startvector
  origin = np.array([x,y,z])
  if radius_list:
    sphere_radii=radius_list
  else:
    sphere_radii = np.arange(r_low,r+increment,increment)
  print 'sphere_radii', sphere_radii
  lowest_radius = sphere_radii[0] # the lowest radius
  total_spherical_milestones = 0
  spheres_in_site = len(sphere_radii)
  for i in range(spheres_in_site): # create a position along each increment
    radius = sphere_radii[i]
    dimensions = {'center_indeces':atomid, 'radius':radius, 'centerx':x, 'centery':y, 'centerz':z}
    diff_rad = find_rad_of_skewed_vectors(startvector*lowest_radius, vector, radius) # radius - lowest_radius # the difference between the current radius and the lowest radius
    anchor =  tuple(origin + startvector*lowest_radius + vector*diff_rad) # choose the location on the next surface.
    milestone = Milestone(anchor=anchor, dimensions=dimensions, index=i, siteid=siteid, center_type='atom')
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    else:
      if not k_off: milestone.end = True

    if i < len(sphere_radii) - 1:
      milestone.neighbors.append(i+1)
    else: # then this is the outermost one, so set BD to true
      milestone.bd = True
      milestone.bd_adjacent = milestones[-1] # make the adjacent milestone the previously created md milestone
      milestone.end = True
      milestone.md = False

    #neighbors = [i-1, i+1] # the milestones to either side
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(spheres_in_site*int(site_index)), i, siteid, anchor[0], anchor[1], anchor[2])
    milestones.append(milestone)
    total_spherical_milestones += 1

  # generate one rotational milestone
  anchor = [1,0,0,0]
  dimensions = {}
  index = len(sphere_radii)
  milestone = Milestone(anchor=anchor, dimensions=dimensions, index=index, siteid=siteid, center_type="irrelevant", shape='rotational')
  milestone.fullname = "rotation_%d" % index
  milestones.append(milestone)
  return milestones

def generate_concentric_spheres_atom_with_rotations(r, atomid, x, y, z, vx, vy, vz, increment, siteid, site_index, k_off, hedron, startvx = 0.0, startvy = 0.0, startvz = 0.0, absolute='False', r_low=1):
  vector = np.array([vx,vy,vz])
  startvector = np.array([startvx, startvy, startvz]) # the direction in which we go up to the first radius
  hedron_quats = positions_orient.get_hedron(hedron) # generates the quaternions themselves from the function in positions_orient
  #print "hedron_quats:", hedron_quats
  milestones = []
  if np.linalg.norm(vector) > 0.0: vector = vector/np.linalg.norm(vector) # normalize the vector
  if np.linalg.norm(startvector) > 0.0: startvector = startvector/np.linalg.norm(startvector) # normalize the startvector
  origin = np.array([x,y,z])
  sphere_radii = np.arange(r_low,r+increment,increment)
  lowest_radius = sphere_radii[0] # the lowest radius
  total_spherical_milestones = 0
  spheres_in_site = len(sphere_radii)
  for i in range(spheres_in_site): # create a position along each increment
    radius = sphere_radii[i]
    dimensions = {'center_indeces':atomid, 'radius':radius, 'centerx':x, 'centery':y, 'centerz':z}
    diff_rad = find_rad_of_skewed_vectors(startvector*lowest_radius, vector, radius) # radius - lowest_radius # the difference between the current radius and the lowest radius
    anchor =  tuple(origin + startvector*lowest_radius + vector*diff_rad) # choose the location on the next surface.
    milestone = Milestone(anchor=anchor, dimensions=dimensions, index=i, siteid=siteid, center_type='atom')
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    else:
      if not k_off: milestone.end = True

    if i < len(sphere_radii) - 1:
      milestone.neighbors.append(i+1)
    else: # then this is the outermost one, so set BD to true
      milestone.bd = True
      milestone.bd_adjacent = milestones[-1] # make the adjacent milestone the previously created md milestone
      milestone.end = True

      milestone.md = False
    #neighbors = [i-1, i+1] # the milestones to either side
    milestones.append(milestone)
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(spheres_in_site*int(site_index)), i, siteid, anchor[0], anchor[1], anchor[2])
    total_spherical_milestones += 1

  # create the rotational milestones
  quat_counter = 0
  for quat in hedron_quats:
    anchor = quat
    dimensions = {}
    index = total_spherical_milestones + quat_counter
    milestone = Milestone(anchor=anchor, dimensions=dimensions, index=index, siteid=siteid, center_type="irrelevant", shape='rotational')
    # add cross and anticross list here...
    #milestone.rotation = quat
    #index = i * len(hedron_quats) + quat_counter # need to construct a unique index for even the rotational milestones
    #new_milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f_%d" % (index, i, siteid, anchor[0], anchor[1], anchor[2], quat_counter) # this is a problem...
    milestone.fullname = "rotation_%d" % index
    milestones.append(milestone)
    quat_counter += 1


  return milestones

def generate_planar_positions_with_rotations(origx,origy,origz,normx, normy, normz, lowest, highest, increment, siteid, site_index, k_off, hedron, absolute='False'):
  '''returns a series of planar anchors and their neighbors'''
  vector = np.array([normx,normy,normz])
  vector = vector/np.linalg.norm(vector) # normalize the vector
  hedron_quats = positions_orient.get_hedron(hedron) # generates the quaternions themselves from the function in positions_orient
  origin = np.array([origx,origy,origz])
  milestones = []
  planes = np.arange(lowest, highest+increment, increment)
  print "planes: ", planes
  total_planar_milestones = 0
  planes_in_site = len(planes)
  for i in range(planes_in_site): # create a position along each increment
    plane = planes[i]
    center = origin + vector * plane
    dimensions = {'centerx':center[0], 'centery':center[1], 'centerz':center[2], 'distance':plane, 'normal':str(tuple(vector)).replace(' ','').replace(',',' ').replace(')','').replace('(','')}
    anchor = tuple(center)
    milestone = Milestone(anchor=anchor, dimensions=dimensions, index=i, siteid=siteid, absolute=absolute, center_type='coord', shape='plane')
    #milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f_%d" % (i, i, siteid, anchor[0], anchor[1], anchor[2], 1)
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(planes_in_site*int(site_index)), i, siteid, anchor[0], anchor[1], anchor[2])
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    else:
      milestone.end = True

    if i < len(planes) - 1:
      milestone.neighbors.append(i+1)
    else:
      milestone.end = True

    milestones.append(milestone)
    total_planar_milestones = 0

  index = len(planes)
  print "number of nonrotational milestones:", index
  # create the rotational milestones
  quat_counter = 0
  for quat in hedron_quats:
    anchor = quat
    dimensions = {}
    index = total_planar_milestones + quat_counter
    milestone = Milestone(anchor=anchor, dimensions=dimensions, index=index, siteid=siteid, center_type="irrelevant", shape='rotational')
    # add cross and anticross list here...
    #milestone.rotation = quat
    #index = i * len(hedron_quats) + quat_counter # need to construct a unique index for even the rotational milestones
    #new_milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f_%d" % (index, i, siteid, anchor[0], anchor[1], anchor[2], quat_counter) # this is a problem...
    milestone.fullname = "rotation_%d" % index
    milestones.append(milestone)
    quat_counter += 1
  return milestones
'''
def generate_hedron_rotations(hedron_type, siteid):
  ''Rotational milestones''
  # NOTE: hedron_neighbors is a dict that contains the indeces of all neighbors within a hedron
  if hedron_type == "simplex":
    hedron_quats = quaternions_simplex
    hedron_neighbors = simplex_neighbors
  elif hedron_type == "tesseract":
    hedron_quats = quaternions_tesseract
    hedron_neighbors = tesseract_neighbors
  #elif hedron_type == ???

  milestones = []
  for i in range(len(hedron_quats)):
    quat = hedron_quats[i]
    anchor = quat
    dimension = {'quat':quat}
    milestone = Milestone(anchor=anchor, dimension=dimension, index=i, sideid=siteid, center_type="quat", shape="rotational")
    milestone.neightbors = hedron_neighbors # append milestone neighbors

    milestones.append(milestone)
  return milestones
'''

def write_milestone_file(site_list, milestone_filename, temperature, md_time_factor, bd_time_factor):
  '''
  given a milestone list and sites, will write an XML file containing information about the milestones and sites
  '''
  ourdoc = Document() # create xml document
  root = ourdoc.createElement("root") # create the xml tree
  ourdoc.appendChild(root)
  #root = ET.Element("root")

  # NOTE: this function needs an overhaul to be inclusive of other milestone types (spherical/planar)
  # Temperature tag
  xmltemp = ourdoc.createElement("temperature") # include temperature
  root.appendChild(xmltemp)
  xmltemptext = ourdoc.createTextNode(str(temperature))
  xmltemp.appendChild(xmltemptext)
  
  # MD time factor
  xmlmd_time_factor = ourdoc.createElement("md_time_factor") # include temperature
  root.appendChild(xmlmd_time_factor)
  xmlmd_time_factor_text = ourdoc.createTextNode(str(md_time_factor))
  xmlmd_time_factor.appendChild(xmlmd_time_factor_text)
  
  # BD time factor
  xmlbd_time_factor = ourdoc.createElement("bd_time_factor") # include temperature
  root.appendChild(xmlbd_time_factor)
  xmlbd_time_factor_text = ourdoc.createTextNode(str(bd_time_factor))
  xmlbd_time_factor.appendChild(xmlbd_time_factor_text)
  
  #root.appendChild(temp_text)
  for site in site_list:
    xmlsite = ourdoc.createElement("site")
    root.appendChild(xmlsite)
    xmlsitename = ourdoc.createElement("name")
    xmlsite.appendChild(xmlsitename)
    xmlsitetext = ourdoc.createTextNode(str(site[0].siteid))
    xmlsitename.appendChild(xmlsitetext)
    
    for milestone in site:
      xmlmilestone = ourdoc.createElement("milestone")
      xmlsite.appendChild(xmlmilestone)
      # write the index
      try:
        xmlattr = ourdoc.createElement("index")
        xmlmilestone.appendChild(xmlattr)
        xmltext = ourdoc.createTextNode(str(milestone.index))
        xmlattr.appendChild(xmltext)
      except NameError:
        continue
      for attribute in ['fullname', 'directory', 'shape', 'anchor', 'end', 'bd', 'md',]: # NOTE: more can be added
        eval_str = '''try:
  xmlattr = ourdoc.createElement("%s")
  xmlmilestone.appendChild(xmlattr)
  xmltext = ourdoc.createTextNode(str(milestone.%s))
  xmlattr.appendChild(xmltext)
except AttributeError:
  pass''' % (attribute, attribute)
        exec(eval_str)
      if milestone.shape == "sphere":
        for dim_attribute in ['radius', ]: # NOTE: more can be added
          eval_str = '''try:
  xmlattr = ourdoc.createElement("%s")
  xmlmilestone.appendChild(xmlattr)
  xmltext = ourdoc.createTextNode(str(milestone.dimensions['%s']))
  xmlattr.appendChild(xmltext)
except AttributeError:
  pass''' % (dim_attribute, dim_attribute)
          exec(eval_str)
      elif milestone.shape == "rotational":
        for dim_attribute in [ ]: # NOTE: more can be added
          eval_str = '''try:
  xmlattr = ourdoc.createElement("%s")
  xmlmilestone.appendChild(xmlattr)
  xmltext = ourdoc.createTextNode(str(milestone.dimensions['%s']))
  xmlattr.appendChild(xmltext)
except AttributeError:
  pass''' % (dim_attribute, dim_attribute)
          exec(eval_str)
      elif milestone.shape == "plane":
        for dim_attribute in [ 'normal', 'distance', ]: # NOTE: more can be added
          eval_str = '''try:
  xmlattr = ourdoc.createElement("%s")
  xmlmilestone.appendChild(xmlattr)
  xmltext = ourdoc.createTextNode(str(milestone.dimensions['%s']))
  xmlattr.appendChild(xmltext)
except AttributeError:
  pass''' % (dim_attribute, dim_attribute)
          exec(eval_str)

  xml_string = ourdoc.toprettyxml(indent="  ")
  #print 'xml_string:', xml_string
  milestone_file = open(milestone_filename,'w')
  milestone_file.write(xml_string)
  milestone_file.close()
  return

def split_milestones_by_site(milestone_list): # splits the given list into a list of sites, each of which contain a list of milestones
  site_list = []
  cur_site = []
  cur_siteid = milestone_list[0].siteid
  for milestone in milestone_list: # for every milestone
    if milestone.siteid != cur_siteid: # if we have a new siteid
      site_list.append(cur_site) # then append off this siteid list and begin anew
      cur_site = []
    cur_site.append(milestone) # at any rate, add this milestone to the current list of site id's
  site_list.append(cur_site)
  return site_list

def main(settings):
  master_temperature = settings['master_temperature']
  milestone_list = []
  site_list = [] # convenient for writing the milestone file
  milestone_filename = settings['milestone_filename']
  site_index=0
  for site in settings['sites']:

    args = ['site_index=%d'%site_index] # start off by defining the function argument for the index of the site
    for arg in site['dimensions'].keys():
      args.append("%s=%s" % (arg, str(site['dimensions'][arg])))

    if site['anchor_function'] == 'ellipsoid':
      command = "generate_elliptic_positions(" + ','.join(args) + ")"
    elif site['anchor_function'] == 'planes':
      command = "generate_planar_positions(" + ','.join(args) + ")"
    elif site['anchor_function'] == 'concentric_spheres':
      command = "generate_concentric_spheres(" + ','.join(args) + ")"
    elif site['anchor_function'] == 'concentric_spheres_atom':
      command = "generate_concentric_spheres_atom(" + ','.join(args) + ")"
    elif site['anchor_function'] == 'concentric_spheres_atom_with_rotations':
      command = "generate_concentric_spheres_atom_with_rotations(" + ','.join(args) + ")"
    elif site['anchor_function'] == 'planes_with_rotations':
      command = "generate_planar_positions_with_rotations(" + ','.join(args) + ")"
    #elif site['anchor_function'] == 'planez':
    #  command = "generate_planar_positions(" + ','.join(args) + ")"
    #elif etc etc
    else:
      raise Exception("Warning: anchor function %s not yet implemented in milestones.py" % site['anchor_function'])

    if verbose: print "Generating anchors/milestones using command:", command

    milestones = eval(command) # execute this string command as if it were code
    milestone_list += milestones
    site_list.append(milestones)
    site_index += 1

  #write_milestone_file(site_list, milestone_filename, master_temperature) # commented out because this will be written later when the filetree is created and we know the names of the directories

  return milestone_list



class Test_milestones(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self): # test this function
    print "WARNING: this module does not have comprehensive unit tests"

  def test_in_ellipsoid(self):
    self.assertTrue(in_ellipsoid(1,1,1,0,0,0,0.5,0,0))
    self.assertFalse(in_ellipsoid(1,1,1,0,0,0,1.5,0,0))
    self.assertFalse(in_ellipsoid(1,1,1,0,0,0,1,1,1))
    self.assertTrue(in_ellipsoid(5,4,3,7,8,9,10,8,9))
    self.assertFalse(in_ellipsoid(5,4,3,7,8,9,10,3,9))
    self.assertFalse(in_ellipsoid(5,4,3,7,8,9,-12,8,9))

  def test_quadratic_formula(self):
    self.assertEqual(quadratic_formula(1.0,4.0,3.0), [-3.0,-1.0]) # test negative roots
    self.assertEqual(quadratic_formula(3.0,-13.0,12.0), [4.0/3.0,3.0]) # test positive roots
    self.assertEqual(quadratic_formula(1.0,-2.0,0.0), [0.0,2.0]) # test zero roots
    self.assertEqual(quadratic_formula(1.0,0.0,0.0), [0.0,0.0]) # test zero roots
    self.assertEqual(quadratic_formula(1.0,0.0,4.0), [complex(0,-2.0),complex(0,2.0)]) # test imaginary roots
    return

  def test_find_rad_of_skewed_vectors(self):
    A = np.array([0.0,1.0,0.0])
    B = np.array([1.0,0.0,0.0])
    B2 = np.array([1.0,1.0,0.0])
    self.assertEqual(find_rad_of_skewed_vectors(A, A, 2.0), 1.0) # check a vector against itself
    self.assertEqual(find_rad_of_skewed_vectors(A, B, 2.0), sqrt(3.0)) # check a vector against right angle
    self.assertEqual(find_rad_of_skewed_vectors(2.0*A, B*2.5, 3.0), sqrt(5.0)) # check a vector against right angle
    self.assertEqual(find_rad_of_skewed_vectors(B2, A, 2.0), sqrt(3.0)-1.0) # check a vector against 45 degree angle
    return

  #def generate_cubic_positions(a, b, c, x, y, z, increment, ellipsoid=False, neighbors_in_ellipsoid_only=False):
  def test_generate_cubic_positions(self):
    print "Warning: no unittests written for function 'generate_cubic_positions' because it's unused in the current version of SEEKR."
    return # temporary until we get this unittest fixed
    #test_milestones = generate_cubic_positions(1.0,1.0,1.0, 0.0, 0.0, 0.0, 1.0)
    #self.assertEqual(len(test_milestones),7)
    #self.assertEqual(len(test_milestones[0].neighbors),1)
    #self.assertEqual(test_milestones[0].anchor,(-1.0, 1))
    #return

  #def generate_elliptic_positions(a,b,c,x,y,z, increment, neighbors_in_ellipsoid_only=False):
  def test_generate_elliptic_positions(self):
    print "Warning: no unittests written for function 'generate_elliptic_positions' because it's unused in the current version of SEEKR."

  #def generate_spherical_positions(r,x,y,z,increment,neighbors_in_ellipsoid_only=False):
  def test_generate_spherical_positions(self):
    print "Warning: no unittests written for function 'generate_spherical_positions' because it's unused in the current version of SEEKR."

  #def generate_planar_positions(origx,origy,origz,normx, normy, normz, lowest, highest, increment, siteid, site_index, k_off, absolute='False'):
  def test_generate_planar_positions(self):
    test_milestones = generate_planar_positions(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -10.0, 10.0, 2.0, "0", 1, False)
    self.assertEqual(len(test_milestones),12)
    self.assertEqual(test_milestones[0].index, 0)
    self.assertEqual(test_milestones[1].end, False)
    self.assertEqual(test_milestones[0].end, True)

  #def generate_concentric_spheres(r,x,y,z,vx,vy,vz,increment, siteid, site_index, k_off, startvx = 0.0, startvy = 0.0, startvz = 0.0, absolute='False', r_low=1):
  def test_generate_concentric_spheres(self):
    #return # temporary until we get this unittest fixed
    test_milestones = generate_concentric_spheres(10.0, 0, 0, 0, 1, 1, 1, 2.0, 1, 1, False, r_low=2.0)
    self.assertEqual(len(test_milestones),6)
    self.assertTrue(test_milestones[0].end)
    self.assertFalse(test_milestones[2].end)
    return

  #def generate_concentric_spheres_atom(r, atomid, x, y, z, vx, vy, vz, increment, siteid, site_index, k_off, startvx = 0.0, startvy = 0.0, startvz = 0.0, absolute='False', r_low=1):
  def test_generate_concentric_spheres_atom(self):
    test_milestones = generate_concentric_spheres_atom(10.0, 234, 0, 0, 0, 1, 1, 1, 2.0, 1, 1, False, r_low=2.0)
    self.assertEqual(len(test_milestones),6)
    self.assertTrue(test_milestones[0].end)
    self.assertFalse(test_milestones[2].end)
    return

  #def generate_concentric_spheres_atom_with_rotations(r, atomid, x, y, z, vx, vy, vz, increment, siteid, site_index, k_off, hedron, startvx = 0.0, startvy = 0.0, startvz = 0.0, absolute='False', r_low=1):
  def test_generate_concentric_spheres_atom_with_rotations(self):
    test_milestones = generate_concentric_spheres_atom_with_rotations(10.0, 234, 0, 0, 0, 1, 1, 1, 2.0, 1, 1, False, 'tesseract', r_low=2.0)
    self.assertEqual(len(test_milestones),13)
    self.assertTrue(test_milestones[0].end)
    self.assertFalse(test_milestones[2].end)
    self.assertEqual(test_milestones[12].shape, 'rotational')
    return

  #def generate_planar_positions_with_rotations(origx,origy,origz,normx, normy, normz, lowest, highest, increment, siteid, site_index, k_off, hedron, absolute='False'):
  def test_generate_planar_positions_with_rotations(self):
    test_milestones = generate_planar_positions_with_rotations(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -10.0, 10.0, 2.0, "0", 1, False, 'tesseract')
    self.assertEqual(len(test_milestones),19)
    self.assertEqual(test_milestones[0].index, 0)
    self.assertEqual(test_milestones[1].end, False)
    self.assertEqual(test_milestones[0].end, True)
    self.assertEqual(test_milestones[15].shape, 'rotational')

  #def write_milestone_file(site_list, milestone_filename, temperature):
  def test_write_milestone_file(self):
    expected_output = '''<?xml version="1.0" ?>\n\n<root>\n\n  <temperature>\n\n    312.0\n\n  </temperature>\n\n  <site>\n\n    <milestone>\n\n      <index>\n\n        0\n\n      </index>\n\n      <fullname>\n\n        2_0_0_0.0_0.0_8.0\n\n      </fullname>\n\n      <directory>\n\n        \n\n      </directory>\n\n      <shape>\n\n        plane\n\n      </shape>\n\n      <anchor>\n\n        (0.0, 0.0, 8.0)\n\n      </anchor>\n\n      <end>\n\n        True\n\n      </end>\n\n      <bd>\n\n        False\n\n      </bd>\n\n      <normal>\n\n        0.0 0.0 1.0\n\n      </normal>\n\n      <distance>\n\n        8.0\n\n      </distance>\n\n    </milestone>\n\n    <milestone>\n\n      <index>\n\n        1\n\n      </index>\n\n      <fullname>\n\n        3_1_0_0.0_0.0_10.0\n\n      </fullname>\n\n      <directory>\n\n        \n\n      </directory>\n\n      <shape>\n\n        plane\n\n      </shape>\n\n      <anchor>\n\n        (0.0, 0.0, 10.0)\n\n      </anchor>\n\n      <end>\n\n        True\n\n      </end>\n\n      <bd>\n\n        False\n\n      </bd>\n\n      <normal>\n\n        0.0 0.0 1.0\n\n      </normal>\n\n      <distance>\n\n        10.0\n\n      </distance>\n\n    </milestone>\n\n    <milestone>\n\n      <index>\n\n        0\n\n      </index>\n\n      <fullname>\n\n        rotation_0\n\n      </fullname>\n\n      <directory>\n\n        \n\n      </directory>\n\n      <shape>\n\n        rotational\n\n      </shape>\n\n      <anchor>\n\n        [1.0, 0.0, 0.0, 0.0]\n\n      </anchor>\n\n      <end>\n\n        False\n\n      </end>\n\n      <bd>\n\n        False\n\n      </bd>\n\n    </milestone>\n\n  </site>\n\n</root>\n'''
    test_filename = '/tmp/test_milestone.xml'
    test_milestones = generate_planar_positions_with_rotations(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 8.0, 10.0, 2.0, "0", 1, False, 'single')
    write_milestone_file([test_milestones], test_filename, 312.0)
    test_file = open(test_filename, 'r')
    test_file_lines = '\n'.join(test_file.readlines())
    test_file.close()
    self.assertEqual(test_file_lines, expected_output)

  #def split_milestones_by_site(milestone_list)
  '''
  def test_split_milestones_by_site(self):
    test_milestones1 = generate_concentric_spheres_atom(10.0, 234, 0, 0, 0, 1, 1, 1, 2.0, 1, 1, False, r_low=2.0)
    test_milestones2 = generate_concentric_spheres_atom(10.0, 234, 10, 10, 10, 1, 1, 1, 2.0, 2, 1, False, r_low=2.0)
    test_milestones_all = test_milestones1 + test_milestones2
    site_list = split_milestones_by_site(test_milestones_all)
    self.assertEqual(site_list, [test_milestones1, test_milestones2])
  '''
  def test_split_milestones_by_site(self):
    print "Warning: no unittests written for function 'split_milestones_by_site' because it is unused in the current version of SEEKR"

if __name__=='__main__':
  print "Now testing milestones.py"
  unittest.main() # run tests of all functions
