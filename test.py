# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 09:39:27 2019

@author: Harry
"""

import raytracer as ray
import numpy as np
import matplotlib.pyplot as plt
import importlib
importlib.reload(ray)

#create rays
r1 = ray.Ray(p=np.array([0.1,0.1,0]), k=np.array([0,0,1]))

#create bundles
b = ray.Bundle(n = 5, r = 5, k = np.array([0,0,1]), p = np.array([0,0,0]))


#create optical elements:
lens = ray.SphericalRefraction(z_0 = 105, r = 50, r_a = 50, n2 = 1, n1 = 1.5168)
plane = ray.PlaneRefraction(z_0 = 100, n = np.array([0,0,-1]), n2 = 1.5168, n1 = 1)
outputplane = ray.OutputPlane(z_0 = 200.86 , n=np.array([0,0,-1]))

#change the order in the 'elements' list to change the order of optical elements that the rays propagate through.
#outputplane should always be last in the list
elements = [plane, lens, outputplane]

#propagate bundles 
for elem in elements:
    b.propagate_bundle(elem)

b.plot2D()     #can only plot 1 of plot2D or plotXY at a time
#b.plotXY()
#b.plot3D()
print('the rms is', b.rms())
ray.focal2(plane, lens)

#propagate rays
for elem in elements:
    elem.propagate_ray(r1)
#ray.plot2D(r1)
