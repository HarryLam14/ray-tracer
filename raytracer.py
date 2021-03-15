# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 09:27:42 2019

@author: Harry
"""
from mpl_toolkits.mplot3d import Axes3D #used to get the 3D projection
import matplotlib.pyplot as plt
import numpy as np

def getx(vertices):
        """
        get x values of a ray
        vertices: 2D array of the vertices of the ray
        """
        x = []
        n = len(vertices)
        for i in range(n):
            x.append(vertices[i][0])
        return x
def gety(vertices):
        """
        get y values of a ray
        vertices: 2D array of the vertices of the ray
        """
        y = []
        n = len(vertices)
        for i in range(n):
            y.append(vertices[i][1])
        return y
def getz(vertices):
        """
        get the z values of a ray
        vertices: 2D array of the vertices of the ray
        """
        z = []
        n=len(vertices)
        for i in range(n):
            z.append(vertices[i][2])
        return z
def plot2D(ray):
        """plot the path of the ray in 2D"""
        x = getx(ray.vertices())
        z = getz(ray.vertices())
        plt.plot(z,x, '-o')
        plt.xlabel('z /mm')
        plt.ylabel('x /mm')
        plt.grid()
        plt.show()
    
def plot3D(ray):
        """plot the path of the ray in 3D"""
        x = getx(ray.vertices())
        y = gety(ray.vertices())
        z = getz(ray.vertices())
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        ax.plot(z,y,x, '-o' )
        ax.set_xlabel('z (optical axis)')
        ax.set_ylabel('y')
        ax.set_zlabel('x')
        plt.show()
def plotXY(ray):
    """plot a spot diagram (cross section) of the ray"""
    x = ray.p()[0]
    y = ray.p()[1]
    plt.plot(x,y, 'o')
    plt.xlabel('x /mm')
    plt.ylabel('y /mm')
    plt.grid()
    plt.axis('equal')
    plt.show()
def plot3D2(ray1, ray2):
        """plot the path of the ray in 3D"""
        x = getx(ray1.vertices())
        y = gety(ray1.vertices())
        z = getz(ray1.vertices())
        x1 = getx(ray2.vertices())
        y1 = gety(ray2.vertices())
        z1 = getz(ray2.vertices())
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        ax.plot(z,y,x, '-o' )
        ax.plot(z1,y1,x1, '-o')
        ax.set_xlabel('z (optical axis)')
        ax.set_ylabel('y')
        ax.set_zlabel('x')
        plt.show()

class Ray:
    """
    Class representing an optical ray
    """
    def __init__(self, p = np.array([0.,0.,0.]), k = np.array([0.,0.,0.])):
        """
        initialise the ray
        p: starting point
        k: direction vector
        """
        self._p = np.array(p)
        self._k = np.array(k)
        self._p_all = [np.array(p)]
        if len(p) != 3:
            raise Exception("parameter p must be a 3-element array")
        if len(k) != 3:
            raise Exception("parameter k must be a 3-element array")
    def p(self):
        """
        return current point of the ray
        """
        return(self._p)
    def k(self):
        """
        return current ray direction
        """
        return(self._k)
    def append(self, p_new, k_new):
        """
        append a new point and direction to the ray
        p_new: array of the new point
        k_new: array of the new k vector
        """
        self._p_all.append(np.array(p_new))        
        self._k = k_new
        self._p = p_new
    def vertices(self):
        """
        return all points along the ray
        """
        return(np.array(self._p_all))
    def __repr__(self):
        return "%s(p = %r, k=%r)" % ("Ray", self._p, self._k)
    def __str__(self):
        return "(%g, %s)" % (self._p, self._k)

class OpticalElement:
    def propagate_ray(self, ray):
        "propagate a ray through the optical element"
        raise NotImplementedError()
        

class SphericalRefraction(OpticalElement):
    def __init__(self, z_0, r, r_a, n1, n2):
        """
        Initialise the spherical refraction surface
        z_0: z intercept of the surface
        r: radius of the surface. Sign determines curvature of the surface
        r_a: aperture radius
        n1: refractive index 1, material before boundary
        n2: refractive index 2, material after boundary
        """
        self._r_c = 1/(r)
        self._z_0 = z_0
        self._r_a = r_a
        self._n1 = n1
        self._n2 = n2
        self._O = np.array([0,0,z_0 + r])
        self._R = r
    
    def intercept(self, ray):
        """
        return the intercept of a ray with a spherical surface
        """
        x = ray.p() - self._O
        k_hat = ray.k()/(np.sqrt(np.dot(ray.k(),ray.k())))
        D = ((np.dot(x,k_hat)*np.dot(x,k_hat))-((np.dot(x,x)) - (self._R*self._R)))
        if D < 0:
            return None
        if self._r_c < 0:
            l = -np.dot(x,k_hat) + np.sqrt(D)
        else:
            l = -np.dot(x,k_hat) - np.sqrt(D)
        if self._r_c is 0:
            u = self._z_0*(ray.k()[0]/ray.k()[2])
            y = self._z_0*(ray.k()[1]/ray.k()[2])
            return (np.array(u,y,self._z_0))
        d = ray.p()+l*k_hat
            
        
        if np.sqrt(d[0]*d[0]+d[1]*d[1]) > self._r_a:
            pass
        else: return np.array(d)        
    
    def Snell(self, ray, n_hat = np.array([0,0,0])):
        """
        implementing Snell's Law
        n_hat: normal to the surface 
        """
        n1 = self._n1
        n2 = self._n2
        k1 = ray.k()
        k1_hat = k1/(np.sqrt(np.dot(k1,k1)))
        r = n1/n2
        c = -np.dot(n_hat, k1_hat)
        p = np.sqrt(1-(np.dot(k1_hat, n_hat))**2)
        if p > 1/r:                                 
            return None
        else:
            k2 = r*k1_hat + (r*c-np.sqrt(1-r*r*(1-c*c)))*n_hat
            k2_hat = k2/(np.sqrt(np.dot(k2,k2)))
            return k2_hat
    
    def propagate_ray(self, ray):
        """propagate the ray for a spherical refraction surface"""
        Q = self.intercept(ray)
        if Q is None:
            pass
        else:
            if self._r_c < 0:
                n = -(Q - self._O)
            else:
                n = Q - self._O
            n1 = np.array(n)
            n2 = n1/(np.sqrt(np.dot(n1,n1)))
            k2 = self.Snell(ray, n_hat=n2)
            ray.append(p_new = Q, k_new = k2)

class PlaneRefraction(OpticalElement):
    def __init__(self, z_0, n, n1, n2):
        """
        z_0: z intercept of the plane
        n: normal of the plane
        n1: refractive index 1
        n2: refractive index 2
        """
        self._z_0 = np.array([0,0,z_0])
        self._n = np.array(n)
        self._n1 = n1
        self._n2 = n2
    def intercept(self, ray):
        """find the intercept of the ray and the plane"""
        l_n = np.dot(ray.k(),self._n)
        k_hat = ray.k()/(np.sqrt(np.dot(ray.k(),ray.k())))
        if l_n is 0:
            return None
        else:
            p_l = self._z_0-ray.p()
            d = (np.dot(p_l,self._n))/(l_n)
            Q = ray.p() + d*k_hat
            return Q     
    def Snell(self, ray, n_hat = np.array([0,0,0])):
        """implement Snell's Law"""
        n1 = self._n1
        n2 = self._n2
        k1 = ray.k()
        k1_hat = k1/(np.sqrt(np.dot(k1,k1)))
        r = n1/n2
        c = -np.dot(n_hat, k1_hat)
        p = np.sqrt(1-(np.dot(k1_hat, n_hat))**2)
        if p > 1/r:
           return None
        else:
            k2 = r*k1_hat + (r*c-np.sqrt(1-r*r*(1-c*c)))*n_hat
            k2_hat = k2/(np.sqrt(np.dot(k2,k2)))
            return k2_hat
        
    def propagate_ray(self, ray):
        """propagate a ray into an output plane"""
        Q = self.intercept(ray)
        if Q is None:
            pass
        else:
            k2 = self.Snell(ray, self._n)
            ray.append(p_new = Q, k_new = k2)
    
    
class OutputPlane(OpticalElement):
    def __init__(self, z_0, n = np.array([0,0,0])):
        """
        z_0: z intercept of output plane
        n: normal of the output plane
        """
        self._z_0 = np.array([0,0,z_0])
        self._n = np.array(n)
    def intercept(self, ray):
        """find the intercept of the ray and the output plane"""
        if ray.k() is None:
            raise Exception("No k value")
        l_n = np.dot(ray.k(), self._n)
        k_hat = ray.k()/(np.sqrt(np.dot(ray.k(),ray.k()))) 
        if l_n is 0:
            return None
        else:
            p_l = self._z_0-ray.p()
            d = (np.dot(p_l,self._n))/(l_n)
            Q = ray.p() + d*k_hat
            return Q
    def propagate_ray(self, ray):
        """propagate a ray into an output plane"""
        Q = self.intercept(ray)
        if Q is None:
            raise Exception ("No Intercept")
        else:
            ray.append(p_new = Q, k_new = None)
        
def distr(n, r):
    """evenly distribute points in a circle
    n: number of rings
    r: radius of the circle
    """
    r_b = r/n
    for i in range(1, n+1):
        theta = 0
        R = r_b*i
        for p in range((i*6)+1):
            if p == 0:
                yield 0., 0.
            else:
                r, t = R, theta
                theta += (2*np.pi)/(i*6)
                yield r, t

class Bundle():       
    def __init__(self, n, r, k = np.array([]), p = np.array([])):
        """
        initalise bundle
        n: number of concentric circles
        r: radius of bundle
        p: starting point of bundle
        k: direction vector of bundle
        """
        self._k = np.array(k)
        self._p = np.array(p)
        self._n = n
        self._r = r
        self._rb = r/n      #radius of each bundle
        self._rays = []     #list containing the rays in the bundle
        
        x = []
        y = []
        
        for r, t in distr(n, r):
            x.append(r*np.cos(t))
            y.append(r*np.sin(t))
        for i in range(len(x)):
            self._rays.append(Ray(p = np.array([self._p[0]+x[i], self._p[1]+y[i], self._p[2]]), k = k))
        
    def propagate_bundle(self, refract):
        """propagate a bundle,
        refract: optical element
        """
        for i in self._rays:
           refract.propagate_ray(i)
    def rays(self):
        """return the list of rays in the bundle"""
        return self._rays
    def plot2D(self):
        """plot the path of all the rays in the bundle in 2D"""
        for i in self._rays:
            plot2D(i)
        plt.show()
    def plot3D(self):
         """plot the path of all the rays in the bundle in 3D"""
         fig = plt.figure()
         ax = fig.add_subplot(111, projection = '3d')
         for i in self._rays:
             x = getx(i.vertices())
             y = gety(i.vertices())
             z = getz(i.vertices())
             ax.plot(z,y,x, '-o' )
         ax.set_xlabel('z (optical axis) /mm')
         ax.set_ylabel('y /mm')
         ax.set_zlabel('x /mm')
         plt.show()
    def plot3D2(self, other):
         """plot the path of all the rays from two bundles in 3D"""
         fig = plt.figure()
         ax = fig.add_subplot(111, projection = '3d')
         for i in self._rays:
             x = getx(i.vertices())
             y = gety(i.vertices())
             z = getz(i.vertices())
             ax.plot(z,y,x, '-o' )
         for i in other._rays:
             x = getx(i.vertices())
             y = gety(i.vertices())
             z = getz(i.vertices())
             ax.plot(z,y,x, '-o' )    
         ax.set_xlabel('z (optical axis)')
         ax.set_ylabel('y')
         ax.set_zlabel('x')
         plt.show()
    def plotXY(self):
        """plot the spot diagram for the bundle"""
        for i in self._rays:
            plotXY(i)
    def rms(self):
        """find the RMS radius of the points in the bundle"""
        r = []
        xn = []
        yn = []
        for i in self._rays:
            xn.append(i.p()[0])
        xm = np.mean(xn)
        for i in self._rays:
            yn.append(i.p()[1])
        ym = np.mean(yn)
        for i in self._rays:
            x = i.p()[0]
            y = i.p()[1]
            r.append(((x*x)-xm)+((y*y)-ym))
        m = np.mean(r)
        return np.sqrt(m)


def focal1(s1):
    """
    function to find the focal length for a single spherical refracting surface
    s1: SphericalRefraction type
    """
    f = s1._r/(s1._n2-1)
    print(f)
def focal2(s1, s2):
    """
    function to find the focal length for a construction of two refracting surfaces,
    assumes one is a plane and the other is a spherical surface.
    Only works if the rays converge to a focal point.
    s1 must be before s2 on the optical axis, i.e. s1._z_0<s2._z_0
    s1: first refractive surface
    s2: second refractive surface
    """
    if type(s1) is PlaneRefraction:
        f = ((s1._z_0[2]+s2._z_0)/2)-(s2._R)/(s2._n1-1)
    if type(s2) is PlaneRefraction:
        f = ((s1._z_0+s2._z_0[2])/2)+(s1._R)/(s1._n2-1)
    print('the focal length is',f)
        
    
      
        
        
        
        
        
        
        
        