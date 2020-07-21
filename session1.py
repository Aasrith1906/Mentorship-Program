from vpython import *
import numpy as np


ground = box(color=color.red)
ground.pos = vector(0,-1.5,0)
ground.size = (22,1,2)
s = sphere(color = color.blue)
s.velocity = vector(10,0,0)
s.pos = vector(-10,0,0)
t = 0
dt = 0.01

while t<2:

    s.pos = s.pos + s.velocity*dt

    t += dt

    rate(1/dt)
