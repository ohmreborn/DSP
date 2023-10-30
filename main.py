import plotly.express as px
import numpy as np
import pandas as pd
import math
class constant:
    def __init__(self):
        self.earth_mass = 5.97219 * 10**24
        self.year = 31556926
        self.distance_earth_sun = 1.5*10**11 # m
        self.G = 6.67*10**-11

        self.k = 1

class Star:
    def __init__(self,
                 mass,
                 position,
                 velocity,
                 name):
        self.m = mass
        self.r = np.array(position)
        self.v = np.array(velocity)
        self.name = name
    def __str__(self):
        return self.name

class Universe:
    def __init__(self,system,time,dt=1e-5):
        self.dt = dt
        self.time = time
        self.system = system

        self.all_pos = []

    def compute(self):
        acclelation = []
        for star1 in self.system:
            a = np.array([0.0,0.0])
            for star2 in self.system:
                if star1.name != star2.name:
                    dis = star2.r-star1.r
                    abs_dis = np.linalg.norm(dis)
                    a += cons.k*star2.m*dis/abs_dis**3
            acclelation.append(a)

        r = []
        for a,star in zip(acclelation,self.system):
            star.v = star.v + a*self.dt
            star.r = star.r + star.v*self.dt
            r.append(star.r)
        return r
        
    def simulate(self):
        t = 0
        while self.time > t:
            pos = self.compute()
            self.all_pos.append(pos)
            t += dt
        
        self.all_pos = np.stack(self.all_pos)
        data = []
        for ind,iter in enumerate(self.all_pos):
            for i in range(len(iter)):
                data.append([iter[i][0],iter[i][1],system[i].name,ind])
        self.all_pos = pd.DataFrame(data,columns=["x","y","star","frame"])
    def animate(self):
        min_range_x = (min(self.all_pos["x"])-1)*3
        max_range_x = (max(self.all_pos["x"])+1)*3
        min_range_y = (min(self.all_pos["y"])-1)*3
        max_range_y = (max(self.all_pos["y"])+1)*3
        range_x = [min_range_x,max_range_x]
        range_y = [min_range_y,max_range_y]
        fig = px.scatter(self.all_pos, x="x", y="y", animation_frame="frame", animation_group="star",range_x=range_x,range_y=range_y,color="star")
        fig.show()

    def show(self):
        fig = px.line(self.all_pos,x="x",y="y",color="star")
        fig.show()

    

cons = constant()
# the input
sun = Star(333000,[0,0],[0,0],"sun")
earth = Star(1,[1,0],[0,(cons.k*sun.m)**0.5],"earth")
earth2 = Star(1000,[-0.5,0],[0,-300],"earth1")

system = [earth,sun]
time = 0.01#2*math.pi/(cons.k*sun.m)**0.5
dt = 5e-5#time/200
universe = Universe(system,time,dt=dt)
universe.simulate()
universe.all_pos.to_csv("python.csv",index=False)
