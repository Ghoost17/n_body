from scipy.interpolate import CubicSpline
import time
from tkinter import *
import math
import numpy as np
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol


G = 6.67*10**(-11)
c = 3*10**8
xc = 600
yc = 300
n = 3
g = 9.8
dots = []
t = 0

wind = Tk()
layer = Canvas(wind, width=2*xc, height=2*yc)
layer.pack()
layer.create_line(6, 0, 6, 2*yc, width=1)
layer.create_line(0, yc, 2*xc, yc, width=1)
layer.create_oval( 3, yc - 3, 9, yc + 3, fill='red', width=0)

class vec:
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def __add__(self, other):
        return vec(self.x + other.x, self.y + other.y)


    def __sub__(self, other):
        return vec(self.x - other.x, self.y - other.y)

    def __abs__(self):
        return math.hypot(self.x, self.y)


class dot:
    def __init__(self):
        self.ball = layer.create_oval( 3, yc - 6, 15, yc + 6, fill='blue', width=0)
        self.m = 1
        self.x = [-30, -20, -10]
        self.y = [3, 2, 5]
        #self.coord = [vec(-30,3), vec(-20,2), vec(-10,5)]



def VerletInt(numDot,t):
    h = 1
    vx_t = (dots[numDot].x[t+h] - dots[numDot].x[t-h])/2*h
    vy_t = (dots[numDot].y[t+h] - dots[numDot].y[t-h])/2*h
    return [vx_t,vy_t]

for a in range(n):
    dot_new = dot()
    dots.append(dot_new)

#Гравитационная сила для двух тел
def AbsForce(d1,d2,t1,t2):
    return G*(dots[d1].m*dots[d2].m/dist(d1,t1,d2,t2)**2)

#Расстояние между двумя точками
def dist(d1,t1,d2,t2):
    distance = math.sqrt(
         (dots[d1].x[t1] - dots[d2].x[t2])**2+
         (dots[d1].y[t1] - dots[d2].y[t2])**2
    )
    return distance
#Сплайны для сегмента
def s(numDot):
    x_s = [dots[numDot].x[0],dots[numDot].x[1],dots[numDot].x[2]]
    y_s = [dots[numDot].y[0],dots[numDot].y[1],dots[numDot].y[2]]
    tck = CubicSpline(x_s,y_s)
    #print(tck.c)
    return tck


def ds(numDot,arg,dx):
    tck = CubicSpline(dots[numDot].x,dots[numDot].y)
    return tck(arg,dx)
#Скорости? ускорения?
#двоичный поиск
#гравитационные силы


def binSerch(numDot1, numDot2, t):

    t0 = math.sqrt((dots[numDot1].x[t])**2+(dots[numDot1].y[t])**2)/c
    return t0


xs = [0.1  * i for i in  range (-1000,500)]
#ys = [f(dots[0].x, dots[0].y, t) for t in xs


for i in xs:
    for j in range(n):
        layer.create_oval((j+1)*i+xc, s(j)(i)+yc, (j+1)*i+xc, s(j)(i)+yc, width=1) #траектория
        layer.moveto(dots[j].ball, (j+1)*i-7+xc,s(j)(i)+yc-7)
        # add X Y in memory
        #if ((j+1)*i - int((j+1)*i)) == 0:
        #  dots[j].y.append(round(int(s(j,i)),3))
        #  dots[j].x.append((j+1)*i)
    wind.update()
    time.sleep(0.001)


fig, ax = plt.subplots()
ax.plot(dots[n-1].x, dots[n-1].y, 'o', label='data')
ax.plot(xs, s(n-1)(xs))
ax.plot(xs, ds(n-1,xs,1))
#plt.show()
#print(dots[n-1].x)
#print(dots[n-1].y)
#print("end")
wind.mainloop()
