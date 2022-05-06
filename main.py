from scipy.interpolate import CubicSpline
import time
from tkinter import *
import math
import numpy as np
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol

dt = 0.1
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
        self.seg = []
        #self.coord = [vec(-30,3), vec(-20,2), vec(-10,5)]

def binSearch(numDot1, numDot2, t):
    delt=t
    while c*delt != l:
        l = math.sqrt((dots[numDot1].x[t]-dots[numDot2].x[t-delt])**2+(dots[numDot1].y[t]-numDot2.seg[0](t-delt)**2))
        if c*delt > l:
            delt -= delt / 2
        elif c*delt < l:
            delt += delt / 2
    return delt


def acc(n1,n2,t):
    t0 = binSearch(n1,n2,t)
    a_x = (G * dots[n2].m * (dots[n2].x[t0] - dots[n1].x[t]))/(abs(dots[n2].x[t0] - dots[n1].x[t])**3)
    a_y = (G * dots[n2].m * (dots[n2].y[t0] - dots[n1].y[t]))/(abs(dots[n2].y[t0] - dots[n1].y[t])**3)
    return [a_x,a_y]


def VerletInt(numDot,t): #xn+1 = 2xn - xn-1 +an*dt^2
    vx_t1 = 2*dots[numDot].x[t] - dots[numDot].x[t-1] + acc(numDot,n)[0]*dt**2
    vy_t1 = 2*dots[numDot].y[t] - dots[numDot].y[t-1] + acc(numDot,n)[1]*dt**2
    return [vx_t1,vy_t1]

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
    dotNum = 3
    x_s = [dots[numDot].x[dotNum-3],dots[numDot].x[dotNum-2],dots[numDot].x[dotNum-1]]
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

xs = [dt  * i for i in  range (-1000,500)]
#ys = [f(dots[0].x, dots[0].y, t) for t in xs
for a in range(n):
    dot_new = dot()
    dots.append(dot_new)
    dots[a].seg.append(s(a))

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
plt.show()
#print(dots[n-1].x)
#print(dots[n-1].y)
#print("end")
wind.mainloop()
