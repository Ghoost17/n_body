from scipy.interpolate import CubicSpline
import time
from tkinter import *
import math
import numpy as np
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol

dt = 0.05
G = 6.67*10**(-3)
c = 3*10**4
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
    def __init__(self,x,y):
        self.ball = layer.create_oval( 3, yc - 6, 15, yc + 6, fill='blue', width=0)
        self.m = 1
        self.time = [-3*dt, -2*dt, -dt, 0]
        self.x = x
        self.y = y
        self.seg_x = []
        self.seg_y = []


def binSearch(numDot1, numDot2, t):
    delt=t
    l = -1
    while round(c*delt) != round(l):

        l = math.sqrt((dots[numDot1].x[t]-sx(numDot2)(t-delt))**2+(dots[numDot1].y[t]-sy(numDot2)(t-delt)**2))
        if c*delt > l:
            delt -= delt / 2
        elif c*delt < l:
            delt += delt / 2
        print('c*dt ',delt*c)
        print('l ',l)
    return delt


def acc(n1,t):
    a_x = 0
    a_y = 0
    for i in range(n):
        if i != n1:
            t0 = binSearch(n1,i,t)
            print(t0)
            a_x += (G * dots[i].m * (sx(i)(t0) - dots[n1].x[t]))/\
                   (math.sqrt((sx(i)(t0) - dots[n1].x[t])**2+(sy(i)(sx(i)(t0)) - dots[n1].y[t])**2)**3)
            a_y += (G * dots[i].m * (sy(i)(sx(i)(t0)) - dots[n1].y[t]))/\
                   (math.sqrt((sx(i)(t0) - dots[n1].x[t])**2+(sy(i)(sx(i)(t0)) - dots[n1].y[t])**2)**3)
    return [a_x,a_y]
print(acc(1,0))

def VerletInt(numDot,t): #xn+1 = 2xn - xn-1 +an*dt^2
    vx_t1 = 2*dots[numDot].x[t] - dots[numDot].x[t-1] + acc(numDot,t)[0]*dt**2
    vy_t1 = 2*dots[numDot].y[t] - dots[numDot].y[t-1] + acc(numDot,t)[1]*dt**2
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
def sx(numDot):
    dotNum = 4
    t_s = dots[numDot].time
    x_s = [dots[numDot].x[dotNum-4],dots[numDot].x[dotNum-3],dots[numDot].x[dotNum-2],dots[numDot].x[dotNum-1]]
    x_tck = CubicSpline(t_s,x_s)
    #print(x_tck(0))
    return x_tck
def sy(numDot):
    dotNum = 3
    t_s = dots[numDot].time
    y_s = [dots[numDot].y[0],dots[numDot].y[1],dots[numDot].y[2],dots[numDot].y[3]]
    y_tck = CubicSpline(t_s,y_s)
    #print(y_tck(0))
    return y_tck

def ds(numDot,arg,dx):
    tck = CubicSpline(dots[numDot].x,dots[numDot].y)
    return tck(arg,dx)
#Скорости? ускорения?
#двоичный поиск
#гравитационные силы


for a in range(n):
    stx = [-0.03+n/2, -0.02+n/3, -0.01+n/5 ,0]
    sty = [0.03+n/2, 0.02+n/3, -0.01+n/5 ,0]
    dot_new = dot(stx,sty)
    dots.append(dot_new)
    dots[a].seg_x.append(sx(a))
    dots[a].seg_y.append(sy(a))


ts = [dt  * i for i in  range (0,500)]
for i in ts:
    for j in range(n):
        layer.create_oval(sx(j)(i)+xc, sy(j)(i)+yc, sx(j)(i)+xc, sy(j)(i)+yc, width=1) #траектория
        layer.moveto(dots[j].ball, sx(j)(i)-7+xc,sy(j)(i)+yc-7)
        # add X Y in memory
        #if ((j+1)*i - int((j+1)*i)) == 0:
        #  dots[j].y.append(round(int(s(j,i)),3))
        #  dots[j].x.append((j+1)*i)
    wind.update()
    time.sleep(0.001)

def ert(numDot):
    dotNum = 3
    x_s = [dots[numDot].x[dotNum-3],dots[numDot].x[dotNum-2],dots[numDot].x[dotNum-1]]
    y_s = [dots[numDot].y[0],dots[numDot].y[1],dots[numDot].y[2]]
    tck = CubicSpline(x_s,y_s)
    #print(tck(0))
    return tck
fig, ax = plt.subplots()
#ax.plot(ts, ert(0)(ts), label='sc')
ax.plot(dots[0].x, dots[0].y, 'o', label='data')
ax.plot(ts, sx(0)(ts))
ax.plot(ts, sy(0)(sx(0)(ts)))
plt.show()
#print(dots[n-1].x)
#print(dots[n-1].y)
#print("end")
wind.mainloop()
