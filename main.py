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
c = 3*10**2
xc = 600
yc = 300
n = 3
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
        self.m = 100
        self.time = [0, dt, 2*dt, 3*dt]
        self.x = x
        self.y = y
        self.seg_x = []
        self.seg_y = []


def binSearch(numDot1, numDot2, t):
    delt=t
    l = -1
    i = int(t/dt)
    while round(c*delt,2) != round(l,2):
        b =  int((t-delt)/dt)
        sgx = dots[numDot2].seg_x[b](t-delt)
        sgy = dots[numDot2].seg_y[b](t-delt)
        x02=(dots[numDot1].x[i]- sgx)**2
        y02=(dots[numDot1].y[i]- sgy)**2
        l = math.sqrt(x02 + y02)
        if c*delt > l:
            delt -= delt / 2
        elif c*delt < l:
            delt += delt / 2
       # print('c*dt ',delt*c)
       #print('l ',l)
    return delt


def acc(n1,t):
    a_x = 0
    a_y = 0
    j = int(t/dt)
    for i in range(n):
        if i != n1:
            t0 = t - binSearch(n1,i,t)
            b = int(t0/dt)
            #print(t0)
            r3 = (math.sqrt((dots[i].seg_x[b](t0) - dots[n1].x[j])**2+(dots[i].seg_y[b](t0) - dots[n1].y[j])**2)**3)
            a_x += (dots[i].m * (dots[i].seg_x[b](t0) - dots[n1].x[j]))/r3
            a_y += (dots[i].m * (dots[i].seg_y[b](t0) - dots[n1].y[j]))/r3

    return G * a_x,G * a_y


def VerletInt(numDot,t): #xn+1 = 2xn - xn-1 +an*dt^2
    i = int(t/dt)
    Wx,Wy = acc(numDot,t)
    x_t1 = 2*dots[numDot].x[i] - dots[numDot].x[i-1] + Wx*dt**2
    y_t1 = 2*dots[numDot].y[i] - dots[numDot].y[i-1] + Wy*dt**2
    return x_t1,y_t1

#Сплайны для сегмента
def sx(numDot,a):
    dotNum = max(a,4)
    t_s = dots[numDot].time[-4:]
    #print(t_s)
    x_s = [dots[numDot].x[dotNum-4],dots[numDot].x[dotNum-3],dots[numDot].x[dotNum-2],dots[numDot].x[dotNum-1]]
    x_tck = CubicSpline(t_s,x_s)
    #print(x_tck(0))
    return x_tck
def sy(numDot,a):
    dotNum = max(a,4)
    t_s = dots[numDot].time[-4:]
    y_s = [dots[numDot].y[dotNum-4],dots[numDot].y[dotNum-3],dots[numDot].y[dotNum-2],dots[numDot].y[dotNum-1]]
    y_tck = CubicSpline(t_s,y_s)
    #print(y_tck(0))
    return y_tck

#Скорости? ускорения?
#двоичный поиск
#гравитационные силы


for a in range(n):

    stx = [-4+5*(-1)**n, -3+5*(-1)**n, -2+5*(-1)**n ,-1+5*(-1)**n]
    sty = [3+n, 2+n, 1+n ,n]
    dot_new = dot(stx,sty)
    dots.append(dot_new)
    for i in range(0,4):
        dots[a].seg_x.append(sx(a,0))
        dots[a].seg_y.append(sy(a,0))
#print('bin s ',binSearch(0,1,dt))
#print('acc ',acc(1,3*dt))
for i in range (3,500):
    t = dt * i
    for j in range(n):
        x1,y1=VerletInt(j,t)
        dots[j].x.append(x1)
        dots[j].y.append(y1)
        dots[j].time.append(t+dt)
        dots[j].seg_x.append(sx(j,i))
        dots[j].seg_y.append(sy(j,i))
        layer.moveto(dots[j].ball, 5*x1-7+xc,5*y1+yc-7)
        layer.create_oval(5*x1 + xc, 5*y1 + yc, 5*x1 + xc, 5*y1 + yc, width=1)  # траектория
    wind.update()
    time.sleep(0.001)

#fig, ax = plt.subplots()
#ax.plot(ts, ert(0)(ts), label='sc')
#ax.plot(dots[0].x, dots[0].y, 'o', label='data')
#ax.plot(ts, sx(0)(ts))
#ax.plot(ts, sy(0)(ts))
#plt.show()
#print(dots[n-1].x)
#print(dots[n-1].y)
#print("end")
wind.mainloop()
