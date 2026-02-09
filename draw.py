import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import sys

class Injector:

    def __init__(self, rho, phi, z0, theta, dtheta, r0, Rinjection, E, Z, M, sign):
        self.rho = rho
        self.phi = phi
        self.z0 = z0
        self.theta = theta
        self.dtheta = dtheta
        self.r0 = r0
        self.Rinjection = Rinjection
        self.E = E
        self.Z = Z
        self.M = M
        self.sign = sign

        self.C = 3e10


    def Omega(self, B):
        return self.Z*B/(self.M*self.C)
    
    def getV(self):
        return np.sqrt(2.*self.E / self.M)
    
    def getVPerp(self):
        return self.getV()*np.sin(self.theta)

    def getVx(self):
        return -self.getVPerp()*np.sin(self.phi) * self.sign
    
    def getVy(self):
        return self.getVPerp()*np.cos(self.phi) * self.sign
    
    def getVz(self):
        return self.getV()*np.cos(self.theta)

    def getRhoLarmor(self, B):
        return self.getVPerp()/self.Omega(B)

    def getDeltaX(self, B):
        return self.getVy() / self.Omega(B)
    
    def getDeltaY(self, B):
        return -self.getVx() / self.Omega(B)

class Draw:

    def __volume(self, iz, ir, iphi):
        
        phi1 = self.phiArray[iphi]
        phi2 = self.phiArray[iphi+1]

        dphi = phi2 - phi1

        z1 = self.zArray[iz]
        z2 = self.zArray[iz+1]

        r1 = self.rArray[ir]
        r2 = self.rArray[ir+1]

        V = dphi * (z2-z1) * (r2*r2 - r1*r1) / 2.

        return V

    def __getIndex(self, iz, ir, iphi):
        return ir + self.nr*iz + self.nz*self.nr*iphi    
    
    def __getColor(self, value : float, min_value, max_value):
        normalized = (value - min_value) / (max_value - min_value) if max_value != min_value else 0.5
        normalized = max(0.0, min(1.0, normalized))
        
        # Получаем цветовую карту
        cmap = cm.get_cmap(self.colormap)
        
        # Получаем цвет
        return cmap(normalized)

    def __readAxis(self, file, data, type):
        line = ''

        while not type + '-axis' in line:
            line = file.readline()

        line = file.readline()

        n = int([i for i in line.split()][-1])
        
        for i in range(n+1):
            line = file.readline()
            value = float([i for i in line.split()][-1])
            data[type + 'Axis'].append(value)

    def __readFileOut(self, filename : str):

        data = {
            'zAxis' : [],
            'rAxis' : [],
            'phiAxis' : [],
            'injectors' : [],
            'intersection' : [],
            'ncapture' : [],
            'ncaptureCenter' : [],
            'intersectionCenter' : [],
            'sarray' : [],
            'nF' : [],
            'nFCenter' : []
        }

        with open(filename, 'r') as file:
            
            line = file.readline()

            line = file.readline()

            self.normaN = float(line.split('=')[-1])

            self.__readAxis(file, data, 'z')
            self.__readAxis(file, data, 'r')
            self.__readAxis(file, data, 'phi')

            nInjector=0

            while not 'count end' in line:

                if 'Bcenter' in line:
                    self.Bcenter = float(line.split("=")[-1])
                if 'injector' in line and not 'injector end' in line:

                    line = file.readline()
                    line = file.readline()
                    theta = float([i for i in line.split('=')][-1]) * np.pi / 180.
                    line = file.readline()
                    line = file.readline()
                    r0 = float([i for i in line.split('=')][-1])
                    line = file.readline()
                    dtheta = float(line.split("=")[-1]) * np.pi / 180.
                    line = file.readline()
                    Rinjection = float(line.split('=')[-1])
                    line = file.readline()
                    line = file.readline()
                    rho = float([i for i in line.split('=')][-1])
                    line = file.readline()
                    z0 = float([i for i in line.split('=')][-1])
                    line = file.readline()
                    phi = float([i for i in line.split('=')][-1]) * np.pi / 180.

                    line = file.readline()
                    line = file.readline()

                    E = float(line.split('=')[-1])*1.6e-12

                    line = file.readline()

                    Z = int(line.split('=')[-1])*4.8e-10

                    line = file.readline()

                    M = int(line.split('=')[-1])*1.67e-24


                    line = file.readline()

                    if '[plus-direction]' in line:
                        sign = 1
                    else:
                        sign = -1

                    rmax = data['rAxis'][-1]
                    print("injection %d l: %f" %(nInjector, 2.*np.sqrt(rmax*rmax-rho*rho) / np.sin(theta)))
                    nInjector+=1

                    data['injectors'].append(Injector(rho, phi, z0, theta, dtheta, r0, Rinjection, E, Z, M, sign))


                line = file.readline()


            nz = len(data['zAxis'])-1
            nr = len(data['rAxis'])-1
            nphi = len(data['phiAxis'])-1

            while not 'nCapture=' in line:
                line = file.readline()

                if 'result injector' in line:
                    line = file.readline()
                    ns = int(line.split("=")[-1])
                    data['sarray'].append(np.zeros(3*ns+1))

                    for i in range(ns):
                        line = file.readline()
                        value = [float(i) for i in line.split()]
                        data['sarray'][-1][1+i] = value[1]
                        data['sarray'][-1][1+i+ns] = value[2]
                        data['sarray'][-1][1+i+2*ns] = value[3]

            self.ncap = float( line.split(("="))[-1].split("%")[0] )

            for iphi in range(nphi):
                line = file.readline()
                for iz in range(nz):
                    
                    line = file.readline()
                    values = [float(i) for i in line.split()]

                    for ir in  range(nr):
                        data['ncapture'].append(float(values[4*ir]))
                        data['ncaptureCenter'].append(values[4*ir+1])
                        data['intersection'].append(bool(values[4*ir+2]))
                        data['intersectionCenter'].append(bool(values[4*ir+3]))
                        if values[4*ir+2]:
                            print("iz:", iz, "ir:", ir, "iphi:", iphi)
            

            while not 'n from r:' in line:
                line = file.readline()

            data['nF'].append(0)
            data['nFCenter'].append(0)
            for ir in range(nr):
                line = file.readline()
                data['nF'].append(float(line.split()[0]))
                data['nFCenter'].append(float(line.split()[1]))
            data['nF'][0] = data['nF'][1]
            data['nFCenter'][0] = data['nFCenter'][1]




        return data

    def __drawLine(self, ax, rho : float, phi : float, theta : float, z0 : float, 
                 Rmax : float, sign : int, injector = None):

        t = np.sqrt(Rmax*Rmax - rho*rho) / np.sin(theta)

        x0 = rho*np.cos(phi)
        y0 = rho*np.sin(phi)

        sx = -np.sin(theta)*np.sin(phi)*sign
        sy = np.sin(theta)*np.cos(phi)*sign
        sz = np.cos(theta)

        ax.plot3D([z0 - 1.5*t*sz, z0 + t*sz], [x0 - 1.5*t * sx, x0 + t*sx], [y0 - 1.5*t * sy, y0 + t*sy], 
                  self.linestyleInjector, linewidth=self.linewidthInjector, color=self.colorInjector)
        
        if not injector == None and self.drawLarmorCenterLine:
            deltaX = injector.getDeltaX(self.Bcenter)
            deltaY = injector.getDeltaY(self.Bcenter)
            ax.plot3D([z0 - 1.5*t*sz, z0 + t*sz], [x0 - 1.5*t * sx + deltaX, x0 + t*sx 
                                                   + deltaX], [y0 - 1.5*t * sy + deltaY, y0 + t*sy + deltaY], 
                  '--', linewidth=self.linewidthInjector, color=self.colorInjector)

        if self.drawPoints:
            ax.scatter3D([z0-t*sz, z0+t*sz], [x0-t*sx, x0+t*sx], [y0-t*sy, y0+t*sy], s=self.pointInjectorSize, color=self.colorInjector)
        
        if self.drawQuiver:
            ax.quiver(z0-0.5*t*sz, 
                    x0 - 0.5*t*sx, 
                    y0 - 0.5*t*sy,
                    0.5*t*sz, 0.5*t*sx, 0.5*t*sy,
                    color=self.colorInjector, arrow_length_ratio=0.3, 
                    linewidth=self.linewidthInjector)

    def __drawSimpleCircleSurface(self, ax, rho, phi, z0, r0, alpha=0.3):
        angles = np.linspace(0, 2*np.pi, 50)
        radii = np.linspace(0, r0, 25)
        
        angles_grid, radii_grid = np.meshgrid(angles, radii)
        
        center_x = rho * np.cos(phi)
        center_y = rho * np.sin(phi)
        
        X = center_x + radii_grid * np.cos(angles_grid)
        Y = center_y + radii_grid * np.sin(angles_grid)
        Z = np.ones_like(X) * z0
        
        ax.plot_surface(Z, X, Y, color=self.colorInjector, alpha=alpha)

    def __drawCircle(self, ax, rho : float, phi : float, z0 : float, r0 : float, linestyle='-', linewidth=0.5):

        phi0 = np.linspace(0, np.pi*2, self.Npoints)
        ax.plot3D(np.ones(self.Npoints)*z0, rho*np.cos(phi)+r0*np.cos(phi0), rho*np.sin(phi) + r0*np.sin(phi0), linestyle, color=self.colorInjector, linewidth=linewidth)

    def __drawMesh(self, ax, value, intersection):

        for iphi in range(self.nphi):
            for iz in range(self.nz):
                for ir in range(self.nr):
                    r1 = self.rArray[ir]
                    r2 = self.rArray[ir+1]
                    z1 = self.zArray[iz]
                    z2 = self.zArray[iz+1]
                    phi1 = self.phiArray[iphi]
                    phi2 = self.phiArray[iphi+1]

                    phi = np.linspace(phi1, phi2, self.Npoints)

                    linewidth = self.linewidthMesh
                    color = self.colorMesh

                    if intersection[self.__getIndex(iz, ir, iphi)]:
                        linewidth = self.linewidthInter
                        color = self.__getColor(value[self.__getIndex(iz, ir, iphi)], 0, max(value))
                    elif not self.drawNotInter:
                        continue

                    ax.plot3D(np.ones(self.Npoints)*self.zArray[iz], r1*np.cos(phi), r1*np.sin(phi), color=color, linewidth=linewidth)
                    ax.plot3D(np.ones(self.Npoints)*self.zArray[iz], r2*np.cos(phi), r2*np.sin(phi), color=color, linewidth=linewidth)

                    ax.plot3D(np.ones(self.Npoints)*self.zArray[iz+1], r1*np.cos(phi), r1*np.sin(phi), color=color, linewidth=linewidth)
                    ax.plot3D(np.ones(self.Npoints)*self.zArray[iz+1], r2*np.cos(phi), r2*np.sin(phi), color=color, linewidth=linewidth)


                    if (r1 != 0):
                        ax.plot3D([z1, z2], [r1*np.cos(phi1), r1*np.cos(phi1)], [r1*np.sin(phi1), r1*np.sin(phi1)], color=color, linewidth=linewidth)
                        ax.plot3D([z1, z2], [r1*np.cos(phi2), r1*np.cos(phi2)], [r1*np.sin(phi2), r1*np.sin(phi2)], color=color, linewidth=linewidth)
                    ax.plot3D([z1, z2], [r2*np.cos(phi1), r2*np.cos(phi1)], [r2*np.sin(phi1), r2*np.sin(phi1)], color=color, linewidth=linewidth)
                    ax.plot3D([z1, z1], [r1*np.cos(phi1), r2*np.cos(phi1)], [r1*np.sin(phi1), r2*np.sin(phi1)], color=color, linewidth=linewidth)
                    ax.plot3D([z1, z2], [r2*np.cos(phi2), r2*np.cos(phi2)], [r2*np.sin(phi2), r2*np.sin(phi2)], color=color, linewidth=linewidth)
                    ax.plot3D([z2, z2], [r1*np.cos(phi1), r2*np.cos(phi1)], [r1*np.sin(phi1), r2*np.sin(phi1)], color=color, linewidth=linewidth)
                    ax.plot3D([z1, z1], [r1*np.cos(phi2), r2*np.cos(phi2)], [r1*np.sin(phi2), r2*np.sin(phi2)], color=color, linewidth=linewidth)
                    ax.plot3D([z2, z2], [r1*np.cos(phi2), r2*np.cos(phi2)], [r1*np.sin(phi2), r2*np.sin(phi2)], color=color, linewidth=linewidth)

    def draw3D(self):
        if len(self.zArray) == 0 or len(self.rArray) == 0 or len(self.phiArray) == 0:
                print("не правильная сетка")
                return

        fig = plt.figure(figsize=(200, 200))
        ax = fig.add_subplot(111, projection='3d')

        if (self.drawMesh):
            #self.__drawMesh(ax, self.concetrationCenter, self.intersectionCenter)
            self.__drawMesh(ax, self.concetration, self.intersection)

        listRho = []
        for injector in self.injectors:

            rho = injector.rho
            phi = injector.phi
            theta = injector.theta
            z0 = injector.z0
            r0 = injector.r0
            sign = injector.sign

            if (self.drawInjectorLine):
                self.__drawLine(ax, rho, phi, theta, z0, injector.Rinjection, sign, injector) # рисуем линию инжекции
            if (self.drawCircleRho and not rho in listRho and rho != 0.):
                self.__drawCircle(ax, 0, 0, z0, rho, linestyle="--", linewidth=1) # рисуем окружность радиусом rho
                listRho.append(rho)
            if r0 != 0 and self.drawInjectorCircle:
                self.__drawCircle(ax, rho, phi, z0, r0) # рисумем разброс начальной точки
                if (self.drawSurfeCircle):
                    self.__drawSimpleCircleSurface(ax, rho, phi, z0, r0)


        if self.drawAxis:
            ax.plot3D([self.zArray[0], self.zArray[-1]], [0, 0], [0, 0], self.linestyleAxis, 
                      color=self.colorAxis, linewidth=self.linewidthAxis) # рисуем ось


        #настройка оси
        ax.set(xlabel='z, см', ylabel='x, см', zlabel='y, см')
        coeff = 1.2
        ax.set_xlim(coeff*self.zArray[0], coeff*self.zArray[-1])
        ax.set_ylim(-coeff*self.rArray[-1], coeff*self.rArray[-1])
        ax.set_zlim(-coeff*self.rArray[-1], coeff*self.rArray[-1])

    def drawCapFromLine(self):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        ax1 = axes[0]
        ax2 = axes[1]

        for sarray in self.sarray:
            ns = int((len(sarray) - 1) / 3)

            x = np.zeros(ns+1)
            y = np.zeros(ns+1)
            z = np.zeros(ns+1)
            sumY = 0

            for i in range(ns):
                x[i+1] = sarray[1+i]
                y[i+1] = sarray[1+i+1*ns]*self.normaN / 1e13
                z[i+1] = sarray[1+i+2*ns] / (sarray[1+i] - sarray[i])
                
                sumY += sarray[1+i+2*ns]

            y[0] = y[1]
            z[0] = z[1]
            ax1.step(x, z)
            ax2.step(x, y)
        ax1.grid()
        ax1.set(xlabel='s, см', ylabel='ncap, см^-1')
        ax2.grid()
        ax2.set(xlabel='s, см', ylabel='ni, 10^13 см^-3')

    def drawCapFromR(self, nf, nFCenter=None):
        fig, axes = plt.subplots(nrows=1, ncols=1)
        r = self.rArray
        ax1 = axes
        ax1.step(r, nf)
        if not nFCenter == None:
            ax1.step(r, nFCenter)
        ax1.grid()
        ax1.set(xlabel='r, см', ylabel='nf, см^-3')

    def __init__(self, fileName, 
                 drawMesh=True, drawInjectorLine=True, drawInjectorCircle=True, 
                 drawCircleRho=True, drawAxis=True, drawNotInter=True,
                 Npoints=1000, colormap='plasma', colorInjector='red', 
                 linestyleInjector='-', drawQuiver=True, drawPoints=True,
                 linewidthInjector=2, colorAxis='black', linewidthAxis=2, 
                 linestyleAxis='--', pointInjectorSize=40, drawSurfeCircle=True,
                 colorMesh='black', linewidthMesh=0.5, linewidthIner=1, drawLarmorCenterLine=True
                 ):
        
        data = self.__readFileOut(fileName)
        self.zArray = data['zAxis']
        self.rArray = data['rAxis']
        self.phiArray = data['phiAxis']
        self.nCapture = data['ncapture']
        self.injectors = data['injectors']
        self.intersection = data['intersection']
        self.intersectionCenter = data['intersectionCenter']
        self.nCaptureCenter = data['ncaptureCenter']
        self.sarray = data['sarray']
        self.nF = data['nF']
        self.nFCenter = data['nFCenter']
        self.nz = len(self.zArray) - 1
        self.nr = len(self.rArray) - 1
        self.nphi = len(self.phiArray) - 1

        # while len(self.intersection) < self.nr*self.nz*self.nphi:
        #     self.intersection.append(False)
        # while len(self.nCapture) < self.nr*self.nz*self.nphi:
        #     self.nCapture.append(0.)


        self.concetration = np.zeros(self.nphi*self.nz*self.nr)
        self.concetrationCenter = np.zeros(self.nphi*self.nz*self.nr)

        for iphi in range(self.nphi):
            for iz in range(self.nz):
                for ir in range(self.nr):
                    self.concetration[self.__getIndex(iz, ir, iphi)] = self.nCapture[self.__getIndex(iz, ir, iphi)] / self.__volume(iz, ir, iphi)
                    self.concetrationCenter[self.__getIndex(iz, ir, iphi)] = self.nCaptureCenter[self.__getIndex(iz, ir, iphi)] / self.__volume(iz, ir, iphi)

        self.drawMesh = drawMesh
        self.drawInjectorLine = drawInjectorLine
        self.drawInjectorCircle = drawInjectorCircle
        self.drawCircleRho = drawCircleRho
        self.drawAxis = drawAxis
        self.drawNotInter = drawNotInter

        self.Npoints = Npoints
        self.colormap = colormap

        self.colorInjector = colorInjector
        self.linewidthInjector = linewidthInjector
        self.linestyleInjector = linestyleInjector

        self.colorAxis = colorAxis
        self.linewidthAxis = linewidthAxis
        self.linestyleAxis = linestyleAxis

        self.drawQuiver = drawQuiver
        self.drawPoints = drawPoints
        self.pointInjectorSize = pointInjectorSize
        self.drawSurfeCircle = drawSurfeCircle

        self.colorMesh = colorMesh
        self.linewidthMesh = linewidthMesh
        self.linewidthInter = linewidthIner
        self.drawLarmorCenterLine = drawLarmorCenterLine

        print("\nncap =", self.ncap, '%')

    def drawPointLine(self, fileName, step=1):

        fig = plt.figure(figsize=(200, 200))
        ax = fig.add_subplot(111, projection='3d')

        rho = []
        z = []
        theta = []
        phi = []
        sign = []

        with open(fileName, 'r') as f:

            data = [i.split() for i in f]

            for i in range(len(data)):
                rho.append(float(data[i][0]))
                theta.append(float(data[i][1]))
                z.append(float(data[i][2]))
                phi.append(float(data[i][3]))
                sign.append(int(data[i][4]))

        index = 0
        print('start')
        for r0, z0, theta0, phi0,sign0 in zip(rho, z, theta, phi, sign):
            if index % step == 0:
                self.__drawLine(ax, r0, phi0, theta0, z0, 15.6, sign0)
            index+=1


        coeff = 1.2
        ax.set(xlabel='z, см', ylabel='x, см', zlabel='y, см')
        ax.set_xlim(coeff*self.zArray[0], coeff*self.zArray[-1])
        ax.set_ylim(-coeff*self.rArray[-1], coeff*self.rArray[-1])
        ax.set_zlim(-coeff*self.rArray[-1], coeff*self.rArray[-1])


    def show(self, draw3d=True, drawFromLine=True, drawCapFromR=True, drawPoints=False):

        if draw3d:
            self.draw3D()
        if drawFromLine:
            self.drawCapFromLine()
        if drawCapFromR:
            self.drawCapFromR(self.nF, self.nFCenter)
        if drawPoints:
            self.drawPointLine("build/points.txt", 5000)
        
        plt.show()


if __name__ == '__main__':

    fileName = sys.argv[1]

    draw = Draw(fileName, drawNotInter=False, drawLarmorCenterLine=True)

    draw.show(drawCapFromR=True, drawFromLine=True, draw3d=False, drawPoints=False)