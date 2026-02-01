import matplotlib.pyplot as plt
import numpy as np


class Injector:

    def __init__(self, rho, phi, z0, theta, r0):
        self.rho = rho
        self.phi = phi
        self.z0 = z0
        self.theta = theta
        self.r0 = r0


def drawGrid(ax, zArray : list, rArray : list, phiArray : list, nz : int, nr : int, nphi : int, color="black", linewidth=0.5, N=1000):

    for iphi in range(nphi):
        for iz in range(nz):
            for ir in range(nr):
                r1 = rArray[ir]
                r2 = rArray[ir+1]
                z1 = zArray[iz]
                z2 = zArray[iz+1]
                phi1 = phiArray[iphi]
                phi2 = phiArray[iphi+1]

                phi = np.linspace(phi1, phi2, N)

                ax.plot3D(np.ones(N)*zArray[iz], r1*np.cos(phi), r1*np.sin(phi), color=color, linewidth=linewidth)
                ax.plot3D(np.ones(N)*zArray[iz], r2*np.cos(phi), r2*np.sin(phi), color=color, linewidth=linewidth)
                if iz == nz-1:
                    ax.plot3D(np.ones(N)*zArray[iz+1], r1*np.cos(phi), r1*np.sin(phi), color=color, linewidth=linewidth)
                    ax.plot3D(np.ones(N)*zArray[iz+1], r2*np.cos(phi), r2*np.sin(phi), color=color, linewidth=linewidth)


                if (r1 != 0):
                    ax.plot3D([z1, z2], [r1*np.cos(phi1), r1*np.cos(phi1)], [r1*np.sin(phi1), r1*np.sin(phi1)], color=color, linewidth=linewidth)
                    ax.plot3D([z1, z2], [r1*np.cos(phi2), r1*np.cos(phi2)], [r1*np.sin(phi2), r1*np.sin(phi2)], color=color, linewidth=linewidth)
                ax.plot3D([z1, z2], [r2*np.cos(phi1), r2*np.cos(phi1)], [r2*np.sin(phi1), r2*np.sin(phi1)], color=color, linewidth=linewidth)
                ax.plot3D([z1, z1], [r1*np.cos(phi1), r2*np.cos(phi1)], [r1*np.sin(phi1), r2*np.sin(phi1)], color=color, linewidth=linewidth)
                ax.plot3D([z1, z2], [r2*np.cos(phi2), r2*np.cos(phi2)], [r2*np.sin(phi2), r2*np.sin(phi2)], color=color, linewidth=linewidth)
                ax.plot3D([z2, z2], [r1*np.cos(phi1), r2*np.cos(phi1)], [r1*np.sin(phi1), r2*np.sin(phi1)], color=color, linewidth=linewidth)
                ax.plot3D([z1, z1], [r1*np.cos(phi2), r2*np.cos(phi2)], [r1*np.sin(phi2), r2*np.sin(phi2)], color=color, linewidth=linewidth)
                ax.plot3D([z2, z2], [r1*np.cos(phi2), r2*np.cos(phi2)], [r1*np.sin(phi2), r2*np.sin(phi2)], color=color, linewidth=linewidth)
                


def drawLine(ax, rho : float, phi : float, theta : float, z0 : float, Rmax : float, color="red", linewidth=2, linestyle="-", drawPoints=True, drawQuiver=True):

    t = np.sqrt(Rmax*Rmax - rho*rho) / np.sin(theta)

    x0 = rho*np.cos(phi)
    y0 = rho*np.sin(phi)

    sx = -np.sin(theta)*np.sin(phi)
    sy = np.sin(theta)*np.cos(phi)
    sz = np.cos(theta)

    ax.plot3D([z0 - 1.5*t*sz, z0 + t*sz], [x0 - 1.5*t * sx, x0 + t*sx], [y0 - 1.5*t * sy, y0 + t*sy], linestyle, linewidth=linewidth, color=color)
    if drawPoints:
        ax.scatter3D([z0-t*sz, z0+t*sz], [x0-t*sx, x0+t*sx], [y0-t*sy, y0+t*sy], s=40, color=color)
    
    if drawQuiver:
        ax.quiver(z0-0.5*t*sz, 
                x0 - 0.5*t*sx, 
                y0 - 0.5*t*sy,
                0.5*t*sz, 0.5*t*sx, 0.5*t*sy,
                color=color, arrow_length_ratio=0.3, linewidth=linewidth)


def drawSimpleCircleSurface(ax, rho, phi, z0, r0, color="red", alpha=0.3):
    angles = np.linspace(0, 2*np.pi, 50)
    radii = np.linspace(0, r0, 25)
    
    angles_grid, radii_grid = np.meshgrid(angles, radii)
    
    center_x = rho * np.cos(phi)
    center_y = rho * np.sin(phi)
    
    X = center_x + radii_grid * np.cos(angles_grid)
    Y = center_y + radii_grid * np.sin(angles_grid)
    Z = np.ones_like(X) * z0
    
    ax.plot_surface(Z, X, Y, color=color, alpha=alpha)

def drawCircle(ax, rho : float, phi : float, z0 : float, r0 : float, color="red", drawSurf=True, linestyle="-", linewidth=0.5, N = 1000):

    phi0 = np.linspace(0, np.pi*2, N)
    ax.plot3D(np.ones(N)*z0, rho*np.cos(phi)+r0*np.cos(phi0), rho*np.sin(phi) + r0*np.sin(phi0), linestyle, color=color, linewidth=linewidth)
    if drawSurf:
        drawSimpleCircleSurface(ax, rho, phi, z0, r0, color)



def draw(zArray : list, rArray : list, phiArray : list, injectorArray : list,
         drawGrid0=True, drawLine0=True, drawInjectionCircle0=True, drawCirle0=True, drawAxis0=True):

    fig = plt.figure(figsize=(200, 200))
    ax = fig.add_subplot(111, projection='3d')

    if (drawGrid0):
        drawGrid(ax, zArray, rArray, phiArray, len(zArray)-1, len(rArray)-1, len(phiArray)-1) # рисуем сетку

    listRho = []
    for injector in injectorArray:

        rho = injector.rho
        phi = injector.phi
        theta = injector.theta
        z0 = injector.z0
        r0 = injector.r0

        if (drawLine0):
            drawLine(ax, rho, phi, theta, z0, rArray[-1], color="red") # рисуем линию инжекции
        if (drawInjectionCircle0 and not rho in listRho):
            drawCircle(ax, 0, 0, z0, rho, drawSurf=False, linestyle="--", linewidth=1) # рисуем окружность радиусом rho
            listRho.append(rho)
        if r0 != 0 and drawCirle0:
            drawCircle(ax, rho, phi, z0, r0, color="red") # рисумем разброс начальной точки

    if drawAxis0:
        ax.plot3D([zArray[0], zArray[-1]], [0, 0], [0, 0], "--", color="black", linewidth=2) # рисуем ось


    #настройка оси
    ax.set(xlabel='z, см', ylabel='x, см', zlabel='y, см')
    ax.set_xlim(1.5*zArray[0], 1.5*zArray[-1])
    ax.set_ylim(-2*rArray[-1], 2*rArray[-1])
    ax.set_zlim(-2*rArray[-1], 2*rArray[-1])


    plt.show()



zArray = [-4, -2, 2, 4]
rArray = [0, 0.5, 1, 2]
phiArray = [0, np.pi/2., np.pi, 3/2*np.pi, 2*np.pi]

rho = 1.
z0 = 0
phi = 0
theta = np.pi / 4
r0 = 0

inArray = [Injector(rho, phi, z0, theta, r0)]


draw(zArray, rArray, phiArray, inArray)