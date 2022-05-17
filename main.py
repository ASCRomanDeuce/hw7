import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import integrate
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
#from celluloid import Camera


class diff:
    def __init__(self, doublea,J, N, deltax, deltat):
        self.doublea,self.N, self.J, self.deltat, self.deltax=doublea,N, J, deltat, deltax
        self.u=np.zeros([N,J], dtype = complex)
        self.u[:,0]=self.InitialCondition(0)
        self.linearAddedMember = np.zeros([N,J])
        for i in range(int(self.N*0.6),int(self.N*0.75)):
            self.linearAddedMember[i,:] = -0.05
        
        
        #for i in range(1,N-1):
        #    self.u[i,0]=np.sin(np.pi*i/(N-1))
        #    print(i, self.u[i,0])

        self.u[0,:]=0
        self.u[self.N-1,:]=0
        print(self.u)

    def InitialCondition(self, n):
        tempArray = np.zeros(self.N, dtype = complex)
        shift = (self.N-1)/2
        for i in range(self.N):
            tempArray[i] = (1.0/self.deltax*(2*np.pi)**(0.5))*np.exp(1j * n * (i-shift) / (self.N - 1) ) * np.exp(-(i-shift)**2/4)#сделал начало в x = 5 изза этого норма не ПИ
        return tempArray
    
    def solver2(self,m):
        ePot = self.linearAddedMember
        z=np.copy(self.u)
        q=self.doublea *self.deltat/(self.deltax)**2
        #q=0.25j
        print(q)
        
        a=np.repeat((-1+0j),self.N-1)
        c=np.copy(a)
        b=np.repeat((2+2/q + 0j),self.N-1)
        d=np.zeros(self.N-1, dtype = complex)
        y=np.copy(d)
        alpha=np.copy(d)
        beta=np.copy(d)

        
        for k in range(m):
            for n in range(1, self.J):
                u=self.u
                z=self.u
                d=np.zeros(self.N-1,dtype = complex)
                y=np.copy(d)
                alpha=np.copy(d)
                beta=np.copy(d)
                for i in range(1, self.N-3, 1):
                    d[i]=u[i,n-1]+(-2+2/q-2*ePot[i+1,n-1]*self.deltat/self.doublea)*u[i+1,n-1]+u[i+2,n-1]
                d[0]=u[0,n]+(-2+2/q-2*ePot[1,n-1]*self.deltat/self.doublea)*u[1,n-1]+u[2,n-1]+u[0,n-1]
                d[self.N-3]=u[self.N-3,n-1]+(-2+2/q-2*ePot[self.N-2,n-1]*self.deltat/self.doublea)*u[self.N-2,n-1]+u[self.N-1,n-1]+u[self.N-1,n]
                y[0]=y[0]+b[0]
                alpha[0]=-c[0]/y[0]
                beta[0]=d[0]/y[0]
                
                for i in range(1,self.N-3):
                    y[i]=b[i]+a[i]*alpha[i-1]
                    alpha[i]=-c[i]/y[i]
                    beta[i]=(d[i]-a[i]*beta[i-1])/y[i]
                y[self.N-3]=b[self.N-3]+a[self.N-3]*alpha[self.N-4]
                beta[self.N-3]=( d[self.N-3]-a[self.N-3]*beta[self.N-1-3])/y[self.N-3]

                u[self.N-2,n]=beta[self.N-3]
                for i in range(self.N-3, 0, -1):
                    u[i,n]=alpha[i-1]*z[i+1,n]+beta[i-1]
                
  

    def plotter(self):
        X = np.arange(-self.J*self.deltax/2,self.J*self.deltax/2,self.deltax)
        T = np.arange(-self.N*self.deltat/2,self.N*self.deltat/2,self.deltat)
        X, T = np.meshgrid(X, T)
        

        fig=plt.figure()
        ax=fig.add_subplot(projection='3d')
        surf = ax.plot_surface(X, T, self.u.real**2 + self.u.imag**2, cmap=cm.coolwarm, alpha = 0.5,
                       linewidth=0, antialiased=False)
       
        surf = ax.plot_surface(X, T, self.linearAddedMember, cmap=cm.coolwarm, alpha = 0.5,
                       linewidth=0, antialiased=False)
        plt.xlabel('t')
        plt.ylabel('x')
        ax.set_title('T(x,t)')
        fig.colorbar(surf, ax=ax, fraction=0.02, pad=0.1, label='Temperature')
        plt.show()

        #fig, ax = plt.subplots()
        #CS = ax.contour(X, T, self.u)
        
        #ax.clabel(CS, inline=True, fontsize=10)
        #ax.set_title('equipotential lines  of temperature')
        #fig.show()
        
    def integrate(self):
        x = np.linspace(-self.N*self.deltax/2, self.N*self.deltax/2, num=self.N)
        self.norm = np.trapz(self.u.real**2 + self.u.imag**2, x=None, dx=self.deltax, axis=0)
        print(self.norm)

   
        
    def update_plot(self,frame_number, zarray, plot):
        ax=self.ax
        x,zarray=self.X,self.zarray
        self.line1.set_ydata(zarray[:,frame_number])
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        plot=ax.plot(x, self.linearAddedMember[:,0], color = 'gray')
        #plot = ax.plot(x, zarray[:,frame_number])

    def animation(self):
    
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot()
        ax=self.ax

        
        X = np.arange(-self.N*self.deltax/2,self.N*self.deltax/2,self.deltax)
        zarray = np.zeros((self.N,self.J))
        self.X,self.zarray = X,zarray


        for i in range(self.J):
            zarray[:,i] = (self.u.real[:,i]**2 + self.u.imag[:,i]**2)/self.norm[i]

        plot = ax.plot(X, zarray[:,0])
        plot=ax.plot(X, self.linearAddedMember[:,0])
        self.line1, =ax.plot(X, zarray[:,0], 'r-')
        animate = animation.FuncAnimation(self.fig, self.update_plot, self.J, fargs=(zarray, plot))
        plt.show()
        
        animate.save('example.gif', writer='imagemagick', fps=30)








x=diff(0.7j, 1000, 100, 10, 10)

x.solver2(1)
#x.plotter()
x.integrate()
x.animation()


