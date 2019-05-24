
from ROOT import *
#import numpy as np
from array import array
import ROOT
import math

gROOT.SetBatch(True)


'''
A = np.array([[1,0],[0,2]])
w, v = np.linalg.eig(A)
print('w=',w)
print('v=',v)

'''


def get_psi0_ij(i,j=0, Nx=500, Ny=500):
    if i<0 or j<0 or i>=Nx or j>=Ny:
        return 0.
    k = 50
    dx = 0.01
    sigma = 0.1
    x0 = 0.2*Nx*dx
    y0 = 0.5*Ny*dx
    x = i*dx
    y = j*dx
    psi0 = math.exp(-(pow(x-x0,2)+pow((y-y0)/2.,2))/pow(sigma,2)/2.)*(math.cos(k*x)+1j*math.sin(k*x))
    #print('psi0=',psi0)
    return psi0

def get_psi0(Nx=500,Ny=500):
    psi = []
    N = 0.
    for i in range(Nx):
        psix = []
        for j in range(Ny):
            y = get_psi0_ij(i,j, Nx, Ny)
            N += pow(abs(y),2)
            psix.append(y)
        psi.append(psix)
    for i in range(Nx):
        for j in range(Ny):
            psi[i][j] *= 1./pow(N,0.5)
    return psi

def get_Vij(i,j, Nx, Ny):
    Vmax = 5100
    C = Nx/2.
    S = 5.
    dx = 0.01
    sigma = 0.1
    slit_width = 0.4*sigma
    slit_distance = 2*sigma
    # potential for tunneling effect
    #Vij = math.exp(-pow((i-C)/S,2)/2.)*Vmax
    # potential for double-slit experiment
    Vij = 0.
    if i == int(0.5*Nx) and (abs(j-0.5*Ny)*dx<slit_distance*0.5 or abs(j-0.5*Ny)*dx>slit_distance*0.5+slit_width):
       Vij = 3000.
    elif abs(i-int(0.5*Nx))==1 and (abs(j-0.5*Ny)*dx<slit_distance*0.4 or abs(j-0.5*Ny)*dx>slit_distance*0.6+slit_width):
       Vij = 3000
    else:
       Vij = 0.
    return Vij

def get_f_ij(i,j=0,PSI=[],Nx=500,Ny=500):

    psi_1_1 = PSI[i][j]

    psi_0_1 = 0.
    psi_0_1 = PSI[i][j]
    if i-1>=0:
        psi_0_1 = PSI[i-1][j]
    psi_2_1 = 0.
    psi_2_1 = PSI[i][j]
    if i+1<Nx:
        psi_2_1 = PSI[i+1][j]

    psi_1_0 = 0.
    psi_1_0 = PSI[i][j]
    if j-1>=0:
        psi_1_0 = PSI[i][j-1]

    psi_1_2 = 0.
    psi_1_2 = PSI[i][j]
    if j+1<Ny:
        psi_1_2 = PSI[i][j+1]

    dx = 0.01
    d2psi_dx2 = (psi_2_1 + psi_0_1 - 2*psi_1_1)/(dx*dx)
    d2psi_dy2 = (psi_1_2 + psi_1_0 - 2*psi_1_1)/(dx*dx)
    Vij = get_Vij(i,j, Nx, Ny)
    Hpsi = -1./2.*(d2psi_dx2 + d2psi_dy2) + Vij*psi_1_1
    f = -1j*Hpsi
    return f
def get_F(PSI, Nx=500, Ny=500):
    F = []
    for i in range(Nx):
        Fx = []
        for j in range(Ny):
            fij = get_f_ij(i,j,PSI,Nx,Ny)
            Fx.append(fij)
        F.append(Fx)
    return F
def update_PSI(PSI,K,c, Nx=500, Ny=500):
    PSI_new = PSI
    #print('PSI len =',len(PSI),len(PSI[0]))
    #print('K len =',len(K),len(K[0]))
    for i in range(Nx):
        for j in range(Ny):
            #continue
            PSI_new[i][j] = PSI[i][j] # +c* K[i][j]
    return PSI_new

def get_psit(PSI,Nx=500,Ny=500):
    dt = 0.000001
    K1 = get_F(PSI,Nx,Ny)
    K2 = get_F(update_PSI(PSI,K1,dt/2., Nx, Ny),Nx,Ny)
    K3 = get_F(update_PSI(PSI,K2,dt/2., Nx, Ny),Nx,Ny)
    K4 = get_F(update_PSI(PSI,K3,dt, Nx, Ny),Nx,Ny)
    PSInew = PSI
    N = 0.
    for i in range(Nx):
        for j in range(Ny):
            PSInew[i][j] = PSI[i][j]+dt/6.*(K1[i][j]+2.*K2[i][j]+2.*K3[i][j]+K4[i][j])
            N += pow(abs(PSInew[i][j]),2)
    for i in range(Nx):
        for j in range(Ny):
            PSInew[i][j] *= 1./pow(N, 0.5)
    return PSInew

def savePSI(PSInew,it):
    Nx = len(PSInew)
    Ny = len(PSInew[0])
    print('Nx=%s, Ny=%s' % (Nx,Ny))
    with open('psi_'+str(it)+'.txt', 'w') as f:
        for i in range(Nx):
            for j in range(Ny):
                f.write(str(PSInew[i][j])+',')
            f.write('\n')
    return
def readPSI(PSIfile):
    PSI = []
    with open(PSIfile, 'r') as f:
        lines = f.readlines()
    Nx = len(lines)
    for line in lines:
        line = line.replace('[', '')
        line = line.replace(']', '')
        line = line.replace('(', '')
        line = line.replace(')', '')
        line = line.split(',')
        psi = []
        for a in line:
            if 'j' not in a:
                continue
            a = a.strip()
            print a, complex(a)
            psi.append(complex(a))
        PSI.append(psi)
    it = PSIfile.replace('psi_','')
    it = it.replace('.txt', '')
    it = int(it)
    return it, PSI

def get_psit_alltimes(Nt=1,Nx=500,Ny=500, PSIfile=''):
    PSIold = get_psi0(Nx,Ny)
    print('Nt =',Nt)
    it0 = 0
    if PSIfile != '':
        it0,PSIold = readPSI(PSIfile)
        it0 += 1
    PSInew = PSIold
    print('it0 =', it0)
    for it in range(it0,it0+Nt):
        if it == 0:
            show_eig(PSInew,Nx,Ny,plotname='t0')
            continue
        PSInew = get_psit(PSIold,Nx,Ny)
        if it%100 == 0:
            show_eig(PSInew,Nx,Ny,plotname='t'+str(it))
        if it == it0+Nt-1 or it%1000 == 0:
            savePSI(PSInew,it)
        PSIold = PSInew
    return


def show_eig(PSI,Nx=500, Ny=500,plotname=''):
    h2 = TH2F('h2', '', Nx, 0, Nx, Ny, 0, Ny)
    Npsisq = 0.
    for i in range(Nx):
        for j in range(Ny):
            psisq = pow(abs(PSI[i][j]),2)
            Npsisq += psisq
            h2.SetBinContent(i+1,j+1, psisq)
    print(plotname,Npsisq)
    h2.Scale(1./Npsisq)
    Cs = TCanvas('Cs', '', 800, 600)
    h2.Draw('colz')
    h2.GetXaxis().SetTitle('x')
    h2.GetYaxis().SetTitle('y')
    h2.GetZaxis().SetTitle('|#psi(x)|^{2}')
    if plotname == '':
        plotname = 'Cs_psix'
    else:
        plotname = 'Cs_psix_'+plotname
    Cs.SaveAs(plotname+'.png')
    #Cs.SaveAs(plotname+'.root')
    return


def main(Nt=1,Nx=100,Ny=100, PSIfile=''):
    get_psit_alltimes(Nt,Nx,Ny, PSIfile)
    return

#readPSI('psi_2999.txt')
main(5000, 100, 100, PSIfile='psi_9999.txt')
