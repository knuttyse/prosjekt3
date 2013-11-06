from scitools.std import *
import numpy as np
oppgave= raw_input("hvilken deloppgave?")

def extract2(filename):
    infile = open(filename,'r')
    i=0
    counter=0
    E_k=[]
    E_pot=[]
    L=[]

    for line in infile:
        words= line.split()
        E_k.append(float(words[0]))
        E_pot.append(float(words[1]))
        L.append(float(words[2]))
    return array(E_k), array(E_pot), array(L)

def extract(filename):
    infile = open(filename,'r')
    i=0
    counter=0

    for line in infile:

        words= line.split()

        if counter==0:
            k=int(words[0])
            n=int(words[1])
            N=int(words[2])
            T=float(words[3])
            matrix=zeros((n+1,k*N))
            t=zeros((n+1))
        if counter==1:
            planets= words[:]

        if counter>1:
            for j in range(k*N):
                t[i]=float(words[0])
                matrix[i, j]=float(words[j+1])
            i+=1
        counter+=1
    t=t[0:i]
    matrix=matrix[0:i,:]
    
    infile.close()
    return k,n,N,T,t,planets, matrix


readfile='data.txt'
readfile2='data2.txt'
k,n,N,T,t,planets,matrix = extract(readfile)
savefile='../../doc/plot%sn%dT%g.png' % (oppgave, n,T)

hold('on')
for i in range(k): 
    plot(matrix[:,i*N],matrix[:,i*N+1],'o',markersize=0.2,legend=planets[i],title=('n=%d steps, T=%7.2f days' % (n,T) ),xlabel='x [AU]',ylabel='y [AU]',hardcopy= savefile)
   
if (oppgave=="c") or (oppgave=="e"):
    theta=linspace(0,2*pi,200)
    x=cos(theta)
    y=sin(theta)
    plot(x,y,'k',hardcopy= savefile)

raw_input()

hold('off')
E_k,E_pot,L= extract2(readfile2)
plot(t,E_k,legend='E_k',xlabel='time [days]',ylabel='energy [Msun*AU^2/day^2]',title=('n=%d steps, T=%7.2f days' % (n,T) ),hardcopy='../../doc/ek%s%d%g.png'% (oppgave, n,T) )
raw_input()
plot(t,E_pot,legend='E_pot',xlabel='time [days]',ylabel='energy [Msun*AU^2/day^2]',title=('n=%d steps, T=%7.2f days' % (n,T) ),hardcopy='../../doc/ep%s%d%g.png' % (oppgave, n,T))
raw_input()
plot(t,L,legend='L',xlabel='time [days]',ylabel='angular momentum [Msun*AU^2/day]',title=('n=%d steps, T=%7.2f days' % (n,T) ), hardcopy='../../doc/L%s%d%g.png' % (oppgave, n,T) )
raw_input()
plot(t,E_pot+E_k,legend='E_tot', xlabel='time [days]', ylabel='energy [Msun*AU^2/day^2]',title=('n=%d steps, T=%7.2f days' % (n,T) ),hardcopy='../../doc/e%s%d%g.png' % (oppgave, n,T))
raw_input()
