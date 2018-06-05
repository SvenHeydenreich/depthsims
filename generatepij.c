#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main(){
  static int num = 400;
  float ds[num],pis[num],xs[num],ys[num];
  int i,j,k,l;
  for(i=0;i<num;i++){
    ds[i]=i*1.5/num;
    pis[i]=M_PI*i*2./num;
    xs[i]=i*1./num;
    ys[i]=i*1./num;
  }
  float result,frac,x1,x2;
  float finalplot[num];
  for(i=0;i<num;i++){
//    printf("%f\n",ds[i]);
    result = 0;
    for(j=0;j<num;j++){
      for(k=0;k<num;k++){
        frac = 0;
        for(l=0;l<num;l++){
          x1 = xs[j]+ds[i]*cos(pis[l]);
          x2 = ys[k]+ds[i]*sin(pis[l]);
          if((x1>=0)&&(x1<=1)&&(x2>=0)&&(x2<=1)){
            frac = frac+1;
          }
        }
        frac = frac/(float)num;
/*        if(frac!=1){
          printf("%f\n",frac);
        }
*/
        result = result+frac;
      }
    }
    result = result/(float)(num*num);
    finalplot[i]=result;
    printf("%f %f\n",ds[i],result);
  }
FILE *fptr;
fptr = fopen("pij.dat","w");
if(fptr == NULL){
  printf("Error!");
  exit(1);
}
for(i=0;i<num;i++){
  fprintf(fptr,"%f %f\n",ds[i],finalplot[i]);
}
fclose(fptr);
}
/*    cdef np.ndarray[np.float_t,ndim=1] ds = np.linspace(0,2,200)
    cdef np.ndarray[np.float_t,ndim=1] pis = np.linspace(0,2,200)*np.pi
    cdef np.ndarray[np.float_t,ndim=1] xs = np.linspace(0,1,200)
    cdef np.ndarray[np.float_t,ndim=1] ys = np.linspace(0,1,200)

    finalplot = np.zeros(np.size(ds))
    cdef int i = 0
    cdef float result
    cdef float frac
    cdef float x1,x2
    for d in ds:
        if(np.abs(i/5.-int(i/5))<0.0001):
            print(i*1./np.size(ds)*100.,'%')
        result = 0
        for x in xs:
            for y in ys:
                frac = 0
                for pi in pis:
                    x1 = x+d*np.cos(pi)
                    x2 = y+d*np.sin(pi)
                    if((x1 >= 0) and (x1 <= 1) and (x2 >= 0) and (x2 <= 1)):
                        frac = frac+1
                frac = frac/np.size(pis)
                result = result+frac
        result = result/(np.size(xs)*np.size(ys))
        finalplot[i] = result
        i = i+1
    eq = finalplot*0.1**2+1
    plt.figure(1)
    plt.plot(ds,finalplot)
    plt.figure(2)
    plt.plot(ds,eq)
    ps = np.fft.rfft(eq)
    plt.figure(3)
    plt.semilogy(np.abs(ps)**2)
    plt.show()
*/
