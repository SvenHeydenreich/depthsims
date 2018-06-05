#include <math.h>
#include <stdlib.h>
#include <stdio.h>

float length(float v[2]);
void analyze(int xposs[],int yposs[],float shearsr[],float shearsi[], float dmax, int dsteps, int galaxynum);
float dot_product(float v[2],float w[2]);
float determinant(float v[2],float w[2]);
float inner_angle(float v[2],float w[2]);
float angle_clockwise(float A[2], float B[2]);
float min(float x,float y);
float max(float x,float y);

int main(){
  int galaxynum = 20000;
  int xpos[galaxynum],ypos[galaxynum];
  float shearr[galaxynum],sheari[galaxynum];
  int i;
  FILE *ptr;
  ptr = fopen("galaxies.txt","r");
  //fscanf(ptr,"%f" "%f" "%f" "%f",&xpos[0] &ypos[0] &shearr[0] &sheari[0]);
  for(i=0;i<galaxynum;i++){
    fscanf(ptr,"%d",&xpos[i]);
    fscanf(ptr,"%d",&ypos[i]);
    fscanf(ptr,"%f",&shearr[i]);
    fscanf(ptr,"%f",&sheari[i]);
    fscanf(ptr,"\n");
  }
  fclose(ptr);
/*
  printf("%d\n",xpos[0]);
  printf("%d\n",ypos[0]);
  printf("%f\n",shearr[0]);
  printf("%f\n",sheari[0]);
  printf("%d\n",xpos[1]);
  printf("%d\n",ypos[1]);
  printf("%f\n",shearr[1]);
*/
  analyze(xpos,ypos,shearr,sheari,1500,750,galaxynum);
}
float min(float x,float y){
  if(x<y){
    return x;
  }
  else{
    return y;
  }
}

float max(float x,float y){
  if(y<x){
    return x;
  }
  else{
    return y;
  }
}

float length(float v[2]){
    return sqrt(pow(v[0],2)+pow(v[1],2));
  }
float dot_product(float v[2],float w[2]){
     return v[0]*w[0]+v[1]*w[1];
   }
float determinant(float v[2],float w[2]){
     return v[0]*w[1]-v[1]*w[0];
   }
float inner_angle(float v[2],float w[2]){
     float cosx = dot_product(v,w) / (length(v)*length(w));
     float rad=acos(min(cosx,1));
     return rad;
  }
float angle_clockwise(float A[2], float B[2]){
      float inner = inner_angle(A,B);
      float det = determinant(A,B);
      if(det<0){
          return inner;
        }
      else{
          return 2*3.14159265359-inner;
        }
      }
void analyze(int xposs[],int yposs[],float shearsr[],float shearsi[], float dmax, int dsteps, int num){
      float xip[dsteps],xim[dsteps],xix[dsteps];
      int mask[dsteps];

      int ii,jj,jj2,jj3;
      int xpos1,xpos2,ypos1,ypos2,point;
      float shear1r,shear1i,shear2r,shear2i,shearx1,shearx2,sheart1,sheart2,point2;
      float phi,dist;
      float x1[2],x2[2];
      for(ii=0;ii<dsteps;ii++){
        xip[ii]=0;
        xim[ii]=0;
        xix[ii]=0;
        mask[ii]=0;
      }
      for(ii=0;ii<num;ii++){
        jj3 = num-ii;
        for(jj2=0;jj2<jj3;jj2++){
              jj = jj2+ii;
              if(ii != jj){
                  xpos1 = xposs[ii];
                  ypos1 = yposs[ii];
                  xpos2 = xposs[jj];
                  ypos2 = yposs[jj];
                  shear1r = shearsr[ii];
                  shear1i = shearsi[ii];
                  shear2r = shearsr[jj];
                  shear2i = shearsi[jj];
                  x1[0] = xpos1;
                  x1[1] = ypos1;
                  x2[0] = xpos2;
                  x2[1] = ypos2;
                  phi = angle_clockwise(x1,x2);
                  sheart1 = -shear1r*cos(2*phi)-shear1i*sin(2*phi);
                  shearx1 = shear1r*sin(2*phi)-shear1i*cos(2*phi);
                  sheart2 = -shear2r*cos(2*phi)-shear2i*sin(2*phi);
                  shearx2 = shear2r*sin(2*phi)-shear2i*cos(2*phi);
                  dist = sqrt(pow((xpos1-xpos2),2)+pow((ypos1-ypos2),2));
        //          dist = max(dist,0.001);
                  if(dist<dmax){
                      point2 = dist/dmax;
                      point2 = point2 * (float) dsteps;
                      point = (int) point2;
                      xip[point] = xip[point]+ sheart1*sheart2+shearx1*shearx2;
                      xim[point] = xim[point]+ sheart1*sheart2-shearx1*shearx2;
                      xix[point] = xix[point]+ sheart1*shearx2;
                      mask[point] = mask[point]+1;
                    }
                  }
                }
              }
          for(ii=0;ii<dsteps;ii++){
            xip[ii] = xip[ii]/mask[ii];
            xim[ii] = xim[ii]/mask[ii];
            xix[ii] = xix[ii]/mask[ii];
            if(mask[ii]<4){
              printf("***************************************");
              printf("%d\n",ii);
            }
          }
          FILE *fptr;
          fptr = fopen("functions.txt","w");
          if(fptr == NULL){
            printf("Error!");
            exit(1);
          }
          for(ii=0;ii<dsteps;ii++){
            fprintf(fptr,"%f %f %f %d\n",xip[ii],xim[ii],xix[ii],mask[ii]);
          }
          fclose(fptr);

            }
  /*    xip = xip/mask;
      xim = xim/mask;
      xix = xix/mask;
      for ii in range(dsteps):
          if mask[ii]==0:
              xip[ii]=np.mean(xip)
              xim[ii]=np.mean(xim)
              xix[ii]=np.mean(xix)

      return xip,xim,xix
    }
*/
