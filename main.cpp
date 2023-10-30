#include "star.h"
#include <fstream>
using namespace std;

double abs(double* a);
void minu(double* a ,double* b,double res[2]);
void plu(double* a ,double* b,double res[2]);
void mul(double a,double b[2],double res[2]);
void slove(Star *system,int n_star,double time,double dt,string filename);



int main(){
    double m1 = 1;
    double m2 = 333000;
    double m3 = 1000;
    double r1[2] = {1,0};
    double v1[2] = {0,sqrt(m2)};
    double r2[2] = {0,0};
    double v2[2] = {0,0};
    double r3[2] = {-0.5,0};
    double v3[2] = {0,-100};
    double time = 1;
    double dt = 5e-5;//time/200;

    Star earth = Star(m1,r1,v1,"earth");
    Star sun = Star(m2,r2,v2,"sun");
    Star earth2 = Star(m3,r3,v3,"earth2");
    Star system[] = {earth,sun,earth2};
    const int n_system = sizeof(system)/sizeof(system[0]);
    slove(system,n_system,time,dt,"results.csv");

    return 0;
}

double abs(double* a){
    double total = 0;
    for (int i=0;i<2;i++){
        total = total + a[i]*a[i];
    }
    return sqrt(total);
}

void minu(double* a ,double* b,double res[2]){
    for (int i=0;i<2;i++){
        res[i] = a[i] - b[i];
    }
}

void plu(double* a ,double* b,double res[2]){
    for (int i=0;i<2;i++){
        res[i] = a[i] + b[i];
    }
}

void mul(double a,double b[2],double res[2]){
    for (int i=0;i<2;i++){
        res[i] = b[i]*a;
    }
}

void slove(Star *system,int n_star,double time,double dt,string filename){
    ofstream myfile;
    myfile.open (filename);
    myfile << "x,y,star,frame" << endl;

    double k = 1;
    double distance;
    double a[n_star][2];
    double acclelation[2];
    double r[n_star][2];
    double dv[2];
    double dr[2];
    double t = 0;
    int frame = 0;



    while (time > t){
        for (int i=0;i<n_star;i++){
            myfile << system[i].r[0] << "," << system[i].r[1] << "," << system[i].name << "," << frame << endl;
            for (int j=0;j<2;j++){
                a[i][j] = 0;
            }
        }

        // conmpute a acclelation;
        for (int i=0;i<n_star;i++){
            for (int j=0;j<n_star;j++){
                if (i != j){
                    minu(system[j].r,system[i].r,r[i]);
                    distance = abs(r[i]);

                    mul(k*system[j].m/pow(distance,3),r[i],acclelation);
                    
                    for (int l=0;l<2;l++){
                        a[i][l] = a[i][l] + acclelation[l];
                    }
                }
            }
        }
    
        for (int i=0;i<n_star;i++){
            // update a velocity

            mul(dt,a[i],dv);
            plu(system[i].v,dv,system[i].v);
            
            // update position
            mul(dt,system[i].v,dr);
            plu(system[i].r,dr,system[i].r);
        }

        frame = frame + 1;
        t = t + dt;
    }

}

// g++ main.cpp -o main && ./main && python3 animation.py
