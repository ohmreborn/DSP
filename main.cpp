#include <iostream>
#include <fstream>
#include <cblas.h>

void count_star(int *n_star,std::string filename){
    *n_star = 0;

    std::fstream newfile;
    newfile.open(filename,std::ios::in);
    std::string line;

    if(newfile.is_open()){
    while(getline(newfile, line)){
        if (line != "mass,x,y,vx,vy,name"){  
            (*n_star)++;
            }
        }
        newfile.close();
    }else{
        std::cout << "error" << std::endl;;
        }

}

void read(std::string filename,int dim,double **v,double **r,double *m,std::string *star_name){
    const char delimeter = ',';
    std::fstream newfile(filename);
    std::string line;
    int length;
    int line_ind = 0;
    int ind;
    std::string element;
    
    if(newfile.is_open()){
        while(getline(newfile, line)){
            if (line != "mass,x,y,vx,vy,name"){  
                length = line.length();
                for (int i=0;i<length;i++){
                    if(line[i] != delimeter){
                        element +=  line[i];
                    }
                    if((line[i] == delimeter) || (i==length-1)){
                        ind = i + 1;
                        m[line_ind] = stod(element);
                        element = "";
                        break;
                    }
                }
                for (int j=0;j<dim;j++){
                    for (int i=ind;i<length;i++){
                        if(line[i] != delimeter){
                            element +=  line[i];
                        }
                        if((line[i] == delimeter) || (i==length-1)){
                            ind = i + 1;
                            r[line_ind][j] = stod(element);
                            element = "";
                            break;
                        }
                    }
                }
                for (int j=0;j<dim;j++){
                    for (int i=ind;i<length;i++){
                        if(line[i] != delimeter){
                            element +=  line[i];
                        }
                        if((line[i] == delimeter) || (i==length-1)){
                            ind = i + 1;
                            v[line_ind][j] = stod(element);
                            element = "";
                            break;
                        }
                    }
                }
                for (int i=ind;i<length;i++){
                    if(line[i] != delimeter){
                        element +=  line[i];
                    }
                    if((line[i] == delimeter) || (i==length-1)){
                        star_name[line_ind] = element;
                        element = "";
                        break;
                    }
                }
                line_ind ++;
            }
        }
        newfile.close();
        }else{
        std::cout << "error";
    }

}

void slove(std::string output,int dim,int n_star,double **v,double **r,double *m,std::string* star_name){
    const double G = 1;
    double time = 1;
    double dt = 2e-5;
    int step = time/dt;

    double distance;
    double *norm_distance = new double[dim];
    double *dr = new double[dim];

    double **dv = new double*[n_star];

    for (int i=0;i<n_star;i++){
        dv[i] = new double[dim];
    }

    std::ofstream myfile;
    myfile.open (output);
    myfile << "x,y,star,frame" << std::endl;

    for (int i=0;i<step;i++){
        for (int ai=0;ai<n_star;ai++){
            for (int aj=0;aj<dim;aj++){
                dv[ai][aj] = 0;
            }
        for (int d=0;d<dim;d++){
            myfile << r[ai][0] << "," << r[ai][1] << ",";
        }
        myfile << star_name[ai] << "," << i << std::endl;
        }
        for (int star1=0;star1<n_star;star1++){
            for (int star2=0;star2<n_star;star2++){
                if (star1 != star2){
                    cblas_dcopy(dim,r[star2],1,norm_distance,1);
                    cblas_daxpy(dim,-1,r[star1],1,norm_distance,1);
                    distance = cblas_dnrm2(dim,norm_distance,1);
                    cblas_dscal(dim,(G * m[star2]/(pow(distance,3))),norm_distance,1);
                    cblas_daxpy(dim,1,norm_distance,1,dv[star1],1);
                }
            }
        }
        for (int star=0;star<n_star;star++){
            cblas_dscal(dim,dt,dv[star],1);
            cblas_daxpy(dim,1,dv[star],1,v[star],1);
            cblas_dcopy(dim,v[star],1,dr,1);
            cblas_dscal(dim,dt,dr,1);
            cblas_daxpy(dim,1,dr,1,r[star],1);
            // std::cout << r[star][0] << ", " << r[star][1];
            // std::cout << std::endl;
        }

    }
    myfile.close();
    for (int i=0;i<n_star;i++){
        delete[] dv[i];
    }
    delete[] norm_distance;
    delete[] dv;
    delete[] dr;
}

int main(){
    int n_star;
    const int dim = 2;
    const char delimeter = ',';
    const std::string filename = "../input.csv";
    const std::string output = "../output.csv";

    count_star(&n_star,filename);
    
    std::cout << n_star << std::endl;
    std::string *star_name = new std::string[n_star];
    double **v = new double*[n_star];
    double **r = new double*[n_star];
    double *m = new double[n_star];
    for (int i=0;i<n_star;i++){
        v[i] = new double[dim];
        r[i] = new double[dim];
    }
    std::cout << "read" << std::endl;
    read(filename,dim,v,r,m,star_name);
    std::cout << "read sucess" << std::endl;

    slove(output,dim,n_star,v,r,m, star_name);

    for (int i=0;i<dim;i++){
        delete[] v[i];
        delete[] r[i];
    }
    delete[] v;
    delete[] r;
    delete[] m;

    return 0;
}