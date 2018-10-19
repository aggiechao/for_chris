#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdio>
#include<cstdlib>
#include<unistd.h>
#include<ctime>
#include<cmath>
#include<math.h>
#include"header.h"
#include"dSFMT.h"
#include"dSFMT.c"
double energy_func(unsigned char *config, double *bonds,double *tail);
void metro(int *neig_index,int *neig, unsigned char **config, double *bonds,double* tail,double* t,double* energy,dsfmt_t* r);
void pt(unsigned char** config,double* energy,double* t,dsfmt_t* r);
void link(double* link_overlap,unsigned char*** config,int* neig_h, int* neig_index_h,int& neig_numb);
void set_ktables(double* kcostable,double*,int l);
double chi_func(unsigned char*** config, double* kcostable,double* ksintable,int& t,int index);
void cluster(int* neig,int* neig_index,const int& index,const int& temp,int& count, char* check_spin,const double& limit,unsigned char*** config);
void  cluster_move(int cluster_or_not,double cluster_ratio,unsigned char*** config,double* bonds,double* tail,double** energy,int* neig_index,int* neig,double* temperature,dsfmt_t* rn); 
long seedgen(void){
        long s, seed, pid;
        time_t seconds;
        pid = getpid();
        s = time ( &seconds ); /* get CPU seconds since 01/01/1970 */
        seed = std::abs(((s*181)*((pid-83)*359))%104729);
        return seed;
}

int main(int argc,char *argv[]){

/*------------read file--------------*/
        double* bonds = new double[spin_numb*spin_numb];
   	for (int i = 0;i <= spin_numb*spin_numb-1;i++){
		bonds[i] = 0.0;
	}
	FILE* fp = fopen(argv[1],"r");
        register int i = 0;
        register int j = 0;
        while (feof(fp)==0){
                fscanf(fp,"%d",&i);
                fscanf(fp,"%d",&j);
                fscanf(fp,"%lf",bonds+i*spin_numb+j);
                *(bonds+j*spin_numb+i) = *(bonds+i*spin_numb+j);
                /*for katz*/
     //           std::fscanf(fp,"%lf",bonds+(i-1)*spin_numb+(j-1));
     //           *(bonds+(j-1)*spin_numb+(i-1)) = *(bonds+(i-1)*spin_numb+(j-1));
        }
        fclose(fp);
        

        /*-------find neighbors------*/
        int* neig = new int[spin_numb];
        int* neig_h = new int[spin_numb];
	for (int i = 0; i <= spin_numb-1;i++){
		neig[i] = 0;
		neig_h[i] = 0;
	}
        int* neig_index = new int[spin_numb*spin_numb];
	int* neig_index_h = new int[spin_numb*spin_numb];
	for (int i = 0; i <= spin_numb*spin_numb-1;i++){
        	neig_index[i] = 0;
		neig_index_h[i] = 0;
        }
        for (int i = 0;i <= spin_numb-1; i++){
                for (int j = 0; j <= spin_numb-1; j++){
                        if (*(bonds+i*spin_numb+j) != 0 && i != j){
				neig_index[i*spin_numb+neig[i]] = j;
                                ++neig[i];                        /*exact the number of neighbors*/
                                if (j > i){
                                        neig_index_h[i*spin_numb+neig_h[i]] = j;
                                	++neig_h[i];
				}
                        }
                }
        }
     
	int neig_numb = 0;
        for (int i = 0;i <= spin_numb-1;i++){
                neig_numb += neig_h[i];
        } 
	printf("neighbor: %d\n",neig_numb);	
/*-------------------initialize energy-------------------*/
        unsigned char*** config = new unsigned char**[replica];
        double* tail = new double[spin_numb];
        for (int rep = 0;rep <= replica-1;rep++){
                config[rep] = new unsigned char*[temp_set];
                for (int t = 0; t <= temp_set-1;++t){
                        *(config[rep]+t) = new unsigned char[spin_numb];
                }       
        }

        /*initialize config randomly*/
        unsigned char base = 1;
        srand(time(0));
        double r =0.0;
        for (int i = 0;i <= replica-1;i++){
                for (int t = 0;t <= temp_set-1;t++){
                        for (int j = 0; j <= spin_numb-1;j++){
                                r = (double)rand()/RAND_MAX;
                                *(*(config[i]+t)+j) = (r <= 0.5)*base;
			}
                }
        }
	
        double**  energy = new double*[replica];
        double** meg = new double*[replica];
	for (int rep = 0;rep <= replica-1;rep++){
                energy[rep] = new double[temp_set];
		meg[rep] = new double[temp_set];
                for (int t = 0; t <= temp_set-1;t++){
                        *(energy[rep]+t) = energy_func(*(config[rep]+t),bonds,tail);
                	*(meg[rep]+t) = 0.0;
		}
        }

	

        /*tail*/
        for (int i = 0; i <= spin_numb-1;i++){
                for (int j = 0; j <= spin_numb-1;j++){
                        tail[i] += (*(bonds+i*spin_numb+j));
                }   
        }        

/*--------------------monte carlo------------------*/
        double* temperature = new double[temp_set]; 
        FILE *fp_t  = fopen("temps.in","r");
	for (int t = 0; t <= temp_set-1; t++){
                fscanf(fp_t,"%lf\n",temperature+t);
//		temperature[t] = 3.0+(7.0/(temp_set-1))*t;
//		printf("%lf\n",temperature[t]);
	}
	fp = fopen(argv[2],"w");
        fprintf(fp,"%s %s  \n","MC","energy","w");
        
	double energy_plot = 0.0;
	
	long seed = 0;
        seed = seedgen();
        dsfmt_t rn;
        dsfmt_init_gen_rand(&rn,seed);
	
	/*set up for sn and cos tables*/
	double* kcostable = new double[spin_numb];
	double* ksintable = new double[spin_numb];
	set_ktables(kcostable, ksintable,4);
	

	/*MC*/
	for (int exp_power = 1;exp_power <= sweep;++exp_power){
                /*---------------measurement-------------*/
                energy_plot = 0.0;
                for (int monte = pow(2,exp_power-1);monte <= pow(2,exp_power);monte++){
               // 	std::cout << monte << std::endl;
                        for (int rep = 0;rep <= replica-1;rep++){
                                /*------------------metro-----------------*/
                                metro(neig_index,neig,config[rep],bonds,tail,temperature,energy[rep],&rn);
                                /*----------------------pt-----------------*/
                                pt(config[rep], energy[rep],temperature,&rn);
                        }
				
	//		cluster_move(1,1.0,config,bonds,tail,energy,neig_index,neig,temperature,&rn);
			energy_plot += (*(energy[0]+ti));
	
		}
                fprintf(fp,"%d %lf \n",exp_power,energy_plot/((double)(pow(2,exp_power-1))+1));
                /*-------------------measurement----------------*/

	}

	/*-------------------------------------------calculaet meg---------------------------------------*/
	double* meg2 = new double[temp_set];
        double* meg4 = new double[temp_set];
	for (int t = 0;t <= temp_set-1;t++){
                meg2[t] = 0.0;
                meg4[t] = 0.0;
        }
	for (int monte = pow(2,sweep-1);monte <=  pow(2,sweep);monte++){
        	for (int t = 0;t<=temp_set-1;t++){
        		double temp_m = 0.0;
               	 	for (int i =0;i<=spin_numb-1;i++){
                		temp_m += ((*(*(config[0]+t)+i))*2-1);
                	}
        		meg2[t] += pow(temp_m/spin_numb,2);
			meg4[t] += pow(temp_m/spin_numb,4);
		}
		/*MC*/
		for (int rep = 0;rep <= replica-1;rep++){
                	metro(neig_index,neig,config[rep],bonds,tail,temperature,energy[rep],&rn);
			/*----------------------pt-----------------*/
                        pt(config[rep], energy[rep],temperature,&rn);
               	}

	}
	/*average*/
	for (int t = 0;t<=temp_set-1;t++){
		meg2[t] = meg2[t]/(pow(2,sweep-1)+1);
		meg4[t] = meg4[t]/(pow(2,sweep-1)+1);
	}
	/*writing files*/
	std::fstream fs_m2;
        std::fstream fs_m4;
	fs_m2.open("m2.dat",std::ios_base::out|std::ios_base::app);
	fs_m4.open("m4.dat",std::ios_base::out|std::ios_base::app);
//	fs_m2 << temperature[0];
//	fs_m4 << temperature[0];
//	for (int t = 0;t <= temp_set-1;t++){
//		fs_m2 << ',' << temperature[t];
//		fs_m4 << ',' << temperature[t];
//	}
//	fs_m2 << std::endl;
 //       fs_m4 << std::endl;
	
	fs_m2 << meg2[0];
        fs_m4 << meg4[0];
	for (int t=1;t<=temp_set-1;t++ ){
		fs_m2 << ','<< meg2[t];
		fs_m4 << ',' << meg4[t]; 		
	}
	fs_m2 << std::endl;
	fs_m4 << std::endl;	
	/*-------------------------------------------calculaet meg---------------------------------------*/


        return 0;
}

double energy_func(unsigned char *config,double *bonds,double *tail){
	double energy = 0.0;    /*-(4ai*aj-2*ai-2aj+1)*Jij+(2*ai-1)*Jii*/
        for (int i = 0; i <= spin_numb-1;i++){
        	for (int j = i+1; j <= spin_numb-1;j++){
                	energy -= (*(bonds+i*spin_numb+j))*(4*config[i]*config[j]-2*config[i]-2*config[j]+1);
                }
                energy+=(2*config[i]-1)*(*(bonds+i*spin_numb+i));
        }
        return energy;
}


void metro(int *neig_index,int *neig, unsigned char **config, double *bonds,double* tail,double *t,double* energy,dsfmt_t* r){
 /*-----------------energy change-----------------*/
	double delta_energy = 0.0;
        unsigned char base = 1;
        int index = 0;
        double randn = 0.0;
        for (int j = 0;j <= temp_set-1;j++){
        	for (int s = 0; s <= spin_numb-1;s++){
                	index = int(dsfmt_genrand_close_open(r)*(spin_numb-1));
                        /*delta energy*/
                        delta_energy= 0.0;
                        for (int i =0; i <= neig[index]-1;++i){
                        	delta_energy += (*(config[j]+*(neig_index+index*spin_numb+i)))*bonds[index*spin_numb+(*(neig_index+index*spin_numb+i))];
                        }
                        delta_energy = (2*(*(config[j]+index))-1)*(4*delta_energy-2*tail[index]);
                        /*--------------------flip-----------------------*/
                        if (delta_energy <= 0.0){
                        	*(config[j]+index) = (*(config[j]+index)) ^ base;
                                energy[j] += delta_energy;
                        }
                        else {
                        	randn = dsfmt_genrand_close_open(r);
                                if (randn<((exp(-delta_energy/t[j])))){
                                	*(config[j]+index) = (*(config[j]+index)) ^ base;
                                        energy[j] += delta_energy;
                                }
                        }

                }
        }
}

void pt(unsigned char** config,double* energy,double* t,dsfmt_t* r){
        double randn = 0.0;
	double beta_diff  = 0.0;
        double delta_energy = 0.0;
        double energy_swap = 0.0;
        unsigned char* config_swap;
	for (int i = 0;i <= temp_set-2;i++){
                beta_diff = 1/t[i]-1/t[i+1];
		delta_energy = energy[i]-energy[i+1];
		if (beta_diff*delta_energy>=0){
                	config_swap = config[i];
                        config[i] = config[i+1];
                        config[i+1] = config_swap;
                        energy_swap = energy[i];
                        energy[i] = energy[i+1];
                        energy[i+1] = energy_swap;
                 }
                 else {
                 	randn = dsfmt_genrand_close_open(r);
                        if ((exp((beta_diff*delta_energy)))>randn){
                        	config_swap = config[i];
                                config[i] = config[i+1];
                                config[i+1] = config_swap;
                                energy_swap = energy[i];
                                energy[i] = energy[i+1];
                                energy[i+1] = energy_swap;
                        }
                }
        }
 }


void cluster(int* neig,int* neig_index,const int& index, const int& temp,int& count, char* check_spin,const double& limit,unsigned char*** config){
	unsigned char overlap = ((*(*(*(config+0)+temp)+index))^(*(*(*(config+1)+temp)+index)));
	switch(overlap){
		case true:
			if ((check_spin[index]!=(temp+1)) && (count <= int(limit))){
				++count;
				check_spin[index] = (temp+1);
				for(int i=0;i <= neig[index]-1;i++){
					cluster(neig,neig_index,neig_index[index*spin_numb+i], temp,count, check_spin,limit,config);
				}				
			}
			break;
		default:
			break;				
	}	
}


void cluster_move(int cluster_or_not,double cluster_ratio,unsigned char*** config,double* bonds,double* tail,double** energy,int* neig_index,int* neig,double* temperature,dsfmt_t* rn){
        unsigned char base = 1;
	if (cluster_or_not){
		double randn = 0.0;
                int count = 0;
                char* check_spin = new char[spin_numb];
                for (int i=0;i<=spin_numb-1;i++){
                	*(check_spin+i) = 0;
                }
                double limit = cluster_ratio*spin_numb;
                for (int temp=0;temp<=temp_set-1;temp++){
               	/*--------------check cluster size---------------*/
                	int cluster_size = 0;
                        for (register int index=0; index<=spin_numb-1; ++index){
				if (((*(*(*(config+0)+temp)+index))^(*(*(*(config+1)+temp)+index)))){
                                	++cluster_size;
                                }
                       	}
                        if(cluster_size > (spin_numb/2)){
                        	for (int index = 0;index <= spin_numb-1;++index){
                                	(*(*(*(config+0)+temp)+index)) = ((*(*(*(config+0)+temp)+index))^base);
                                }
                                *(energy[0]+temp) = energy_func(*(config[0]+temp),bonds,tail);
                        }
                        /*-------------------------------*/
                        count = 0;
                        int index = int(dsfmt_genrand_close_open(rn)*(spin_numb-1));
                        cluster(neig,neig_index,index,temp,count,check_spin,limit,config);
			/*flip*/
                       	double energy_diff = 0;
                        unsigned char** temp_spin  = new unsigned char*[2];
                        temp_spin[0] = new unsigned char[spin_numb];
                        temp_spin[1] = new unsigned char[spin_numb];
                        for(int j=0;j<=spin_numb-1;j++){
                        	if (check_spin[j] == (temp+1)){
                                	*(*(temp_spin+0)+j) = *(*(*(config+1)+temp)+j);
                                        *(*(temp_spin+1)+j) = *(*(*(config+0)+temp)+j);
				}
                               	else{
                                	*(*(temp_spin+0)+j) = *(*(*(config+0)+temp)+j);
                                       	*(*(temp_spin+1)+j) = *(*(*(config+1)+temp)+j);
                               	}
                       	}
                     //   energy_diff = (energy_func(*(temp_spin+0),bonds,tail)+ energy_func(*(temp_spin+1),bonds,tail)) - (energy_func(*(*(config+0)+temp),bonds,tail)+ energy_func(*(*(config+1)+temp),bonds,tail));
	      //                  printf("%lf\n",energy_diff);
 		//	randn = dsfmt_genrand_close_open(rn);
	        //        if (randn<(exp(-energy_diff/temperature[temp]))){
                      if (randn<2.0){   // flip without condition
               			unsigned char base = 1;
				for (int j=0;j<=spin_numb-1;j++){
					if (check_spin[j] == (temp+1)){
						*(*(*(config+0)+temp)+j) = (*(*(*(config+0)+temp)+j))^ base;
						*(*(*(config+1)+temp)+j) = (*(*(*(config+1)+temp)+j))^ base;
					}
				}
			}
			delete[] temp_spin[0];
			delete[] temp_spin[1];
			delete[] temp_spin;
			/*update energy for cluster*/
			for (int rep = 0;rep <= replica-1;rep++){
				*(energy[rep]+temp) = energy_func(*(config[rep]+temp),bonds,tail);
			}
		}
	}

}
void link(double* link_overlap,unsigned char*** config,int* neig_h, int* neig_index_h,int& neig_numb){
 	double temp_link;
        for (int t = 0;t <= temp_set-1; t++){
        	link_overlap[t] = 0.0;
                for (int i = 0; i <= spin_numb-1;++i){
                	for (int j = 0;j <= neig_h[i]-1; j++){
                        	temp_link = (2*(*(*(config[0]+t)+i))-1)*(2*(*(*(config[0]+t)+neig_index_h[i*spin_numb+j]))-1);
                                for (int rep = 1; rep <= replica-1; rep++){
                                	temp_link *= (2*(*(*(config[rep]+t)+i))-1)*(2*(*(*(config[rep]+t)+neig_index_h[i*spin_numb+j]))-1);
                                }
                                link_overlap[t] += temp_link;
                        }
                }
                link_overlap[t] =  link_overlap[t]/neig_numb;
        }
}

void set_ktables(double* kcostable,double* ksintable,int l){
	double pi = 0.0, q = 0.0;
	pi = 4*atan(1.0);
	q = 2*pi/((double) l);

	for (int i = 0; i <= spin_numb-1;i++){
		*(kcostable+i) = cos((double) q*((i)%l));
		*(ksintable+i) = sin ((double) q*((i)%l));
	}
}   

double chi_func(unsigned char*** config, double* kcostable,double* ksintable,int& t,int index){
	double qkx = 0.0, qky = 0.0, qk2 = 0.0;
	if (index != 0){
		for (int j =0;j <= spin_numb-1;j++){
			qkx += ((2*(*(*(config[0]+t)+j))-1)*(2*(*(*(config[1]+t)+j))-1)*kcostable[j]);
			qky += ((2*(*(*(config[0]+t)+j))-1)*(2*(*(*(config[1]+t)+j))-1)*ksintable[j]);
		}
	}
	else if (index == 0){
		for (int j =0;j <= spin_numb-1;j++){
                        qkx += (2*(*(*(config[0]+t)+j))-1)*(2*(*(*(config[1]+t)+j))-1);
                }
		qky = 0.0;
	}
	qk2 = (qkx*qkx) + (qky*qky);	
	return qk2;
}



