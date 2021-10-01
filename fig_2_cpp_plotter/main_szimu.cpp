#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>  /* time */
#include <math.h> //log iylesmi
#include <random>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <sys/time.h>

using namespace std;

struct szint //Ebben tárolom az egy szinthez tartozó rátákat és a szint várt méretét
{
    long double rscdplusd,rscd,racd,delta,pk,qk;
    long long int Nk,N_w;
    vector<long long int>N_m;
    vector<long long int>N_m_prev;
    vector<long long int>N_m_dt;
};



long double dt;//egy idõlépés hossza
long double ido; //amennyi ideig fut a szim.
long double legyartando;
vector<szint> szintek; //vektor amiben minden szinthez minden ráta meg van adva
long double mutrata;
long double szamlalo; // ebbe fogom összeadni a rátákat
long long int szintszam;
long long int instanciak;
long long int instancia;
long long int gen;
long long int maxMut;//maximálisan begyűjthető mutációk száma
long double randdriver;// mutrata*drivergének száma/összes   ~ 0.84*70/20000
long double s_acd;
long double s_scd;
long double s_scdd; //fittnessz erősség
long long int sumnacd;  // összesített driver mutációk az adott szinten acd,scd,scdplusdből
long long int sumnscd;
long long int sumnscdplusd;
long double utolsodelta;
long double idosum;
long double Nsuminstanc=0;
long double Nnegyzetatlag_Ref=0;
long double Nnegyzetatlag_Varied=0;
long double supernegyzetart_Varied;
long double Natlag_Ref=0;
long double Natlag_Varied;
long long int vartosszes;
long long int Nrak = 0;
long double B;
long long int rakosodas;
long double STOP;
long double Dszum=0;
long double Dszumnegyzet;
long double Dnegyzetszum=0;
long double Dnegyzetatlag=0;
long double Dn_atlagsum=0;
long double D_natlag=0;
long double legyartottatlag=0;
long long int rendszermeret=0;
long long int boolkapcs = 0;
long double differencialt_sum =0;
long double mert_tranziensido;
long double N_idoatlagsum=0;
long long int N_idokanter = 0;
long double gamma_manual = 0;

long double superatl_Ref=0;
long double superatl_Varied=0;
long double summ;
long double summNegyzet;
long double sosszeg;
long long int insrendszmeret;
long long int systemmutnumber=0;
long double atlsystemmutnumber=0;
long double atlmut_scdd = 0;
long double atlmut_scd = 0;
long long int mutsum_scdd=0;
long long int mutsum_scd=0;
long double stemMutTime = 0;
long long int mutansKimos = 0;
long long int Dload = 0;


long long int Nmutmax=0;
char flag_optimal; //opt gamma?
char flag_stemMut; //őssejt mut?
char flag_gnuplot;
char flag_idKiir;
char flag_kimenoKiir;
char flag_pkqk;
char flag_idoplotok;
long long int flag_tranziens=0;
char flag_legyartando_ido_Futas;
char flag_neutr_s;

string outputDir;
long long int mutCounterStrike;
long long int mutCounterStrikeszum = 0;
long double D_n=0;
long long int D_max_lineage = 0;
long long int D_max_lineage_sum = 0;
long double D_n_negyzet=0;
long long int D_n_kilepes = 0;
long long int D_n_sum=0;
long long int D_n_max=0;
long long int differencialt_sejtek=0;
long double ettol=0;
long double eddig=0;
long long int Renmut=0;
long long int maxMutszubpop=0;
long long int maxMutKimosszubpop=0;
long long int maxMutelement=0;
long long int maxMutelementKimos=0;
long double S = 0;
long long int S_rak = 0;
long long int ossejt_kritmut = 0;
long long int ossejt_kritmut_szam=0;
long long int idCounter = 0;
long long int Nsum = 0;

long long int scd_counter = 0.0;
long long int scdd_counter = 0.0;
long long int acd_counter = 0.0;
long long int scd_Mut_counter = 0.0;
long long int scdd_Mut_counter = 0.0;
long long int acd_Mut_counter = 0.0;

long long int scd_mut_selfrep=0.0;
long long int scdd_mut_selfrep=0.0;
long long int acd_mut_selfrep=0.0;


long double stem_scd_counter = 0.0;
long double stem_scdd_counter = 0.0;
long double stem_acd_counter = 0.0;
long long int elseCounter = 0;
long double tranziens_ido=0;
long long int kihalas = 0;
long double legyartott_sejtek=0.0;
long long int ossejt_kezdomutszam = 0;
long long int stopMut = 0;
long long int D_n_summa = 0;
long double NsumperTime = 0;
long double eltelt_ido_atl = 0;

long double eltelt_ido=0;
long double eltelt_ido_atlag=0;
double Kappa = 0;




void txtsorcsere(string ezt, string erre, ifstream& in,const char *inStreamFileName)
{
    ofstream out("tmp.txt");
    string tmp;
    while(true){
        getline(in,tmp);
        if(tmp==ezt){
            out<<erre<<endl;
        }else{
            out<<tmp<<endl;
        }
        if(in.eof()) break;
    }
    rename("tmp.txt",inStreamFileName);
    out.close();


}
void fancyRead(ifstream & ifs, string & tmp){

    while(true){
        getline(ifs,tmp);
        if(ifs.eof()){
            break;
        }else if(tmp.length()!=0){
            if(tmp.at(0)!='#') break;
        }
    }
}


ifstream f;
void txtsorcsere(string ezt, string erre, ifstream& in, const char *inStreamFileName);
void fancyRead(ifstream & ifs, string & tmp);






int main()
{
	
	uniform_real_distribution<long double> dis(0.0, 1.0); //random generálás


    f.open("adatok.txt");
    ofstream g; //kimenodadatok

    string tempRead;

    szintszam=100;

    while(szintszam>=100){
	    while(!f.eof()){

	        f>>tempRead;
	        if(tempRead=="instanciak="){
	            f>>instanciak;
	        }else if(tempRead=="legyartando="){
	            f>>legyartando;
	        }else if(tempRead=="ido="){
	            f>>ido;
	        }else if(tempRead=="rakosodas="){
	            f>>rakosodas;
	        }else if(tempRead=="STOP="){
	            f>>STOP;
	        }else if(tempRead=="szintszam="){
	            f>>szintszam;
	        }else if(tempRead=="mutrata="){
	            f>>mutrata;
	        }else if(tempRead=="s_acd="){
	            f>>s_acd;
	        }else if(tempRead=="s_scd="){
	            f>>s_scd;
	        }else if(tempRead=="s_scdd="){
	            f>>s_scdd;
	        }else if(tempRead=="Kappa="){
	            f>>Kappa;
	        }else if(tempRead=="maxMut="){
	            f>>maxMut;
	        }else if(tempRead=="flag_optimal="){
	            f>>flag_optimal;
	        }else if(tempRead=="gamma_manual="){
	            f>>gamma_manual;
	        }else if(tempRead=="flag_stemMut="){
	            f>>flag_stemMut;
	        }else if(tempRead=="outputDir="){
	            f>>outputDir;
	        }else if(tempRead=="flag_gnuplot="){
	            f>>flag_gnuplot;
	        }else if(tempRead=="ossejt_kezdomutszam="){
	            f>>ossejt_kezdomutszam;
	        }else if(tempRead=="stopMut="){
	            f>>stopMut;
	        }else if(tempRead=="flag_idKiir="){
	            f>>flag_idKiir;
	        }else if(tempRead=="flag_pkqk="){
	            f>>flag_pkqk;
	        }else if(tempRead=="flag_legyartando_ido_futas="){
	            f>>flag_legyartando_ido_Futas;
	        }else if(tempRead=="flag_neutr_s="){
	            f>>flag_neutr_s;
	        }else if(tempRead=="flag_idoplotok="){
	            f>>flag_idoplotok;
	        }else if(tempRead=="matrix="){
	            szintszam++;
	            for(int i=0; i<szintszam;i++)
	            {
	                szint ujszint;
	                f>>ujszint.delta;
	                f>>ujszint.Nk;
	                
			if(i==0){
	                    f>>ujszint.pk;
	                    ujszint.qk = 0;
	                }
	                if(flag_pkqk=='p'&& i>0) {
	                    f>>ujszint.pk;
	                    ujszint.qk = 0;
	                }

	                if(flag_pkqk=='q' && i>0){
	                    f>>ujszint.qk;
	                    ujszint.pk = 0;
	                }




	                ujszint.racd=0;
	                ujszint.rscd=0;
	                ujszint.rscdplusd=0;
	                szintek.push_back(ujszint);
	            }
	        }else{
	            getline(f,tempRead);
	        }
	    }
	    if(outputDir=="default"){
	        outputDir="";
	    }
	    f.close();
	}

	long double gamma_optimalis;
	gamma_optimalis=gamma_manual;
	utolsodelta=1;
	cout<<"legy: "<<legyartando<<endl;

	 if(gamma_manual<2.0)
	 {
	    cout<<"Nem fog elmenni a megadott szintszamig ekkora N-el, mert az optimális gamma tul nagy lesz. Ha ez zavar, noveld meg N-t"<<endl;
	 }

	 szintek[szintszam-1].delta=2.0;
	 szintek[szintszam-2].delta=1.0;
	 for(int i=szintszam-2;i>=0;i--){
	    szintek[i-1].delta=szintek[i].delta/gamma_optimalis;
	 }

	 for(int i=0;i<szintszam;i++)
	 {
	    if (i>0){
	        
	        if(flag_pkqk=='p'){
	        szintek[i].qk=2*(szintek[i-1].delta/szintek[i].delta)/szintek[i].pk;
		    }

		    else if(flag_pkqk=='q'){
		        szintek[i].pk=2*(szintek[i-1].delta/szintek[i].delta)/szintek[i].qk;
		    }

		    else{
		        cout<<"Nem sikerult valasztani, irj pk, vagy qk-t"<<endl;
		        return 1;
		    }
	    }


	    szintek[i].racd=szintek[i].delta*(1-szintek[i].pk)/szintek[i].Nk;
	    szintek[i].rscd=0.5*szintek[i].delta*(1-szintek[i].qk)*szintek[i].pk/szintek[i].Nk;
	    szintek[i].rscdplusd=0.5*szintek[i].delta*szintek[i].pk/szintek[i].Nk;


		if(szintek[i].pk>1 || szintek[i].qk>1){
		    cout<<"baj van, a pk vagy qk nagyobb lett 1-nel, valassz jobb pk-t ezen a szinten:"<<i<<endl;
		    cout<<szintek[i].pk<<endl;
		    cout<<szintek[i].qk<<endl;

		    if(gamma_optimalis<2.0)
		    {
		        cout<<"Az opt. gamma it túl kicsi (<2.0). Próbáld meg nagyobb N-re ha több szintet akarsz (az optimális görbe csak idáig tart)";
		    }

		    return 1;
		}

		if((szintek[i+1].delta/szintek[i].delta)<2.0 && i!=szintszam-1){
		    cout<<"baj van, a delta(k+1)/delta(k) kisebb lett 2.0-nel ezen a szinten:"<<i<<endl;
		    return 3;
		}

	}

	int kritmut = 0;
	double kritmut_float = 0;
	if(s_scd>0){
			kritmut = ceil((1.0/(gamma_manual-1))/s_scd) ;


			kritmut_float = (1.0/(gamma_manual-1))/s_scd;

	}   
	if(s_scdd>0){

			kritmut = ceil((1.0/(gamma_manual-1))/s_scdd) ;


			kritmut_float = (1.0/(gamma_manual-1))/s_scdd;

	} 

	if(s_scdd==0 && s_scd==0){

			kritmut = ceil((1.0/(gamma_manual-1))/s_acd) ;


			kritmut_float = (1.0/(gamma_manual-1))/s_acd;


	}

	for(int k=0;k<szintszam; k++){ //insitialize mutation vector

		for(int i=0;i<maxMut+1; i++){	
			int null = 0;
			szintek[k].N_m.push_back(null);
			szintek[k].N_m_prev.push_back(null);
			szintek[k].N_m_dt.push_back(null);
		}	
	}

	cout<<kritmut<<endl;

	vector<int> mutans_Kapcsolo;
	vector<int> talaltam_Mutanst;
	vector<long long int> talaltam_Mutanst_ossz;
	vector<long double> atlag_sejtszamok;

	for(int j=0; j<maxMut+1;++j){
	    int zero =0;
	    talaltam_Mutanst.push_back(zero);
	    talaltam_Mutanst_ossz.push_back(zero);
	    mutans_Kapcsolo.push_back(zero);
	    atlag_sejtszamok.push_back(zero);
	}

	long double transient_time = 0;
	for(int l=1;l<szintszam; l++){
		transient_time+=szintek[l].Nk/szintek[l-1].delta;	
	}
				


	for(int instancia = 0; instancia<instanciak; instancia++){

		random_device rd;
  		mt19937_64 generator (rd());  	

		cout<<"Instancia: "<<instancia<<endl;

		for(int j=0; j<maxMut+1;++j){
		    mutans_Kapcsolo[j]=0;
		}

		for(int k=0;k<szintszam; k++){ //insitialize mutation vector
			for(int i=0;i<maxMut+1; i++){	
				szintek[k].N_m[i]=0;
				szintek[k].N_m_prev[i]=0;
				szintek[k].N_m_dt[i]=0;
			}	
		}

		if(ido==0){
			ido =legyartando*szintek[0].Nk;

		}

		eltelt_ido = 0;

		legyartott_sejtek = 0;
		D_n=(szintszam-2)*(gamma_manual-1)+1;

		for(int k= 0; k<szintszam; k++){
			szintek[k].N_m[0] = szintek[k].Nk;
			szintek[k].N_m_prev[0] = szintek[k].Nk;	
		}
		
		int breakKapcs = 0;
		while(eltelt_ido<=ido){


			dt = (1.0/szintek[szintszam-1].rscdplusd)*2;
			
			//cout<<eltelt_ido<<'\t'<<ido<<endl<<endl;
	
			Nsum = 0;
			
			for(int k= 0; k<szintszam; k++){
				/*
				for(int l=0;l<maxMut+1;l++){
					szintek[k].N_m_prev[l] = szintek[k].N_m[l];
					szintek[k].N_m_dt[l] = (long long int) 0 ;	
				}
				*/
				szintek[k].N_m_prev = szintek[k].N_m; 
				fill(szintek[k].N_m_dt.begin(), szintek[k].N_m_dt.end(), 0); ;								
			}			

			for(int k=0;k<szintszam;k++){

				for(int i=0;i<maxMut; i++){
					
					
					if(k==szintszam-1){

						long double mean = dt*szintek[k].rscdplusd*szintek[k].N_m_prev[i];
						
						if(mean<10 and mean>0){

					    	poisson_distribution<long long int> distribution_scdd (mean);;
					    	scdd_mut_selfrep = distribution_scdd(generator);
					    }
					    else if(mean==0){
					    	scdd_mut_selfrep = 0;	
					    }
					    else{
					    	scdd_mut_selfrep =  mean;
					    }
						
						
						
						szintek[k].N_m_dt[i] -= scdd_mut_selfrep;
						legyartott_sejtek+=scdd_mut_selfrep;
						legyartottatlag+=scdd_mut_selfrep;
					}

					else{
			
						//mutated cells self reproduce
						if(k>0){

							long double scd_mut = dt*min(szintek[k].rscd+( pow((double)(i/kritmut_float),Kappa)*kritmut_float*s_scd*(szintek[k].delta-szintek[k-1].delta))/szintek[k].Nk,(long double)szintek[k].rscdplusd)*szintek[k].N_m_prev[i];		
							
							long double mean = scd_mut;
							if(mean<10 and mean>0){

						    	poisson_distribution<long long int> distribution_scd(mean);
						    	scd_mut_selfrep = distribution_scd(generator);
						    }
						    else if(mean==0){
						    	scd_mut_selfrep = 0;	
						    }

						    else{
						    	scd_mut_selfrep =  mean;
						    }

							long double scdd_mut = dt*max(szintek[k].rscdplusd - ( pow((double)(i/kritmut_float),Kappa)*kritmut_float*s_scdd*(szintek[k].delta-szintek[k-1].delta))/szintek[k].Nk,(long double) szintek[k].rscd)*szintek[k].N_m_prev[i];		
						
							
							mean =  scdd_mut;
							if(mean<10 and mean>0){

						    	poisson_distribution<long long int> distribution_scdd(mean);
						    	scdd_mut_selfrep = distribution_scdd(generator);
						    }

						    else if(mean==0){
						    	scdd_mut_selfrep = 0;	
						    }

						    else{
						    	scdd_mut_selfrep =  mean;
						    }						

						}

						else {
						
							scd_mut_selfrep = 0;
							scdd_mut_selfrep = 0;

						}
						

						
						long double mean = dt*szintek[k].racd*szintek[k].N_m_prev[i];
						if(mean<10 and mean>0){

					    	poisson_distribution<long long int> distribution_acd(mean);
					    	acd_mut_selfrep = distribution_acd(generator);
					    }
					    else if(mean==0){
					    	acd_mut_selfrep = 0;	
					    }
					    else{
					    	acd_mut_selfrep =  mean;
					    }

						szintek[k].N_m_dt[i] += scd_mut_selfrep - scdd_mut_selfrep;
						szintek[k+1].N_m_dt[i] += acd_mut_selfrep + 2*scdd_mut_selfrep;

						// cells acquire additional mutations
						//scd


		
						mean = scd_mut_selfrep*mutrata;
						long int mutscd1 = mean;
						long int mutscd2 = mean;
						if(mean<10 and mean>0){

					    	poisson_distribution<long long int> distribution_scd_Mut(mean);
					    	mutscd1 = distribution_scd_Mut(generator);
					    	mutscd2 = distribution_scd_Mut(generator);
					    }
					    else if(mean==0){
					    	mutscd1 = 0;
					    	mutscd2 = 0;	
					    }
					    else{
					    	mutscd1 =  mean;
					    	mutscd2 =  mean;
					    }


						szintek[k].N_m_dt[i+1] += mutscd1;
						szintek[k].N_m_dt[i] -= mutscd1;


						szintek[k].N_m_dt[i+1] += mutscd2;
						
						//mutations with scdd
						mean = scdd_mut_selfrep*mutrata;
						long int mutscdd1 = mean;
						long int mutscdd2 = mean;
						if(mean<10 and mean>0){

					    	poisson_distribution<long long int> distribution_scdd_Mut(mean);
					    	mutscdd1 = distribution_scdd_Mut(generator);
					    	mutscdd2 = distribution_scdd_Mut(generator);
					    }
					    else if(mean==0){
					    	mutscdd1 = 0;
					    	mutscdd2 = 0;	
					    }
					    else{
					    	mutscdd1 =  mean;
					    	mutscdd2 =  mean;

					    }

						
						szintek[k+1].N_m_dt[i+1] += mutscdd1;
						szintek[k].N_m_dt[i] -= mutscdd1;


						szintek[k+1].N_m_dt[i+1] += mutscdd2;
						

						//mutations with acd
						mean = acd_mut_selfrep*mutrata;
						long int mutacd1 = mean;
						long int mutacd2 = mean;
						if(mean<10 and mean>0){

					    	poisson_distribution<long long int> distribution_acd_Mut(mean);
					    	mutacd1 = distribution_acd_Mut(generator);
					    	mutacd2 = distribution_acd_Mut(generator);
					    }
					  	else if(mean==0){
					    	mutacd1 = 0;						    	
					    	mutacd2 = 0;
					    }
					    else{
					    	mutacd1 =  mean;
					    	mutacd2 =  mean;
					    }

						szintek[k+1].N_m_dt[i+1] += mutacd1;

						szintek[k].N_m_dt[i+1] += mutacd2;
						szintek[k].N_m_dt[i] -= mutacd2;						


					}

					szintek[k].N_m[i]=max(szintek[k].N_m_prev[i]+szintek[k].N_m_dt[i], (long long int) 0);

					if(szintek[k].N_m[i]>0 and mutans_Kapcsolo[i]==0){
									
						mutans_Kapcsolo[i]=1;
						talaltam_Mutanst[i]+=1;
						if(i==kritmut){
							breakKapcs = 1;
							break;
						}
					}
					
				}

				if(breakKapcs==1){
					break;
				}
				

			} // szintszam for vege
			


			//cout<<endl;
			eltelt_ido+=dt;
			eltelt_ido_atlag+=dt;
			if(breakKapcs==1){
				break;
			}

		} // while vege


	} //instanciak vege



g.open(outputDir+"kimenoadatok.txt" );

int w=17;
g<<"Az egyes szintek adatai:"<<endl<<endl;
g<<setw(5)<<"szint";
g<<setw(8)<<"Nk";
g<<setw(w)<<"racd";
g<<setw(w)<<"rscd";
g<<setw(w)<<"rscdplusd";
g<<setw(w)<<"pk";
g<<setw(w)<<"qk";
g<<setw(w)<<"delta"<<endl;
rendszermeret=0;

for(int i=0;i<szintszam;i++)
{
    if(i == szintszam-1){
        g<<setw(5)<<"TERM";
    }
    else{
        g<<setw(5)<<i;
    }

    g<<setw(8)<<szintek[i].Nk;
    g<<setw(w)<<szintek[i].racd;
    g<<setw(w)<<szintek[i].rscd;
    g<<setw(w)<<szintek[i].rscdplusd;
    g<<setw(w)<<szintek[i].pk;
    g<<setw(w)<<szintek[i].qk;
    g<<setw(w)<<szintek[i].delta<<endl;
}
g<<endl;
g<<"A szimulacio soran legyartott sejtek szama: "<<legyartottatlag/(long double)instanciak<<endl;
g<<"A szimulaciok soran eltelt ido atlaga: "<<eltelt_ido_atlag/(long double) instanciak<<endl;
g<<"D_n: "<<D_n_summa/(long double)instanciak<<endl;
g<<"Kritmut valsz: "<<Nrak/(long double)instanciak<<endl;
g<<"Mutaciok varhato erteke az adott szinten: "<<endl; 


for(int i = 0; i<maxMut; i++){
	g<<"az "<<i<<".mutans sejtszam/kezdeti sejtszam: "<<atlag_sejtszamok[i]/(long double)szintek[0].Nk/(long double) instanciak/(long double) ido / (long double) szintszam<<endl;
}
	
g<<"mc-1->mc: "<<talaltam_Mutanst[kritmut]/(long double)instanciak<<endl;	

for(int i=1; i<maxMut;i++){
	g<<"szimzu: "<<i<<" "<<talaltam_Mutanst[i]/(long double)instanciak<<endl;
	long double tmp=0;
	tmp=talaltam_Mutanst_ossz[i];
	
	g<<i<<"_simu_mossz: "<<tmp/(long double)instanciak/(long double)szintek[0].Nk<<endl;
}

g.close();
return 0;
}

