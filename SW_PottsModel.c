#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include "rand.c"
// 2D Ising: z=0.35 with sw, z=2.215 with metropolis
// 3D Ising: z=0.75 with sw, z=2.0 with metropolis        
// tau ~ L**z
// MCDIS should be 10*tau (according to prof.Wangjs)

/* number of Potts states, spin = 0,1...Q */                                             
int Q; 
/* dimension */
int D;       
/* lattice linear size */                                                   
int L;    
/* total number of spins */                                            
int N; 
/* total Monte Carlo sweeps */                                       
int MCTOT; 
/* sweeps discarded in the beginning */                                      
int MCDIS;
/* In PottsModel, Tr = Tc/2 in IsingModel: D=2 Tc=2.2692;D=3 Tc= 4.5115;D=4 Tc=D*/                              
double Tr;
/* Pc: at Tr, Probability of generate a bond between two nodes with same spin*/
double P; 
/* Potts spin, 0 to Q-1 */                                                  
int *s;
/* bond[D*i+d] at site i, direction d */                                              
int *bond;
/* site i has label list[i] */                             
int *list;                                       

void sw(void);
void init(double T);
int energy(void);
int magnet(void);
char* get_name(char* name, int L, int D, double T, const char bc[4]);

int main(int argc, char **argv)
{
   int e,m;
   int mc,d;
   int ii,jj,kk;
   double T;
   int gen_slices = 0;
   int is_fBC = 0;
   char *e_fname, *m_fname, *log_fname, *spin_fname, *slices_fname, *cube_fname;
   const char any[30]; // to be discarded from infile
   const char BC[4]; // fBC/pBC in infile
   FILE *fin, *e_file, *m_file, *log_file, *spin_file, *slices_file, *cube_file;
   fin = stdin;
   log_file = stdout;

   /* "r"read, "a"append, "w"write */
   /* get parameters*/
   fin = fopen(argv[1], "r");
   assert(fin != NULL);
   fscanf(fin, "%s", &(*any));
   fscanf(fin, "%lf", &T);
   fscanf(fin, "%s", &(*any));
   fscanf(fin, "%d", &Q);
   fscanf(fin, "%s", &(*any));
   fscanf(fin, "%d", &D);
   fscanf(fin, "%s", &(*any));
   fscanf(fin, "%d", &L);
   fscanf(fin, "%s", &(*any));
   fscanf(fin, "%d", &MCTOT);
   fscanf(fin, "%s", &(*any));
   fscanf(fin, "%s", &(*BC));
   // fscanf(fin, "%s", &(*any));
   // fscanf(fin, "%d", &MCDIS);
   fclose(fin);

   assert(D>0 && L>0);
   is_fBC = (BC[0] == 'f');
   assert(argc > 2); gen_slices = (atoi(argv[2])); 
   /* N = L^D is total number of sites */
   N = 1; for(d = 0; d < D; ++d) N *= L;

   /* get file names */
   e_fname = get_name("e_dist", L, D, T, BC);
   m_fname = get_name("m_dist", L, D, T, BC);
   log_fname = get_name("logger", L, D, T, BC);
   spin_fname = get_name("spin_lattice", L, D, T, BC);
   slices_fname = get_name("slices", L, D, T, BC);
   cube_fname = get_name("cube", L, D, T, BC);

   /* initialize s[], list[] */
   init(T);   
   /* record parameters*/
   log_file = fopen(log_fname, "w"); assert(log_file != NULL);
   fprintf(log_file,"\n T = %.4f\n Q = %d\n D = %d\n L = %d\n MCTOT = %d\n is_fBC = %d\n gen_slices = %s\n",Tr*2,Q,D,L,MCTOT,is_fBC,argv[2]);
   fclose(log_file);

   if(gen_slices != 0){
      spin_file = fopen(spin_fname, "r"); assert(spin_file != NULL);
      for(ii = 0; ii < N; ++ii)
         fscanf(spin_file, " %d", &s[ii]);
      fclose(spin_file);
   }

   /* MonteCarlo Using Swendsen_Wang Algorithm*/
   for(mc = 1; mc <= MCTOT; ++mc) {
      sw();
      /* do statistics */  
      e_file = fopen(e_fname, "a"); assert(e_file != NULL); 
      e = energy(); fprintf(e_file, "%d\n", e);
      fclose(e_file);
      m_file = fopen(m_fname, "a"); assert(m_file != NULL); 
      m = magnet(); fprintf(m_file, "%d\n", m);
      fclose(m_file);

      // slice is just the row in the middle of 2D lattice, or the plane in the middle of 3D lattice
      if(L == 1000){
         slices_file = fopen(slices_fname, "a"); assert(slices_file != NULL);
         if(D == 2)
            for(ii = (int)(L*L/2+0); ii < (int)(L*L/2+L); ++ii)
               fprintf(slices_file, " %d", s[ii]);
         else if(D == 3)
            for(ii = (int)(L*L*L/2+L*L/2+0); ii < (int)(L*L*L/2+L*L/2+L); ++ii)
               fprintf(slices_file, " %d", s[ii]);
         fprintf(slices_file, "\n");
         fclose(slices_file);
      }

      // cube is the central part of the whole lattice
      // if(L == 1000){
      //    cube_file = fopen(cube_fname, "a"); assert(cube_file != NULL);
      //    if(D == 3)
      //       for(ii = (int)(L/2-50); ii < (int)(L/2+50); ++ii)
      //          for(jj = (int)(L/2-50); jj < (int)(L/2+50); ++jj)
      //             for(kk = (int)(L/2-50); kk < (int)(L/2+50); ++kk)
      //                fprintf(cube_file, " %d", s[ii*L*L+jj*L+kk]);
      //    else if(D == 2)
      //          for(jj = (int)(L/2-50); jj < (int)(L/2+50); ++jj)
      //             for(kk = (int)(L/2-50); kk < (int)(L/2+50); ++kk)
      //                fprintf(cube_file, " %d", s[jj*L+kk]);
      //    fprintf(cube_file, "\n");
      //    fclose(cube_file);
      }

      if(mc%100 == 0){
         log_file = fopen(log_fname, "a"); assert(log_file != NULL);
         fprintf(log_file, "mc=%d, ", mc);
         fclose(log_file);
      }
   }

   /* record final spin_lattice */
   spin_file = fopen(spin_fname, "w+"); assert(spin_file != NULL);
   for(ii = 0; ii < N; ++ii)
      fprintf(spin_file, " %d",s[ii]);
   fprintf(spin_file, "\n");
   fclose(spin_file);

   return 0;
}

char* get_name(char* name, int L, int D, double T, const char bc[4]){
   char Lc[10];
   char Dc[2];
   char Tc1[2];
   char Tc2[5];
   int T1 = (int)floor(T);
   int T2 = (int)floor(T*10000)-T1*10000;
   sprintf(Lc, "%d", L);
   sprintf(Dc, "%d", D);
   sprintf(Tc1, "%d", T1);
   sprintf(Tc2, "%d", T2);
   
   char* final_name;
   final_name = malloc(strlen(name)+2+strlen(Lc)+2+strlen(Dc)+2+strlen(Tc1)+1+strlen(Tc2)+1+strlen(bc)+1);
   strcpy(final_name, name);
   strcat(final_name, "_L");
   strcat(final_name, Lc);
   strcat(final_name, "_D");
   strcat(final_name, Dc);
   strcat(final_name, "_T");
   strcat(final_name, Tc1);
   strcat(final_name, ".");
   strcat(final_name, Tc2);
   strcat(final_name, "_");
   strcat(final_name, bc);
   return final_name;
}

/* Do initialization of spins, coupling, mat[], etc */
void init(double T)
{
   int i, k;

   list = (int *) malloc(N*sizeof(int));
   assert(list != NULL);

   s = (int *) malloc(N*sizeof(int));
   assert(s != NULL);

   if(D == 2) {
      Tr = (1.0/(log(sqrt((double)Q)+1.0))); //1.1346 = (Tc/2=2.2692/2)
      Tr = T/2;
      P = 1.0 - exp(-1.0/Tr);
      //P = (1.0-1.0/(sqrt((double)Q)+1.0)); //0.58
   } else if (D == 3 && Q == 2) {
      Tr = 1.0/(0.221657*2.0); 
      Tr = T/2;            //2.25 = (Tc/2=4.5/2)            /* using Potts scaled */
      P = 1.0 - exp(-1.0/Tr);
   } else {
      Tr = D/2.0;
      Tr = T/2;
      P = 1.0 - exp(-1.0/Tr);
   }

   srand64(time(NULL));
   
   for(i = 0; i < N; ++i) 
      s[i] = Q * drand64();
      
}


int energy(void)
{
   int b, k;
   int i, ip, ie, si, j;
   int r, p, q;


   ie = 0;
   for(i = 0; i < N; ++i) {
      si = s[i];
      assert(si >= 0 && si < Q);
      r = i;
      p = 1 - L;
      q = 1;
      for(j = 0; j < D; ++j) {
         ip = (r + 1) % L == 0 ? i + p : i + q;
         ie += (si==s[ip]);
         r = r/L;
         p *= L;
         q *= L;
      }
   }
   return  -ie;
}

int magnet(void)
{
   int i;
   int m = 0;
   for(i = 0; i < N; ++i)
      m += s[i];
   return m;
}

void sw(void)     /* Swendsen-Wang algorithm, perform one Swendsen-Wang step */
{
   int i, ip, j, cnt, inc, a, b, min, max;
   int r, p, q;

   /* initially each site is a cluster by itself */
   for(i = 0; i < N; ++i)                   
      list[i] = i;               
      
   for(i = 0; i < N; ++i) {                      /* set the bond with prob P */
      r = i;
      p = 1 - L;
      q = 1;
      for(cnt = 0; cnt < D; ++cnt) {
         ip = (r + 1) % L == 0 ? i + p : i + q;
         /* implement Hoshen-Kopelman */
         if(s[i] == s[ip] && drand64() < P) {      
            a = list[i];
            /* run through until a == list[a] */
            while (a > list[a]) {             
               a = list[a];
            } 
            b = list[ip];
            while (b > list[b]) {
               b = list[b];
            } 
            /* find min and max of two labels */
            if (a > b) {                      
               min = b;
               max = a;
            } else {
               min = a;
               max = b;
            }
            list[max] = min;
            list[i] = min;
         }   
         r = r/L;
         p *= L;
         q *= L;      
      }
   }

   /* last sweep to make list pointing to the final label */
   inc = 0;           
   for(i = 0; i < N; ++i) {
      if(i == list[i]) {
         /* spin value over-written by label */
         s[i] = inc;  
         /* a new cluster find */                   
         ++inc;                                        
      } 
      else {
        j = list[i];
        while (j > list[j])
           j = list[j];
        assert(j < i);
        /* one of old cluster */
        s[i] = s[j];                                   
      }
   }

   assert(inc <= N);
   for(j = 0; j < inc; ++j) 
      list[j] = Q*drand64(); 
   /* new spin for each cluster */  
   /* old s[i] is cluster number */                 
   for(i = 0; i < N; ++i) {
      assert(s[i] < inc);
      s[i] = list[s[i]];                       
   }

}

