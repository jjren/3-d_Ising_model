#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int Energy( int a[20][20] );
int Random();
double DRandom();
int main()
{
  int i, j, c, b, E1, E2, mesh_a, mesh_b, net_spin = 0, a[20][20]; 
  double tmp,T;

  srand(time(0));
  for( i = 0; i < 20; i++ )
  { for( j = 0; j < 20; j++ ) {
    if(Random()<=9) a[i][j] = 1; else a[i][j] = -1 ; }
  }
  printf( "Please set the temperature (number between 1 to 10, crital temperature set as 5):\n");
  scanf("%lf",&T);
  printf( "Please set the iteration number:\n");
  scanf("%d",&b);
  
  for( i = 0; i < 20; i++){ 
    for( j = 0; j < 20; j++) { printf("%2d", a[i][j]); }
    printf("\n");
 }

  for ( c = 1; c <= b; c++){
  E1 = Energy( a );
 // printf("E1=%d\n",E1);
  mesh_a = Random(); mesh_b = Random(); 
 // printf("%d,%d\n", mesh_a,mesh_b);
 // printf("inital%d\n",a[mesh_a][mesh_b]);
  if( a[mesh_a][mesh_b] == 1 ) {
	  a[mesh_a][mesh_b] = -1 ;
	  E2 = Energy( a ); 
      //  printf("E2=%d\n",E2);
	  if( E2 - E1 <= 0 ) {
		  a[mesh_a][mesh_b] = -1; }
	  else {
		  tmp=DRandom();
            //  printf("tmp%f\n",tmp);
            //  printf("%f\n",pow(e,-(E2-E1)/(2.203*T)));
	        if ( exp(-(E2-E1)/T) > tmp) {
		      a[mesh_a][mesh_b] = -1; }
	        else a[mesh_a][mesh_b] = 1;
	       }
   }
  else { 
	  a[mesh_a][mesh_b] = 1; 
	  E2 = Energy( a ); 
     //   printf("E2=%d\n",E2);
	  if( E2 - E1 <= 0 ) a[mesh_a][mesh_b] = 1; 
	  else {
		  tmp=DRandom();
            //  printf("tmp%f\n",tmp);
             // printf("%f\n",pow(e,-(E2-E1)/(2.203*T)));
	        if ( exp(-(E2-E1)/T) > tmp) a[mesh_a][mesh_b] = 1; 
	  else a[mesh_a][mesh_b] = -1;}
	 }
  
//  printf("after%d\n",a[mesh_a][mesh_b]);
//  printf("%d\n",E2-E1);
 }  
  
  printf( "After iteration, the spin configuration is shown here:\n" );
  for( i = 0; i < 20; i++){ 
    for( j = 0; j < 20; j++) { printf("%2d", a[i][j]); net_spin += a[i][j];}
    printf("\n");
 }
  printf( "The average of spin is %d",net_spin );
  E2 = Energy(a) ;
  printf("%d\n",E2);
 //  system("pause");
  return 0;
}

int Energy( int a[20][20] )
{
  int i, j, E;
  E=0;
  for( i = 0; i < 20; i++){ 
    for( j = 0; j < 20; j++) {
	    if(i<19) {
	      if( a[i][j] * a[i+1][j] == 1 ) E -= 1; else E +=1; }
	    if(j<19) {
            if( a[i][j] * a[i][j+1] == 1 ) E -= 1; else E +=1; }
    }
  }
  for( i = 0; i < 20; i++){ 
      if( a[i][19] * a[i][0] == 1 ) E -= 1; else E +=1;}
   for( j = 0; j < 20; j++) {if( a[19][j] * a[0][j] == 1 ) E -= 1; else E +=1; }
   return E; 
 }
    
int Random()
{
  double _d;
  int m;
  _d = (double)rand() / ((double)RAND_MAX + 1.0);
  m=(int)(_d * 20);
  return (m);
}

double DRandom()
{
  double _d;
  _d = (double)rand() / ((double)RAND_MAX + 1.0);
  // printf("_d=%f\n",_d);
  return (_d);
}
