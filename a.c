main(){
double m, sigma;
int i, PML_LAYERS;
double a;
m = 3.5;
PML_LAYERS = 10;	
	for(i=0;i<10;i++){
		a =(double) i/(PML_LAYERS-1);
		sigma=pow(a,m);
		printf("i %d sigma %e a %e m %e\n", i, sigma, a, m);
		
	}// end for


}
