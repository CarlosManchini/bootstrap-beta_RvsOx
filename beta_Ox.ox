#include <oxstd.oxh>
#include <oxprob.oxh>

#import <maximize>

const decl nobs=50;
static decl x;

loglik(const vP, const adFunc, const avScore, const amHess){

	//decl vone = ones(1,nobs);
	//adFunc[0] = double( nobs*loggamma(vP[0]+vP[1]) - nobs*loggamma(vP[0]) - nobs*loggamma(vP[1])
	//				    + (vP[0]-1)*(vone *log(s_vx)) + (vP[1]-2)*(vone *log(1-s_vx)) );


	//adFunc[0] = double( sumc( loggamma(vP[0]+vP[1]) - loggamma(vP[0]) - loggamma(vP[1])
	//			+ vP[0]*log(s_vx) + vP[1]*log(1-s_vx) - log(s_vx) - log(1-s_vx) ));

	adFunc[0] = double( sumc(log(densbeta(x, vP[0],vP[1]))) );
	
	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0;
	else
		return 1;
}

main(){	  
	
	decl vp, dfunc, ir, nu, falhas=0, i,emv;
	
	ranseed("MWC_52"); //ranseed("MWC_52");
	
	nu = <3.0;4.0>;
	
	vp = <2.0;3.0>; //chute inicial	zeros(2,1)

	x = ranbeta(nobs,1,nu[0],nu[1]);
				
	ir = MaxBFGS(loglik, &vp, &dfunc, 0, TRUE);
	if(ir == MAX_CONV || ir == MAX_WEAK_CONV){
		print("SUCESSO!\n Log-vessom maximizada: ", dfunc,
			  "\n Valores verdadeiro:",nu',
			  "\n EMVs:", vp');
	} 	
	else{
		falhas++; print("falha na converg");
    }														  
}