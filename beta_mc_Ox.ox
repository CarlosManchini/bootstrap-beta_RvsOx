#include <oxstd.oxh>
#include <oxprob.h>
#import <maximize>
					   
static decl s_vx;
const decl nobs=500;

loglik(const vP, const adFunc, const avScore, const amHess){

 	//decl vone = ones(1,nobs);
	//adFunc[0] = double( nobs*loggamma(vP[0]+vP[1]) - nobs*loggamma(vP[0]) - nobs*loggamma(vP[1])
	//				    + (vP[0]-1)*(vone *log(s_vx)) + (vP[1]-2)*(vone *log(1-s_vx)) );


	adFunc[0] = double( sumc( loggamma(vP[0]+vP[1]) - loggamma(vP[0]) - loggamma(vP[1])
				+ vP[0]*log(s_vx) + vP[1]*log(1-s_vx) - log(s_vx) - log(1-s_vx) ));

	//adFunc[0] = double( sumc(log(densbeta(s_vx, vP[0],vP[1]))) );

	if(avScore){
	 	(avScore[0])[0] = nobs*polygamma(vP[0]+vP[1],0)-nobs*polygamma(vP[0],0)+sumc(log(s_vx));
		(avScore[0])[1] = nobs*polygamma(vP[0]+vP[1],0)-nobs*polygamma(vP[1],0)+sumc(log(1-s_vx));
	}
	
	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0;
	else
		return 1;
}

main(){
	decl exectime;
	exectime = timer();
	decl x, vp, dfunc, ir, nu, falhas=0;
	decl R=1000, i, emv, x_bar, sigma2_hat, alpha_hat, beta_hat, emm;
	
	ranseed("MWC_52"); //ranseed("MWC_52");
	ranseed(369);
	
	nu = <5.0;3.0>;
							
	emv = zeros(R,2);
	emm = zeros(R,2);

	for(i=0 ; i<R ; i++){
		s_vx = ranbeta(nobs,1,nu[0],nu[1]);

		x_bar = meanc(s_vx);
		sigma2_hat = (1/(nobs-1))*sumc((s_vx-x_bar).^2);
		alpha_hat =	x_bar*(((x_bar*(1-x_bar))/sigma2_hat)-1);
		beta_hat = (1-x_bar)*(((x_bar*(1-x_bar))/sigma2_hat)-1);
		vp = alpha_hat|beta_hat;
		emm[i][0] = vp[0];
		emm[i][1] = vp[1];
		
		ir = MaxBFGS(loglik, &vp, &dfunc, 0, FALSE);
		if(ir == MAX_CONV || ir == MAX_WEAK_CONV){
			emv[i][0] = vp[0];
			emv[i][1] = vp[1];
		} 	
		else{
			falhas++; print("falha na converg");
    	}
	 }

print("\nReplicas: ",R,"\nn = ", nobs, "\nFalhas Converg: ",falhas);
print("\nMedias EMV: ", meanc(emv[][0]), meanc(emv[][1]) );//double(meanc(emv)));
print("\nMedias EMM: ", meanc(emm[][0]), meanc(emm[][1]) );

print("\nViés EMV\n", (double(meanc(emv[][0]))-nu[0]),"\n",(double(meanc(emv[][1]))-nu[1]),"\n" );
print("\nViés EMM\n", (double(meanc(emm[][0]))-nu[0]),"\n",(double(meanc(emm[][1]))-nu[1]),"\n" );

print("\nVR (%) EMV\n", (double(meanc(emv[][0]))-nu[0])/nu[0]*100,
							"\n",(double(meanc(emv[][1]))-nu[1])/nu[1]*100,"\n" );
print("\nVR (%) EMM\n", (double(meanc(emm[][0]))-nu[0])/nu[0]*100,
							"\n",(double(meanc(emm[][1]))-nu[1])/nu[1]*100,"\n" );
							
print("\nErro padrão EMV", sqrt(varc(emv[][0])),sqrt(varc(emv[][1]))  );
print("\nErro padrão EMM", sqrt(varc(emm[][0])),sqrt(varc(emm[][1]))  );

print("\nEQM EMV", varc(emv[][0])+((meanc(emv[][0]))-nu[0])^2, varc(emv[][1])+((meanc(emv[][1]))-nu[1])^2  );
print("\nEQM EMM", varc(emm[][0])+((meanc(emm[][0]))-nu[0])^2, varc(emm[][1])+((meanc(emm[][1]))-nu[1])^2  );


println( "\nEXECUTION TIME: ", timespan(exectime) );

}
