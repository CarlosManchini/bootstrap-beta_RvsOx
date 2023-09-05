#include <oxstd.oxh>
#include <oxprob.h>
#import <maximize>
					   
static decl s_vx;
static decl x_star;
const decl nobs=500;

loglik(const vP, const adFunc, const avScore, const amHess){
	//adFunc[0] = double( sumc(log(densbeta(s_vx, vP[0],vP[1]))) );
	adFunc[0] = double( sumc( loggamma(vP[0]+vP[1]) - loggamma(vP[0]) - loggamma(vP[1])
				+ vP[0]*log(s_vx) + vP[1]*log(1-s_vx) - log(s_vx) - log(1-s_vx) ));
	if(avScore){
	 	(avScore[0])[0] = nobs*polygamma(vP[0]+vP[1],0)-nobs*polygamma(vP[0],0)+sumc(log(s_vx));
		(avScore[0])[1] = nobs*polygamma(vP[0]+vP[1],0)-nobs*polygamma(vP[1],0)+sumc(log(1-s_vx));
	}
	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0;
	else
		return 1;
}

loglikB(const vP, const adFunc, const avScore, const amHess){
	//adFunc[0] = double( sumc(log(densbeta(x_star, vP[0],vP[1]))) );
		adFunc[0] = double( sumc( loggamma(vP[0]+vP[1]) - loggamma(vP[0]) - loggamma(vP[1])
				+ vP[0]*log(x_star) + vP[1]*log(1-x_star) - log(x_star) - log(1-x_star) ));
	if(avScore){
	 	(avScore[0])[0] = nobs*polygamma(vP[0]+vP[1],0)-nobs*polygamma(vP[0],0)+sumc(log(x_star));
		(avScore[0])[1] = nobs*polygamma(vP[0]+vP[1],0)-nobs*polygamma(vP[1],0)+sumc(log(1-x_star));
	}
	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
		return 0;
	else
		return 1;
}

main(){
	decl exectime;
	decl x, vp, dfunc, ir, nu, falhas=0;
	decl R=1000, i, emv, x_bar, sigma2_hat, alpha_hat, beta_hat, emm;
	decl B=500, emv_boot, j, irB, vpB, falhasB=0, dfuncB, emv_star, emm_star;
	decl x_bar_star, sigma2_hat_star, alpha_hat_star, beta_hat_star, emm_boot;
	exectime = timer();
	ranseed("MWC_52"); 
	ranseed(369);
	
	nu = <5.0;3.0>;	//<2.0;8.0>;
							
	emv = zeros(R,2);
	emm = zeros(R,2);
	emv_boot = zeros(B,2);
	emm_boot = zeros(B,2);
	emm_star = zeros(R,2);
	emv_star = zeros(R,2);

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
			   
		for(j=0; j<B; j++){
			x_star = ranbeta(nobs,1,emv[i][0],emv[i][1]); //vp[0],vp[1]

			x_bar_star = meanc(x_star);
			sigma2_hat_star = (1/(nobs-1))*sumc((x_star-x_bar_star).^2);
			alpha_hat_star = x_bar_star*(((x_bar_star*(1-x_bar_star))/sigma2_hat_star)-1);
			beta_hat_star = (1-x_bar_star)*(((x_bar_star*(1-x_bar_star))/sigma2_hat_star)-1);
			vpB =  alpha_hat_star|beta_hat_star;
			emm_boot[j][0] = vpB[0];
			emm_boot[j][1] = vpB[1];
			
			irB = MaxBFGS(loglikB, &vpB, &dfuncB, 0, FALSE);
			if(irB == MAX_CONV || irB == MAX_WEAK_CONV){
				emv_boot[j][0] = vpB[0];
				emv_boot[j][1] = vpB[1];
			} 	
			else{
				falhasB++; print("falha na converg BOOT");
    		}
		}
		emv_star[i][0] = 2*emv[i][0] - meanc(emv_boot[][0]);
		emv_star[i][1] = 2*emv[i][1] - meanc(emv_boot[][1]);

		emm_star[i][0] = 2*emm[i][0] - meanc(emm_boot[][0]);
		emm_star[i][1] = 2*emm[i][1] - meanc(emm_boot[][1]);
	 }

print("\nReplicas: ",R,"\nReplicas boot: ",B,"\nn = ", nobs, "\nFalhas Converg: ",falhas, "\nFalhas Converg BOOT: ",falhasB);

print("\nMedias EMV: ", meanc(emv)); //  meanc(emv[][0]), meanc(emv[][1])
print("\nMedias EMM: ", meanc(emm));	//meanc(emm[][0]), meanc(emm[][1])
print("\nMedias EMV_boot: ", meanc(emv_star)); //meanc(emv_star[][0]), meanc(emv_star[][1]) );
print("\nMedias EMM_boot: ", meanc(emm_star));
												 
print("\nViés EMV\n", (double(meanc(emv[][0]))-nu[0]),"\n",(double(meanc(emv[][1]))-nu[1]),"\n");
print("\nViés EMM\n", (double(meanc(emm[][0]))-nu[0]),"\n",(double(meanc(emm[][1]))-nu[1]),"\n");
print("\nViés EMV_boot\n", (double(meanc(emv_star[][0]))-nu[0]),"\n",(double(meanc(emv_star[][1]))-nu[1]),"\n");
print("\nViés EMM_boot\n", (double(meanc(emm_star[][0]))-nu[0]),"\n",(double(meanc(emm_star[][1]))-nu[1]),"\n");

print("\nVR (%) EMV\n", (double(meanc(emv[][0]))-nu[0])/nu[0]*100,
							"\n",(double(meanc(emv[][1]))-nu[1])/nu[1]*100,"\n");
print("\nVR (%) EMM\n", (double(meanc(emm[][0]))-nu[0])/nu[0]*100,
							"\n",(double(meanc(emm[][1]))-nu[1])/nu[1]*100,"\n");
print("\nVR (%) EMV_boot\n", (double(meanc(emv_star[][0]))-nu[0])/nu[0]*100,
							"\n",(double(meanc(emv_star[][1]))-nu[1])/nu[1]*100,"\n");
print("\nVR (%) EMM_boot\n", (double(meanc(emm_star[][0]))-nu[0])/nu[0]*100,
							"\n",(double(meanc(emm_star[][1]))-nu[1])/nu[1]*100,"\n");
							
print("\nErro padrão EMV", sqrt(varc(emv[][0])),sqrt(varc(emv[][1])));
print("\nErro padrão EMM", sqrt(varc(emm[][0])),sqrt(varc(emm[][1])));
print("\nErro padrão EMV_boot", sqrt(varc(emv_star[][0])),sqrt(varc(emv_star[][1])));
print("\nErro padrão EMM_boot", sqrt(varc(emm_star[][0])),sqrt(varc(emm_star[][1])));

print("\nEQM EMV", varc(emv[][0])+((meanc(emv[][0]))-nu[0])^2, varc(emv[][1])+((meanc(emv[][1]))-nu[1])^2);
print("\nEQM EMM", varc(emm[][0])+((meanc(emm[][0]))-nu[0])^2, varc(emm[][1])+((meanc(emm[][1]))-nu[1])^2);
print("\nEQM EMV_boot", varc(emv_star[][0])+((meanc(emv_star[][0]))-nu[0])^2, varc(emv_star[][1])+((meanc(emv_star[][1]))-nu[1])^2);
print("\nEQM EMM_boot", varc(emm_star[][0])+((meanc(emm_star[][0]))-nu[0])^2, varc(emm_star[][1])+((meanc(emm_star[][1]))-nu[1])^2);

println("\nEXECUTION TIME: ", timespan(exectime));	
}
