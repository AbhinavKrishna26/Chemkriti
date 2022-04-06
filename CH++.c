#include <stdio.h> 
#include <math.h> 
#include <stdlib.h> 

#define Tc 647.096 /*K*/ 
#define pc 22064000 /*pa*/ 
#define rho_c 322 /*kg/cu m*/ 
#define alpha_0 1000 /*J/kg*/  
double phi_0=1000/647.096;
double a[6]={-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719,1.80122502};
double b[6]={1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45}; 
double c[6]={-2.03150240,-2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063}; 
double d[5]={-0.0000000565134998, 2690.66631, 127.287297, -135.003439, 0.981825814};
double d_alpha=-1135.905627715;
double d_phi=2319.5246;

double sat_pres(double sat_temp)
{
    double ans, temp_var, theta, tau; 
    theta=sat_temp/Tc; 
    tau=1-theta; 
    temp_var=(Tc/sat_temp)*a[0]*tau+a[1]*pow(tau,1.5)+a[2]*pow(tau,3) +a[3]*pow(tau,3.5)+a[4]*pow(tau,4)+a[5]*pow(tau,7.5);
    ans=pc*exp(temp_var); 
    return ans;
}  

double rho_l(double sat_temp)
{
    double ans, temp_var, theta, tau, lntau; 
    theta=sat_temp/Tc; 
    tau=1-theta; 
    lntau=log(tau); 
    temp_var=1+b[0]*exp(lntau/3.0)+ b[1]*exp(2.0*lntau/3.0+b[2]*exp(5.0*lntau/3.0)+b[3]*exp(16.0*lntau/3.0)+b[4]*exp(43.0*lntau/3.0)+b[5]*exp(110.0*lntau/3.0));  
    ans=rho_c*(temp_var*1.0); 
    return ans; 
} 
double rho_v(double sat_temp)
{
    double ans, temp_var, theta, tau, lntau;
    theta=sat_temp/Tc;  
    tau=1-theta; 
    lntau=log(tau); 
    temp_var= c[0]*exp(lntau/3.0)+c[1]*exp(2.0*lntau/3.0)+c[2]*exp(4.0*lntau/3.0)+c[3]*exp(3.0*lntau)+c[4]*exp(37.0*lntau/6.0)+c[5]*exp(71.0*lntau/6.0);  
    ans=rho_c*exp(temp_var*1.0); 
    return ans; 
}
double sp_vol_l(double sat_temp)
{
    double ans; 
    ans=1/rho_l(sat_temp); 
    return ans;
} 
 
double sp_vol_v(double sat_temp)
{
    double ans; 
    ans=1/rho_v(sat_temp);  
    return ans;
}  

double alpha(double sat_temp) 
{
    double ans, temp_var, theta; 
    theta=sat_temp/Tc; 
    temp_var=d_alpha+ d[0]*pow(theta,-19)+d[1]*pow(theta,1)+d[2]*pow(theta,4.5)+d[3]*pow(theta,5)+d[4]*pow(theta,-54.5);
    ans=alpha_0*temp_var;
    double phi (double sat_temp);
    return ans;
}

double phi(double sat_temp)
{
    double ans, temp_var, theta, lntheta;
    theta=sat_temp/Tc;
    lntheta=log(theta); 
    temp_var=d_phi+((19.0/20.0)*d[0]*exp(-20.0*lntheta))+(d[1]*lntheta)+((9.0/7.0)*d[02]*exp(3.5*lntheta))+((5.0/4.0)*d[3]*exp(4.0*lntheta))+((109.0/107.0)*d[4]*exp(53.5*lntheta));
    ans=phi_0*(temp_var*1.0); 
    return ans;
} 

double dp_dT(double pres, double sat_temp)
{
    double ans, temp_var, theta, tau, ln_tau0, ln_tau1, ln_tau2;
    theta=sat_temp/Tc; 
    tau=1-theta; 
    ln_tau0=6.5*log(tau);
    ln_tau1=2.5*log(tau);
    ln_tau2=0.5*log(tau);  
    temp_var=(7.5*a[5]*exp(ln_tau0))+(4*a[4]*pow(tau,3))+(3.5*a[3]*exp(ln_tau1))+(3*a[2]*pow(tau,2))+(1.5*a[1]*exp(ln_tau2))+a[0]+(log(pres/pc));
    ans=-1*(pres/sat_temp)*temp_var; 
    return ans; 
}

double enthalpy_l(double sat_temp) 
{
    double ans, pres; 
    pres=sat_pres(sat_temp) ; 
    ans=alpha(sat_temp)+((sat_temp/rho_l(sat_temp))*dp_dT(pres,sat_temp)); 
    return ans;
} 
double enthalpy_v(double sat_temp)
{
    double ans, pres; 
    pres=sat_pres(sat_temp); 
    ans=alpha(sat_temp/rho_v(sat_temp))*dp_dT(pres,sat_temp); 
    return ans;
}

double entropy_1(double sat_temp) 
{
    double ans, pres; 
    pres=sat_pres(sat_temp); 
    ans=phi(sat_temp)+((1.0/rho_l(sat_temp))*dp_dT(pres,sat_temp)); 
    return ans; 
}

double entropy_v(double sat_temp) 
{
    double ans, pres; 
    pres=sat_press(sat_temp); 
    ans=phi(sat_temp)+((1.0/rho_v(sat_temp))*dp_dT(pres,sat_temp)); 
    return ans;
} 

int main ()
{
    double temp, pres, dens_l, dens_v, vol_l, vol_v, h_l, h_v, s_l, s_v; 
    printf ("Steam Property Program for 273.15K > T < 647.14 K\n"); 
    printf ("Enter temperature (K); ");
    scanf ("&If",&temp); 
    pres=sat_pres(temp);
    dens_l=rho_1(temp); 
    dens_v=rho_v(temp);  
    vol_l=sp_vol_l(temp); 
    vol_v= sp_vol_v(temp); 
    h_l=enthalpy_l(temp); 
    h_v=enthalpy_v(temp); 
    s_l=entropy_l(temp); 
    s_v=entropy_v(temp); 
    printf("\n The saturated pressure is %lf Pa\n",pres); 
    printf("The density (liquid) is %lf kg/cu m\n",dens_l); 
    printf("The density (vapor) is %lf kg/ cu m\n",dens_v); 
    printf("The specific volume (liquid) is %lf cu m/ kg \n",vol_l); 
    printf("The specific volume (vapor) is %lf cu m/kg\n",vol_v); 
    printf("The enthalpy (liquid) is %lf J/kg\n",h_l); 
    printf("The enthalpy (vapor) is %lf J/kg\n",h_v); 
    printf("The entropy (liquid) is %lf J/kg-K\n",s_l);
    printf("The entropy (vapor) is %lf J/kg-K\n",s_v); 
    system("pause");
}
