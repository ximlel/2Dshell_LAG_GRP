#include<stdio.h>
#include<math.h>
#include<malloc.h>
#define sec (0.5) //sec=0 1st order godunov,sec=0.5 godunov MUSCL
#define pi (4.*atan(1.0))
#define Md 420 // max vector dimension
#define Mt 25  // max theta dimension
#define GAMMAL (1.4)
#define GAMMAR (3.0) //Ratio of special heats Gamma=1.4 or 3.0
#define Ncell 410  // Number of computing cells in r direction
#define Ndis 400
#define Tcell 20 // Number of computing cells in theta direction
#define Diaph1 (10.)
#define Diaph2 (10.2)
#define Domlen (10.25) // Domain length
#define Timeout (0.01)   // Output time
#define CFL (0.45)// CFL condition
#define m (3.)// m=1 planar;m=2 cylindrical;m=3 spherical
#define Epsilon (1.5) //GRP limiter parameter
#define DL0 (0.00129)
#define DR0 (19.237)
#define UL0 (0.)
#define UR0 (-200)
#define PL0 (1.01325)
#define PR0 (1.01325)
#define KK 0 //ALE number
#include"./inp.h"
#include"./VIPLimiter.h"
int main()
{	//parameters
	double Gamma, GammaL=GAMMAL, GammaR=GAMMAR; 
	double DL=DL0,DR=DR0,UL=UL0,UR=UR0,PL=PL0,PR=PR0;//D:Density;U:Velocity;P:Pressure
	double CL,CR;//Sound speed
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);
	if(2.0*CL/(GammaL-1.)+2.0*CR/(GammaR-1.)<=UR-UL)
		{
			printf("Error:Vacuum is generated by initial data!\n");
			return 0;
		}
	double dt; //delta_t
	double F1[Md],F2[Md],F3[Md],E1[Md],E2[Md],E3[Md],Speed1[Md],Speed2[Md];
	//flux, conservative variable and wave speed
	double M[Md];//(derivative)centers distance
	double PM,UM,DML,DMR,Smax,time=0.;//P_star, U_star, rho_starL, roh_starR, max wave speed
	double dr,r,S,D,U,P;//initial d_raidus
	dr=(double)Domlen/Ncell;
	double dtheta;//initial d_angle
	dtheta=0.5*pi/Tcell;
	double RR[Md],DD[Md],UU[Md],PP[Md],CC[Md],GammaGamma[Md];//radius and variable in cells
	RR[0]=0.;
	RR[Ncell+1]=Domlen*cos(0.5*dtheta); 
	double Ddr[Md],deltar;//distance of boundary in a cell
	Ddr[0]=0.5*dr*cos(0.5*dtheta);
	double Rb[Md],Ub[Md];//radius and moving velocity of outer cell boundary
	double mass[Md],vol[Md],Rbh[Md];//Rbh:Rb half time step
	int i,k;
	for(i=1;i<=Ncell;i++)
		{
			RR[i]=(i-0.5)*dr*cos(0.5*dtheta);
			Ddr[i]=dr*cos(0.5*dtheta);
			Rb[i]=i*dr*cos(0.5*dtheta);
			if(RR[i]<=Diaph1)
				{
					DD[i]=DL;
					UU[i]=UL;
					PP[i]=PL;
				}
			else if(RR[i]<=Diaph2)
				{
					DD[i]=DR;
					UU[i]=UL;
					PP[i]=PR;
				}
			else
				{
					DD[i]=DR;
					UU[i]=UR;
					PP[i]=PR;
				}				
			if(i<=Ndis)
				GammaGamma[i]=GammaL;
			else
				GammaGamma[i]=GammaR;						
			CC[i]=sqrt(GammaGamma[i]*PP[i]/DD[i]);
			M[i] =Ddr[i];
			E1[i]=DD[i];//*pow(RR[i],m-1.);
			E2[i]=UU[i];//*pow(RR[i],m-1.);
			E3[i]=(0.5*UU[i]*UU[i]+PP[i]/(DD[i]*(GammaGamma[i]-1.)));//*pow(RR[i],m-1.);
		}//initial value
	Rb[0]=0.;
	Rbh[0]=0.;
	for(i=1;i<=Ncell;i++)
		{
			vol[i]=RR[i]*RR[i]*(Rb[i]-Rb[i-1]);//m=3.
			mass[i]=DD[i]*vol[i];
		}
	deltar=dr*cos(0.5*dtheta);
	E3[0]=E3[1]*DD[1];
	E1[0]=E1[1];
	for(i=1;i<=Ncell;i++)
		{
			deltar=(deltar<Ddr[i]?deltar:Ddr[i]);
			if(deltar<0)
				printf("error\n");
		}
	Smax=0.;
	Smax=(Smax>(fabs(UL)+CL)?Smax:(fabs(UL)+CL));
	Smax=(Smax>(fabs(UR)+CR)?Smax:(fabs(UR)+CR));
	if(Smax>0.)
		{
			dt=CFL*deltar/Smax;
			printf("Smax=%lf,deltar=%lf\n",Smax,deltar);
		}
	else		
		dt=CFL*0.1*dr;
	printf("dt=%lf\n",dt);
	M[0]=M[1];
	M[Ncell+1]=M[Ncell];
	RR[0]=0.;
	RR[Ncell+1]=Domlen*cos(0.5*dtheta);
	Rb[0]=0.;
	Rb[Ncell+1]=Domlen*cos(0.5*dtheta);
	Ddr[0]=Ddr[1];
	Ddr[Ncell+1]=Ddr[Ncell];
    DD[0]=DD[1];
	UU[0]=UU[1];
	PP[0]=PP[1];
	DD[Ncell+1]=DD[Ncell];
	UU[Ncell+1]=UU[Ncell];
	PP[Ncell+1]=PP[Ncell];//boundary
	
	double DmD[Md],DmU[Md],DmP[Md],Dmtau[Md],slopeL,slopeR,DDL,DDR,DUL,DUR,DPL,DPR,tau_star,C_star,Dttau,DtU,DtP,DtD,TDSL,TDSR,DpsiL,DphiR,//DtauL,DtauR,
		Us,Ps,Ds;//derivative
	double Umin[Md],Pmin[Md],Dmin[Md],sD,stau,sU,sP,C_starL,C_starR,F2P[Md],F1P[Md],F3P[Md];
	for(i=1;i<=Ncell;i++)
		{
			DmD[i]=(DD[i]-DD[i-1])/M[i];
			DmU[i]=(UU[i]-UU[i-1])/M[i];
			DmP[i]=(PP[i]-PP[i-1])/M[i];
			Dmtau[i]=(1./DD[i]-1./DD[i-1])/M[i];
		}
	for(i=1;i<=Ncell;i++)
		{
			DmD[i]=0.;
			DmU[i]=0.;
			DmP[i]=0.;
			Dmtau[i]=0.;
		}
	DmD[0]=DmD[1];
	DmU[0]=DmU[1];
	DmP[0]=DmP[1];
	DmD[Ncell+1]=DmD[Ncell];
	DmU[Ncell+1]=DmU[Ncell];
	DmP[Ncell+1]=DmP[Ncell];

	Dmtau[0]=Dmtau[1];
	Dmtau[Ncell+1]=Dmtau[Ncell];
	dt=1e-4;
	DD[0]=DD[1];
	UU[0]=UU[1];
	PP[0]=PP[1];
	DD[Ncell+1]=DD[Ncell];
	UU[Ncell+1]=UU[Ncell];
	PP[Ncell+1]=PP[Ncell];
	for(k=1;k<=1000000;k++)
		{
		
			for(i=1;i<=Ncell;i++)
				{
					GammaL = GammaGamma[i];
					GammaR = GammaGamma[i+1];
					slopeL=DmD[i];
					slopeR=DmD[i+1];
					DL=DD[i]+0.5*Ddr[i]*slopeL;
					DR=DD[i+1]-0.5*Ddr[i+1]*slopeR;
					slopeL=DmU[i];
					slopeR=DmU[i+1];
					UL=UU[i]+0.5*Ddr[i]*slopeL;
					UR=UU[i+1]-0.5*Ddr[i+1]*slopeR;
					slopeL=DmP[i];
					slopeR=DmP[i+1];
					PL=PP[i]+0.5*Ddr[i]*slopeL;
					PR=PP[i+1]-0.5*Ddr[i+1]*slopeR;
					CL=sqrt(GammaL*PL/DL);
					CR=sqrt(GammaR*PR/DR);
					DUL=DmU[i];
					DUR=DmU[i+1];
					DPL=DmP[i];
					DPR=DmP[i+1];
					DDL=DmD[i];
					DDR=DmD[i+1];
					StarPU(PM,UM,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					if(PM>PL)//left shock
						{
							Speed1[i]=UL-CL*sqrt(PM/PL*(GammaL+1)/(2.*GammaL)+(GammaL-1.)/(2.*GammaL));//Euler
							DML=DL*(Speed1[i]-UL)/(Speed1[i]-UM);
						}
					else//left fan
						{
							Speed1[i]=UL-CL;
							DML=DL*pow(PM/PL,1./GammaL);
						}
					if(PM>PR)//right shock
						{
							Speed2[i]=UR+CR*sqrt(PM/PR*(GammaR+1)/(2.*GammaR)+(GammaR-1.)/(2.*GammaR));
							DMR=DR*(Speed2[i]-UR)/(Speed2[i]-UM);
						}
					else//right fan
						{
							Speed2[i]=UR+CR;
							DMR=DR*pow(PM/PR,1./GammaR);
						}
				}//end for 1 round
			for(i=1;i<Ncell;i++)
				Smax=(Smax>fabs(Speed1[i])?Smax:fabs(Speed1[i]));
			for(i=1;i<Ncell;i++)
				Smax=(Smax>fabs(Speed2[i])?Smax:fabs(Speed2[i]));
			for(i=1;i<=Ncell;i++)
				{
					deltar=(deltar<Ddr[i]?deltar:Ddr[i]);
					if(deltar<0)
						printf("error\n");
				}	
			if(Smax>0.)
				dt=CFL*deltar/Smax;			
			else
				dt=CFL*0.1*dr;
			if(time<Timeout&&(time+dt)>Timeout)
				dt=Timeout-time;//compute for time step
	
			for(i=1;i<=Ncell;i++)
				{
					GammaL = GammaGamma[i];
					GammaR = GammaGamma[i+1];
					slopeL=DmD[i];
					slopeR=DmD[i+1];
					DL=DD[i]+0.5*Ddr[i]*slopeL;
					DR=DD[i+1]-0.5*Ddr[i+1]*slopeR;
					slopeL=DmU[i];
					slopeR=DmU[i+1];
					UL=UU[i]+0.5*Ddr[i]*slopeL;
					UR=UU[i+1]-0.5*Ddr[i+1]*slopeR;
					slopeL=DmP[i];
					slopeR=DmP[i+1];
					PL=PP[i]+0.5*Ddr[i]*slopeL;
					PR=PP[i+1]-0.5*Ddr[i+1]*slopeR;
					CL=sqrt(GammaL*PL/DL);
					CR=sqrt(GammaR*PR/DR);
					DUL=DmU[i];
					DUR=DmU[i+1];
					DPL=DmP[i];
					DPR=DmP[i+1];
					DDL=DmD[i];
					DDR=DmD[i+1];
					StarPU(PM,UM,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					if(PM>PL)//left shock
						{
							Speed1[i]=UL-CL*sqrt(PM/PL*(GammaL+1)/(2.*GammaL)+(GammaL-1.)/(2.*GammaL));//Euler
							DML=DL*(Speed1[i]-UL)/(Speed1[i]-UM);
						}
					else//left fan
						{
							Speed1[i]=UL-CL;
							DML=DL*pow(PM/PL,1./GammaL);
						}
					if(PM>PR)//right shock
						{
							Speed2[i]=UR+CR*sqrt(PM/PR*(GammaR+1)/(2.*GammaR)+(GammaR-1.)/(2.*GammaR));
							DMR=DR*(Speed2[i]-UR)/(Speed2[i]-UM);
						}
					else//right fan
						{
							Speed2[i]=UR+CR;
							DMR=DR*pow(PM/PR,1./GammaR);
						}
					if(i<=KK)
						S=0.;
					else
						S=UM;
					r=Rb[i]; 
					Sample(D,U,P,Gamma,PM,UM,S,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					if((DL-DR)*(DL-DR)+(UL-UR)*(UL-UR)+(PL-PR)*(PL-PR)<1e-6)//Acoustic case
						{
							tau_star=1./D;
							C_star=sqrt(Gamma*P/D);
							AcousticS(DtD,DtU,DtP,DDL,DDR,DUL,DUR,DPL,DPR,D,U,C_star,S,r);
						}
					else// non Acoustic case
						{
							TDSL=-CL*CL/(DL*(GammaL-1.))*DDL+1./(DL*(GammaL-1.))*DPL;
							TDSR=-CR*CR/(DR*(GammaR-1.))*DDR+1./(DR*(GammaR-1.))*DPR;
							DpsiL=DUL+GammaL/((GammaL-1.)*CL*DL)*DPL-CL/(DL*(GammaL-1.))*DDL;
							DphiR=DUR-GammaR/((GammaR-1.)*CR*DR)*DPR+CR/(DR*(GammaR-1.))*DDR;
							C_starL=CL*pow(PM/PL,(GammaL-1.)/(2.*GammaL));
							C_starR=CR*pow(PM/PR,(GammaR-1.)/(2.*GammaR));
							if(S>=(UM-C_starL)&&S<=(UM+C_starR))
								{//nonsonic case
									GRPsolverES(DtD,DtU,DtP,UM,PM,DL,DR,UL,UR,PL,PR,DDL,DDR,DUL,DUR,DPL,DPR,TDSL,TDSR,DpsiL,DphiR,S,r,GammaL,GammaR);
								}
							else
								{//sonic case
									GRPsolverSS(DtD,DtU,DtP,UM,PM,DL,DR,UL,UR,PL,PR,DDR,DUR,DPR,TDSL,TDSR,DpsiL,DphiR,S,r,GammaL,GammaR);
								}
							if(S<=(UL-CL)||S>=(UR+CR))
								GRPsolverSS1(DtD,DtU,DtP,DL,DR,UL,UR,PL,PR,DDL,DDR,DUL,DUR,DPL,DPR,S,r,GammaL,GammaR);								
						}
					Ub[i]=S+0.5*dt*DtU;
					if(i<=KK)
						Ub[i]=0.;
				}//end for 2 round
			Rb[Ncell+1]=Rb[Ncell];
			for(i=1;i<=Ncell;i++)
				{
					GammaL = GammaGamma[i];
					GammaR = GammaGamma[i+1];
					slopeL=DmD[i];
					slopeR=DmD[i+1];
					DL=DD[i]+0.5*Ddr[i]*slopeL;
					DR=DD[i+1]-0.5*Ddr[i+1]*slopeR;
					slopeL=DmU[i];
					slopeR=DmU[i+1];
					UL=UU[i]+0.5*Ddr[i]*slopeL;
					UR=UU[i+1]-0.5*Ddr[i+1]*slopeR;
					slopeL=DmP[i];
					slopeR=DmP[i+1];
					PL=PP[i]+0.5*Ddr[i]*slopeL;
					PR=PP[i+1]-0.5*Ddr[i+1]*slopeR;
					CL=sqrt(GammaL*PL/DL);
					CR=sqrt(GammaR*PR/DR);
					DUL=DmU[i];
					DUR=DmU[i+1];
					DPL=DmP[i];
					DPR=DmP[i+1];
					DDL=DmD[i];
					DDR=DmD[i+1];
					StarPU(PM,UM,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					if(PM>PL)//left shock
						{
							Speed1[i]=UL-CL*sqrt(PM/PL*(GammaL+1)/(2.*GammaL)+(GammaL-1.)/(2.*GammaL));//Euler
							DML=DL*(Speed1[i]-UL)/(Speed1[i]-UM);
						}
					else//left fan
						{
							Speed1[i]=UL-CL;
							DML=DL*pow(PM/PL,1./GammaL);
						}
					if(PM>PR)//right shock
						{
							Speed2[i]=UR+CR*sqrt(PM/PR*(GammaR+1)/(2.*GammaR)+(GammaR-1.)/(2.*GammaR));
							DMR=DR*(Speed2[i]-UR)/(Speed2[i]-UM);
						}
					else//right fan
						{
							Speed2[i]=UR+CR;
							DMR=DR*pow(PM/PR,1./GammaR);
						}
					if(i<=KK)
						S=0.;
					else
						S=Ub[i];
					r=Rb[i]; 
					Sample(D,U,P,Gamma,PM,UM,S,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					if((DL-DR)*(DL-DR)+(UL-UR)*(UL-UR)+(PL-PR)*(PL-PR)<1e-6)//Acoustic case
						{
							tau_star=1./D;
							C_star=sqrt(Gamma*P/D);
							AcousticS(DtD,DtU,DtP,DDL,DDR,DUL,DUR,DPL,DPR,D,U,C_star,S,r);
						}
					else// non Acoustic case
						{
							TDSL=-CL*CL/(DL*(GammaL-1.))*DDL+1./(DL*(GammaL-1.))*DPL;
							TDSR=-CR*CR/(DR*(GammaR-1.))*DDR+1./(DR*(GammaR-1.))*DPR;
							DpsiL=DUL+GammaL/((GammaL-1.)*CL*DL)*DPL-CL/(DL*(GammaL-1.))*DDL;
							DphiR=DUR-GammaR/((GammaR-1.)*CR*DR)*DPR+CR/(DR*(GammaR-1.))*DDR;
							C_starL=CL*pow(PM/PL,(GammaL-1.)/(2.*GammaL));
							C_starR=CR*pow(PM/PR,(GammaR-1.)/(2.*GammaR));
							if(S>=(UM-C_starL)&&S<=(UM+C_starR))
								{//nonsonic case
									GRPsolverES(DtD,DtU,DtP,UM,PM,DL,DR,UL,UR,PL,PR,DDL,DDR,DUL,DUR,DPL,DPR,TDSL,TDSR,DpsiL,DphiR,S,r,GammaL,GammaR);
								}
							else
								{//sonic case
									GRPsolverSS(DtD,DtU,DtP,UM,PM,DL,DR,UL,UR,PL,PR,DDR,DUR,DPR,TDSL,TDSR,DpsiL,DphiR,S,r,GammaL,GammaR);
								}
							if(S<=(UL-CL)||S>=(UR+CR))							
								GRPsolverSS1(DtD,DtU,DtP,DL,DR,UL,UR,PL,PR,DDL,DDR,DUL,DUR,DPL,DPR,S,r,GammaL,GammaR);
						}
					Rbh[i]=Rb[i]+Ub[i]*dt*0.5;
					Rb[i]=Rb[i]+Ub[i]*dt;
					Us=U+0.5*dt*DtU;
					Ps=P+0.5*dt*DtP;
					Ds=D+0.5*dt*DtD;
					F1[i]=Ds;
					F2[i]=Ps;
					F3[i]=Ps*Us;
					if(m==1.)
						{
							F2P[i]=0.;
							F1P[i]=0.;
							F3P[i]=0.;
						}
					else
						{
							F2P[i]=-(m-1.)/Rb[i]*Ds*Us*Us;
							F1P[i]=-(m-1.)/Rb[i]*Ds*Us;
							F3P[i]=-(m-1.)/Rb[i]*Us*(Ps+Ds*(0.5*Us*Us+Ps/(Ds*(Gamma-1.))));
						}
					Umin[i]=U+dt*DtU;
					Pmin[i]=P+dt*DtP;
					Dmin[i]=D+dt*DtD;
				}//end for 3 round
			DD[0]=DD[0]-dt*2.*F1[1]*F3[1]/F2[1]/dr;
			E3[0]=E3[0]-dt*2.*(F3[1]/dr+F3[1]/F2[1]*F1[1]*(0.5*F3[1]/F2[1]*F3[1]/F2[1]+F2[1]/((GammaGamma[1]-1)*F1[1]))/dr);
			UU[0]=E3[0]*(GammaGamma[0]-1);
			F2[0]=F2[1]-DmP[1]*(Rbh[1]-0.);
			F3[0]=0.;
			PP[0]=UU[0];
			UU[0]=0.;		
			for(i=1;i<=Ncell;i++)
				{
					vol[i]=(Rb[i]-Rb[i-1])*0.5*(Rb[i]+Rb[i-1])*0.5*(Rb[i]+Rb[i-1]);
					DD[i]=mass[i]/vol[i];
					M[i]=Rb[i]-Rb[i-1];
					Ddr[i]=M[i];
				}		
			M[0]=M[1];
			M[Ncell+1]=M[Ncell];
			Ddr[0]=Ddr[1];
			Ddr[Ncell+1]=Ddr[Ncell];
			Rb[0]=0.;
			Rbh[0]=0.;
			for(i=1;i<Ncell;i++)//m=3
				{
					E2[i]=E2[i]-dt/mass[i]*(F2[i]*Rbh[i]*Rbh[i]-F2[i-1]*Rbh[i-1]*Rbh[i-1])+0.5*dt*(Rbh[i]-Rbh[i-1])*2.*0.5*(Rbh[i]+Rbh[i-1])*(F2[i]+F2[i-1])/mass[i];
					E3[i]=E3[i]-dt/mass[i]*(F3[i]*Rbh[i]*Rbh[i]-F3[i-1]*Rbh[i-1]*Rbh[i-1]);
				}
			E2[1]=0.;
			//Decoding to get physical variables
			for(i=1;i<Ncell;i++)
				{
					RR[i]=0.5*(Rb[i]+Rb[i-1]);
					UU[i]=E2[i];
					PP[i]=(E3[i]-0.5*UU[i]*UU[i])*(GammaGamma[i]-1.)*DD[i];
				}
			UU[1]=0.;
			UU[0]=0.;
			DD[Ncell]=DD[Ncell-1];
			UU[Ncell]=UU[Ncell-1];
			PP[Ncell]=PP[Ncell-1];	
			DD[Ncell+1]=DD[Ncell];
			UU[Ncell+1]=UU[Ncell];
			PP[Ncell+1]=PP[Ncell];
			for(i=2;i<Ncell;i++)
				{
					sU=(Umin[i]-Umin[i-1])/M[i];sP=(Pmin[i]-Pmin[i-1])/M[i];sD=(Dmin[i]-Dmin[i-1])/M[i];
					DmD[i]=minmod(Epsilon*(DD[i]-DD[i-1])/(RR[i]-RR[i-1]),sD,Epsilon*(DD[i+1]-DD[i])/(RR[i+1]-RR[i]));
					DmU[i]=minmod(Epsilon*(UU[i]-UU[i-1])/(RR[i]-RR[i-1]),sU,Epsilon*(UU[i+1]-UU[i])/(RR[i+1]-RR[i]));
					DmP[i]=minmod(Epsilon*(PP[i]-PP[i-1])/(RR[i]-RR[i-1]),sP,Epsilon*(PP[i+1]-PP[i])/(RR[i+1]-RR[i]));
				}
			DmD[1]=minmod2((Dmin[1]-DD[1])/(0.5*M[1]),DmD[2]);
			DmU[1]=minmod2((Umin[1]-UU[1])/(0.5*M[1]),DmU[2]);
			DmP[1]=minmod2((Pmin[1]-PP[1])/(0.5*M[1]),DmP[2]);
			DmU[Ncell]=DmU[Ncell-1];
			DmP[Ncell]=DmP[Ncell-1];
			DmD[Ncell]=DmD[Ncell-1];
			DmU[Ncell+1]=DmU[Ncell];
			DmP[Ncell+1]=DmP[Ncell];
			DmD[Ncell+1]=DmD[Ncell];

			time=time+dt;
			for(i=0;i<=(Ncell+1);i++){Ddr[i]=M[i];}
			printf("Time[%d]=%e,dt=%e\n",k,time,dt);
			if(fabs(time-Timeout)<1e-6){break;}
		}//end k
	DD[0]=DD[1];
	PP[0]=PP[1];
	DD[Ncell+1]=DD[Ncell];
	UU[Ncell+1]=UU[Ncell];
	PP[Ncell+1]=PP[Ncell];
	RR[0]=0.;
	RR[Ncell+1]=Domlen*cos(0.5*dtheta);
	FILE *out,*outs;
	outs=fopen("../data_out/datas.m","w");
	fprintf(outs,"RR=[");
	Write(outs,RR,Ncell+1);
	fprintf(outs,"DD=[");
	Write(outs,DD,Ncell+1);
	fprintf(outs,"UU=[");
	Write(outs,UU,Ncell+1);
	fprintf(outs,"PP=[");
	Write(outs,PP,Ncell+1);	
	DL=DL0;
	DR=DR0;
	UL=UL0;
	UR=UR0;
	PL=PL0;
	GammaL=GAMMAL;
	GammaR=GAMMAR;
	PR=PR0;
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);
	printf("OK!\n");
	StarPU(PM,UM,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
	printf("P=%lf,U=%lf\n",PM,UM);
	int j;double rb[Md][Mt],zb[Md][Mt],DD2[Md][Mt],UUxi2[Md][Mt],PP2[Md][Mt];
	for(i=0;i<=Ncell;i++)
		for(j=0;j<=Tcell;j++)
			{
				rb[i][j]=Rb[i]/cos(0.5*dtheta)*sin(j*dtheta);
				zb[i][j]=Rb[i]/cos(0.5*dtheta)*cos(j*dtheta);
			}
	for(i=1;i<Ncell;i++)
		for(j=1;j<Tcell;j++)
			{
				DD2[i][j]=0.5*(DD[i]+DD[i+1]);
				UUxi2[i][j]=0.5*(UU[i]+UU[i+1]);
				PP2[i][j]=0.5*(PP[i]+PP[i+1]);
			}
	for(i=1;i<Ncell;i++)
		{
			DD2[i][0]=DD2[i][1];
			UUxi2[i][0]=UUxi2[i][1];
			PP2[i][0]=PP2[i][1];
			DD2[i][Tcell]=DD2[i][Tcell-1];
			UUxi2[i][Tcell]=UUxi2[i][Tcell-1];
			PP2[i][Tcell]=PP2[i][Tcell-1];
		}
	for(j=0;j<=Tcell;j++)
		{
			DD2[0][j]=DD[0];
			UUxi2[0][j]=UU[0];
			PP2[0][j]=PP[0];
			DD2[Ncell][j]=DD2[Ncell-1][j];
			UUxi2[Ncell][j]=UUxi2[Ncell-1][j];
			PP2[Ncell][j]=PP2[Ncell-1][j];
		}
	out=fopen("../data_out/datad2d.dat","w");
	wrin2s(out,rb,zb,DD2);
	out=fopen("../data_out/datap2d.dat","w");
	wrin2s(out,rb,zb,PP2);
	out=fopen("../data_out/datau2d.dat","w");
	wrin2s(out,rb,zb,UUxi2);
	return 1;
}