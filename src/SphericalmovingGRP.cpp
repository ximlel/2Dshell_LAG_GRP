#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#define pi (4.*atan(1.0))
#define EPS (1e-8)
#define CFL (0.45) // CFL condition
#define m (2.)    // m=1 planar; m=2 cylindrical; m=3 spherical
#define Epsilon (1.) // r_0=Epsilon*dr
#include "./initdata.h"
#define Md Ncell+5 // max vector dimension
#define Mt Tcell+5  // max theta dimension
#include "./inp.h"
#include "./Riemann.h"
#include "./VIPLimiter.h"
int main()
{	//parameters
	double GammaL=GAMMAL, GammaR=GAMMAR; 
	double DL=DL0,DR=DR0,UL=UL0,UR=UR0,VL,VR,PL=PL0,PR=PR0;//D:Density;U,V:Velocity;P:Pressure
	double CL,CR;//Sound speed
	CL=sqrt(GammaL*PL/DL);
	CR=sqrt(GammaR*PR/DR);
	if(2.0*CL/(GammaL-1.)+2.0*CR/(GammaR-1.)<=UR-UL)
		{
			printf("Error:Vacuum is generated by initial data!\n");
			return 0;
		}
	double dt; //delta_t
	double F2[Md],F3[Md],E[Md],Speed1[Md],Speed2[Md];
	//flux, conservative variable and wave speed
	double dRc[Md];//(derivative)centers distance
	double PM,UM,DML,DMR,Smax_deltar,time=0.;//P_star, U_star, rho_starL, roh_starR, max wave speed
	double dr,r;//initial d_raidus
	dr=(double)Domlen/Ncell;
	double dtheta,dtheta_plot;//initial d_angle
	dtheta=0.5*pi/Tcell;
	dtheta_plot=0.5*pi/Tcell_plot;
	double RR[Md],DD[Md],UU[Md],PP[Md],CC[Md],GammaGamma[Md];//centroidal radius and variable in cells
	double DdrL[Md],DdrR[Md],Ddr[Md];//distance from boundary to center in a cell
	double Rb[Md],Lb[Md];//radius and length of outer cell boundary
	double Rbh[Md],Lbh[Md];//h: half time step
	double Rb_NStep,Lb_NStep;
	//double Rb_side[Md],Lb_side[Md],Rbh_side[Md],Lbh_side[Md],Sh[Md];
	double mass[Md],vol[Md];
	FILE *out,*outs;
	mkdir(DATAOUT, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	char file_data[FILENAME_MAX];
	double plot_t=D_PLOT_T;
	int i,j,k;
	
	for(i=1;i<=Ncell;i++)//center cell is cell 0
		{			
			Rb[i]=(Epsilon+i-1.)*dr*cos(0.5*dtheta);//outer cell boundary
			Lb[i]=2.*(Epsilon+i-1.)*dr*sin(0.5*dtheta);
			RR[i]=(Epsilon+i-(3.*Epsilon+3.*i-2.)/(6.*Epsilon+6.*i-3.))*dr*cos(0.5*dtheta);
			//centroid of triangular and trapezoid
			if(RR[i]<=Diaph1)
				{
					DD[i]=DL;
					UU[i]=UL;
					PP[i]=PL;
					GammaGamma[i]=GammaL;
				}
			else if(RR[i]<=Diaph2)
				{
					DD[i]=DR;
					UU[i]=UL;
					PP[i]=PR;
					GammaGamma[i]=GammaR;						
				}
			else if(RR[i]<=Diaph3)
				{
					DD[i]=DR;
					UU[i]=UR;
					PP[i]=PR;
					GammaGamma[i]=GammaR;						
				}
			else
				{
					DD[i]=DL1;
					UU[i]=UL1;
					PP[i]=PL1;
					GammaGamma[i]=GammaL;						
				}										
			CC[i]=sqrt(GammaGamma[i]*PP[i]/DD[i]);
			E[i]=0.5*UU[i]*UU[i]+PP[i]/(DD[i]*(GammaGamma[i]-1.));
		}//initial value
	Rb[0]=0.;
	Lb[0]=0.;
	RR[0]=(2./3.*Epsilon)*dr*cos(0.5*dtheta);
	DD[0]=DD[1];
	UU[0]=0.;	
	PP[0]=PP[1];
	GammaGamma[0]=GammaGamma[1];
	CC[0]=sqrt(GammaGamma[0]*PP[0]/DD[0]);
	E[0]=E[1]-0.5*UU[1]*UU[1];
	for(i=1;i<Ncell;i++)
		{
			DdrL[i]=Rb[i+1]-RR[i];
			DdrR[i]=RR[i]-Rb[i];
			Ddr[i] =DdrL[i]+DdrR[i];
			dRc[i] =RR[i]-RR[i-1];
			vol[i] =RR[i]*0.5*(Lb[i]+Lb[i+1])*(Rb[i+1]-Rb[i]);//m=2.
			mass[i]=DD[i]*vol[i];
		}
	DdrL[0]=Rb[1]-RR[0];
	Ddr[0] =Rb[1];
	dRc[0] =DdrL[0];
	vol[0] =RR[0]*0.5*Lb[1]*Rb[1];//m=2.
	mass[0]=DD[0]*vol[0];
	DdrR[Ncell]=DdrR[Ncell-1];
	Ddr[Ncell] =Ddr[Ncell-1];
	dRc[Ncell] =DdrL[Ncell-1];//boundary condition
	/*
	  DdrR[Ncell]=RR[Ncell]-Rb[Ncell];
	  Ddr[Ncell] =Ddr[Ncell-1];
	  dRc[Ncell] =RR[Ncell]-RR[Ncell-1];//boundary condition
	*/
	Rbh[0]=0.;
	Lbh[0]=0.;
	
	Smax_deltar=0.;
	Smax_deltar=(Smax_deltar>(fabs(UL)+CL)/dr?Smax_deltar:(fabs(UL)+CL)/dr);
	Smax_deltar=(Smax_deltar>(fabs(UR)+CR)/dr?Smax_deltar:(fabs(UR)+CR)/dr);
	dt=CFL/Smax_deltar;

	double DmD[Md],DmU[Md],DmP[Md],TmV[Md],slopeL,slopeR,DDL,DDR,DUL,DUR,DPL,DPR,TVL,TVR;//spacial derivative
	double C_star,DtU,DtP,DtDL,DtDR,TDSL,TDSR,DpsiL,DphiR,Us,Ps,Ds;//GRP variables
	double Umin[Md],VLmin[Md],Pmin[Md],DLmin[Md],DRmin[Md],sD,sU,sP,sV,C_starL,C_starR,F2P[Md];
	double rb[Md][Mt],zb[Md][Mt],DD2[Md][Mt],UUxi2[Md][Mt],PP2[Md][Mt],GammaGamma2[Md][Mt];
	double VIP_lim, Vave[4][2], V0[2], Vp1[2], Vp2[2], Vp3[2];//VIP limiter
	int wrong_idx = 0;
	/*
	  for(i=1;i<=Ncell;i++)
	  {
	  DmD[i]=(DD[i]-DD[i-1])/dRc[i];
	  DmU[i]=(UU[i]-UU[i-1])/dRc[i];
	  DmP[i]=(PP[i]-PP[i-1])/dRc[i];
	  }
	*/
	for(i=0;i<=Ncell;i++)
		{
			DmD[i]=0.;
			DmU[i]=0.;
			DmP[i]=0.;
			TmV[i]=0.;
		}
	for(k=1;k<=1e10;k++)
		{		
			for(i=0;i<Ncell;i++)
				{
					GammaL = GammaGamma[i];
					GammaR = GammaGamma[i+1];
					DDL=DmD[i];
					DDR=DmD[i+1];
					DUL=DmU[i];
					DUR=DmU[i+1];
					DPL=DmP[i];
					DPR=DmP[i+1];
					DL=DD[i]  +DdrL[i]  *DDL;
					DR=DD[i+1]-DdrR[i+1]*DDR;
					UL=UU[i]  +DdrL[i]  *DUL;
					UR=UU[i+1]-DdrR[i+1]*DUR;
					PL=PP[i]  +DdrL[i]  *DPL;
					PR=PP[i+1]-DdrR[i+1]*DPR;
					CL=sqrt(GammaL*PL/DL);
					CR=sqrt(GammaR*PR/DR);
						
					StarPU(PM,UM,DML,DMR,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					//Riemann_solver_exact(PM,UM,DML,DMR,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					//					if(i>=396 & i<=403)
					//						printf("%lf,%lf,%lf,%lf,%d\n",PM,UM,DML,DMR,i);
					
					if(PM>PL)//left shock
						Speed1[i+1]=UL-CL*sqrt(PM/PL*(GammaL+1)/(2.*GammaL)+(GammaL-1.)/(2.*GammaL));
					else//left fan
						Speed1[i+1]=UL-CL;
					if(PM>PR)//right shock
						Speed2[i+1]=UR+CR*sqrt(PM/PR*(GammaR+1)/(2.*GammaR)+(GammaR-1.)/(2.*GammaR));
					else//right fan
						Speed2[i+1]=UR+CR;
				}//end for 1 round
			Smax_deltar=0.;
			for(i=0;i<Ncell;i++)
				{
					Smax_deltar=(Smax_deltar>fabs(Speed1[i+1])/Ddr[i] ? Smax_deltar:fabs(Speed1[i+1])/Ddr[i]);
					Smax_deltar=(Smax_deltar>fabs(Speed2[i+1])/Ddr[i] ? Smax_deltar:fabs(Speed2[i+1])/Ddr[i]);
					Smax_deltar=(Smax_deltar>fabs(Speed1[i+1])/Ddr[i+1]?Smax_deltar:fabs(Speed1[i+1])/Ddr[i+1]);
					Smax_deltar=(Smax_deltar>fabs(Speed2[i+1])/Ddr[i+1]?Smax_deltar:fabs(Speed2[i+1])/Ddr[i+1]);
				}
			dt=CFL/Smax_deltar;
			if(time<Timeout&&(time+dt)>Timeout)
				dt=Timeout-time;//compute for time step

			for(i=0;i<Ncell;i++)
				{
					GammaL = GammaGamma[i];
					DDL=DmD[i];
					DUL=DmU[i];
					DPL=DmP[i];
					TVL=TmV[i];
					DL=DD[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*DDL;
					UL=UU[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*DUL;
					PL=PP[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*DPL;
					VL=0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)*TVL;
					CL=sqrt(GammaL*PL/DL);
					StarPU(PM,UM,DML,DMR,DL,DL,-UL*sin(0.5*dtheta)+VL*cos(0.5*dtheta),UL*sin(0.5*dtheta)-VL*cos(0.5*dtheta),PL,PL,CL,CL,GammaL,GammaL);
					//Riemann_solver_exact(PM,UM,DML,DMR,DL,DL,-UL*sin(0.5*dtheta),UL*sin(0.5*dtheta),PL,PL,CL,CL,GammaL,GammaL);
					r=0.5*(Rb[i]+Rb[i+1])/cos(0.5*dtheta);
					C_star=sqrt(GammaL*PM/DML);
					AcousticSLagTangent(DtP,DtU,DUL+TVL,DPL*cos(0.5*dtheta),DML,UL*cos(0.5*dtheta)+VL*sin(0.5*dtheta),C_star,r);
					F2P[i]=PM+dt*DtP;
					VLmin[i]=(UL*cos(0.5*dtheta)+VL*sin(0.5*dtheta)+dt*DtU)*sin(0.5*dtheta);
				}//end for 2 round
			
			for(i=0;i<Ncell;i++)
				{
					GammaL = GammaGamma[i];
					GammaR = GammaGamma[i+1];
					DDL=DmD[i];
					DDR=DmD[i+1];
					DUL=DmU[i];
					DUR=DmU[i+1];
					DPL=DmP[i];
					DPR=DmP[i+1];
					DL=DD[i]  +DdrL[i]  *DDL;
					DR=DD[i+1]-DdrR[i+1]*DDR;
					UL=UU[i]  +DdrL[i]  *DUL;
					UR=UU[i+1]-DdrR[i+1]*DUR;
					PL=PP[i]  +DdrL[i]  *DPL;
					PR=PP[i+1]-DdrR[i+1]*DPR;
					CL=sqrt(GammaL*PL/DL);
					CR=sqrt(GammaR*PR/DR);
					StarPU(PM,UM,DML,DMR,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					//Riemann_solver_exact(PM,UM,DML,DMR,DL,DR,UL,UR,PL,PR,CL,CR,GammaL,GammaR);
					r=Rb[i+1];
					if((DL-DR)*(DL-DR)+(UL-UR)*(UL-UR)+(PL-PR)*(PL-PR)+(CL-CR)*(CL-CR)<EPS)//Acoustic case
						{
							C_star=0.5*(CL+CR);
							AcousticSLag(DtDL,DtDR,DtU,DtP,DUL,DUR,DPL,DPR,0.5*(DL+DR),UM,C_star,r);
						}
					else// non-Acoustic case
						{
							TDSL=-CL*CL/(DL*(GammaL-1.))*DDL+1./(DL*(GammaL-1.))*DPL;
							TDSR=-CR*CR/(DR*(GammaR-1.))*DDR+1./(DR*(GammaR-1.))*DPR;
							DpsiL=DUL+GammaL/((GammaL-1.)*CL*DL)*DPL-CL/(DL*(GammaL-1.))*DDL;
							DphiR=DUR-GammaR/((GammaR-1.)*CR*DR)*DPR+CR/(DR*(GammaR-1.))*DDR;
							GRPsolverSLag(DtDL,DtDR,DtU,DtP,UM,PM,DL,DR,UL,UR,PL,PR,DDL,DDR,DUL,DUR,DPL,DPR,TDSL,TDSR,DpsiL,DphiR,r,GammaL,GammaR); 
						}
					Us=UM+0.5*dt*DtU;
					Ps=PM+0.5*dt*DtP;
					F2[i+1]=Ps;
					F3[i+1]=Ps*Us;			
					Umin[i+1] =UM +dt*DtU;
					Pmin[i+1] =PM +dt*DtP;
					DLmin[i+1]=DML+dt*DtDL;
					DRmin[i+1]=DMR+dt*DtDR;
					
					Rb_NStep=Rb[i+1]+Us*dt;
					Lb_NStep=2.*Rb_NStep*sin(0.5*dtheta)/cos(0.5*dtheta);
					Lbh[i+1]=0.5*(Lb[i+1]+Lb_NStep);
					Rbh[i+1]=(Rb[i+1]*(2.*Lb[i+1]+Lb_NStep)+Rb_NStep*(Lb[i+1]+2.*Lb_NStep))/(3.*(Lb[i+1]+Lb_NStep));
					Rb[i+1]=Rb_NStep;
					Lb[i+1]=Lb_NStep;
					RR[i]=Rb[i+1]-(2.*Lb[i]+Lb[i+1])/(3.*(Lb[i]+Lb[i+1]))*(Rb[i+1]-Rb[i]);
					/*
					  Rb_side[i]=0.25*(Lb[i]+Lb[i+1])/sin(0.5*dtheta);
					  Lb_side[i]=0.5*(Lb[i+1]-Lb[i])/sin(0.5*dtheta);
					  Rbh_side[i]=0.25*(Lbh[i]+Lbh[i+1])/sin(0.5*dtheta);
					  Lbh_side[i]=0.5*(Lbh[i+1]-Lbh[i])/sin(0.5*dtheta);
					  Sh[i]=Rbh[i+1]*Lbh[i+1]-Rbh[i]*Lbh[i]-sin(dtheta)*Rbh_side[i]*Lbh_side[i];
					*/
				}//end for 3 round			
			for(i=1;i<Ncell;i++)
				{
					DdrL[i]=Rb[i+1]-RR[i];
					DdrR[i]=RR[i]-Rb[i];
					Ddr[i] =DdrL[i]+DdrR[i];
					dRc[i] =RR[i]-RR[i-1];
					vol[i] =RR[i]*0.5*(Lb[i]+Lb[i+1])*(Rb[i+1]-Rb[i]);//m=2.
					DD[i]=mass[i]/vol[i];
					if(Ddr[i]<0.)
						{									
							printf("deltar<0,error!\n");
							return 0;
						}
				}
			DdrL[0]=Rb[1]-RR[0];
			Ddr[0] =Rb[1];
			dRc[0] =DdrL[0];
			vol[0] =RR[0]*0.5*Lb[1]*Rb[1];//m=2.
			DD[0]=mass[0]/vol[0];
			DdrR[Ncell]=DdrR[Ncell-1];
			Ddr[Ncell] =Ddr[Ncell-1];
			dRc[Ncell] =DdrL[Ncell-1];//boundary condition
			
			for(i=1;i<Ncell;i++)//m=2
				{
					UU[i]=UU[i]-dt/mass[i]*((F2[i+1]-F2P[i])*Rbh[i+1]*Lbh[i+1]-(F2[i]-F2P[i])*Rbh[i]*Lbh[i]);
					E[i] =E[i] -dt/mass[i]*(F3[i+1]*Rbh[i+1]*Lbh[i+1]-F3[i]*Rbh[i]*Lbh[i]);
				}
			UU[0]=0.;
			E[0]=E[0]-dt/mass[0]*(F3[1]*Rbh[1]*Lbh[1]);			
			//Decoding to get physical variables
			for(i=0;i<Ncell;i++)
				{				   				
					PP[i]=(E[i]-0.5*UU[i]*UU[i])*(GammaGamma[i]-1.)*DD[i];
					if(isnan(PP[i])||isnan(UU[i])||isnan(DD[i]))
						{									
							printf("variable is nan,error!\n");
							wrong_idx = 1;
						}
					else if (PP[i]<0)
						{									
							printf("p<0,error!\n");
							wrong_idx = 1;
						}
				}
			DD[Ncell]=DD[Ncell-1];
			UU[Ncell]=UU[Ncell-1];
			PP[Ncell]=PP[Ncell-1];
			if (wrong_idx)
				break;
			/*									
			for(i=1;i<Ncell;i++)
				{
					sU=(Umin[i+1] -Umin[i]) /Ddr[i];
					sP=(Pmin[i+1] -Pmin[i]) /Ddr[i];
					sD=(DLmin[i+1]-DRmin[i])/Ddr[i];
					DmD[i]=minmod(Alpha*(DD[i]-DD[i-1])/dRc[i],sD,Alpha*(DD[i+1]-DD[i])/dRc[i+1]);
					DmU[i]=minmod(Alpha*(UU[i]-UU[i-1])/dRc[i],sU,Alpha*(UU[i+1]-UU[i])/dRc[i+1]);
					DmP[i]=minmod(Alpha*(PP[i]-PP[i-1])/dRc[i],sP,Alpha*(PP[i+1]-PP[i])/dRc[i+1]);
				}
			DmD[0]=minmod2((DLmin[1]-DD[0])/dRc[0],DmD[1]);
			DmU[0]=minmod2((Umin[1]-UU[0]) /dRc[0],DmU[1]);
			DmP[0]=minmod2((Pmin[1]-PP[0]) /dRc[0],DmP[1]);
			DmD[Ncell]=minmod2((DLmin[Ncell]-DD[Ncell-1])/dRc[Ncell],DmD[Ncell-1]);
			DmU[Ncell]=minmod2((Umin[Ncell]-UU[Ncell-1]) /dRc[Ncell],DmU[Ncell-1]);
			DmP[Ncell]=minmod2((Pmin[Ncell]-PP[Ncell-1]) /dRc[Ncell],DmP[Ncell-1]);
			*/
			//VIP limiter update
			for(i=1;i<Ncell;i++)
				{
					sU=(Umin[i+1] -Umin[i]) /Ddr[i];
					sP=(Pmin[i+1] -Pmin[i]) /Ddr[i];
					sD=(DLmin[i+1]-DRmin[i])/Ddr[i];
					//sV=0.;
					//sV=VLmin[i]/(0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta));
					sV=UU[i]/(0.5*(Rb[i]+Rb[i+1]));
					Vave[0][0] = UU[i+1];
					Vave[0][1] = 0.;
					Vave[1][0] = UU[i-1];
					Vave[1][1] = 0.;
					Vave[2][0] = UU[i]*cos(dtheta);
					Vave[2][1] = UU[i]*sin(dtheta);
					Vave[3][0] = UU[i]*cos(dtheta);
					Vave[3][1] =-UU[i]*sin(dtheta);
					V0[0] = UU[i];
					V0[1] = 0.;
					Vp1[0] = UU[i]+(0.5*(Rb[i]+Rb[i+1])-RR[i])*sU;
					Vp1[1] = 0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)*sV;
					Vp2[0] = UU[i]+DdrL[i]*sU;
					Vp2[1] = 0.0;
					Vp3[0] = UU[i]-DdrR[i]*sU;
					Vp3[1] = 0.0;
					VIP_lim = fmin(1.0,     useVIPLimiter(4, Vave, V0, Vp1));
					VIP_lim = fmin(VIP_lim, useVIPLimiter(4, Vave, V0, Vp2));
					VIP_lim = fmin(VIP_lim, useVIPLimiter(4, Vave, V0, Vp3));
					DmU[i]=VIP_lim*sU;
					TmV[i]=VIP_lim*sV;
					if (abs(LIMITER_CONF)==1)
						{
							DmD[i]=minmod(Alpha*(DD[i]-DD[i-1])/dRc[i],sD,Alpha*(DD[i+1]-DD[i])/dRc[i+1]);
							DmP[i]=minmod(Alpha*(PP[i]-PP[i-1])/dRc[i],sP,Alpha*(PP[i+1]-PP[i])/dRc[i+1]);
							if(LIMITER_CONF>0)
								DmU[i]=minmod(Alpha*(UU[i]-UU[i-1])/dRc[i],sU,Alpha*(UU[i+1]-UU[i])/dRc[i+1]);
						}
					else if (abs(LIMITER_CONF)==2)
						{
							DmD[i]=minmod(Alpha*(DD[i]-DD[i-1])/2./DdrR[i],sD,Alpha*(DD[i+1]-DD[i])/2./DdrL[i]);
							DmP[i]=minmod(Alpha*(PP[i]-PP[i-1])/2./DdrR[i],sP,Alpha*(PP[i+1]-PP[i])/2./DdrL[i]);
							if(LIMITER_CONF>0)
								DmU[i]=minmod(Alpha*(UU[i]-UU[i-1])/2./DdrR[i],sU,Alpha*(UU[i+1]-UU[i])/2./DdrL[i]);
						}
					if (LIMITER_CONF>0)						
						TmV[i]=minmod2(Alpha*(UU[i]*sin(dtheta))/2./(0.5*(Rb[i]+Rb[i+1])*tan(0.5*dtheta)),sV);
				}
			// i = 0
			sU=(Umin[1]-UU[0]) /dRc[0];
			sP=(Pmin[1]-PP[0]) /dRc[0];
			sD=(DLmin[1]-DD[0])/dRc[0];
			//sV=0.;
			//sV=VLmin[1]/(0.5*Rb[1]*tan(0.5*dtheta));
			sV=UU[0]/Rb[1];
			Vave[0][0] = UU[1];
			Vave[0][1] = 0.;
			Vave[1][0] = UU[0]*cos(dtheta);
			Vave[1][1] = UU[0]*sin(dtheta);
			Vave[2][0] = UU[0]*cos(dtheta);
			Vave[2][1] =-UU[0]*sin(dtheta);
			V0[0] = UU[0];
			V0[1] = 0.;
			Vp1[0] = UU[0]+(0.5*Rb[1]-RR[0])*sU;
			Vp1[1] = 0.5*Rb[1]*tan(0.5*dtheta)*sV;
			Vp2[0] = UU[0]+DdrL[0]*sU;
			Vp2[1] = 0.;
			VIP_lim = fmin(1.0,     useVIPLimiter(3, Vave, V0, Vp1));
			VIP_lim = fmin(VIP_lim, useVIPLimiter(3, Vave, V0, Vp2));
			DmU[0]=VIP_lim*sU;
			TmV[0]=VIP_lim*sV;
			DmD[0]=minmod2(sD,DmD[1]);
			DmP[0]=minmod2(sP,DmP[1]);
			if (LIMITER_CONF>0)
				{									
					DmU[0]=minmod2(sU,DmU[1]);
					TmV[0]=minmod2(sV,TmV[1]);
				}
			// i = Ncell
			sU=(Umin[Ncell]-UU[Ncell-1]) /dRc[Ncell];
			sP=(Pmin[Ncell]-PP[Ncell-1]) /dRc[Ncell];
			sD=(DLmin[Ncell]-DD[Ncell-1])/dRc[Ncell];
			DmD[Ncell]=minmod2(sD,DmD[Ncell-1]);
			DmP[Ncell]=minmod2(sP,DmP[Ncell-1]);
			DmU[Ncell]=minmod2(sU,DmU[Ncell-1]);

			time=time+dt;
			printf("Time[%10d]=%e,dt=%e\n",k,time,dt);
			if (time>plot_t)
				{
					for(i=0;i<=Ncell;i++)
						for(j=0;j<=Tcell_plot;j++)
							{
								rb[i][j]=Rb[i]/cos(0.5*dtheta_plot)*sin(j*dtheta_plot);
								zb[i][j]=Rb[i]/cos(0.5*dtheta_plot)*cos(j*dtheta_plot);
							}
					for(i=1;i<=Ncell;i++)
						for(j=0;j<=Tcell_plot;j++)
							{
								DD2[i][j]=0.5*(DD[i-1]+DD[i]);
								UUxi2[i][j]=0.5*(UU[i-1]+UU[i]);
								PP2[i][j]=0.5*(PP[i-1]+PP[i]);
								GammaGamma2[i][j]=GammaGamma[i];
							}
					for(j=0;j<=Tcell_plot;j++)
						{
							DD2[0][j]=DD[0];
							UUxi2[0][j]=UU[0];
							PP2[0][j]=PP[0];
							GammaGamma2[0][j]=GammaGamma[0];
						}	
					sprintf(file_data, "%s/FLU_VAR_%.5g.dat", DATAOUT,plot_t);
					out=fopen(file_data,"w");
					wrin2s(out,rb,zb,DD2,UUxi2,PP2,GammaGamma2,time);	
					fclose(out);
					plot_t+=D_PLOT_T;
				}
			
			if(time-Timeout>EPS)
				break;
		}//end k

	outs=fopen("../datas_fin.m","w");
	fprintf(outs,"RR=[");
	Write(outs,RR,Ncell);
	fprintf(outs,"];\n");
	fprintf(outs,"DD=[");
	Write(outs,DD,Ncell);
	fprintf(outs,"];\n");
	fprintf(outs,"UU=[");
	Write(outs,UU,Ncell);
	fprintf(outs,"];\n");
	fprintf(outs,"PP=[");
	Write(outs,PP,Ncell);
	fprintf(outs,"];\n");
	fclose(outs);
	
	for(i=0;i<=Ncell;i++)
		for(j=0;j<=Tcell_plot;j++)
			{
				rb[i][j]=Rb[i]/cos(0.5*dtheta_plot)*sin(j*dtheta_plot);
				zb[i][j]=Rb[i]/cos(0.5*dtheta_plot)*cos(j*dtheta_plot);
			}
	for(i=1;i<=Ncell;i++)
		for(j=0;j<=Tcell_plot;j++)
			{
				DD2[i][j]=0.5*(DD[i-1]+DD[i]);
				UUxi2[i][j]=0.5*(UU[i-1]+UU[i]);
				PP2[i][j]=0.5*(PP[i-1]+PP[i]);
				GammaGamma2[i][j]=0.5*(GammaGamma[i-1]+GammaGamma[i]);
			}
	for(j=0;j<=Tcell_plot;j++)
		{
			DD2[0][j]=DD[0];
			UUxi2[0][j]=UU[0];
			PP2[0][j]=PP[0];
			GammaGamma2[0][j]=GammaGamma[0];
		}	
	sprintf(file_data, "%s/FLU_VAR_%.5g.dat", DATAOUT, time);
	out=fopen(file_data,"w");
	wrin2s(out,rb,zb,DD2,UUxi2,PP2,GammaGamma2,time);	
	fclose(out);
					
	return 1;
}
