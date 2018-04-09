#ifndef _INP_H
#define _INP_H
#define OK 1
#define ERROR 0
typedef int Status;
Status wrin2s(FILE *out,double R[Md][Mt],double Z[Md][Mt],double U[Md][Mt])//for mesh
{
	fprintf(out,"Solution For Euler Equation\n");
	fprintf(out,"variables=Z,R,U\n");
	fprintf(out,"zone I=%d,J=%d,F=POINT\n",Tcell+1,Ncell/3+1);
	int it,jr;
	for(jr=0;jr<=Ncell;jr++)
	{
	for(it=0;it<=(Tcell);it++)
	{
		if(jr%3==0){
			fprintf(out,"%lf %lf %lf\n",Z[jr][it],R[jr][it],U[jr][it]);}
	}
	}
	fclose(out);
return OK;
}
Status Write(FILE *out,double XX[Md],int N)// write data in file
{
int i;
for(i=1;i<=N;i++)
{
	if(i<N){fprintf(out,"%lf ",XX[i]);}
	else{fprintf(out,"%lf];\n",XX[i]);}
}
return OK;
}
double interp1(double x,double x1,double y1,double x2,double y2)//linear interpolate
{
double y;
y=(x-x2)*(y2-y1)/(x2-x1)+y2;
return(y);
}
double minmod2(double a,double b)// return minmod value
{
double s=0.;
if(a*b>=0)
{
if(a>=0.)
{ 
	s=(fabs(a)<fabs(b)?fabs(a):fabs(b));
	//s=(fabs(s)<fabs(c)?fabs(s):fabs(c));
}
else
{
s=-(fabs(a)<fabs(b)?fabs(a):fabs(b));
//s=-(fabs(s)<fabs(c)?fabs(s):fabs(c));
}
}
else
{s=0.;}

return(s);
}
double minmod(double a,double b,double c)// return minmod value
{
double s=0.;
if(a*b>=0&&a*c>=0.)
{
if(a>=0.)
{ 
	s=(fabs(a)<fabs(b)?fabs(a):fabs(b));
	s=(fabs(s)<fabs(c)?fabs(s):fabs(c));
}
else
{
s=-(fabs(a)<fabs(b)?fabs(a):fabs(b));
s=-(fabs(s)<fabs(c)?fabs(s):fabs(c));
}
}
else
{s=0.;}

return(s);
}
Status GuessP(double &PM,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double Gamma)
//provide a guess value for pressure PM in the star region.PVRS approximate Riemann solvers.
{
	double CUP,GEL,GER,PMAX,PMIN,PPV,PNU,PDE,QMAX,QUSER=2.0;
	CUP=0.25*(DL+DR)*(CL+CR);
	PPV=0.5*(PL+PR)+0.5*(UL-UR)*CUP;
	PPV=(1e-10>PPV?1e-10:PPV);
	PMIN=(PL<PR?PL:PR);
	PMAX=(PL>PR?PL:PR);
	QMAX=PMAX/PMIN;
	if(QMAX<=QUSER&&PMIN<=PPV&&PPV<=PMAX)
	{PM=PPV;
	//printf("Qmax=%lf,PPV=%lf\n",QMAX,PPV);
	}//PVRS solution
	else
	{
		if(PPV<PMIN){//Two rarefaction
			//PQ=pow(PL/PR,(Gamma-1.)/(2.*Gamma));
		//UM=(PQ*UL/CL+UR/CR+2./(Gamma-1.)*(PQ-1.))/(PQ/CL+1./CR);
		//PTL=1.+(Gamma-1.)/2.*(UL-UM)/CL;
		//PTR=1.+(Gamma-1.)/2.*(UM-UR)/CR;
		//	PM=0.5*(pow(PL*PTL,2.*Gamma/(Gamma-1.))+pow(PR*PTR,2.*Gamma/(Gamma-1.)));
			PNU=CL+CR-0.5*(Gamma-1.)*(UR-UL);
			PDE=CL/(pow(PL,0.5*(Gamma-1.)/Gamma))+CR/(pow(PR,0.5*(Gamma-1.)/Gamma));
		PM=pow(PNU/PDE,2.*Gamma/(Gamma-1.));
		//printf("PTR=%lf,PDE=%lf\n",PM,PDE);
		}
		else{//Two shock
			GEL=sqrt((2./((Gamma+1.)*DL))/((Gamma-1.)/(Gamma+1.)*PL+PPV));
		    GER=sqrt((2./((Gamma+1.)*DR))/((Gamma-1.)/(Gamma+1.)*PR+PPV));
			PM=(GEL*PL+GER*PR-(UR-UL))/(GEL+GER);
			PM=(1e-10>PM?1e-10:PM);
		//	printf("PTS=%lf\n",PM);
		}
	
	}

	return OK;
}


Status PreFun(double &F,double &FD,double P,double DK,double PK,double CK,double Gamma)
// evaluate the pressure functions FL and FR in exact Riemann solver
{
	double AK,BK,PRAT,QRT;
	if(P<=PK)//Rarefaction wave
	{
		PRAT=P/PK;
		F=2./(Gamma-1.)*CK*(pow(PRAT,(Gamma-1.)/(2.*Gamma))-1.);
		FD=1./(DK*CK)*pow(PRAT,-(Gamma+1.)/(2.*Gamma));
	
	}
	else// Shock wave
	{
		AK=2./((Gamma+1.)*DK);
		BK=PK*(Gamma-1.)/(Gamma+1.);
		QRT=sqrt(AK/(BK+P));
		F=(P-PK)*QRT;
		FD=(1.-0.5*(P-PK)/(BK+P))*QRT;
	}

	return OK;
}

Status StarPU(double &P,double &U,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double Gamma)
// compute the solution for pressure and velocity in the star region
{
	double change,FL,FR,FLD,FRD,POLD,PSTART,TOLPRE=1e-6,UDIFF;
	int i,NRITER=20;
	GuessP(PSTART,DL,DR,UL,UR,PL,PR,CL,CR,Gamma);
	POLD=PSTART;
    UDIFF=UR-UL;
	for(i=1;i<=NRITER;i++)
	{
		PreFun(FL,FLD,POLD,DL,PL,CL,Gamma);//GammaL
		PreFun(FR,FRD,POLD,DR,PR,CR,Gamma);//GammaR
		P=POLD-(FL+FR+UDIFF)/(FLD+FRD);
		change=2.*fabs((P-POLD)/(P+POLD));
		//printf("iteration number=%d,change=%lf\n",i,change);
		if(change<=TOLPRE)//compute velocity in star region
		{
			U=0.5*(UL+UR+FR-FL);
			//printf("iteration number=%d,Pressure=%lf,Velocity=%lf\n",i,P,U);
			break;
		}
		if(P<0)
		{
			P=TOLPRE;
		}
		POLD=P;
	
	
	}
	return OK;
}
Status Sample(double &D,double &U,double &P,double PM,double UM,double S,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double Gamma)
//Sampling is performed in terms of the 'speed' S=X/T. Sampled values are D,U,P. PM and UM in the star region are known.
{
	double C,CML,CMR,PML,PMR,SHL,SHR,SL,SR,STL,STR;
	if(S<=UM)//(4)
	{// Sampling points lies to the left of contact discontinuity
		if(PM<=PL)//(3)
		{// Left rarefaction
			SHL=UL-CL;
			//printf("SHL=%lf\n",SHL);
			if(S<=SHL)//(2)
			{
				D=DL;U=UL;P=PL;
			}
			else//(2)
			{
				CML=CL*pow(PM/PL,(Gamma-1.)/(2.*Gamma));
				STL=UM-CML;
				if(S>STL)//(1)
				{//Smapled point is star left state
					D=DL*pow(PM/PL,1./Gamma);
					U=UM;
					P=PM;
				}
				else//(1)
				{//Sampled point is inside left fan
					U=2./(Gamma+1.)*(CL+(Gamma-1.)/2.*UL+S);
					C=2./(Gamma+1.)*(CL+0.5*(Gamma-1.)*(UL-S));
					D=DL*pow(C/CL,2./(Gamma-1.));
					P=PL*pow(C/CL,2.*Gamma/(Gamma-1.));
				}
			}
		}
		else//(3)
		{//Left shock
			PML=PM/PL;
			SL=UL-CL*sqrt((Gamma+1.)/(2.*Gamma)*PML+(Gamma-1.)/(2.*Gamma));
		//	printf("SL=%lf\n",SL);
			if(S<=SL)//(3)1
			{//Sampled point is left data state
				D=DL;U=UL;P=PL;
			}
			else//(3)1
			{//Sampled point is star left state
				D=DL*(PML+(Gamma-1.)/(Gamma+1.))/(PML*(Gamma-1.)/(Gamma+1.)+1.);
				U=UM;P=PM;
			}
		}
	}
	else//(4)
	{//Sampling point lies to the right os the contact discontinuity
		if(PM>PR)//i
		{//Right Shock
			PMR=PM/PR;SR=UR+CR*sqrt((Gamma+1.)/(2.*Gamma)*PMR+(Gamma-1.)/(2.*Gamma));
			if(S>=SR)//ii
			{//Sampled point is right data state
				D=DR;U=UR;P=PR;
			}
			else//ii
			{//Sampled point is star right state
				D=DR*(PMR+(Gamma-1.)/(Gamma+1.))/((Gamma-1.)/(Gamma+1.)*PMR+1.);
				U=UM;P=PM;
			}
		}
		else//i
		{//Right rarefaction
			SHR=UR+CR;
			if(S>=SHR)//iii
			{//Sampled point is right data state
				D=DR;U=UR;P=PR;
			}
			else//iii
			{
				CMR=CR*pow(PM/PR,(Gamma-1.)/(2.*Gamma));STR=UM+CMR;
				if(S<=STR)//iv
				{//Sampled point is star right state
					D=DR*pow(PM/PR,1./Gamma);
					U=UM;P=PM;
				}
				else//iv
				{//sampled point is inside left fan
					U=2./(Gamma+1.)*(-CR+0.5*(Gamma-1.)*UR+S);
					C=2./(Gamma+1.)*(CR-0.5*(Gamma-1.)*(UR-S));
					D=DR*pow(C/CR,2./(Gamma-1.));
					P=PR*pow(C/CR,(2.*Gamma)/(Gamma-1.));
				}
			}
		}

	}
return OK;
}
Status SampleL(double &D,double &U,double &P,double PM,double UM,double S,double DL,double DR,double UL,double UR,double PL,double PR,double CL,double CR,double DML,double DMR,double Gamma)
//Sampling is performed in terms of the 'speed' S=X/T. Sampled values are D,U,P. PM and UM in the star region are known.
{
	double CML,CMR,PML,PMR,SHL,SHR,SL,SR,STL,STR,DM;
	DM=DMR;
	if(S<=0.)//(4)
	{// Sampling points lies to the left of contact discontinuity
		if(PM<=PL)//(3)
		{// Left rarefaction
			SHL=-CL*DL;
			//printf("SHL=%lf\n",SHL);
			if(S<=SHL)//(2)
			{
				D=DL;U=UL;P=PL;
			}
			else//(2)
			{
				CML=CL*pow(PM/PL,(Gamma-1.)/(2.*Gamma));
				STL=-CML*DML;//*pow(PM/PL,1./Gamma);
				if(S>STL)//(1)
				{//Smapled point is star left state
					D=DL*pow(PM/PL,1./Gamma);
					U=UM;
					P=PM;
				}
				else//(1)
				{//Sampled point is inside left fan
					P=PL*pow(S/(-DL*CL),2.*Gamma/(Gamma+1.));
					if(P<0){P=1e-6;}
					D=DL*pow(P/PL,1./Gamma);
					//U=2./(Gamma+1.)*(CL+(Gamma-1.)/2.*UL+S);
					U=UL+2./(Gamma-1.)*(CL+S/(D+1e-10));
				
				}
			}
		}
		else//(3)
		{//Left shock
			PML=PM/PL;
			SL=-CL*DL*sqrt((Gamma+1.)/(2.*Gamma)*PML+(Gamma-1.)/(2.*Gamma));
		//	printf("SL=%lf\n",SL);
			if(S<=SL)//(3)1
			{//Sampled point is left data state
				D=DL;U=UL;P=PL;
			}
			else//(3)1
			{//Sampled point is star left state
				D=DL*(PML+(Gamma-1.)/(Gamma+1.))/(PML*(Gamma-1.)/(Gamma+1.)+1.);
				U=UM;P=PM;
			}
		}
	}
	else//(4)
	{//Sampling point lies to the right os the contact discontinuity
		if(PM>PR)//i
		{//Right Shock
			PMR=PM/PR;SR=CR*DR*sqrt((Gamma+1.)/(2.*Gamma)*PMR+(Gamma-1.)/(2.*Gamma));
			if(S>=SR)//ii
			{//Sampled point is right data state
				D=DR;U=UR;P=PR;
			}
			else//ii
			{//Sampled point is star right state
				D=DR*(PMR+(Gamma-1.)/(Gamma+1.))/((Gamma-1.)/(Gamma+1.)*PMR+1.);
				U=UM;P=PM;
			}
		}
		else//i
		{//Right rarefaction
			SHR=CR*DR;
			if(S>=SHR)//iii
			{//Sampled point is right data state
				D=DR;U=UR;P=PR;
			}
			else//iii
			{
				CMR=CR*pow(PM/PR,(Gamma-1.)/(2.*Gamma));
				STR=CMR*DMR;//*pow(PM/PR,1./Gamma);
				if(S<=STR)//iv
				{//Sampled point is star right state
					D=DR*pow(PM/PR,1./Gamma);
					U=UM;P=PM;
				}
				else//iv
				{//sampled point is inside left fan
					P=PR*pow(S/(DR*CR),2.*Gamma/(Gamma+1.));
					if(P<0){P=1e-6;}
					D=DR*pow(P/PR,1./Gamma);
					U=UR+2./(Gamma-1.)*(-CR+S/(D+1e-10));
				}
			}
		}

	}
return OK;
}
Status Acoustic(double &Dttau,double &DtU,double &DtP,//double DL,double DR,double UL,double UR,double PL,double PR,
				double DUL,double DUR,double DPL,double DPR,double tau_star,double C_star)
// Acoustic case for U_star=UL=UR, caculate for the time derivative (dU/dt)_star, GRP slover,output:Dttau,Dtu,Dtp;Lagrangian version
{
	//if((DL-DR)*(DL-DR)+(UL-UR)*(UL-UR)+(PL-PR)*(PL-PR)<1e-12)//Acoustic case
	//{
	DtU=-0.5*(DPL+DPR+(DUL-DUR)*C_star/tau_star);
	DtP=-0.5*C_star/tau_star*(DPL-DPR+(DUL+DUR)*C_star/tau_star);
	Dttau=0.5*tau_star/C_star*(DPL-DPR+(DUL+DUR)*C_star/tau_star);
	//Dttau=-tau_star*tau_star/(C_star*C_star)*DtP;
	// }

	return OK;
}
Status AcousticE(double &DtD,double &DtU,double &DtP,double DDL,double DDR,double DUL,double DUR,double DPL,
				 double DPR,double D,double U,double C_star,double S)//Euler version(moving mesh)
{
	DtP=-0.5*D*C_star*((U-S+C_star)*(DUL+DPL/(D*C_star))-(U-S-C_star)*(DUR-DPR/(D*C_star)));
	DtU=-0.5*((U-S+C_star)*(DUL+DPL/(D*C_star))+(U-S-C_star)*(DUR-DPR/(D*C_star)));
	if(U>=0.)
	{
	DtD=(-0.5*D*C_star*((U-S+C_star)*(DUL+DPL/(D*C_star))-(U-S-C_star)*(DUR-DPR/(D*C_star)))+(U-S)*(DPL-C_star*C_star*DDL))/(C_star*C_star);
	}
	else
	{
    DtD=(-0.5*D*C_star*((U-S+C_star)*(DUL+DPL/(D*C_star))-(U-S-C_star)*(DUR-DPR/(D*C_star)))+(U-S)*(DPR-C_star*C_star*DDR))/(C_star*C_star);
	}
	return OK;
}
Status AcousticS(double &DtD,double &DtU,double &DtP,double DDL,double DDR,double DUL,double DUR,double DPL,
				 double DPR,double D,double U,double C_star,double S,double r)//Euler version(moving mesh) for spherical symmetry(r!=0)
{
	DtP=-0.5*D*C_star*((U-S+C_star)*(DUL+DPL/(D*C_star))-(U-S-C_star)*(DUR-DPR/(D*C_star)))-(m-1)/r*D*(U-S)*C_star*C_star;
	DtU=-0.5*((U-S+C_star)*(DUL+DPL/(D*C_star))+(U-S-C_star)*(DUR-DPR/(D*C_star)));
	if(U>=0.)
	{
	DtD=(-0.5*D*C_star*((U-S+C_star)*(DUL+DPL/(D*C_star))-(U-S-C_star)*(DUR-DPR/(D*C_star)))-(m-1)/r*D*(U-S)*C_star*C_star+(U-S)*(DPL-C_star*C_star*DDL))/(C_star*C_star);
	}
	else
	{
    DtD=(-0.5*D*C_star*((U-S+C_star)*(DUL+DPL/(D*C_star))-(U-S-C_star)*(DUR-DPR/(D*C_star)))-(m-1)/r*D*(U-S)*C_star*C_star+(U-S)*(DPR-C_star*C_star*DDR))/(C_star*C_star);
	}
	return OK;
}
Status GRPsolver(double &Dttau,double &DtU,double &DtP,double UM,double PM,double DL,double DR,double UL,double UR,double PL,
double PR,double DtauL,double DtauR,double DUL,double DUR,double DPL,double DPR,double TDSL,double TDSR,double DpsiL,double DphiR,double Gamma)
//GRP solver for non acoustic case.Lagrangian version
{
	double aL,bL,dL,aR,bR,dR,CL,CR,C_star,D,U,P,theta,SR,SL,phi1,phi2,phi3,sigmaL,sigmaR,mus=(Gamma-1.)/(Gamma+1.);
	CL=sqrt(Gamma*PL/DL);CR=sqrt(Gamma*PR/DR);
	if(PM<=PL)//Left rarefaction wave
	{
// middle left state
		D=DL*pow(PM/PL,1./Gamma);U=UM;P=PM;
		C_star=sqrt(Gamma*P/D);
		aL=1.;bL=1./(D*C_star);theta=C_star/CL;
		dL=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSL*DL
			-pow(theta,0.5/mus)*CL*DL*DpsiL;
        if(PM<=PR)//Right rarefaction wave
		{
			//middle right state
			D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			aR=1.;bR=-1./(D*C_star);theta=C_star/CR;
            dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR*DR
			+pow(theta,0.5/mus)*CR*DR*DphiR;
            DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
            Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
    //        printf("RR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
		else//Right shock wave
		{
		//middle Right state (behind shock)
          SR=DR*CR*sqrt(PM/PR*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DR*(PM/PR+(Gamma-1.)/(Gamma+1.))/(PM/PR*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM+PR*(1.+2.*mus))/(PM+mus*PR);
		  phi2=-0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM*(2.+mus)+mus*PR)/(PM+mus*PR);
          phi3=0.5*(PM-PR)*DR*sqrt((1.-mus)/(DR*(PM+mus*PR)));
		  sigmaR=-(UM-UR)/(1./D-1./DR);
          aR=1.+sigmaR*phi1;
		  bR=-phi1-sigmaR/(D*D*C_star*C_star);
		  dR=(sigmaR*phi2-1.)*DPR+(sigmaR+phi3-CR*CR*DR*DR*phi2)*DUR+sigmaR*phi3*DtauR;
		  DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
          Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
	//	  printf("RS\n");
		}//end Right shock wave
	}//end Left rarefaction wave
	else//Left shock wave
	{
		//middle Left state(behind shock)
		  SL=-CL*DL*(PM/PL*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DL*(PM/PL+(Gamma-1.)/(Gamma+1.))/(PM/PL*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DL*(PM+mus*PL)))*(PM+PL*(1.+2.*mus))/(PM+mus*PL);
		  phi2=-0.5*sqrt((1.-mus)/(DL*(PM+mus*PL)))*(PM*(2.+mus)+mus*PL)/(PM+mus*PL);
          phi3=0.5*(PM-PL)*DL*sqrt((1.-mus)/(DL*(PM+mus*PL)));
		  sigmaL=-(UM-UL)/(1./D-1./DL);
          aL=1.-sigmaL*phi1;
		  bL=phi1-sigmaL/(D*D*C_star*C_star);
		  dL=(-sigmaL*phi2-1.)*DPL+(sigmaL-phi3+CL*CL*DL*DL*phi2)*DUL-sigmaL*phi3*DtauL;
		  if(PM<=PR)//Right rarefaction wave
		{
			//middle right state
			D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			aR=1.;bR=-1./(D*C_star);theta=C_star/CR;
            dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR*DR
			+pow(theta,0.5/mus)*CR*DR*DphiR;
            DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
            Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
 //           printf("SR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
		else//Right shock wave
		{
		//middle Right state (behind shock)
          SR=DR*CR*sqrt(PM/PR*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DR*(PM/PR+(Gamma-1.)/(Gamma+1.))/(PM/PR*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM+PR*(1.+2.*mus))/(PM+mus*PR);
		  phi2=-0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM*(2.+mus)+mus*PR)/(PM+mus*PR);
          phi3=0.5*(PM-PR)*DR*sqrt((1.-mus)/(DR*(PM+mus*PR)));
		  sigmaR=-(UM-UR)/(1./D-1./DR);
          aR=1.+sigmaR*phi1;
		  bR=-phi1-sigmaR/(D*D*C_star*C_star);
		  dR=(sigmaR*phi2-1.)*DPR+(sigmaR+phi3-CR*CR*DR*DR*phi2)*DUR+sigmaR*phi3*DtauR;
		  DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
          Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
//		  printf("SS\n");
		}//end Right shock wave 
	}//end Left shock wave
	

	return OK;
}
Status GRPsolverE(double &DtD,double &DtU,double &DtP,double UM,double PM,double DL,double DR,double UL,double UR,double PL,
double PR,double DDL,double DDR,double DUL,double DUR,double DPL,double DPR,double TDSL,double TDSR,double DpsiL,double DphiR,double S,double Gamma)
//GRP solver for non acoustic case.Euler version(moving mesh)
{
	double aL,bL,dL,aR,bR,dR,CL,CR,C_star,D,U,P,theta,SR,SL,phi1,phi2,phi3,sigmaL,sigmaR,mus=(Gamma-1.)/(Gamma+1.);
	CL=sqrt(Gamma*PL/DL);CR=sqrt(Gamma*PR/DR);
double gD,gU,gP,H1,H2,H3,f;
	if(PM<=PL)//Left rarefaction wave
	{
// middle left state
		D=DL*pow(PM/PL,1./Gamma);U=UM;P=PM;
		C_star=sqrt(Gamma*P/D);theta=C_star/CL;
		aL=1.;bL=1./(D*C_star);
		dL=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSL
			-pow(theta,0.5/mus)*CL*DpsiL;
        if(PM<=PR)//Right rarefaction wave
		{
			//middle right state
			D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			aR=1.;bR=-1./(D*C_star);theta=C_star/CR;
            dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR
			+pow(theta,0.5/mus)*CR*DphiR;
            DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
            //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
			if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*((dL*aR-dR*aL)/(bL*aR-bR*aL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=(UM-S)-sigmaR;gU=(UM-S)*D*(sigmaR-(UM-S))*H1;gP=sigmaR/(C_star*C_star)-(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*DUR;
			DtD=((UM-S)*f-gU*(dL*bR-dR*bL)/(aL*bR-aR*bL)-gP*(dL*aR-dR*aL)/(bL*aR-bR*aL))/gD;
			}
      //      printf("RR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
		else//Right shock wave
		{
		//middle Right state (behind shock)
          SR=DR*CR*sqrt(PM/PR*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DR*(PM/PR+(Gamma-1.)/(Gamma+1.))/(PM/PR*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM+PR*(1.+2.*mus))/(PM+mus*PR);
		  phi2=-0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM*(2.+mus)+mus*PR)/(PM+mus*PR);
          phi3=-0.5*(PM-PR)/DR*sqrt((1.-mus)/(DR*(PM+mus*PR)));
		  sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
          aR=1.+D*(sigmaR-(UM-S))*phi1;
		  bR=-phi1-(sigmaR-(UM-S))/(D*C_star*C_star);
		  dR=((sigmaR-(UR-S))*phi2-1./DR)*DPR+(sigmaR-(UR-S)-DR*phi3-CR*CR*DR*phi2)*DUR+(sigmaR-(UR-S))*phi3*DDR;
		  DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
          //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
          if(UM>=0.)
			{
			  theta=C_star/CL;
			DtD=1./(C_star*C_star)*((dL*aR-dR*aL)/(bL*aR-bR*aL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=(UM-S)-sigmaR;gU=(UM-S)*D*(sigmaR-(UM-S))*H1;gP=sigmaR/(C_star*C_star)-(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*DUR;
			DtD=((UM-S)*f-gU*(dL*bR-dR*bL)/(aL*bR-aR*bL)-gP*(dL*aR-dR*aL)/(bL*aR-bR*aL))/gD;
			}
	//	  printf("RS\n");
		}//end Right shock wave
	}//end Left rarefaction wave
	else//Left shock wave
	{
		//middle Left state(behind shock)
		  SL=-CL*DL*(PM/PL*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DL*(PM/PL+(Gamma-1.)/(Gamma+1.))/(PM/PL*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DL*(PM+mus*PL)))*(PM+PL*(1.+2.*mus))/(PM+mus*PL);
		  phi2=-0.5*sqrt((1.-mus)/(DL*(PM+mus*PL)))*(PM*(2.+mus)+mus*PL)/(PM+mus*PL);
          phi3=-0.5*(PM-PL)/DL*sqrt((1.-mus)/(DL*(PM+mus*PL)));
		  sigmaL=(D*(U-S)-DL*(UL-S))/(D-DL);
          aL=1.-D*(sigmaL-(UM-S))*phi1;
		  bL=phi1-(sigmaL-(UM-S))/(D*C_star*C_star);
		  dL=(-(sigmaL-(UL-S))*phi2-1./DL)*DPL+(sigmaL-(UL-S)+DL*phi3+CL*CL*DL*phi2)*DUL-(sigmaL-(UL-S))*phi3*DDL;
		  if(PM<=PR)//Right rarefaction wave
		{
			//middle right state
			D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			aR=1.;bR=-1./(D*C_star);theta=C_star/CR;
            dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR
			+pow(theta,0.5/mus)*CR*DphiR;
            DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
            //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
            if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*((dL*aR-dR*aL)/(bL*aR-bR*aL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=(UM-S)-sigmaR;gU=(UM-S)*D*(sigmaR-(UM-S))*H1;gP=sigmaR/(C_star*C_star)-(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*DUR;
			DtD=((UM-S)*f-gU*(dL*bR-dR*bL)/(aL*bR-aR*bL)-gP*(dL*aR-dR*aL)/(bL*aR-bR*aL))/gD;
			}
        //    printf("SR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
		else//Right shock wave
		{
		//middle Right state (behind shock)
          SR=DR*CR*sqrt(PM/PR*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DR*(PM/PR+(Gamma-1.)/(Gamma+1.))/(PM/PR*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM+PR*(1.+2.*mus))/(PM+mus*PR);
		  phi2=-0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM*(2.+mus)+mus*PR)/(PM+mus*PR);
          phi3=-0.5*(PM-PR)/DR*sqrt((1.-mus)/(DR*(PM+mus*PR)));
		  sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
          aR=1.+D*(sigmaR-(UM-S))*phi1;
		  bR=-phi1-(sigmaR-(UM-S))/(D*C_star*C_star);
		  dR=((sigmaR-(UR-S))*phi2-1./DR)*DPR+(sigmaR-(UR-S)-DR*phi3-CR*CR*DR*phi2)*DUR+(sigmaR-(UR-S))*phi3*DDR;
		  DtU=(dL*bR-dR*bL)/(aL*bR-aR*bL);DtP=(dL*aR-dR*aL)/(bL*aR-bR*aL);
          //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
		  if(UM>=0.)
			{
			  theta=C_star/CL;
			DtD=1./(C_star*C_star)*((dL*aR-dR*aL)/(bL*aR-bR*aL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=(UM-S)-sigmaR;gU=(UM-S)*D*(sigmaR-(UM-S))*H1;gP=sigmaR/(C_star*C_star)-(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*DUR;
			DtD=((UM-S)*f-gU*(dL*bR-dR*bL)/(aL*bR-aR*bL)-gP*(dL*aR-dR*aL)/(bL*aR-bR*aL))/gD;
			}
		//  printf("SS\n");
		}//end Right shock wave 
	}//end Left shock wave
	return OK;
}
Status GRPsolverS(double &DtD,double &DtU,double &DtP,double UM,double PM,double DL,double DR,double UL,double UR,double PL,
double PR,double DDR,double DUR,double DPR,double TDSL,double TDSR,double DpsiL,double DphiR,double S,double Gamma)
//sonic case
{
		double dL,dR,C,CL,CR,C_star,C_starL,C_starR,D,U,P,theta,sigmaR,mus;mus=(Gamma-1.)/(Gamma+1.);
	CL=sqrt(Gamma*PL/DL);CR=sqrt(Gamma*PR/DR);C_starL=CL*pow(PM/PL,(Gamma-1.)/(2.*Gamma));C_starR=CR*pow(PM/PR,(Gamma-1.)/(2.*Gamma));
double gD,gU,gP,H1,H2,H3,f,SHL,STL,SHR,STR;SHL=UL-CL;STL=UM-C_starL;SHR=UR+CR;STR=UM+C_starR;
	if(S>SHL&&S<=STL)//Left rarefaction wave
	{
// inside left fan
		U=2./(Gamma-1.)*(CL+0.5*(Gamma-1.)*UL+S);
		C=2./(Gamma-1.)*(CL+0.5*(Gamma-1.)*(UL-S));
		D=DL*pow(C/CL,2./(Gamma-1.));
		P=PL*pow(C/CL,2.*Gamma/(Gamma-1.));
		//D=DL*pow(PM/PL,1./Gamma);U=UM;P=PM;
		C_star=sqrt(Gamma*P/D);theta=C_star/CL;
		//aL=1.;bL=1./(D*C_star);
		dL=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSL
			-pow(theta,0.5/mus)*CL*DpsiL;
		DtU=dL;DtP=D*(UM-S)*dL;
		if(UM>=0.)
			{
			theta=C_star/CL;
			DtD=1./(C_star*C_star)*(D*(UM-S)*dL+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=(UM-S)-sigmaR;gU=(UM-S)*D*(sigmaR-(UM-S))*H1;gP=sigmaR/(C_star*C_star)-(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*DUR;
			DtD=((UM-S)*f-gU*dL-gP*D*(UM-S)*dL)/gD;
			}

	}//end Left rarefaction
        if(S>STR&&S<SHR)//Right rarefaction wave
		{
			//inside right fan
        U=2./(Gamma-1.)*(-CR+0.5*(Gamma-1.)*UL+S);
		C=2./(Gamma-1.)*(CR-0.5*(Gamma-1.)*(UL-S));
		D=DL*pow(C/CR,2./(Gamma-1.));
		P=PL*pow(C/CR,2.*Gamma/(Gamma-1.));
			//D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			//aR=1.;bR=-1./(D*C_star);
			theta=C_star/CR;
            dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR
			+pow(theta,0.5/mus)*CR*DphiR;
         	DtU=dR;DtP=D*(UM-S)*dR;
			if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*(D*(UM-S)*dR+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=(UM-S)-sigmaR;gU=(UM-S)*D*(sigmaR-(UM-S))*H1;gP=sigmaR/(C_star*C_star)-(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*DUR;
			DtD=((UM-S)*f-gU*dR-gP*D*(UM-S)*dR)/gD;
			}
    //        printf("RR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
	
return OK;
}
Status GRPsolverES(double &DtD,double &DtU,double &DtP,double UM,double PM,double DL,double DR,double UL,double UR,double PL,
double PR,double DDL,double DDR,double DUL,double DUR,double DPL,double DPR,double TDSL,double TDSR,double DpsiL,double DphiR,double S,double r,double Gamma)
//GRP solver for non acoustic case.Euler version(moving mesh) spherical case
{
	double aL,bL,dL,aR,bR,dR,CL,CR,C_star,D,U,P,theta,SR,SL,phi1,phi2,phi3,sigmaL,sigmaR,mus=(Gamma-1.)/(Gamma+1.);
	CL=sqrt(Gamma*PL/DL);CR=sqrt(Gamma*PR/DR);
double gD,gU,gP,H1,H2,H3,f,phic,hL,qL,kL,hR,qR,kR;
	if(PM<=PL)//Left rarefaction wave
	{
// middle left state
		D=DL*pow(PM/PL,1./Gamma);U=UM;P=PM;
		C_star=sqrt(Gamma*P/D);theta=C_star/CL;
		aL=1.;bL=1./(D*C_star);
		if(Gamma==(5./3.))
		{
		phic=-2.*(3.*C_star*log(theta)-(-1.)*(UL+2.*CL/(Gamma-1.))*(1.-theta));
		}
		else
		{
		phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))+(UL+2.*CL/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		}
		dL=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSL
			-pow(theta,0.5/mus)*CL*(DpsiL+(m-1.)/(2.*r)*(UL-S))+(m-1.)/(2.*r)*C_star*(phic-(U-S));
		hL=aL-D*(U-S)*bL;qL=bL-(U-S)/(D*C_star*C_star)*aL;
		kL=(m-1.)/r*(U-S)*(U-S)*hL+(C_star*C_star-(U-S)*(U-S))/(C_star*C_star)*dL;
        if(PM<=PR)//Right rarefaction wave
		{
			//middle right state
			D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			aR=1.;bR=-1./(D*C_star);theta=C_star/CR;
			if(Gamma==(5./3.))
		{
		phic=-2.*(3.*C_star*log(theta)-(1.)*(UR-2.*CR/(Gamma-1.))*(1.-theta));
		}
		else
		{
		//phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))+(UL+2.*CL/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))-(UR-2.*CR/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		}
			//phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))-(UR-2.*CR/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
            dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR
			+pow(theta,0.5/mus)*CR*(DphiR+(m-1.)/(2.*r)*(UR-S))+(m-1.)/(2.*r)*C_star*(phic+(U-S));
			hR=aR-D*(U-S)*bR;qR=bR-(U-S)/(D*C_star*C_star)*aR;
		    kR=(m-1.)/r*(U-S)*(U-S)*hR+(C_star*C_star-(U-S)*(U-S))/(C_star*C_star)*dR;
            DtU=(kL*qR-kR*qL)/(hL*qR-hR*qL);DtP=(kL*hR-kR*hL)/(qL*hR-qR*hL);
            //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
			if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*((kL*hR-kR*hL)/(qL*hR-qR*hL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=((UM-S)-sigmaR)*(C_star*C_star-(UM-S)*(UM-S));gU=(UM-S)*D*(sigmaR)*(H1*C_star*C_star-1.);
			gP=sigmaR-((UM-S)*(sigmaR-(UM-S))+C_star*C_star)*(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*(DUR+UR*(m-1.)/r);
			f=(UM-S)*(C_star*C_star-(UM-S)*(UM-S))*f+(UM-S)*(UM-S)*gU*(m-1.)/r;
			DtD=(f-gU*(kL*qR-kR*qL)/(hL*qR-hR*qL)-gP*(kL*hR-kR*hL)/(qL*hR-qR*hL))/gD;
			}
      //      printf("RR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
		else//Right shock wave
		{
		//middle Right state (behind shock)
          SR=DR*CR*sqrt(PM/PR*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DR*(PM/PR+(Gamma-1.)/(Gamma+1.))/(PM/PR*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM+PR*(1.+2.*mus))/(PM+mus*PR);
		  phi2=-0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM*(2.+mus)+mus*PR)/(PM+mus*PR);
          phi3=-0.5*(PM-PR)/DR*sqrt((1.-mus)/(DR*(PM+mus*PR)));
		  sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
          aR=1.+D*(sigmaR-(UM-S))*phi1;
		  bR=-phi1-(sigmaR-(UM-S))/(D*C_star*C_star);
		  dR=((sigmaR-(UR-S))*phi2-1./DR)*DPR+(sigmaR-(UR-S)-DR*phi3-CR*CR*DR*phi2)*DUR+(sigmaR-(UR-S))*phi3*DDR-(DR*(UR-S)*CR*CR*phi2+DR*(UR-S)*phi3
			  -(sigmaR-(UM-S))*(UM-S))*(m-1.)/r;
          hR=aR-D*(U-S)*bR;qR=bR-(U-S)/(D*C_star*C_star)*aR;
		  kR=(m-1.)/r*(U-S)*(U-S)*hR+(C_star*C_star-(U-S)*(U-S))/(C_star*C_star)*dR;
		  DtU=(kL*qR-kR*qL)/(hL*qR-hR*qL);DtP=(kL*hR-kR*hL)/(qL*hR-qR*hL);
          //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
          if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*((kL*hR-kR*hL)/(qL*hR-qR*hL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=((UM-S)-sigmaR)*(C_star*C_star-(UM-S)*(UM-S));gU=(UM-S)*D*(sigmaR)*(H1*C_star*C_star-1.);
			gP=sigmaR-((UM-S)*(sigmaR-(UM-S))+C_star*C_star)*(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*(DUR+UR*(m-1.)/r);
			f=(UM-S)*(C_star*C_star-(UM-S)*(UM-S))*f+(UM-S)*(UM-S)*gU*(m-1.)/r;
			DtD=(f-gU*(kL*qR-kR*qL)/(hL*qR-hR*qL)-gP*(kL*hR-kR*hL)/(qL*hR-qR*hL))/gD;
			}
	//	  printf("RS\n");
		}//end Right shock wave
	}//end Left rarefaction wave
	else//Left shock wave
	{
		//middle Left state(behind shock)
		  SL=-CL*DL*(PM/PL*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DL*(PM/PL+(Gamma-1.)/(Gamma+1.))/(PM/PL*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DL*(PM+mus*PL)))*(PM+PL*(1.+2.*mus))/(PM+mus*PL);
		  phi2=-0.5*sqrt((1.-mus)/(DL*(PM+mus*PL)))*(PM*(2.+mus)+mus*PL)/(PM+mus*PL);
          phi3=-0.5*(PM-PL)/DL*sqrt((1.-mus)/(DL*(PM+mus*PL)));
		  sigmaL=(D*(U-S)-DL*(UL-S))/(D-DL);
          aL=1.-D*(sigmaL-(UM-S))*phi1;
		  bL=phi1-(sigmaL-(UM-S))/(D*C_star*C_star);
		  dL=(-(sigmaL-(UL-S))*phi2-1./DL)*DPL+(sigmaL-(UL-S)+DL*phi3+CL*CL*DL*phi2)*DUL-(sigmaL-(UL-S))*phi3*DDL+(DL*(UL-S)*CL*CL*phi2+DL*(UL-S)*phi3
			  -(sigmaL-(UM-S))*(UM-S))*(m-1.)/r;
	    hL=aL-D*(U-S)*bL;qL=bL-(U-S)/(D*C_star*C_star)*aL;
		kL=(m-1.)/r*(U-S)*(U-S)*hL+(C_star*C_star-(U-S)*(U-S))/(C_star*C_star)*dL;
		  if(PM<=PR)//Right rarefaction wave
		{
			//middle right state
			D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			aR=1.;bR=-1./(D*C_star);theta=C_star/CR;
            //dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR
			//+pow(theta,0.5/mus)*CR*DphiR;
if(Gamma==(5./3.))
		{
		phic=-2.*(3.*C_star*log(theta)-(1.)*(UR-2.*CR/(Gamma-1.))*(1.-theta));
		}
		else
		{
		//phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))+(UL+2.*CL/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))-(UR-2.*CR/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		}
		//	phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))-(UR-2.*CR/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSL
			+pow(theta,0.5/mus)*CR*(DphiR+(m-1.)/(2.*r)*(UR-S))+(m-1.)/(2.*r)*C_star*(phic+(U-S));
		hR=aR-D*(U-S)*bR;qR=bR-(U-S)/(D*C_star*C_star)*aR;
		kR=(m-1.)/r*(U-S)*(U-S)*hR+(C_star*C_star-(U-S)*(U-S))/(C_star*C_star)*dR;
         DtU=(kL*qR-kR*qL)/(hL*qR-hR*qL);DtP=(kL*hR-kR*hL)/(qL*hR-qR*hL);
            //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
          if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*((kL*hR-kR*hL)/(qL*hR-qR*hL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=((UM-S)-sigmaR)*(C_star*C_star-(UM-S)*(UM-S));gU=(UM-S)*D*(sigmaR)*(H1*C_star*C_star-1.);
			gP=sigmaR-((UM-S)*(sigmaR-(UM-S))+C_star*C_star)*(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*(DUR+UR*(m-1.)/r);
			f=(UM-S)*(C_star*C_star-(UM-S)*(UM-S))*f+(UM-S)*(UM-S)*gU*(m-1.)/r;
			DtD=(f-gU*(kL*qR-kR*qL)/(hL*qR-hR*qL)-gP*(kL*hR-kR*hL)/(qL*hR-qR*hL))/gD;
			}  
        //    printf("SR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
		else//Right shock wave
		{
		//middle Right state (behind shock)
          SR=DR*CR*sqrt(PM/PR*(Gamma+1.)/(2.*Gamma)+(Gamma-1.)/(2.*Gamma));
		  D=DR*(PM/PR+(Gamma-1.)/(Gamma+1.))/(PM/PR*(Gamma-1.)/(Gamma+1.)+1.);
		  U=UM;P=PM;C_star=sqrt(Gamma*P/D);
		  phi1=0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM+PR*(1.+2.*mus))/(PM+mus*PR);
		  phi2=-0.5*sqrt((1.-mus)/(DR*(PM+mus*PR)))*(PM*(2.+mus)+mus*PR)/(PM+mus*PR);
          phi3=-0.5*(PM-PR)/DR*sqrt((1.-mus)/(DR*(PM+mus*PR)));
		  sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
		  aR=1.+D*(sigmaR-(UM-S))*phi1;
		  bR=-phi1-(sigmaR-(UM-S))/(D*C_star*C_star);
          dR=((sigmaR-(UR-S))*phi2-1./DR)*DPR+(sigmaR-(UR-S)-DR*phi3-CR*CR*DR*phi2)*DUR+(sigmaR-(UR-S))*phi3*DDR-(DR*(UR-S)*CR*CR*phi2+DR*(UR-S)*phi3
			  -(sigmaR-(UM-S))*(UM-S))*(m-1.)/r;
          hR=aR-D*(U-S)*bR;qR=bR-(U-S)/(D*C_star*C_star)*aR;
		  kR=(m-1.)/r*(U-S)*(U-S)*hR+(C_star*C_star-(U-S)*(U-S))/(C_star*C_star)*dR;
		  DtU=(kL*qR-kR*qL)/(hL*qR-hR*qL);DtP=(kL*hR-kR*hL)/(qL*hR-qR*hL);
          //Dttau=-1./(D*D*C_star*C_star)*(dL*aR-dR*aL)/(bL*aR-bR*aL);
          if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*((kL*hR-kR*hL)/(qL*hR-qR*hL)+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=((UM-S)-sigmaR)*(C_star*C_star-(UM-S)*(UM-S));gU=(UM-S)*D*(sigmaR)*(H1*C_star*C_star-1.);
			gP=sigmaR-((UM-S)*(sigmaR-(UM-S))+C_star*C_star)*(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*(DUR+UR*(m-1.)/r);
			f=(UM-S)*(C_star*C_star-(UM-S)*(UM-S))*f+(UM-S)*(UM-S)*gU*(m-1.)/r;
			DtD=(f-gU*(kL*qR-kR*qL)/(hL*qR-hR*qL)-gP*(kL*hR-kR*hL)/(qL*hR-qR*hL))/gD;
			}
		//  printf("SS\n");
		}//end Right shock wave 
	}//end Left shock wave
	return OK;
}
Status GRPsolverSS(double &DtD,double &DtU,double &DtP,double UM,double PM,double DL,double DR,double UL,double UR,double PL,
double PR,double DDR,double DUR,double DPR,double TDSL,double TDSR,double DpsiL,double DphiR,double S,double r,double Gamma)
//sonic case spherical case
{
		double dL,dR,C,CL,CR,C_star,C_starL,C_starR,D,U,P,theta,sigmaR,mus;mus=(Gamma-1.)/(Gamma+1.);
	CL=sqrt(Gamma*PL/DL);CR=sqrt(Gamma*PR/DR);C_starL=CL*pow(PM/PL,(Gamma-1.)/(2.*Gamma));C_starR=CR*pow(PM/PR,(Gamma-1.)/(2.*Gamma));
double gD,gU,gP,H1,H2,H3,f,SHL,STL,SHR,STR,phic;SHL=UL-CL;STL=UM-C_starL;SHR=UR+CR;STR=UM+C_starR;
	if(S>SHL&&S<=STL)//Left rarefaction wave
	{
// inside left fan
		U=2./(Gamma-1.)*(CL+0.5*(Gamma-1.)*UL+S);
		C=2./(Gamma-1.)*(CL+0.5*(Gamma-1.)*(UL-S));
		D=DL*pow(C/CL,2./(Gamma-1.));
		P=PL*pow(C/CL,2.*Gamma/(Gamma-1.));
		D=DL*pow(PM/PL,1./Gamma);U=UM;P=PM;
		C_star=sqrt(Gamma*P/D);theta=C_star/CL;
		if(Gamma==(5./3.))
		{
		phic=-2.*(3.*C_star*log(theta)-(-1.)*(UL+2.*CL/(Gamma-1.))*(1.-theta));
		}
		else
		{
		phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))+(UL+2.*CL/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		//phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))-(UR-2.*CR/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		}
		//phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))+(UL+2.*CL/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		//aL=1.;bL=1./(D*C_star);
		//dL=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSL
		//	-pow(theta,0.5/mus)*CL*DpsiL;
		dL=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSL
			-pow(theta,0.5/mus)*CL*(DpsiL+(m-1.)/(2.*r)*(UL-S))+(m-1.)/(2.*r)*C_star*(phic-(U-S));
		DtU=dL+(UM-S)*(UM-S)*(m-1.)/r;DtP=D*(UM-S)*dL;
	    if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*(D*(UM-S)*dL+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
			else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=((UM-S)-sigmaR)*(C_star*C_star-(UM-S)*(UM-S));gU=(UM-S)*D*(sigmaR)*(H1*C_star*C_star-1.);
			gP=sigmaR-((UM-S)*(sigmaR-(UM-S))+C_star*C_star)*(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*(DUR+UR*(m-1.)/r);
			f=(UM-S)*(C_star*C_star-(UM-S)*(UM-S))*f+(UM-S)*(UM-S)*gU*(m-1.)/r;
			DtD=(f-gU*(dL+(UM-S)*(UM-S)*(m-1.)/r)-gP*D*(UM-S)*dL)/gD;
			}
			//printf("dL=%lf,D=%lf,U=%lf,DL=%lf,UM=%lf,LR\n",dL,D,U,DL,UM);

	}//end Left rarefaction
        if(S>STR&&S<SHR)//Right rarefaction wave
		{
			//inside right fan
        U=2./(Gamma-1.)*(-CR+0.5*(Gamma-1.)*UL+S);
		C=2./(Gamma-1.)*(CR-0.5*(Gamma-1.)*(UL-S));
		D=DL*pow(C/CR,2./(Gamma-1.));
		P=PL*pow(C/CR,2.*Gamma/(Gamma-1.));
			D=DR*pow(PM/PR,1./Gamma);U=UM;P=PM;
			C_star=sqrt(Gamma*P/D);
			//aR=1.;bR=-1./(D*C_star);
			theta=C_star/CR;
			if(Gamma==(5./3.))
		{
		phic=-2.*(3.*C_star*log(theta)-(1.)*(UR-2.*CR/(Gamma-1.))*(1.-theta));
		}
		else
		{
		//phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))+(UL+2.*CL/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))-(UR-2.*CR/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
		}
			//phic=(mus-1.)*C_star/(mus*(4.*mus-1.))*(1.-pow(theta,(1.-4.*mus)/(2.*mus)))-(UR-2.*CR/(Gamma-1.))/(2.*mus-1.)*(1.-pow(theta,(1.-2*mus)/(2.*mus)));
            dR=((1.+mus)/(1.+2.*mus)*pow(theta,0.5/mus)+mus/(1.+2.*mus)*pow(theta,(1.+mus)/mus))*TDSR
			+pow(theta,0.5/mus)*CR*(DphiR+(m-1.)/(2.*r)*(UR-S))+(m-1.)/(2.*r)*C_star*(phic+(U-S));
         	DtU=dR+(UM-S)*(UM-S)*(m-1.)/r;DtP=D*(UM-S)*dR;
			if(UM>=0.)
			{
				theta=C_star/CL;
			DtD=1./(C_star*C_star)*(D*(UM-S)*dR+(Gamma-1)*D*(UM-S)*pow(theta,(1.+mus)/mus)*TDSL);
			}
		else
			{
			H1=DR*PR*(1-mus*mus)/((PR+mus*PM)*(PR+mus*PM));
            H2=DR*PM*(-1+mus*mus)/((PR+mus*PM)*(PR+mus*PM));
			H3=(PM+mus*PR)/(PR+mus*PM);
			sigmaR=(D*(U-S)-DR*(UR-S))/(D-DR);
			gD=((UM-S)-sigmaR)*(C_star*C_star-(UM-S)*(UM-S));gU=(UM-S)*D*(sigmaR)*(H1*C_star*C_star-1.);
			gP=sigmaR-((UM-S)*(sigmaR-(UM-S))+C_star*C_star)*(UM-S)*H1;
			f=(sigmaR-UR)*H2*DPR+(sigmaR-UR)*H3*DDR-DR*(H2*CR*CR+H3)*(DUR+UR*(m-1.)/r);
			f=(UM-S)*(C_star*C_star-(UM-S)*(UM-S))*f+(UM-S)*(UM-S)*gU*(m-1.)/r;
			DtD=(f-gU*(dR+(UM-S)*(UM-S)*(m-1.)/r)-gP*D*(UM-S)*dR)/gD;
			}
           printf("RR\n");
//Dttau=-1./(D*D*C_star*C_star)*DtP;
		}//end Right rarefaction wave
	
return OK;
}
Status GRPsolverSS1(double &DtD,double &DtU,double &DtP,double DL,double DR,double UL,double UR,double PL,
double PR,double DDL,double DDR,double DUL,double DUR,double DPL,double DPR,double S,double r,double Gamma)//smooth
{
	double SHL,SHR,CL,CR;CL=sqrt(Gamma*PL/DL);CR=sqrt(Gamma*PR/DR);
SHL=UL-CL;SHR=UR+CR;
if(S<=SHL)//Left region
{
DtD=-(DL*DUL+(UL-S)*DDL+(UL-S)*DL*(m-1.)/r);
DtU=-(DPL/DL+(UL-S)*DUL);
DtP=-((UL-S)*DPL+Gamma*PL*DUL+Gamma*PL*(UL-S)*(m-1.)/r);
}
if(S>=SHR)//Right region
{
DtD=-(DR*DUR+(UR-S)*DDR+(UR-S)*DR*(m-1.)/r);
DtU=-(DPR/DR+(UR-S)*DUR);
DtP=-((UR-S)*DPR+Gamma*PR*DUR+Gamma*PR*(UR-S)*(m-1.)/r);
}
return OK;
}
#endif