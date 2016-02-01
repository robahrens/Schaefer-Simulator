//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: UBC ADMB May Workshop													 
//Date:	May 5-9, 2013														 
//Purpose:pm for hake data set											 
//Notes: 				 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	int sim;
	int rseed;
	LOCAL_CALCS
		sim=0;
		rseed=0;
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
		}
	END_CALCS
	init_adstring datafile;
	init_adstring ctlfile;
	//change to the new data file
	!!ad_comm::change_datafile_name(datafile);
	init_int syr;
	init_int eyr; //read in the max age from the data file
	init_ivector iyrs(syr,eyr);
	init_vector ct(syr,eyr); // vector for ages
	init_vector yt(syr,eyr); //vector for lengths
	init_int eof;
	int iter;
	!!iter=0;
	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	!!ad_comm::change_datafile_name(ctlfile);
	init_number ir;
	init_number ik;
	init_number iqo;
	init_number iao; //read in the max age from the data file
	init_number ipsig; //read in the max age from the data file
	init_number iptau; //read in the max age from the data file
	init_vector ieffort(syr,eyr);
	init_int eofc;
	vector fBt(1,20);
	LOCAL_CALCS
		if(eofc!=999)
		{
			cout<<"Error reading control file.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	//!!cout<<ieffort<<endl;
	//!!ad_exit(1);

PARAMETER_SECTION
	init_bounded_number r(0.001,1);//use bounded number to prevent  instability
	init_number log_k;
	init_number log_q;
	init_number log_psig(3);
	init_number log_ptau(3);
	init_bounded_vector xt(syr,eyr,-15.,15.,2);	
	//random_effects_vector xt(syr,eyr,2);
	
	number fpen;
	!!r=ir;// set initial values can also use the .pin file for this. Need to add this later
	!!log_k=log(ik);
	!!log_q=log(iqo);
	!!log_psig=log(ipsig);
	!!log_ptau=log(iptau);

	objective_function_value nll;//must define the objective function
	
	sdreport_number k;
	number q;
	number sig;
	number tau;
	vector Bt(syr,eyr+1);// we are not searching on this vector but do need to keep track of the derivatives of these values
	vector it(syr,eyr);// we are not searching on this vector but do need to keep track of the derivatives of these values
	vector nu(syr,eyr);// we are not searching on this vector but do need to keep track of the derivatives of these values

PRELIMINARY_CALCS_SECTION
  if(sim)
  {
  	run_data_simulation();
  }
	//cout<<"HEY!"<<endl;
	//ad_exit(1);

PROCEDURE_SECTION
	statedynamics();
	observationmodel();
	objectivefunction();

	if(mceval_phase())
	{ 
		forecast();
		mcmc_output();
	}
	if(last_phase()) forecast();

FUNCTION statedynamics
	k=mfexp(log_k);//transfor the log parameter that are being searched on to regular space.
	Bt(syr) = k;
	fpen = 0.0;
	for(int i=syr; i<=eyr; i++)
	{
		Bt(i+1)=posfun(Bt(i)+r*Bt(i)*(1.-Bt(i)/k)*mfexp(xt(i))-ct(i),0.001,fpen);
	}

FUNCTION observationmodel
	q=mfexp(log_q);
	it=q*Bt(syr,eyr);
	nu=log(it)-log(yt);


FUNCTION objectivefunction 
	//calculate objective function
	dvar_vector likevec(1,8);
	//likelihoods
	tau=sqrt(1./mfexp(log_ptau));
	likevec(1)=dnorm(nu,tau);
	sig=sqrt(1./mfexp(log_psig));
	likevec(2)=dnorm(xt,sig);
	//penalties
	likevec(3)=1.e5*fpen;
	//priors
	likevec(4)=dlnorm(k,log(2000),1);
	likevec(5)=dlnorm(r,log(0.35),0.2);
	likevec(6)=-log(q);

	likevec(7)=dgamma(1.0/(sig*sig),5.,0.2);
	likevec(8)=dgamma(1.0/(tau*tau),5.,0.2);

	//cout<<fpen<<endl;
	nll=sum(likevec);
	if(fpen>0)cout<<"fpen= "<<fpen<<endl;

FUNCTION mcmc_output
	
	if(iter==0)
	{
		ofstream ofs("refpar.mcmc");
		ofs<<"fmsy\t bmsy\t msy\t b/bmsy\t f/fmsy"<<endl;
		ofstream ofs1("Bt.mcmc");
		ofstream ofs2("xt.mcmc");
		ofstream ofs3("fBt.mcmc");
		//ofs<<"r\t k\t q\t sig\t"<<endl;
	}
	iter++;
	double fratio=value(-log(1-ct[eyr]/Bt[eyr])/(r/2));
	double bratio=value(Bt[eyr]/(k/2));
	ofstream ofs("refpar.mcmc",ios::app);
	ofs<<r/2<<"\t"<<k/2<<"\t"<<r*k/4<<"\t"<<bratio<<"\t"<<fratio<<endl;
	ofstream ofs1("Bt.mcmc",ios::app);
	ofs1<<Bt<<endl;
	ofstream ofs2("xt.mcmc",ios::app);
	ofs2<<xt<<endl;
	ofstream ofs3("fBt.mcmc",ios::app);
	ofs3<<fBt<<endl;

FUNCTION forecast
	
	//dvector fBt(1,20);
	fBt(1)=value(Bt(eyr+1));
	for(int i=2;i<=20;i++)
	{
		fBt(i)=fBt(i-1)+value(r)*fBt(i-1)*(1.-fBt(i-1)/value(k))-value(r)/2*fBt(i-1);
	}

FUNCTION run_data_simulation
	random_number_generator rng(rseed);
	dvector tmp_xt(syr,eyr);
	dvector tmp_nu(syr,eyr);
	dvector qt(syr,eyr);
	dvector tmp_Bt(syr,eyr+1);

	tmp_xt.fill_randn(rng);
	tmp_nu.fill_randn(rng);

	tmp_xt*=sqrt(1./ipsig);
	tmp_nu*=sqrt(1./iptau);

	k=mfexp(log_k);
	tmp_Bt(syr)=value(k);
	for(int i=syr;i<=eyr;i++)
	{
		qt(i)=iqo/(iao+(1-iao)*tmp_Bt(i)/value(k));
		ct(i)=(1.-mfexp(-qt(i)*mfexp(tmp_nu(i))*ieffort(i)))*tmp_Bt(i);
		tmp_Bt(i+1)=tmp_Bt(i)+value(r)*tmp_Bt(i)*(1.-tmp_Bt(i)/value(k))*mfexp(tmp_xt(i))-ct(i);
	}
	yt=elem_div(ct,ieffort);


REPORT_SECTION
	//if(last_phase()) forecast();
	REPORT(k);// write to the report file using the REPORT function
	REPORT(r);// write to the report file using the REPORT function
	REPORT(q);// write to the report file using the REPORT function
	REPORT(sig);// write to the report file using the REPORT function
	REPORT(tau);// write to the report file using the REPORT function
	REPORT(ct);// write to the report file using the REPORT function
	REPORT(yt);// write to the report file using the REPORT function
	REPORT(Bt);// write to the report file using the REPORT function
	REPORT(it);// write to the report file using the REPORT function
	REPORT(fBt);

TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include<contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;


