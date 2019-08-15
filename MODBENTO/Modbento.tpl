GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");

TOP_OF_MAIN_SECTION
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5955); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 


DATA_SECTION
 init_int ntime  
 init_int nedades
 init_number minedad
 init_int ntallas

 init_matrix mdatos(1,ntime,1,6)
 init_vector Tallas(1,ntallas)

 init_matrix Ctot(1,ntime,1,ntallas)
 init_vector msex(1,ntallas)
 init_vector Wmed(1,ntallas)

//!! ad_comm::change_datafile_name("modbento_opt.ctl");
 init_number sigmaR
 init_vector dt(1,2)
 init_vector Par_bio(1,5)
 init_number hprior
 init_number bprior


  number log_Loprior
  number log_cva_prior
  number log_h_prior
  number log_b_prior

  
  !! log_Loprior = log(Par_bio(3));
  !! log_cva_prior = log(Par_bio(4));
  !! log_b_prior = log(bprior);

 init_number L50prior
 init_number s1prior

 number log_L50fprior
 number log_s1prior

 !! log_L50fprior = log(L50prior);
 !! log_s1prior = log(s1prior);

 init_int    nbloques1
 init_vector ybloques1(1,nbloques1)

 init_int    nqbloques
 init_vector yqbloques(1,nqbloques)

 init_int    opt_qf
 init_int    opt_bpow
 init_int    opt1_fase
 init_int    opt_tiposel
 init_int    opt_Lo
 init_int    opt_cva
 init_int    opt_F
 init_int    opt_devRt
 init_int    opt_devNo

// Proyeccion de poblaci�n 
 init_int    npbr
 init_vector pbr(1,npbr)
 init_int ntime_sim
 init_int opt_Frms//*********agregar
 init_number Frms//*********agregar

INITIALIZATION_SECTION

  log_Lo        log_Loprior
  log_cv_edad   log_cva_prior
  log_L50        log_L50fprior 
  log_sigma1     log_s1prior 
  log_sigma2     9.2
  log_b          log_b_prior 
  log_F          -1.6



  
PARAMETER_SECTION

// selectividad param�trica a la talla com�n
// init_bounded_vector log_L50f(1,nbloques1,-5,8,opt1_fase)  
 
 init_vector log_L50(1,nbloques1,opt1_fase)  
 init_vector log_sigma1(1,nbloques1,opt1_fase)
 init_vector log_sigma2(1,nbloques1,opt_tiposel)

// parametros reclutamientos y mortalidades)
 init_number log_Ro(1)
 init_bounded_dev_vector dev_log_Ro(1,ntime,-10,10,opt_devRt)
 init_bounded_vector dev_log_No(1,nedades,-10,10,opt_devNo)
 init_bounded_vector log_F(1,ntime,-20,0.7,opt_F) // log  mortalidad por pesca por flota

// capturabilidades
 init_vector log_qflo(1,nqbloques,opt_qf)
 init_number log_b(opt_bpow)

// Crecim
 init_number log_Lo(opt_Lo)
 init_number log_cv_edad(opt_cva)

//---------------------------------------------------------------------------------
//Defino las variables de estado 
 vector BMflo(1,ntime)
 vector Brec(1,ntime)
 vector pred_CPUE(1,ntime);
 vector pred_Desemb(1,ntime);
 vector likeval(1,5);
 vector Neq(1,nedades);

 vector Rpred(1,ntime);
 sdreport_vector Rest(1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_anos(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 vector BDo(1,ntime);
 vector No(1,nedades)
 vector prior(1,7)

 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 matrix cv_index(1,3,1,ntime)
 vector nm(1,ntime)
 vector Lmed_obs(1,ntime)
 vector Lmed_pred(1,ntime)
 vector edades(1,nedades)

 matrix S1(1,nbloques1,1,nedades)
 matrix S2(1,nbloques1,1,nedades)

 matrix Sel(1,ntime,1,nedades)
 matrix F(1,ntime,1,nedades)
 matrix Z(1,ntime,1,nedades)
 matrix S(1,ntime,1,nedades)


 matrix N(1,ntime,1,nedades)

 matrix NM(1,ntime,1,nedades)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)

 matrix pred_Ctot(1,ntime,1,ntallas)
 matrix pred_Ctot_a(1,ntime,1,nedades)

 matrix pobs(1,ntime,1,ntallas)
 matrix ppred(1,ntime,1,ntallas)

 matrix Prob_talla(1,nedades,1,ntallas)

 matrix P1(1,nedades,1,ntallas)
 matrix P2(1,nedades,1,ntallas)
 matrix P3(1,nedades,1,ntallas)
 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,nedades)

 number suma1
 number suma2
 number suma3
 number suma4

 number penalty

 number So
 number alfa
 number beta

 number Linf
 number k
 number Linfh
 number M
 number h

 number BDp
 number Npplus
 number Bp_anch 

 number nm1;
 number cuenta1;

 vector Np(1,nedades)
 vector Zpbr(1,nedades)
 vector Fpbr(1,nedades)
 vector Sp(1,nedades)

 matrix Bp(1,npbr,1,ntime_sim)
 vector CTPp(1,nedades)
 matrix Yp(1,npbr,1,ntime_sim)

 
 objective_function_value f
  
 sdreport_vector BD(1,ntime) // 
 sdreport_vector BT(1,ntime) // 
 sdreport_vector RPR(1,ntime) // 
 vector RPRlp(1,ntime) // 
 
 sdreport_number SSBo
 sdreport_vector YTPp(1,npbr)//*********afregar
// sdreport_vector RPRp(1,npbr) // RPR proyectado en la simulacion


PRELIMINARY_CALCS_SECTION

 yrs=column(mdatos,1);
 Desemb=column(mdatos,2);
 CPUE=column(mdatos,4);
 nm=column(mdatos,6);

 edades.fill_seqadd(minedad,1);

 cv_index(1)=column(mdatos,3);
 cv_index(2)=column(mdatos,5);


 Linf=Par_bio(1);
 k=Par_bio(2);
 M=Par_bio(5);
 

 Unos_edad=1;// lo uso en  operaciones matriciales con la edad
 Unos_anos=1;// lo uso en operaciones matriciales con el a�o
 Unos_tallas=1;// lo uso en operaciones matriciales con el a�o


RUNTIME_SECTION
 maximum_function_evaluations 500,2000,5000
 convergence_criteria  1e-2,1e-5,1e-5


PROCEDURE_SECTION
// se listan las funciones que contienen los calculos
 Eval_prob_talla_edad();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_deinteres();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_logverosim();
 Eval_funcion_objetivo();

 if(last_phase()){Eval_CTP();}

FUNCTION Eval_prob_talla_edad


 int i, j;

// genero una clave edad-talla para otros calculos. Se modela desde L(1)
 mu_edad(1)=exp(log_Lo);
 for (i=2;i<=nedades;i++)
  {
  mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
  }

 sigma_edad=exp(log_cv_edad)*mu_edad;

  Prob_talla = ALK( mu_edad, sigma_edad, Tallas);

//----------------------------------------------------------------------
FUNCTION dvar_matrix ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);


//----------------------------------------------------------------------

FUNCTION Eval_selectividad
 int i,j;

 for (j=1;j<=nbloques1;j++){

 S1(j)=exp(-0.5*square(edades-exp(log_L50(j)))/square(exp(log_sigma1(j))));


    for (i=1;i<=nedades;i++){

      if(edades(i)>=exp(log_L50(j))){
      S1(j,i)= exp(-0.5*square(edades(i)-exp(log_L50(j)))/square(exp(log_sigma2(j))));
      }

 }}

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sel(i)=S1(j);}
       }
   }



FUNCTION Eval_mortalidades

 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_edad));

 Z=F+M;

 S=mfexp(-1.0*Z);


FUNCTION Eval_abundancia
 int i, j;


 h=hprior;

 // Biomasa desovante virgen de largo plazo
 No(1)=exp(log_Ro+0.5*square(sigmaR)); //

 for (int j=2;j<=nedades;j++)
 {   No(j)=No(j-1)*exp(-1.*M);}
     No(nedades)=No(nedades)/(1-exp(-1.*M));

 SSBo=sum(elem_prod(No*exp(-dt(1)*M)*Prob_talla,elem_prod(msex,Wmed)));

 alfa=4*h*exp(log_Ro+0.5*square(sigmaR))/(5*h-1);//
 beta=(1-h)*SSBo/(5*h-1);// Reclutamiento

// Abundancia inicial

 N(1)=elem_prod(No,exp(dev_log_No));

 BD(1)=sum(elem_prod(elem_prod(N(1),exp(-dt(1)*Z(1)))*Prob_talla,elem_prod(msex,Wmed)));

 Rpred(1)=exp(log_Ro+0.5*square(sigmaR));//


// se estima la sobrevivencia por edad(a+1) y a�o(t+1)
 for (i=1;i<ntime;i++)
 {
    Rpred(i+1)=exp(log_Ro+0.5*square(sigmaR)); // 

    if(i>minedad){
     Rpred(i+1)=(alfa*BD(i-minedad)/(beta+BD(i-minedad)));
    } // 

     N(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i));  // 

     N(i+1)(2,nedades)=++elem_prod(N(i)(1,nedades-1),S(i)(1,nedades-1));
     N(i+1,nedades)+=N(i,nedades)*S(i,nedades);// grupo plus

     BD(i+1)=sum(elem_prod(elem_prod(N(i+1),exp(-dt(1)*Z(i+1)))*Prob_talla,elem_prod(msex,Wmed)));
 }

   Rest=column(N,1);

FUNCTION Eval_deinteres

// Rutina para calcular RPR
 Nv=N;// solo para empezar los calculos

// se estima la sobrevivencia por edad(a+1) y a�o(t+1)
 for (int i=1;i<ntime;i++)
 {
     Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*M);
     Nv(i+1,nedades)=Nv(i+1,nedades)+Nv(i,nedades)*exp(-1.0*M);// grupo plus
 }


 NDv=elem_prod((Nv*exp(-dt(1)*M))*Prob_talla,outer_prod(Unos_anos,msex));
 BDo=NDv*Wmed;
 RPR=elem_div(BD,BDo);


 RPRlp=BD/SSBo;


FUNCTION Eval_biomasas
 
 NMD=elem_prod(N,mfexp(-dt(1)*Z))*Prob_talla;
 NMD=elem_prod(NMD,outer_prod(Unos_anos,msex));
 
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel)*Prob_talla;

// vectores de biomasas derivadas
 BD=NMD*Wmed;
 BMflo=NVflo*Wmed;
 BT=(N*Prob_talla)*Wmed;


FUNCTION Eval_capturas_predichas

// matrices de capturas predichas por edad y a�o
 pred_Ctot_a=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
 pred_Ctot=pred_Ctot_a*Prob_talla;


// vectores de desembarques predichos por a�o
 pred_Desemb=pred_Ctot*Wmed;

// matrices de proporcion de capturas por talla y a�o
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));

 
 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);

FUNCTION Eval_indices
 

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*pow(BMflo(i),exp(log_b));}
       }
   }


FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0; suma2=0; penalty=0;

 for (i=1;i<=ntime;i++)
 {
  if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}
 }

FUNCTION Eval_funcion_objetivo

 suma3=0; suma4=0; penalty=0;

 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb

 for (int i=1;i<=ntime;i++){
 suma3+=-nm(i)*sum(elem_prod(pobs(i),log(ppred(i))));
 }

 likeval(3)=suma3;//

// lognormal Ninicial y Reclutas
 if(active(dev_log_Ro)){
 likeval(4)=1./(2*square(sigmaR))*norm2(dev_log_Ro);}

 if(active(dev_log_No)){
 likeval(5)=1./(2*square(sigmaR))*norm2(dev_log_No);}

 if (active(log_F)){
 penalty+=1000*norm2(log_F-mean(log_F));}

 f=(sum(likeval)+penalty);

 if(last_phase){
 f=sum(likeval);}

FUNCTION  Eval_CTP
//-----------------------------------------------------------------

  for (int i=1;i<=npbr;i++){ // ciclo de PBR

  Np=N(ntime);
  Sp=S(ntime);

  for (int j=1;j<=ntime_sim;j++){ // ciclo de a�os

  if(j==1){
  Np(1)=(alfa*BD(ntime)/(beta+BD(ntime)));}
  if(j>1){
  Np(1)=(alfa*Bp(i,j-1)/(beta+Bp(i,j-1)));}
 
  Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
  Np(nedades)+=Np(nedades)*Sp(nedades);

  Fpbr=F(ntime)*pbr(i);//
  
 if(opt_Frms>0)//agregada 1 Activado Frms/ -1 Activa F ultimo a�o
  {
   Fpbr=Frms*pbr(i);//agregada Activar o desactivar -1
  }

  Zpbr=Fpbr+M;

  Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dt(1)*Zpbr))*Prob_talla,elem_prod(msex,Wmed)));
  CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  Yp(i,j)=sum(elem_prod(CTPp*Prob_talla,Wmed));
  Sp=exp(-1.*Zpbr);
  }
 YTPp(i)=Yp(i,1);//agregada. Ver en STD

 }



REPORT_SECTION

 report << "Resultados Modelo Captura a la Talla (Modbento)" << endl;
 report << "cristian.canales@ifop.cl, 2014" << endl;
 report << "-----------------------------------------------------------"<<endl;

 report << "A�os" << endl;
 report << yrs << endl;
 report << "CPUE_obs & pred" << endl;
 report << CPUE << endl;
 report << pred_CPUE << endl;
 report << "Desemb_obs & pred" << endl;
 report << Desemb << endl;
 report << pred_Desemb << endl;
 report << "Lmed_obs & pred" << endl;
 report << Lmed_obs << endl;
 report << Lmed_pred << endl;
 report << "Biomasa_desovante" << endl;
 report << BD << endl;
 report << "Biomasa_total" << endl;
 report << BT << endl;
 report << "Biomasa_explotable" << endl;
 report << BMflo << endl;
 report << "Reclutamiento  predicho & Estimado" << endl;
 report << Rpred<< endl;
 report << column(N,1)<< endl;
 report << "F " << endl;
 report << exp(log_F) << endl;
 report << "Reducci�n del stock (BD/BDo) de LP" << endl;
 report << RPRlp << endl;
 report << "Edades"<< endl;
 report << edades<< endl;
 report <<"Abundancia a la edad por a�o"<<endl;
 report <<N<< endl;
 report <<"Selectividad a la edad por a�o"<<endl;
 report <<Sel<< endl;
 report <<"Mort por pesca a la edad por a�o"<<endl;
 report <<F<< endl;
 report <<"CAPTURAS"<< endl;
 report <<"pobs"<< endl;
 report <<pobs<< endl;
 report <<"ppred"<< endl;
 report <<ppred<< endl;
 
// report << "-----------------------------------------------" << endl;
 report << "Tallas"<< endl;
 report << Tallas<< endl;
// report << "-----------------------------------------------" << endl;
 report << "Abundancia_talla_pobl" << endl;
 report << N*Prob_talla<< endl;
 report << "Proptalla_obs" << endl;
 report << pobs<< endl;
 report << "Proptalla_pred" << endl;
 report << ppred<< endl;
 report << "Prob_talla" << endl;
 report << Prob_talla << endl;
// report << "-----------------------------------------------" << endl;
 report << "BD_virginal" << endl;
 report << SSBo << endl;
 report << "Talla_media_edad" << endl;
 report << edades<< endl;
 report << mu_edad<< endl;
// report << "-----------------------------------------------" << endl;


 //report << "BD virgen anual" << endl;
 //report << BDo << endl;


 //report << "Reducci�n del stock (BD/BDo) din�mico" << endl;
 //report << RPR << endl;

// report << "-----------------------------------------------" << endl;
// report << "Lo    CV    " << endl;
// report << exp(log_Lo) <<" "<<exp(log_cv_edad)<<endl;
// report << "-----------------------------------------------" << endl;

 

 //report << "-----------------------------------------------" << endl;
 //report << "M   #edades   h    bCPUE   sigmaR" << endl;
 //report << M <<" "<<nedades<<" "<<exp(log_h)<<" "<<exp(log_b)<<" "<<sigmaR<< endl;

//-------------------------------------------------------------------
// ESTIMA nm y CV

  suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppred(i),1-ppred(i)));
      suma2=norm2(pobs(i)-ppred(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}
 report << "Tama�o muestra ideal" <<endl;
 report <<pow(nm1,1/cuenta1)<< endl;

 
 //report << "-----------------------------------------------" << endl;
 report << "Biomasa_desovante_proyectada_para_cada_mF" << endl;
 report << Bp << endl;
 report << "Capturas proyectadas para cada mF" << endl;
 report << Yp << endl;
 report << "-----------------------------------------------" << endl;
 report << "Componentes de log-verosimilitud" << endl;
 report << "CPUE_Desemb   Prop    dev_Rt   dev_No"<<endl;
 report << likeval << endl;






FINAL_SECTION

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;


