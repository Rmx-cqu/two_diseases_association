### new file 
### EAS liver to Ps
exposures='bbj-a-158'
outcome='ebi-a-GCST90018687'
psoriasis=extract_instruments(outcomes = exposures,
                              clump = TRUE,
                              r2 = 0.01)
chd_out_dat<-extract_outcome_data(snps=psoriasis$SNP,
                                  outcomes = outcome,
                                  proxies = T)
dat <- harmonise_data(
  exposure_dat = psoriasis, 
  outcome_dat = chd_out_dat)
single <- mr_leaveoneout(dat)
###check the outlier 
mr_leaveoneout_plot(single)
dat=dat[-3,]
calculate_F=function(beta,se){
  F_value=beta^2/se^2
  return(F_value)
}
ps1=psoriasis[-3,]
calculate_F(ps1$beta.exposure,ps1$se.exposure) |> mean()
res=mr(dat)
generate_odds_ratios(res)
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
#get_IV_F=function(exposure_data1,Ncases,Ncontrols,exposure=T){
  ### calculate the Rsqure
  IV=dim(exposure_data1)[1]
  r_2=function(eaf,b,se){
    if(IV<=5){
      return(2*eaf*(1-eaf)*(b^2))}
    else{
      N=2/(1/Ncases+1/Ncontrols)
      x=2*eaf*(1-eaf)*(b^2)
      y=2*eaf*(1-eaf)*N*(se^2)
      return(x/(x+y))
    }
  }
  if (exposure){
    R2=r_2(exposure_data1$eaf.exposure,
           exposure_data1$beta.exposure,
           exposure_data1$se.exposure)
  }
  else(
    R2=r_2(exposure_data1$eaf.outcome,
           exposure_data1$beta.outcome,
           exposure_data1$se.outcome)
  )
  print(R2)
  print(sum(R2))
  
  IV=length(R2)
  ### calculate the F value
  calculate_F=function(Ncase=Ncases,Ncontrol=Ncontrols,k=IV,R_2){
    N=2/(1/Ncase+1/Ncontrol)
    return((N-k-1)/k*R_2/(1-R_2))
  }
  return(calculate_F(R_2 = R2))
}
#get_IV_F(dat1,1866,195745,exposure=T)

#####EUR
##exposure liver cancer
### outcome psoriasis
library(data.table)
### read outcome
df=fread('C:/Users/rmx/Downloads/summary_stats_finngen_R7_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC.gz',
         header = T)
df1=df[df$pval<5e-08,] 

df1$phenotype_col='liver'
df1=as.data.frame(df1)
exposure=format_data(
  df1,
  type='exposure',
  snp_col='rsids',
  phenotype_col = 'phenotype_col',
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col="af_alt",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval")
exposure_data=clump_data(exposure,
                         clump_r2 = 0.01)
calculate_F(exposure_data$beta.exposure,exposure_data$se.exposure) |> mean()
chd_out_dat<-extract_outcome_data(snps=exposure_data$SNP,
                                  outcomes = 'finn-b-L12_PSORIASIS',
                                  proxies =T)
dat=harmonise_data(exposure_dat=exposure_data,
                   outcome_dat = chd_out_dat)
res=mr(dat)
res
generate_odds_ratios(res)
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

check_SNP=phenoscanner(snpquery = dat$SNP)
check_SNP$results$trait
get_IV_F(dat,518,238678,exposure=T)

#### 
get_data=function(filename,type='exposure',clump_R2=0.01,phenotype='liver'){
  filename=file.path('C:/Users/rmx/Downloads/',filename)
  df=fread(filename)
  if(type=='exposure'){
    df1=df[df$pval<5e-08,] 
    df1$phenotype_col=phenotype
    df1=as.data.frame(df1)
    exposure=format_data(
      df1,
      type='exposure',
      snp_col='rsids',
      phenotype_col = 'phenotype_col',
      beta_col = "beta",
      se_col = "sebeta",
      eaf_col="af_alt",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval")
    data=clump_data(exposure,
                    clump_r2 = clump_R2)
  }else{
    df1=df[df$pval>1e-05,]
    df1$phenotype_col=phenotype
    df1=as.data.frame(df1)
    outcome=format_data(
      df1,
      type='outcome',
      snp_col='rsids',
      phenotype_col = 'phenotype_col',
      beta_col = "beta",
      se_col = "sebeta",
      eaf_col="af_alt",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval")
    data=outcome
  }
  return(data)
}
liver_exp=get_data(filename='summary_stats_finngen_R7_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC.gz',
                   type='exposure',
                   clump_R2=0.01,
                   phenotype='liver')
liver_out=get_data(filename='summary_stats_finngen_R7_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC.gz',
                   type='outcome',
                   clump_R2=0.01,
                   phenotype='liver')
psoriasis_exp=get_data(filename='finngen_R7_L12_PSORIASIS.gz',
                       type='exposure',
                       clump_R2=0.01,
                       phenotype='psoriasis')
psoriasis_out=get_data(filename='finngen_R7_L12_PSORIASIS.gz',
                       type='outcome',
                       clump_R2=0.01,
                       phenotype='psoriasis')

dat=harmonise_data(liver_exp,
                   psoriasis_out)


### revesre  EUR,Ps to liver
exposures='finn-b-L12_PSORIASIS'
ps_intruc=extract_instruments(outcomes = exposures,
                              clump = TRUE,
                              r2 = 0.01)
dat=harmonise_data(ps_intruc,liver_out)
mr(dat)
#get_IV_F(dat,4510,212242,exposure=T)
calculate_F(ps_intruc$beta.exposure,ps_intruc$se.exposure) |> mean()
check_snp=phenoscanner(dat$SNP)
check_snp$results[232]
check_snp$results[232,]
index=grep('B',check_snp$results$trait)
check_snp$results[index,]


### reverse EAS,Ps to liver
outcome='bbj-a-158'
exposures='ebi-a-GCST90018687'
psoriasis=extract_instruments(outcomes = exposures,
                              clump = TRUE,
                              r2 = 0.01)

psoriasis
df=phenoscanner(snpquery = psoriasis$SNP)
df$results$ancestry

calculate_F=function(beta,se){
  F_value=beta^2/se^2
  return(F_value)
}

calculate_F(psoriasis$beta.exposure,psoriasis$se.exposure) |> mean()
chd_out_dat<-extract_outcome_data(snps=psoriasis$SNP,
                                  outcomes = outcome,
                                  proxies = T)

dat <- harmonise_data(
  exposure_dat = psoriasis, 
  outcome_dat = chd_out_dat)
res=mr(dat)
#get_IV_F(dat,206,172289, exposure = T)
phenoscanner(snpquery = dat$SNP)
