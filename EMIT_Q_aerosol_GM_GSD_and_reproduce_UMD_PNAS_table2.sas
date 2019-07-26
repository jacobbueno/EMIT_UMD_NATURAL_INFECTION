* Script: reproduce_PNAS_table2.SAS;
* Author: Jacob Bueno de Mesquita;
* Date: March 2019;
* Summary: This script imports dataset (produced with R script in the EMIT_QUARANTINE project),
and produces a tobit mixed model with outcome of G-II measured aerosol shedding (per half hour).
* I also do a tobit model for the Quarantine study - but it is simplified given the nature of the data.


* Import the data;
/* Generated Code (IMPORT) */
/* Source File: pos_samples_tobit_table2.csv */
/* Source Path: /home/jbueno/sasuser.v94/Milton Lab/Dissertation */
/* Code generated on: 3/5/19, 3:28 PM */

%web_drop_table(WORK.table2);


FILENAME REFFILE '/home/jbueno/sasuser.v94/Milton Lab/Dissertation/pos_samples_tobit_table2.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.table2;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=WORK.table2; RUN;


%web_open_table(WORK.table2);



data table2;
	set table2;

	if final.copies ~=. and Fine=1 then
		logfinalcopies=log10(final.copies);

	if final.copies=. and Fine=1 and type='A' then
		logfinalcopies=log10(2000);

	if final.copies=. and Fine=1 and type='B' then
		logfinalcopies=log10(9000);

	if final.copies ~=. and Coarse=1 then
		logfinalcopies=log10(final.copies);

	if final.copies=. and Coarse=1 and type='A' then
		logfinalcopies=log10(2000);

	if final.copies=. and Coarse=1 and type='B' then
		logfinalcopies=log10(9000);
run;



* Manipulate the study.day variable to get dummy variables. ;
data table2;
	set table2;

	if dpo = 2 then
		X1 = 1;
	else
		X1 = 0;

	if dpo = 3 then
		X2 = 1;
	else
		X2 = 0;

run;
quit;


data table2; 
	set table2;
	subject_id_char = put(subject.id, 3.);
*	rename subject.id = subject_id;
RUN;


data finetotal;
	set table2;

	if Fine=0 then
		delete;
run;
quit;

data coarsetotal;
	set table2;

	if Coarse=0 then
		delete;
run;
quit;



* UMD_EMIT_Reproducing Table 2 GM and GSD.
**********************************************
**********************************************
* Run proc genmod and get out parameter estimates for each parameter. 
* GEMOD with marginal effects for FINE aerosol as outcome (all table 3 adjusted parameters);
proc genmod data = finetotal descending;
	class subject_id_char sex vax_bothyear X1 X2 ;
	model logfinalcopies = age sex BMI vax_bothyear X1 X2 cough_number sex*cough_number /dist = poisson;
	repeated subject = subject_id_char/ printmle;
	ods OUTPUT parameterestimates = params_fine;
run;


data params_fine;
	set params_fine (keep = Parameter Estimate);
run;

proc transpose data = params_fine out = paramst_fine;
run;

data _null_;
	set paramst_fine;
	call symput('unfBeta0', COL1);
	call symput('unfBeta1', COL2);
	call symput('unfBeta2', COL4);
	call symput('unfBeta3', COL5);
	call symput('unfBeta4', COL7);
	call symput('unfBeta5', COL9);
	call symput('unfBeta6', COL11);
	call symput('unfBeta7', COL12);
	call symput('unfBeta8', COL14);
run;


* According to the UCLA tutorial, the proc mixed is another way of doing this, 
* However it looks like the proc mixed uses a different transformation of the data and we will stick with the genmod for now.;
*proc mixed data = quar_alpha;
*  class subject_id_char;
*  model logfinalcopies = cough / solution;
*  random int / subject = subject_id_char;
*run;


* Apply Tobit regression with NLMixed for fine aerosol;
proc nlmixed data = finetotal XTOL = 1E-12 method = GAUSS qpoints = 100;
	Title "Tobit Regression of fine Particles";
	parms sigma2_u = 1 sigma2 = 1 beta0 = &unfBeta0 beta1 = &unfBeta1 beta2 = &unfBeta2 beta3 = &unfBeta3 
	beta4 = &unfBeta4 beta5 = &unfBeta5 beta6 = &unfBeta6 beta7 = &unfBeta7 beta8 = &unfBeta8;
	bounds sigma2_u sigma2 >= 0;
	pi = constant('pi');
	mu = beta0 + b_0j + beta1*age + beta2*sex + beta3*BMI + beta4*vax_bothyear + beta5*X1 + beta6*X2 +
	beta7*cough_number + beta8*cough_number*sex;
	if final.copies ne . then
		ll = (1 / (sqrt(2*pi*sigma2))) * exp(-(logfinalcopies - mu)**2 / (2*sigma2) );
	if final.copies = .  then
		ll = probnorm((logfinalcopies - mu) / sqrt(sigma2) );
	L = log(ll);
	model logfinalcopies ~ general(L);
	random b_0j ~ normal(0, sigma2_u) subject = subject_id_char;
	predict mu out = output_mu_fine;
run;
quit;


**********************************************
*** New model -- using only the significant main effects (From original submission PNAS table 3 adjusted model) ;

* GEMOD with marginal effects for FINE aerosol as outcome (all table 3 adjusted parameters);
proc genmod data = finetotal descending;
	class subject_id_char vax_bothyear X1 X2 sex ;
	model logfinalcopies = BMI vax_bothyear X1 X2 cough_number sex*cough_number /dist = poisson;
	repeated subject = subject_id_char/ printmle;
	ods OUTPUT parameterestimates = params_fine;
run;


data params_fine;
	set params_fine (keep = Parameter Estimate);
run;

proc transpose data = params_fine out = paramst_fine;
run;

data _null_;
	set paramst_fine;
	call symput('unfBeta0', COL1);
	call symput('unfBeta1', COL2);
	call symput('unfBeta2', COL4);
	call symput('unfBeta3', COL6);
	call symput('unfBeta4', COL8);
	call symput('unfBeta5', COL9);
	call symput('unfBeta6', COL11);
run;

* Apply Tobit regression with NLMixed for fine aerosol;
proc nlmixed data = finetotal XTOL = 1E-12 method = GAUSS qpoints = 100;
	Title "Tobit Regression of fine Particles";
	parms sigma2_u = 1 sigma2 = 1 beta0 = &unfBeta0 beta1 = &unfBeta1 beta2 = &unfBeta2 beta3 = &unfBeta3 
	beta4 = &unfBeta4 beta5 = &unfBeta5 beta6 = &unfBeta6;
	bounds sigma2_u sigma2 >= 0;
	pi = constant('pi');
	mu = beta0 + b_0j + beta1*BMI + beta2*vax_bothyear + beta3*X1 + beta4*X2 + beta5*cough_number + 
	beta6*cough_number*sex;
	if final.copies ne . then
		ll = (1 / (sqrt(2*pi*sigma2))) * exp(-(logfinalcopies - mu)**2 / (2*sigma2) );
	if final.copies = .  then
		ll = probnorm((logfinalcopies - mu) / sqrt(sigma2) );
	L = log(ll);
	model logfinalcopies ~ general(L);
	random b_0j ~ normal(0, sigma2_u) subject = subject_id_char;
	predict mu out = output_mu_fine_UMD_positives;
run;
quit;



**********************************************
**********************************************
* GEMOD with marginal effects for COARSE aerosol as outcome;
proc genmod data = coarsetotal descending;
	class subject_id_char sex vax_bothyear X1 X2 ;
	model logfinalcopies = age sex BMI vax_bothyear X1 X2 cough_number sex*cough_number /dist = poisson;
	repeated subject = subject_id_char/ printmle;
	ods OUTPUT parameterestimates = params_coarse;
run;


data params_coarse;
	set params_coarse (keep = Parameter Estimate);
run;

proc transpose data = params_coarse out = paramst_coarse;
run;

data _null_;
	set paramst_coarse;
	call symput('unfBeta0', COL1);
	call symput('unfBeta1', COL2);
	call symput('unfBeta2', COL3);
	call symput('unfBeta3', COL5);
	call symput('unfBeta4', COL6);
	call symput('unfBeta5', COL8);
	call symput('unfBeta6', COL10);
	call symput('unfBeta7', COL12);
	call symput('unfBeta8', COL13);
run;

* Apply Tobit regression with NLMixed for coarse aerosol;
proc nlmixed data = coarsetotal XTOL = 1E-12 method = GAUSS qpoints = 100;
	Title "Tobit Regression of coarse Particles";
	parms sigma2_u = 1 sigma2 = 1 beta0 = &unfBeta0 beta1 = &unfBeta1 beta2 = &unfBeta2 beta3 = &unfBeta3 
	beta4 = &unfBeta4 beta5 = &unfBeta5 beta6 = &unfBeta6 beta7 = &unfBeta7 beta8 = &unfBeta8;
	bounds sigma2_u sigma2 >= 0;
	pi = constant('pi');
	mu = beta0 + b_0j + beta1*age + beta2*sex + beta3*BMI + beta4*vax_bothyear + beta5*X1 + beta6*X2 +
	beta7*cough_number + beta8*cough_number*sex;
	if final.copies ne . then
		ll = (1 / (sqrt(2*pi*sigma2))) * exp(-(logfinalcopies - mu)**2 / (2*sigma2) );
	if final.copies = .  then
		ll = probnorm((logfinalcopies - mu) / sqrt(sigma2) );
	L = log(ll);
	model logfinalcopies ~ general(L);
	random b_0j ~ normal(0, sigma2_u) subject = subject_id_char;
	predict mu out = output_mu_coarse;
run;
quit;


**********************************************
*** New model -- using only the significant main effects (From original submission PNAS table 3 adjusted model) ;

* GEMOD with marginal effects for COARSE aerosol as outcome;
proc genmod data = coarsetotal descending;
	class subject_id_char X1 X2 ;
	model logfinalcopies = BMI X1 X2 cough_number /dist = poisson;
	repeated subject = subject_id_char/ printmle;
	ods OUTPUT parameterestimates = params_coarse;
run;


data params_coarse;
	set params_coarse (keep = Parameter Estimate);
run;

proc transpose data = params_coarse out = paramst_coarse;
run;

data _null_;
	set paramst_coarse;
	call symput('unfBeta0', COL1);
	call symput('unfBeta1', COL2);
	call symput('unfBeta2', COL3);
	call symput('unfBeta3', COL5);
	call symput('unfBeta4', COL7);
run;

* Apply Tobit regression with NLMixed for coarse aerosol;
proc nlmixed data = coarsetotal XTOL = 1E-12 method = GAUSS qpoints = 100;
	Title "Tobit Regression of coarse Particles";
	parms sigma2_u = 1 sigma2 = 1 beta0 = &unfBeta0 beta1 = &unfBeta1 beta2 = &unfBeta2 beta3 = &unfBeta3 
	beta4 = &unfBeta4 ;
	bounds sigma2_u sigma2 >= 0;
	pi = constant('pi');
	mu = beta0 + b_0j + beta1*BMI + beta2*X1 + beta3*X2 + beta4*cough_number ;
	if final.copies ne . then
		ll = (1 / (sqrt(2*pi*sigma2))) * exp(-(logfinalcopies - mu)**2 / (2*sigma2) );
	if final.copies = .  then
		ll = probnorm((logfinalcopies - mu) / sqrt(sigma2) );
	L = log(ll);
	model logfinalcopies ~ general(L);
	random b_0j ~ normal(0, sigma2_u) subject = subject_id_char;
	predict mu out = output_mu_coarse_UMD_positives;
run;
quit;



**********************************************
**********************************************
**********************************************
**********************************************

EMIT_UK Computing GM and GSD for fine and coarse aerosols. ;


*** Quarantine virus shedding into fine particle aerosols. ;

/* Generated Code (IMPORT) */
/* Source File: Q_positive_fine_samples.csv */
/* Source Path: /home/jbueno/sasuser.v94/Milton Lab/Dissertation */
/* Code generated on: 3/6/19, 9:42 AM */

%web_drop_table(WORK.quarantine_fine_pos);


FILENAME REFFILE '/home/jbueno/sasuser.v94/Milton Lab/Dissertation/Q_positive_fine_samples.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.quarantine_fine_pos;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=WORK.quarantine_fine_pos; RUN;


%web_open_table(WORK.quarantine_fine_pos);


data quarantine_fine_pos; 
	set quarantine_fine_pos;
	cough_number_numeric = input(cough_number, 3.);
	Ct_numeric = input(Ct, 5.2);
	subject_id_char = put(subject.id, 3.);
	study_day_char = put(study.day, 1.);
	drop cough_number Ct;
	rename cough_number_numeric = cough_number Ct_numeric = Ct;
*	rename subject.id = subject_id study.day = study_day;
RUN;


data quarantine_fine_pos;
	set quarantine_fine_pos;

	if final.copies ~= . and Fine = 1 then
		logfinalcopies = log10(final.copies);

	if final.copies = . and Fine = 1 then
		logfinalcopies = log10(2000);
		
run;


data quarantine_fine_pos;
	set quarantine_fine_pos;

	if study.day = 3 then
		X1 = 1;
	else
		X1 = 0;

	if study.day = 4 then
		X2 = 1;
	else
		X2 = 0;
run;
quit;


* GEMOD with marginal effects of cough and study day;
proc genmod data = quarantine_fine_pos descending;
	class subject_id_char X1 X2 ;
	model logfinalcopies = cough X1 X2 /dist = poisson;
	repeated subject = subject_id_char/ printmle;
	ods OUTPUT parameterestimates = params_fine_Q;
run;

data params_fine_Q;
	set params_fine_Q (keep = Parameter Estimate);
run;

proc transpose data = params_fine_Q out = paramst_fine_Q;
run;

data _null_;
	set paramst_fine_Q;
	call symput('unfBeta0', COL1);
	call symput('unfBeta1', COL2);
	call symput('unfBeta2', COL4);
	call symput('unfBeta3', COL6);
run;


* Apply Tobit regression with NLMixed, with fixed effect cough and study day and random effect of person ;
proc nlmixed data = quarantine_fine_pos XTOL = 1E-12 method = GAUSS qpoints = 100;
	Title "Tobit Regression of fine Particles Quarantine: Main effects of cough and study day with random effect of person";
	parms sigma2_u = 1 sigma2 = 1 beta0 = &unfBeta0 beta1 = &unfBeta1 beta2 = &unfBeta2 beta3 = &unfBeta3 ;
	bounds sigma2_u sigma2 >= 0;
	pi = constant('pi');
	mu = beta0 + b_0j + beta1*cough + beta2*X1 + beta3*X2;
	if final.copies ne . then
		ll = (1 / (sqrt(2*pi*sigma2))) * exp(-(logfinalcopies - mu)**2 / (2*sigma2) );
	if final.copies = .  then
		ll = probnorm((logfinalcopies - mu) / sqrt(sigma2) );
	L = log(ll);
	model logfinalcopies ~ general(L);
	random b_0j ~ normal(0, sigma2_u) subject = subject_id_char out = output_random_effect;
	predict mu out = output_mu_Q_fine;
run;
quit;




*** Quarantine virus shedding into coarse particle aerosols. ;

%web_drop_table(WORK.quarantine_coarse_pos);


FILENAME REFFILE '/home/jbueno/sasuser.v94/Milton Lab/Dissertation/Q_positive_coarse_samples.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.quarantine_coarse_pos;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=WORK.quarantine_coarse_pos; RUN;


%web_open_table(WORK.quarantine_coarse_pos);


data quarantine_coarse_pos; 
	set quarantine_coarse_pos;
	cough_number_numeric = input(cough_number, 3.);
	Ct_numeric = input(Ct, 5.2);
	subject_id_char = put(subject.id, 3.);
	study_day_char = put(study.day, 1.);
	drop cough_number Ct;
	rename cough_number_numeric = cough_number Ct_numeric = Ct;
*	rename subject.id = subject_id study.day = study_day;
RUN;


data quarantine_coarse_pos;
	set quarantine_coarse_pos;

	if final.copies ~= . and Coarse = 1 then
		logfinalcopies = log10(final.copies);

	if final.copies = . and Coarse = 1 then
		logfinalcopies = log10(2000);
		
run;


data quarantine_coarse_pos;
	set quarantine_coarse_pos;

	if study.day = 2 then
		X1 = 1;
	else
		X1 = 0;

	if study.day = 3 then
		X2 = 1;
	else
		X2 = 0;
run;
quit;


* GEMOD with marginal effects;
proc genmod data = quarantine_coarse_pos descending;
	class subject_id_char ;
	model logfinalcopies = age /dist = poisson;
	repeated subject = subject_id_char/ printmle;
	ods OUTPUT parameterestimates = params_coarse_Q;
run;

data params_coarse_Q;
	set params_coarse_Q (keep = Parameter Estimate);
run;

proc transpose data = params_coarse_Q out = paramst_coarse_Q;
run;

data _null_;
	set paramst_coarse_Q;
	call symput('unfBeta0', COL1);
	call symput('unfBeta1', COL2);
run;


* Apply Tobit regression with NLMixed;
proc nlmixed data = quarantine_coarse_pos XTOL = 1E-12 method = GAUSS qpoints = 100;
	Title "Tobit Regression of coarse particles from Quarantine";
	parms sigma2_u = 1 sigma2 = 1 beta0 = &unfBeta0 beta1 = &unfBeta1 ;
	bounds sigma2_u sigma2 >= 0;
	pi = constant('pi');
	mu = beta0 + b_0j + beta1*age ;
	if final.copies ne . then
		ll = (1 / (sqrt(2*pi*sigma2))) * exp(-(logfinalcopies - mu)**2 / (2*sigma2) );
	if final.copies = .  then
		ll = probnorm((logfinalcopies - mu) / sqrt(sigma2) );
	L = log(ll);
	model logfinalcopies ~ general(L);
	random b_0j ~ normal(0, sigma2_u) subject = subject_id_char out = output_random_effect;
	predict mu out = output_mu_Q_coarse;
run;
quit;
* It appears that there isn't enough data here to produce any predicted mus from this dataset. 
* I tried to run this model with just cough and study days 2 and 3 (there was no data on study day 4), just cough, just age, none of them produced output predicted mus. ;
* My conclusion for now is that with just 8 observations with data (And 4 that we are trying to impute), there is just not enough data to make this model converge. 




















