PROC IMPORT OUT= WORK.fuckYou 
            DATAFILE= "C:\Users\kreidles\Documents\Git\bitbucket\Kenward
RogerPowerFixedUnbalanced\validationExperiment\data\longitudinalEmpirica
l.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
