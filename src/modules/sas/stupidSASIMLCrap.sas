
options mprint macrogen;
%macro deb(number);
%do j = 1 %to &number;
proc iml;
  do i = 1 to 10;
    array = array//i;
  end;

name = {'mama'};

datasetName = "foo&j";
call execute("create ", datasetName, " from array[colname = name];");
append from array;
close one&j;

free array;
quit;

proc print data = foo&j;
run;

%end;
%mend;

%deb(3);
run;

