function ics = nrlmass_hist(t)

global N numddevars myrand

ics = zeros(numddevars*N,1);

for i=1:N
    ics((i-1)*3+1:i*3,1) = [(myrand(i)-0.5)*0.8-0.2; (myrand(i)-0.5)*0.6+0.3; (myrand(i)-0.5)*0.16+0.05];
end;

