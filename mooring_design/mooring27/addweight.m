clear all
load mdcodes
[n,m]=size(floats);
if m>50, 'ALREADY DONE!', clear all; return; end 
% code to add weight column to buoyancy for float weight (all others=0)
text=floats;
[n,m]=size(text);
text(:,36:59)=text(:,27:50); % shift over make room for weight
fill='     0.0 '; % set the fill text string
for r=1:n, 
    text(r,27:35)=fill;
end
floats=text;
% 
text=acrels;
[n,m]=size(text);
text(:,36:59)=text(:,27:50); % shift over make room for weight
for r=1:n, 
    text(r,27:35)=fill;
end
acrels=text;
%
text=anchors;
[n,m]=size(text);
text(:,36:59)=text(:,27:50); % shift over make room for weight
for r=1:n, 
    text(r,27:35)=fill;
end
anchors=text;
%
text=chains;
[n,m]=size(text);
text(:,36:59)=text(:,27:50); % shift over make room for weight
for r=1:n, 
    text(r,27:35)=fill;
end
chains=text;
%
text=cms;
[n,m]=size(text);
text(:,36:59)=text(:,27:50); % shift over make room for weight
for r=1:n, 
    text(r,27:35)=fill;
end
cms=text;
text=miscs;
[n,m]=size(text);
text(:,36:59)=text(:,27:50); % shift over make room for weight
for r=1:n, 
    text(r,27:35)=fill;
end
miscs=text;
%
text=wires;
[n,m]=size(text);
text(:,36:59)=text(:,27:50); % shift over make room for weight
for r=1:n, 
    text(r,27:35)=fill;
end
wires=text;
clear r fill text n m
save mdcodes

