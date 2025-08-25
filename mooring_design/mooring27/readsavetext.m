% script to read and save help text 
clear all
fid=fopen('addhelp.txt');
for l=1:22,
    line=fgets(fid);
    text(l,:)=line;
end
fclose(fid);
addhelp=text;
clear fid l line text ans