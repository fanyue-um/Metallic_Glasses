clear;clc;
N = 2000;
Nsamples = 200;

pathFile = split(pwd, '/');
coolingRate = pathFile{end - 2};
Temperature = pathFile{end - 1};
Rexp = coolingRate(8:9);
idx = coolingRate(end);

for ii=1: Nsamples
    disp(ii);
    fid1=fopen(['../ML.', num2str(ii), '.txt'],'rt');
    fid2=fopen(['ML', Rexp, '_', idx, '_', Temperature, '_', num2str(ii), '.txt'],'w');
    while feof(fid1)~=1
        file=fgetl(fid1);
        if strcmp(file,'2 atom types')==1
            file=fgetl(fid1);
            break;
        end
    end
    bound1=split(fgetl(fid1));
    x1=str2double(bound1(1));x2=str2double(bound1(2));
    bound2=split(fgetl(fid1));
    y1=str2double(bound2(1));y2=str2double(bound2(2));
    bound3=split(fgetl(fid1));
    z1=str2double(bound3(1));z2=str2double(bound3(2));

    while feof(fid1)~=1
        file=fgetl(fid1);
        if strcmp(file,'Atoms # atomic')==1
            file=fgetl(fid1);
            break;
        end
    end

    atom=zeros(6,N);
    for tmp=1:N
        file=fgetl(fid1);
        data=split(file);
        atom(1,tmp)=str2double(data(1));%%ID
        atom(2,tmp)=str2double(data(2));%%Type
        atom(3,tmp)=str2double(data(3));%%X
        atom(4,tmp)=str2double(data(4));%%Y
        atom(5,tmp)=str2double(data(5));%%Z
    end

    atom_sort = sortrows(atom', [2, 5])';
    atom = atom_sort;

    fprintf(fid2,'\n');
    fprintf(fid2,[num2str(N),' atoms\n']);
    fprintf(fid2,'2 atom types\n');

    fprintf(fid2,'\n');
    fprintf(fid2,[num2str(x1,12),' ',num2str(x2,12),' xlo xhi\n']);
    fprintf(fid2,[num2str(y1,12),' ',num2str(y2,12),' ylo yhi\n']);
    fprintf(fid2,[num2str(z1,12),' ',num2str(z2,12),' zlo zhi\n']);
    fprintf(fid2,'\n');
    fprintf(fid2,'Masses\n');
    fprintf(fid2,'\n');
    fprintf(fid2,'1 91.22\n');% Zr
    fprintf(fid2,'2 63.55\n');% Cu
    fprintf(fid2,'\n');
    fprintf(fid2,'Atoms # atomic\n');
    fprintf(fid2,'\n');

    for i=1:N
        fprintf(fid2,[num2str(i,20),' ',num2str(atom(2,i),20),' ',num2str(atom(3,i),20),' ',num2str(atom(4,i),20),' ',num2str(atom(5,i),20),'\n']);
    end
    fclose('all');
end



%%%%%% Energy Part %%%%%%%
fid1=fopen('../minimizeenergy.txt','rt');
fid2=fopen(['ML', Rexp, '_', idx, '_Energy.txt'],'w');
for i = 1: Nsamples*3
    file=fgetl(fid1);
    if mod(i, 3) == 0
        data = split(file, ',');
        energy = str2double(data(2));
        fprintf(fid2, [num2str(energy,20), '\n']);
    end
end
fclose('all');







