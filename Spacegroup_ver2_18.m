%% make atom
clear all;clear all;tic;
X=[ 0 0 -1; -0.5 0 -0.5; 0 -0.5 -0.5; 0 0.5 -0.5; ...
    0.5 0 -0.5; -1 0 0; -0.5 -0.5 0; -0.5 0.5 0; ...
    0 -1 0; 0 1 0; 0.5 -0.5 0; 0.5 0.5 0;1 0 0;...
    -0.5 0 0.5;0 -0.5 0.5;0 0.5 0.5;0.5 0 0.5;0 0 1];
X=unique(X,'rows');
atom=length(X);
S=[];
% for kks=1:length(atom)
soll={};
%% group computing
for z=1:24%% group counter
    sollt=[];
for i=1:atom
    a=X(i,:);
    grou={[a(1),a(2),a(3)],[-a(1),a(2),-a(3)],[a(1),-a(2),-a(3)],[-a(3),-a(2),-a(1)],[a(3),-a(2),a(1)]... % group C2
           ,[-a(1),-a(3),-a(2)],[-a(1),a(3),a(2)],[a(2),a(1),-a(3)],[-a(2),-a(1),-a(3)],[-a(1),-a(2),a(3)]...
           ,[a(2),a(3),a(1)],[-a(2),a(3),-a(1)],[-a(2),-a(3),a(1)],[a(2),-a(3),-a(1)]...% C3
           ,[a(3),a(1),a(2)],[-a(3),-a(1),a(2)],[a(3),-a(1),-a(2)],[-a(3),a(1),-a(2)]...
           ,[-a(2),a(1),a(3)],[a(1),-a(3),a(2)],[a(3),a(2),-a(1)],[a(2),-a(1),a(3)],[a(1),a(3),-a(2)],[-a(3),a(2),a(1)]};% C4
    nX(i,:)=grou{z};
end
%% group compare%%%%
%%%%%%%%%%%%%%%%%%%%
for i=1:length(X)
    for j=1:length(X)
        if nX(i,:)==X(j,:)
            sollt=[sollt;[i,j]];
        end
    end
end
soll{end+1}=sollt;
end
%% find group orbits%%%
%%%%%%%%%%%%%%%%%%%%%%%
gdn={}; 
for i=1:length(soll)
    gdt=[];
    for j=1:length(soll{i})
        gd=[];
        gd(end+1)=soll{i}(j,1);
        k=j;
        while 1
        if soll{i}(j,1)==soll{i}(soll{i}(k,2),2)
            gd(end+1)=soll{i}(soll{i}(k,2),1);
            break
        else
            gd(end+1)=soll{i}(soll{i}(k,2),1);
            k=soll{i}(k,2);
        end
        end
        gdt{end+1}=gd;
    end
    gdn{end+1}=gdt;
end
%% group orbit number%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:length(gdn)
     for j=1:length(gdn{i})
         if gdn{i}{j}~=0
            for k=2:length(gdn{i}{j})
                gdn{i}{gdn{i}{j}(k)}=0;
            end
        end
    end
 end
 for i=1:length(gdn)
     for j=1:length(gdn{i})
        if gdn{i}{j}==0
            gdn{i}{j}=[];
        end
     end
 end
gdt={};
same=[];
for i=1:length(gdn)
    gd=[];
    for j=1:length(gdn{i})
            gd(end+1)=length(gdn{i}{j});
    end
    gdt{end+1}=gd; %number of orbits
    same(end+1)=atom-sum(gdt{i}); %how many self-cycle orbit
end
%%
all=[];
alls=[];
ps=[2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61];
while sum(alls) ~= prod(1:18)/prod(1:6)/prod(1:6)/prod(1:6)%/prod(1:3)/prod(1:3)/prod(1:3)
disp(sum(alls))
allt=[];
%%%%%%%%%%%random%%%%%%%%%%
na=randperm(18);
atom(na(1:6))=1;
atom(na(7:12))=-1;
atom(na(13:18))=2;
% atom(na(10:12))=-2;
% atom(na(9:10))=3;
% atom(na(11:12))=-3;
atomt=atom;
%%%%%%%%%%%group build%%%%%%%%%%%%%%%%
for i=1:length(gdn)
    atom=atomt;
    x=cell2mat(gdn{i}');
    y=size(x);
	if y(2)==2
        t=atom(x(:,1));
        atom(x(:,1))=atom(x(:,2));
        atom(x(:,2))=t;
    elseif y(2)==3
        t=atom(x(:,3));
        atom(x(:,3))=atom(x(:,2));
        atom(x(:,2))=atom(x(:,1));
        atom(x(:,1))=t;
    elseif y(2)==4
        t=atom(x(:,4));
        atom(x(:,4))=atom(x(:,3));
        atom(x(:,3))=atom(x(:,2));
        atom(x(:,2))=atom(x(:,1));
        atom(x(:,1))=t;
    end
     allt(end+1)=prod(ps(find(atom==1)))*prod(ps(find(atom==-1)).^-1)*prod(ps(find(atom==2)).^2);...
%          *prod(ps(find(atom==-2)).^-2);%*prod(ps(find(atom==3)).^3)*prod(ps(find(atom==-3)).^-3);
%%%%%%%%%%%%change pls
end
allt=unique(allt);
alltemp=length(allt);
allt=max(allt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(all)
all(end+1)=allt;
alls(end+1)=alltemp;
elseif sum(ismember(all,allt)) == 0
all(end+1)=allt;
alls(end+1)=alltemp;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
su=sum(alls);

for i=1:length(alls)
    allss(i)=alls(i)/su;
end
S=0;
for i=1:length(allss)
    S=S-allss(i)*log(allss(i));
end
toc;