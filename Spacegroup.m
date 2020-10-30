%% make atom
clear all;clear all;tic;
% Xt1=[];
% Xt2=[];
% Xt3=[];
X=[ -0.5 0 0.5; 0.0 -0.5 0.5; 0.5 0.0 0.5; 0 0.5 0.5; ...
    -0.5 -0.5 0; 0.5 -0.5 0; -0.5 0.5 0; 0.5 0.5 0; ...
    0 -0.5 -0.5; -0.5 0 -0.5; 0 0.5 -0.5; 0.5 0 -0.5];
% for k=0:2
%       Xt1=[Xt1;[X1(:,1),X1(:,2),X1(:,3)+k]];
% end
% for k=0:2
%       Xt2=[Xt2;[Xt1(:,1),Xt1(:,2)+k,Xt1(:,3)]];
% end
% for k=0:2
%       Xt3=[Xt3;[Xt2(:,1)+k,Xt2(:,2),Xt2(:,3)]];
% end
% X=[Xt1;Xt2;Xt3];
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
%%%%%%%%%%%%%%%%%%%%%%%%%
%% find atom location%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:length(atom) 
%     at=nX(i,:)-[3,0,0];
%     bt=nX(i,:)-[0,3,0];
%     ct=nX(i,:)-[0,0,3];
%     att=[];
%     btt=[];
%     ctt=[];
%     for j=1:length(atom)
%         att(end+1)=sum(nX(j,:)==at)==3;
%         btt(end+1)=sum(nX(j,:)==bt)==3;
%         ctt(end+1)=sum(nX(j,:)==ct)==3;
%     end
%     if sum(att)==1&sum(btt)==1&sum(ctt)==1
%         miNx=i;
%     end
% end
% v=X(length(X))-nX(miNx,:);
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
all={};
total=0;
ps=[2 3 5 7 11 13 17 19 23 29 31 37];
nnn=[1 -1 2 -2 3 -3 4 -4 5 -5 6 -6];
while sum(total)~=prod(1:12)/prod(1:4)/prod(1:4)/prod(1:4)
disp(sum(total))
total=[];
allt=[];
ck=0;
%%%%%%%%%%%random%%%%%%%%%%
na=randperm(12);
for i=1:4:12
atom(na(i:i+3))=ps(i);
end
%%%%%%%%%%%group build%%%%%%%%%%%%%%%%
atomt=atom;
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
    allts=0;
    for j=1:12
    allts=allts+atom(j)^nnn(j);
    end
    allt(end+1)=allts;
%%%%%%%%%%%%change pls
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(all)
all{end+1}=allt;
end
for kk=length(all):-1:1
        if sum(ismember(all{kk},allt)) > 0
            ck=kk;
            break
        end
end
if ck==0
    all{end+1}=allt;
    ck=length(all);
else
         all{ck}=[all{ck},allt];
end
all{ck}=unique(all{ck});
for ka=1:length(all)
    total(end+1)=length(all{ka});
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
su=sum(total);
for i=1:length(total)
    total(i)=total(i)/su;
end
S=0;
for i=1:length(total)
    S=S-total(i)*log(total(i));
end
toc;




