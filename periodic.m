%% make atom
clear all;clc;
Xt1=[];
Xt2=[];
Xt3=[];
coe=1;
X1=[ 0 0 0; 1 0 0; 0.5 0.5 0; 0 1 0; 1 1 0; 0.5 0 0.5; 0 0.5 0.5; 1 0.5 0.5; 0.5 1 0.5; 0 0 1; 1 0 1; 0.5 0.5 1; 0 1 1; 1 1 1];
for k=0:coe
      Xt1=[Xt1;[X1(:,1),X1(:,2),X1(:,3)+k]];
end
for k=0:coe
      Xt2=[Xt2;[Xt1(:,1),Xt1(:,2)+k,Xt1(:,3)]];
end
for k=0:coe
      Xt3=[Xt3;[Xt2(:,1)+k,Xt2(:,2),Xt2(:,3)]];
end
X=[Xt1;Xt2;Xt3];
X=unique(X,'rows');
for i=1:length(X)
    if sum(X(i,:)==0)>0
        X(i,:)=[5 5 5];
    end
end
X=unique(X,'rows');
X(end,:)=[];
atom=ones(1,length(X));
atom(16:25)=-1;
b=[];
soll={};
%% group computing
for y=1:4
for z=1:23  %% group counter
    aa={[0 0 0],[0.5 0 0.5],[0 0.5 0.5],[0.5 0.5 0]};
    sollt=[];
for i=1:length(atom)
    a=X(i,:);
    grou={[-a(1),a(2),-a(3)]+aa{y},[a(1),-a(2),-a(3)]+aa{y},[-a(3),-a(2),-a(1)]+aa{y},[a(3),-a(2),a(1)]+aa{y}... % group C2
           ,[-a(1),-a(3),-a(2)]+aa{y},[-a(1),a(3),a(2)]+aa{y},[a(2),a(1),-a(3)]+aa{y},[-a(2),-a(1),-a(3)]+aa{y},[-a(1),-a(2),a(3)]+aa{y}...
           ,[a(3),a(1),a(2)]+aa{y},[-a(3),-a(1),a(2)]+aa{y},[a(3),-a(1),-a(2)]+aa{y},[-a(3),a(1),-a(2)]+aa{y}...% C3
           ,[a(2),a(3),a(1)]+aa{y},[-a(2),a(3),-a(1)]+aa{y},[a(2),-a(3),-a(1)]+aa{y},[-a(2),-a(3),a(1)]+aa{y}... 
           ,[a(2),-a(1),a(3)]+aa{y},[-a(2),a(1),a(3)]+aa{y},[a(1),a(3),-a(2)]+aa{y}... % C4
           ,[a(1),-a(3),a(2)]+aa{y},[-a(3),a(2),a(1)]+aa{y},[a(3),a(2),-a(1)]+aa{y}}; 
    nX(i,:)=grou{z};
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%% find atom location%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(atom) 
    if nX(i,1) > coe+1
        nX(i,1) = nX(i,1)-coe-1;
    elseif nX(i,1) < 0
        nX(i,1) = nX(i,1)+coe+1;
    end
    if nX(i,2) > coe+1 
        nX(i,2) = nX(i,2)-coe-1;
    elseif nX(i,2) < 0
        nX(i,2) = nX(i,2)+coe+1;
    end
    if nX(i,3) > coe+1
        nX(i,3) = nX(i,3)-coe-1;
    elseif nX(i,3) < 0 
        nX(i,3) = nX(i,3)+coe+1;
    end
    if nX(i,1) == 0
        nX(i,1) = nX(i,1)+coe+1;
        end
    if nX(i,2) == 0 
        nX(i,2) = nX(i,2)+coe+1;
    end
    if nX(i,3) == 0
        nX(i,3) = nX(i,3)+coe+1;
    end
end
%% group compare%%%%
%%%%%%%%%%%%%%%%%%%%
for i=1:length(nX)
    for j=1:length(X)
        if nX(i,:)==X(j,:)
            sollt=[sollt;[i,j]];
        end
    end
end
soll{end+1}=sollt;
end
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
    same(end+1)=length(atom)-sum(gdt{i}); %how many self-cycle orbit
end
 E={};
 for i=1:length(gdt)
     e=zeros(1,length(grou));
     a=nonzeros(gdt{i});
    for j=1:length(a)
        e(a(j))=e(a(j))+1;
    end
        e(1)=same(i);
        E{end+1}=e;
 end
 all=prod(1:length(atom))/prod(1:length(find(atom==1)))/prod(1:length(find(atom==-1)));
 P=[];
 for i=1
    p=[];
    n=[];
    coef=[];
    de1=[];
    de2=[];
    part=0;
    p=find(E{i}~=0);
    n=E{i}(p);
    lll={};
    kkk={};
    for j=1:length(p)
        coef(end+1)=n(j); % how many atom needed in different orbits
        de1(end+1)=floor(length(find(atom==1))/p(j));
        de2(end+1)=floor(length(find(atom==-1))/p(j));
    end
    for ll=1:length(de1)-1
        lll{end+1}=[1:de1(ll)]*p(ll);
    end
        lll{end+1}=[0:de1(end)]*p(end);
    for kk=1:length(de2)-1
        kkk{end+1}=[1:de2(kk)]*p(kk);
    end
        kkk{end+1}=[0:de2(end)]*p(end);
    d1=length(lll)-1;
    d2=length(kkk)-1;
    ss=[];
    sk=[];
    sss={};
    ssk={};
    for t=1:length(lll{end})
        ss(end+1)=lll{end}(t);
    end
    for t=1:length(kkk{end})
        sk(end+1)=kkk{end}(t);
    end
    while d1>0;
        tt=length(ss);
        for s=1:length(lll{d1})
            for t=1:tt
                ss(end+1)=lll{d1}(s)+ss(t);
            end
        end
        d1=d1-1;
    end
    while d2>0;
        tt=length(sk);
        for s=1:length(kkk{d2})
            for t=1:tt
                sk(end+1)=kkk{d2}(s)+sk(t);
            end
        end
        d2=d2-1;
    end
    ss=find(ss==de1(1));
    sk=find(sk==de2(1));
    for st=1:length(ss)
    if length(de1)==2
        sss{end+1}=[mod(ss(st),de1(2)+1),mod(ceil(ss(st)/(de1(2)+1)),de1(1)+1)];
    elseif length(de1)==3
        sss{end+1}=[mod(ss(st),de1(3)+1),mod(ceil(ss(st)/(de1(3)+1)),de1(2)+1),mod(ceil(ss(st)/(de1(2)+1)/(de1(3)+1)),de1(1)+1)];
    end
    end
    
    for st=1:length(sk)
    if length(de2)==2
        ssk{end+1}=[mod(sk(st),de2(2)+1),mod(ceil(sk(st)/(de2(2)+1)),de2(1)+1)];
    elseif length(de2)==3
        ssk{end+1}=[mod(sk(st),de2(3)+1),mod(ceil(sk(st)/(de2(3)+1)),de2(2)+1),mod(ceil(sk(st)/(de2(2)+1)/(de2(3)+1)),de2(1)+1)];
    end
    end
    sc1={};
    sc2={};
    for st=1:length(sss)
        for stt=1:length(ssk)
            sss{st}(find(sss{st}==0))=de1(end)+1;
            ssk{stt}(find(ssk{stt}==0))=de2(end)+1;
            sss{st}=sss{st}-1;
            ssk{stt}=ssk{stt}-1;
            if sss{st}+ssk{stt}==coef(end:-1:1)
                sc1{end+1}=sss{st};
                sc2{end+1}=ssk{stt};
            end
        end
    end
    scoef=coef(end:-1:1);
    if length(sc1)==0
        part=0;
    else
        part=1;
    end
    for st=1:length(sc1)
        for stt=1:length(sc1{st})
            part=part*prod(1:scoef(stt))/prod(1:sc1{st}(stt))/prod(1:sc2{st}(stt));
                        P(end+1)=part;
        end
    end
end  



