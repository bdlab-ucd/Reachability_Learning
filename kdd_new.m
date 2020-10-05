function [ X, Yt ] = kdd_new( fullM, Mdata,k ,limit)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%fully known data as ground truth
%input partially observed data
% rank k

%%%% Initialize variables
visM=full(spones(Mdata));
numadded=0;
infectedC=[];
infectedR=[];

[nr nc]=size(visM);
colsums=sum(visM,1);
rowsums=sum(visM,2);

rqueue=[];
cqueue=[];

X = zeros(nr,k);
Yt = zeros(k,nc);

newk=k+1;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 1: Create an order on the rows and columns %%%%%
structRC=[[colsums';rowsums],[ones(nc,1);zeros(nr,1)]];
[vals inds]=sortrows(structRC,-1);
rqueue=inds-nc;
cqueue=inds;
rqueue(find(vals(:,2)==1))=0;
cqueue(find(vals(:,2)==0))=0;
orderR=nonzeros(rqueue);
orderC=nonzeros(cqueue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2: Pick the first kxk submatrix A to initialize the infection and solving%%%%%
infectedR=rqueue(1:k);
infectedC=cqueue(1:k);
i=k;
while(nnz(rqueue(1:i))<k)
    i=i+1;
end
infectedR=nonzeros(rqueue(1:i));
j=k;
while(nnz(cqueue(1:j))<k)
    j=j+1;
end
infectedC=nonzeros(cqueue(1:j));


mclique=visM(infectedR,infectedC);

if((numadded+numel(mclique)-nnz(mclique))>limit)
    return
end
numadded=numadded+numel(mclique)-nnz(mclique);
visM(infectedR,infectedC)=1;
Mdata(infectedR,infectedC)=fullM(infectedR,infectedC);
for q=1:k
    indF=find(cqueue==infectedC(q));
    cqueue(indF)=0;
end
for q=1:k
    indF=find(rqueue==infectedR(q));
    rqueue(indF)=0;
end
ij=min(i,j);
rqueue(1:ij)=[];
cqueue(1:ij)=[];

if(size(infectedC,1)==1)
    infectedC=infectedC';
end
if(size(infectedR,1)==1)
    infectedR=infectedR';
end
orderC=[infectedC;nonzeros(cqueue)];
orderR=[infectedR;nonzeros(rqueue)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step 3: Postprocessing step -- alterations on the order %%%%

 
%     for i=(newk+1):length(orderR)
%         node=orderR(i);
%         row=visM(node,:);
%         if(sum(row)==0)
%             continue
%         end
%         if(sum(row)<newk )
%             loc=[];
%             neigh=find(row);
%             
%             if(length(setdiff(neigh,infectedC))==0)
%                 insertL=numel(rqueue)-1;
%             else
%                 neigh=setdiff(neigh,infectedC);
%                 for(r=1:length(neigh))
%                     loc(r)=find(cqueue==neigh(r));
%                 end
%                 insertL=max(loc);
%             end
%             
%             rqueue=[rqueue(1:insertL), node, ...
%                 rqueue((insertL+1):length(rqueue))];
%             cqueue=[cqueue(1:insertL), 0, ...
%                 cqueue((insertL+1):length(cqueue))];
%             whereinQ=find(rqueue==node);
%             whereinQ=setdiff(whereinQ,insertL+1);
%             rqueue(whereinQ)=[];
%             if(cqueue(whereinQ)>0)
%                 return
%             end
%             cqueue(whereinQ)=[];
%             nnz(cqueue) ;
%             
%         else
%             loc=[];
%             neigh=find(row);
%             for r=1:length(neigh)
%                 l=find(cqueue==neigh(r));
%                 if(length(l)==0)
%                     loc(r)=0;
%                 else
%                     loc(r)=l;
%                 end
%             end
%             sloc=sort(loc);
%             insertL=sloc(newk);
%             
%             rqueue=[rqueue(1:insertL), node, ...
%                 rqueue((insertL+1):length(rqueue))];
%             cqueue=[cqueue(1:insertL), 0, ...
%                 cqueue((insertL+1):length(cqueue))];
%             whereinQ=find(rqueue==node);
%             whereinQ=setdiff(whereinQ,insertL+1);
%             rqueue(whereinQ)=[];
%             cqueue(whereinQ)=[];
%         end
%     end
%     %length(rqueue)
%     %length(cqueue)
%     for j=(newk+1):length(orderC)
%         node=orderC(j);
%         col=visM(:,node);
%         if(sum(col)==0)
%             continue
%         end
%         if(sum(col)<newk && sum(col)>0)
%             loc=[];
%             neigh=find(col);
%             
%             if length(setdiff(neigh,infectedR))==0
%                 insertL=0;
%             else
%                 neigh=setdiff(neigh,infectedR);
%                 for r=1:length(neigh)
%                     loc(r)=find(rqueue==neigh(r));
%                 end
%                 insertL=max(loc);
%             end
%             % insertL
%             cqueue=[cqueue(1:insertL), node, ...
%                 cqueue((insertL+1):length(cqueue))];
%             rqueue=[rqueue(1:insertL), 0, ...
%                 rqueue((insertL+1):length(rqueue))];
%             whereinQ=find(cqueue==node);
%             whereinQ=setdiff(whereinQ,insertL+1);
%             rqueue(whereinQ)=[];
%             cqueue(whereinQ)=[];
%         else
%             loc=[];
%             neigh=find(col);
%             %neigh=setdiff(neigh,infectedR);
%             for r=1:length(neigh)
%                 l=find(rqueue==neigh(r));
%                 if(length(l)==0)
%                     loc(r)=0;
%                 else
%                     loc(r)=find(rqueue==neigh(r));
%                 end
%             end
%             sloc=sort(loc);
%             insertL=sloc(newk);
%             
%             cqueue=[cqueue(1:insertL), node, ...
%                 cqueue((insertL+1):length(cqueue))];
%             rqueue=[rqueue(1:insertL), 0, ...
%                 rqueue((insertL+1):length(rqueue))];
%             whereinQ=find(cqueue==node);
%             whereinQ=setdiff(whereinQ,insertL+1);
%             rqueue(whereinQ)=[];
%             cqueue(whereinQ)=[];
%         end
%     end
%     
%     orderC=[infectedC;nonzeros(cqueue)];
%     orderR=[infectedR;nonzeros(rqueue)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


infectedR=sort(infectedR);
infectedC=sort(infectedC);

initR=infectedR;
initC=infectedC;

cqueue=reshape(cqueue,1,numel(cqueue));
rqueue=reshape(rqueue,1,numel(rqueue));

X(infectedR(1:k),1:k)=eye(k);

%%%%% Solve the initial columns of Y %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:length(infectedC)
    col=infectedC(g);
    s = intersect(find(visM(:,col)), infectedR);
    %gensUsed=[gensUsed,zeros(1,length(s))];
    if det(X(s,:)' * X(s,:)) == 0
        Yt(:,col) = pinv(X(s,:)) * nonzeros(Mdata(infectedR,j));
    else
        Yt(:,col) = X(s,:) \Mdata(infectedR,col) ;
        [~, wid] = lastwarn;
        if (strcmp(wid, 'MATLAB:singularMatrix') || strcmp(wid, 'MATLAB:nearlySingularMatrix')|| strcmp(wid, 'MATLAB:rankDeficientMatrix'))
            fprintf('INIT ENTRIES NOT GOOD INIT\n')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


z=1;
while(rqueue(z)==0)
    z=z+1;
end
w=1;
while(cqueue(w)==0)
    w=w+1;
end

next=infectedR(1);
row=next;

visM(infectedR,infectedC)=1;
Mdata(infectedR,infectedC)=fullM(infectedR,infectedC);

orderC=[infectedC;nonzeros(cqueue)];
orderR=[infectedR;nonzeros(rqueue)];

pointer=2*newk+1;
if(numel(find(fullM(infectedR,infectedC)==0))>0)
    full(fullM(infectedR,infectedC))
    full(visM(infectedR,infectedC))
    fprintf('err in begining, zeros')
    return
end
if(numel(find(isnan(X)))>0 || numel(find(isnan(Yt)))>0)
    full(fullM(infectedR,infectedC))
    X(infectedR,:)
    Yt(:,infectedC)
    fprintf('err in begining')
    return
end

rows=1;
cols=1;

condsC(initC)=1;
condsR(initR)=1;

if(length(infectedC)==nc)
    cols=0;
end
if(length(infectedR)==nr)
    rows=0;
end
next=-1;

addedR=zeros(1,nr);
addedC=zeros(1,nc);

initCond=cond(Yt(:,infectedC))






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Infection Step!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while(rows || cols)
    orderC=[reshape(infectedC,numel(infectedC),1);nonzeros(cqueue)];
    orderR=[reshape(infectedR,numel(infectedR),1);nonzeros(rqueue)];
    
    %%% if the next element in the order is a row
    if(mod(pointer,2)==1)
        if(rows)
            next=rqueue(1);
            if(next==0)
                pointer=pointer+1;
                rqueue(1)=[];
                if(isempty(rqueue))
                    rows=0;
                end
                continue;
            end
            
            knowns=find(visM);
            
            %%% how many edges need to be added
            toAdd=max(0,(k-indF));
            
            if(toAdd>0)
                % if we've gone over the limit, try adding to the end
                if((numadded+toAdd)>limit && ~ismember(rqueue(1),movedR))
                    fprintf('limit reached')
                    rqueue=[rqueue,rqueue(1)];
                    movedR=[movedR,rqueue(1)];
                    rqueue(1)=[];
                    cqueue=[cqueue,0];
                    pointer=pointer+1;
                    continue
                elseif((numadded+toAdd)>limit)
                    rqueue(1)=[];
                    pointer=pointer+1;
                    if(isempty(rqueue))
                        rows=0;
                    end
                    continue
                end
                
            end
            
            neigh=find(visM(next, infectedC));
            cand=infectedC(find(~visM(next, infectedC)));
            colIs=[];
            
            if(toAdd>0)
                %%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%  from the candidates pick the most orthogonal
                if(numel(cand)>toAdd)
                    left=toAdd;
                    colIs=[];
                    bk=[];
                    bk=cand(b);
                    tempC=sort([infectedC(neigh),colIs]);
                    tempC=reshape(tempC,1,numel(tempC));
                    
                    mm=fullM(next,tempC);
                    e=mm/Yt(:,tempC);
                    me=e*Yt(:,tempC);
                    while(left>0)
                        mm=fullM(next,tempC);
                        mmn=[mm,fullM(knowns(randi(numel(knowns))))];
                        for i=1:min(numel(bk),100)
                            nvs=[nvs,norm(Yt(:,tempC)'*Yt(:,bk(i)),'fro')];
                        end
                        [v w2]=min(nvs);
                        colIs=[colIs,bk(w2)];
                        tempC=sort([tempC,bk(w2)]);
                        bk(w2)=[];
                        left=left-1;
                    end
                    
                else
                    colIs=cand;
                end
            end
            
            infectedR=reshape(infectedR,1,numel(infectedR));
            infectedR=union(infectedR,next);
            rqueue(1)=[];
            
            if(isempty(rqueue) || nnz(rqueue)==0)
                rows=0;
            end
            
            %%%%%%%%solve%%%%%%%%%%
            colIs=reshape(colIs,1,numel(colIs));
            infectedC=reshape(infectedC,1,numel(infectedC));
            s=union(infectedC(find(visM(next, infectedC))),sort(colIs));
            %genR(next)=max(genC(s))+1;
            
            row=next;
            s=sort(s);
            cand=setdiff(infectedC,s);
            
            X(row,:) = Mdata(row,s)/(Yt(:,s));
            
            neigh=find(visM(next,infectedC));
            cand=setdiff(infectedC,infectedC(neigh));
            
            if(numel(cand)>0 && numadded+1<=limit)
                
                ad=0;
                nEnd=0;
                tempM=fullM(next,s);
                e=(tempM)/Yt(:,s);
                e1=(tempM+randn(size(tempM))*std(tempM)*.1)/Yt(:,s);
                %
                tempMe=e1*Yt(:,s);
                nEnd=norm(e-e1,'fro')/norm(e,'fro');
                
                normsL=[];nvs=[];
                tempS=sort(reshape(s,1,numel(s)));
                cand=setdiff(infectedC,tempS);
                ad=0;
                noise=randn(1,numel(tempS)+1)'*std(tempM)*.1;
                maxVV=100;
                x=e';
                
                while(numel(cand)>0 && (maxVV>10 || nEnd>.01) && ad<10)
                    ns=[];
                    normsL=[];nvs=[];normsL2=[];normsL3=[];normsL1=[];
                    
                    A=Yt(:,tempS)';
                    B=inv(A'*A);
                    
                    if(numel(tempS)+1>numel(noise))
                        noise=[noise;randn()*std(tempM)*.1];
                    end
                    knowns=find(visM);
                    rp=randperm(numel(cand));
                    cand=cand(rp);
                    for i=1:min(100,numel(cand))
                        tempC=sort([tempS,cand(i)]);
                        v=Yt(:,cand(i))';
                        %b2=sherman2(v,A);
                        b2=ones(4,1);
                        A2=[A',v'];
                        mm=[fullM(next,s)';fullM(knowns(randi(numel(knowns),ad+1,1)))];
                        
                        est=(b2*A2*mm)';
                        t=est*A2;
                        mmn=mm+noise;
                        est2=(b2*A2*mmn)';
                        
                        normsL1(i)=norm(est-est2,'fro')/norm(noise,'fro');
                        
                    end
                    
                    [nEnd c]=min(normsL1);
                    numadded=numadded+toAdd;
                    tempS=sort([tempS,cand(c)]);
                    ad=ad+1;
                    cand=setdiff(infectedC,tempS);
                    visM(row,tempS)=1;
                    
                    X(next,:)=(fullM(row,tempS))/(Yt(:,tempS));
                end
                
                tempS=sort(tempS);
                Mdata(row,tempS)=fullM(row,tempS);
                X(next,:)=(fullM(row,tempS))/(Yt(:,tempS));
                
                numadded=numadded+max(1,ad);
                
                pointer=pointer+1;
                continue
                
            else
                if( ~addedR(row) )
                    rqueue(end+1)=row;
                    cqueue(end+1)=0;
                    infectedR(find(infectedR==next))=[];
                    addedR(row)=max(numel(s),1);
                    pointer=pointer+1;
                    %removed=removed+1;
                    continue
                end
            end
            
        end
        pointer=pointer+1;
        
    end
    
    
    if(mod(pointer,2)==0)
        if(cols)
            next=cqueue(1);
            
            if(next==0)
                pointer=pointer+1;
                cqueue(1)=[];
                if(isempty(cqueue))
                    cols=0;
                end
                continue;
            end
        end
        pointer=pointer+1;
    end
end
end
