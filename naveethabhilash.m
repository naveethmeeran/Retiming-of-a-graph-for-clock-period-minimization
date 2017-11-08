%This script retimes a graph for clock period minimization
clc
clear all;
%reading inputs from file
a=load('input.txt');
%assigning inputs to suitable variables
v=a(1,1);
e=a(1,2);
b=[0];
for i=2:v+1
    b=[b,a(i,2)];
end
t=b(1,2:v+1);
tmax=max(t);
T=t;
for i=1:v-1
    t=[t;T];
end
M=tmax*v;
c=a(v+2:v+e+1,:);
graph=Inf*ones(v,v);
for i=1:e
    graph(c(i,1),c(i,2))=c(i,3);
end
%Computing graphdash from graph
graphdash=zeros(v,v);
for i=1:v
    for j=1:v
        graphdash(i,j)=M*graph(i,j)-T(i);
    end
end
%Apply floyd wallace to graphdash
for i = 1:v
	for j=1:v
  	dist(i,j) = graphdash(i,j);
end
end

for k = 1:v
	for i=1:v
		for j=1:v
                if  (dist(i,k) + dist(k,j) < dist(i,j))
                    dist(i,j) = dist(i,k) + dist(k,j);
            end
end
end
end
c=zeros(1,v);
%computing w and d for graphdash
W=diag(c);
W=ceil(dist./M);
D=diag(T);
D=W.*M-dist+t;
%Computing retiming equations
c=[0,0,0];
for i=1:v
    for j=1:v
        if graphdash(i,j)~=Inf
          c=[c;i j W(i,j)];
        end
    end
end
d=[0 0 0];
for i=1:v
    for j=1:v
        if D(i,j)>tmax
            d=[d;i j W(i,j)-1];
        end
    end
end
%Eliminating unrequired retiming equations
[f,~]=size(d);
for i=2:e+1
    var=0;
     for j=2:f
        if ((c(i,1)==d(j,1))&&(c(i,2)==d(j,2)))
            var=var+1;
            
        end
     end
    if(var==0)
        d=[d;c(i,:)];
    end
end
%computing cnstraint graph
[f,~]=size(d); 
congraph=Inf*ones(v,v);
for i=2:f
    congraph(d(i,2),d(i,1))=d(i,3);
end

r1=Inf*ones(v+1,v+1);
r1(1:v,1:v)=congraph;
z=zeros(1,v);
r1(5,:)=[z,Inf];
%computing retiming values using floyd wallace 
for i = 1:v+1
	for j=1:v+1
  	R(i,j)=r1(i,j);
end
end

for k = 1:v+1
	for i=1:v+1
		for j=1:v+1
                if  (R(i,k) + R(k,j) < R(i,j))
                    R(i,j) = R(i,k) + R(k,j);
                end
        end
    end
end
%retime values are stored in variable retime
retime=R(v+1,1:v);
%computing the retimed graph
retimegraph=graph;
for i=1:v
    if retime(i)<0
       retimegraph(i,:)=retimegraph(i,:)-retime(i);
        retimegraph(:,i)=retimegraph(:,i)+retime(i);
        
    else 
        retimegraph(:,i)=retimegraph(:,i)-retime(i);
        retimegraph(i,:)=retimegraph(i,:)+retime(i);
        end
end
%computing retimegraphdash
retimegraphdash=zeros(v,v);
for i=1:v
    for j=1:v
        retimegraphdash(i,j)=M*retimegraph(i,j)-T(i);
    end
end

%computing w and d matrix for retimed graph using floyd wallace
for i = 1:v
	for j=1:v
  	NR(i,j) = retimegraphdash(i,j);
end
end

for k = 1:v
	for i=1:v
		for j=1:v
                if  (NR(i,k) + NR(k,j) < NR(i,j))
                    NR(i,j) = NR(i,k) + NR(k,j);
            end
end
end
end

W1=ceil(NR./M);
for i=1:v
    W1(i,i)=0;
end
D1=W1.*M-NR+t;
for i=1:v
    D1(i,i)=t(i,i);
end
%Finding the critical path of both input graph and retimed graph
q=[0];
r=[0];
for i=1:v
    for j=1:v
        if W(i,j)==0
        q=[q,D(i,j)];
        end
        if W1(i,j)==0
           r=[r,D1(i,j)];
        end
    end
end
s=max(q);
u=max(r);
%displaying critical path of graph and retimed graph for comparison
graph
retime
retimegraph
fprintf('The critical path delay of the given graph is %d\n',s);
fprintf('The critical path delay of the retimed graph is %d\n',u);

    