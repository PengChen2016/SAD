function [ b ] = check_Nagaoka( a )
%根据螺线管线圈2R/l的值,查表得长冈系数k，用来矫正理想螺线管电感。
%0<2R/l<10
A=textread('Nagaoka.txt');
[n,m]=size(A);
A1=A(:,1);
A2=A(:,2);

%按区间查表
if(a<A1(1)|a>=A1(n))
   b=0;
else
    for i=1:n
        if(a>=A1(i)&a<A1(i+1))
            b=A2(i)+(A2(i+1)-A2(i))/(A1(i+1)-A1(i))*(a-A1(i));
            break;
        end
    end
end

