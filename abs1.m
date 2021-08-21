function output=abs1(arg1,arg2)
if nargin==1
    output=abs(arg1);
else
    if arg2~=0
        output=sign(arg2);
    else
        output=1;
    end
end