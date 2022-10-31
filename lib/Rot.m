function bRi=Rot(R,F)
% inputs:
%    R   - the rotation matrix
%    F   - the coordinates before the rotation
% Output:
%    bRi  - the rotated coordinates by R
 N=length(F);
 bRi=F*0;
 for i=1:N
     bRi(i,:)=R*F(i,:)'; % dot product on each vertex
 end
end