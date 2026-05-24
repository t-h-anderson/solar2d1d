

function [G,Isort]=sortmat(Gu,Du)
%
% sort the columns of Gu and Du in order of
% the imaginary part of Du
%
[tmp,Isort]=sort(imag(Du),'descend');
%D=Du(Isort);
G=Gu(:,Isort);

end

