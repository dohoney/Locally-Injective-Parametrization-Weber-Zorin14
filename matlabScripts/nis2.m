%mn_border =  (m_border-min(min(m_border)))/(max(max(m_border))-min(min(m_border)));
%m_border = [ m_border(:,1) m_border(:,2) ];
stage2;
uiwait(stage2);
%E = sort([m_faces(:,1) m_faces(:,2); m_faces(:,2) m_faces(:,3); m_faces(:,3) m_faces(:,1)]')';
%[u,m,n] = unique(E,'rows');
%counts = accumarray(n(:), 1);
%O = u(counts==1,:);
%%%%%%%
%OO=O;
%bp=zeros(length(OO),2);
%bp(1,1) = m_points(OO(1,1),1);
%bp(1,2) = m_points(OO(1,1),2);
%bp(2,1) = m_points(OO(1,2),1);
%bp(2,2) = m_points(OO(1,2),2);
%index=3;
%first = OO(1,1);
%cur = OO(1,2);
%OO(1,2) = -1;
%while cur~=first
%    for i=1:length(OO)
%        if cur==OO(i,1)
%            bp(index,1) = m_points(OO(i,2),1);
%            bp(index,2) = m_points(OO(i,2),2);
%            index=index+1;
%            cur = OO(i,2);
%            OO(i,2)= -1;
%            break;
%        end
%    end
%    for i=1:length(OO)
%        if cur==OO(i,2)
%            bp(index,1) = m_points(OO(i,1),1);
%            bp(index,2) = m_points(OO(i,1),2);
%            index=index+1;
%            cur = OO(i,1);
%            OO(i,1)= -1;
%            break;
%        end
%    end
%end
    
%%%%%%%%%
stage3;
uiwait(stage3);