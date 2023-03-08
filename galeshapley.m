%---------GALE - SHAPLEY ALGORITHM (Men propose) --------%
function stablematch = galeshapley(N, M, men_pref, women_pref)

men_free = zeros(M,1);
women_suitor = zeros(N,M);
women_partner = zeros(N,1);
rank = zeros(N,M);

% zero_count=zeros(N,1);
% for i=1:N
%     for j=1:M
%         if women_pref(i,j)==0
%             zero_count(i,1)=zero_count(i,1)+1;
%         end
%     end
% end

for i = 1:N
    for j = 1:M
        for k = 1:M
            if(women_pref(i,k) == j)
                rank(i,j) = k;
            end
        end
    end
end
for i=1:M
    if (max(men_pref(i,:))==0)
        men_free(i,1)=1;
    end
end
for i=1:N
    if (max(women_pref(i,:))==0)
        women_partner(i,1)=9999;
    end
end
y=size(nonzeros(women_partner(:,1)),1);
x=size(nonzeros(men_free(:,1)),1);
counter=1;
while (min(women_partner) == 0 &&  counter<=(N-y)*(M-x))
    for i = 1:M
        if (men_free(i,1) == 0)
            next = find(men_pref(i,:) > 0, 1);
            if rank(men_pref(i,next),i)~=0
                women_suitor(men_pref(i,next),i) = i;
            end
            men_pref(i,next) = 0;
        end
    end
    for i = 1:N
        for j = 1:M
            if(women_suitor(i,j) ~= 0)
                if(women_partner(i,1) == 0 && rank(i,j) ~= 0)
                    women_partner(i,1) = women_suitor(i,j);
                    men_free(j,1) = 1;
                end
                if(women_partner(i,1) ~= 0)
                    if(rank(i,women_suitor(i,j)) < rank(i,women_partner(i,1)))
                        men_free(women_partner(i,1),1) = 0;
                        women_partner(i,1) = women_suitor(i,j);
                        men_free(j,1) = 1;
                    end
                end
            end
        end
    end
    if (min(men_free(:,1)==1))
        break;
    end
    counter=counter+1;
end

stablematch = women_partner;

