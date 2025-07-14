
%%% U 2G %% 



function ret = algoF(gi, agi)

  global W
  global allGroups
  global theta % the lower bound on the mutual trust among agents
  global gamma_tau
  global groups
  global tau

  global SUM_RANKING;
  global SUM_VOTING;

  %printf("algoF(%d, %d)... \n", gi, agi);

  % to store user --> trust(gi, user)
  allAgents = cell2mat(groups(1,gi));

  if(numel(allAgents)>0)  

    g2t = zeros(columns(allAgents)+1, 2);
    g2t(:, 1) = union(allAgents, [agi])';

    for i=1:rows(g2t)
      g2t(i,2) = g2u_tau(gi, g2t(i,1));
    endfor

    g2t_sorted = sortrows(g2t, -2);

    % to store user -- > trust(user, agi)
    u2t = zeros(columns(allAgents), 2);
    u2t(:, 1) = allAgents;

    for i=1:rows(u2t)
      u2t(i,2) = tau(u2t(i,1),agi);
    endfor

    %compute voting
    voting=0;
    nv=0;
    for i=1:rows(u2t)
      if(u2t(i,2)>=gamma_tau)
	nv++;
      endif
    endfor

    if(nv>=floor(columns(allAgents))/2)
      voting =1;
    else 
      voting = 0;
    endif

    % now decide whether to admit the agent 
    result = 0;	  
    if(columns(allAgents)>=W) % ranking
      SUM_RANKING+=1;
      if(g2t_sorted(rows(g2t_sorted),1)==agi) % agi is the worst in the ranking
	result = 0; % reject the request
	newAgents = allAgents;
      else % 
	agentOut = g2t_sorted(rows(g2t_sorted),1);
	newAgents = setdiff(allAgents, [agentOut]);
	newAgents = union(newAgents, [agi]);
	groups(1,gi) = mat2cell(newAgents, 1, columns(newAgents)); % put the new set of agents into the correspondent entry
	result = 1; % accept the request
      endif
    else % voting
      SUM_VOTING+=1;
      if(voting==1)
	result =1; % accept the request
	newAgents = union(allAgents, [agi]);
	groups(1,gi) = mat2cell(newAgents, 1, columns(newAgents)); % put the new set of agents into the correspondent entry
   %   printf("algoF(%d,%d), voting has accepted the request of %d", agi, agi, agi);
      else
	result =0;
	newAgents = allAgents;
    %  printf("algoF(%d,%d), voting has rejected the request of %d", gi, agi, agi);
      endif
    endif
  
  else
    printf("WARN: empty group %d", gi);
    result = 1;
    groups(1,gi) = mat2cell([agi], 1, 1); % put the new set of agents into the correspondent entry
  endif

 % printf("\t algoF(%d, %d): Result: %d, new set of agents: %s \n", gi, agi, result, num2str(newAgents));
  ret = result;
endfunction
