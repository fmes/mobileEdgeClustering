
function algoA(agi)
  global W
  global M
  global allGroups
  global theta % the lower bound on the mutual trust among agents
  global groups
  global alpha

  try

  %printf("algoA(%d)..\n", agi);

  allgr = grps(agi);
  newgr = setdiff(allGroups, allgr); % groups not containing agi

  % select W new groups from the set newgr  
  if(numel(newgr)>0)
    rndindeces = unique(randi([1,columns(newgr)], [1,W]));

    % all the candidate groups
    if(numel(allgr)==0)
      allcand=newgr(rndindeces); % column or row vector?
    else  
      allcand=union(newgr(rndindeces), allgr); % column or row vector?
    endif
    %printf("algoA(%d): All the candidates: %s\n", agi, num2str(allcand));

    % vector of agent and trust vs other groups
    % first column the id of the single candidate group
    % second column the trust of the agent vs the group
    a2t = zeros(rows(allcand'), 2);
    % first column
    a2t(:,1) = allcand';

    % second column
    for i=1:rows(a2t)
      a2t(i,2) = u2g_tau(agi, a2t(i,1));
    endfor

    % sort in decreasing order
    sorted_a2t = sortrows(a2t, -2);

    %printf("algoA(%d): All the other groups ordered by trust: %s\n", agi, num2str(sorted_a2t));

    % then select the best W groups
    nmax = 0;
    %goodCand = setdiff(sorted_a2t(:,1)', allgr);
    goodCand = sorted_a2t(:,1);
    for i = 1:rows(goodCand)
      if(sorted_a2t(i,2)<theta)
	break;
      else
	nmax = i;
      endif
    endfor

    goodCand = goodCand(1:(min(M, nmax)),1); % at most W candidates which are not already in allgr and for which tau(agi, cand)>=theta, a row vector

    c=0;
    joined = [];
    %printf("algoA(%d): Best candidates (MAX W): %s\n", agi, num2str(goodCand))
    for j=1:rows(goodCand)
      selectedG = goodCand(j,1); %request to join with the group
      if(!any(allgr == selectedG))
	  if(algoF(selectedG, agi)==1)% ok to join
	    c=c+1;
	    joined(1,c) = selectedG;
	  endif
      else
	c = c+1;
	joined(1,c) = selectedG;
      endif
    endfor

    newsize = numel(union(joined,allgr));
    c = newsize-M;

    % printf("\t algoA(%d), joined %d new groups: %s \n", agi, c, num2str(joined));

    if(c>0)
      % groups to leave
      toLeave = setdiff(allgr, joined);
      for j=1:columns(toLeave)
	if(c==0)
	  break;
	else
	  c--;
	endif
	grcell = groups(1,toLeave(1,j));
	gr = cell2mat(grcell);
	newgr = setdiff(gr, [agi]);
    %  printf("algoA(%d): leaving group %d, new set of agents: %s\n", agi, toLeave(1,j), num2str(newgr));
	groups(1,toLeave(1,j)) = mat2cell(newgr, rows(newgr), columns(newgr));
      endfor
    endif
  endif

  catch err
    error(err.message, err.identifier);
    printf("ERROR-LINE: %d\n", err.stack.line);
  end_try_catch

endfunction
