% END functions from Information sciences

%% NEW FUNCTIONS %%

pkg load parallel
pkg load statistics
%pkg load u2g.m

global DEBUG=0

% algoritm to run 
global U2G=1
global MERIT=2
global ALGO=MERIT

%poisson, expected value of the number of feedbacks per user per tick
global lambda_F=20;
%poisson, expected value of the number of recommendations per user
%per tick
global lambda_R=20;
%stdev recommendation values
global delta_R = 0.05;

%Interactions
global F_LP_NU=0.2; % Low performance
global F_LP_SIGMA=0.1;
global F_MP_NU=0.6; % Medium performance
global F_MP_SIGMA=0.1;
global F_HP_NU=0.9; % High performance 
global F_HP_SIGMA=0.1;

global P_high_perf=0.2; 

%parameters relater to the disadvantage, best contrbutors etc
global BC_ts=0.5
global BR_ts=0.5 

% ``standard'' mean and standard deviation, if any 
global FEED_NU=F_HP_NU;
global FEED_SIGMA=F_HP_SIGMA;

% standard deviation for generated recc
global stdev_R=0.2;

global N=100; % number of agents 
global K=10; %max number of groups
global W=12; % max number of agents per groups
global M=3; % max number of groups to analyze for a single agent during the execution of the algorithm 
global NF=5 % max number of friends 

% behaviour of the N agents 
global CH

global dVector = zeros(1, N); 
global dPerformance = zeros(1, N);

global NI=5; % number of interactions of agents 

global alpha = 0.5; % ?? 
global theta = 0.2; % ?? 
global gamma_tau = 0.5; % ?? 

global groups = cell(1,K); % set of groups / a tuple containing the group composition
% NEW 
global friends = zeros(N, N); 
global allGroups = [1:K]; % array of indeces of groups 
global tau = zeros(N, N); % mutual trust matrix 

global NSTEPS=10; % steps of simulations 
global ALGON=100; % ??? 

global SUM_VOTING=0; % ??? 
global SUM_RANKING=0; % ??? 
global CURRENT_STEP=0; % ???

%% INIT FUNCTIONS MERITOCRATIC / Mobile edge computing 

% NF RANDOM FRIENDS FOR ALL THE AGENTS 
function initFriends(agi)
	global friends % friends matrix
	global N
	global NF
	global DEBUG
	for i=1:N
		findex=unique(randi([1,N], 1, NF));
		findex(find(findex==i))=[];
		friends(i,findex)=1;
		for j=1:length(findex)
			friends(findex(j), i)=1; 
		endfor 
		if(DEBUG)
			printf("\n Friends of agent %d: ", i);
			friends(i,:)
		endif
		if(find(findex==i))
			printf("\n ERROR! agent %d friend of itself! ", i); 
			exit
		endif
		%friends(1,i) 
	endfor
	% printf("DEBUG: friends matrix")
	% friends
endfunction

function addFriend(agi, agj)
	global friends
	friends(agi, agj)=1
	friends(agj, agi)=1
endfunction 

function removeFriend(agi, agj)
	global friends
	friends(agi, agj)=0
	friends(agj, agi)=0
endfunction 

% RETURN THE SET OF FRIENDS FOR A GIVEN AGENT AGI 
function f=F(agi)
	global friends
	f=find(friends(agi,:)==1);
endfunction

function ret=agentsInGroup(gi)
		global groups
		ret=cell2mat(groups(1,gi)); 
endfunction 

% UNION OF ALL AGENTS FROM ALL AGENTS' GROUPS
function allAg=allAgentsFromGroups(agi)
	% union of all the agents in the groups
	% of the agent agi
	% looking for all the agents becoming from all the groups of agent agi
	global groups
	ret=[];
	allGroupsAgi=grps(agi);
	for i=1:length(allGroupsAgi)
		gindex=allGroupsAgi(1,i);
		allAg=cell2mat(groups(1,gindex)); 
		ret=union(ret, allAg); 
	endfor
	allAg=ret;
endfunction 

% BEST CONTRIBUTORS 
function bestCon=BC(agi)
	global tau
	global BC_ts
	global DEBUG
	if(DEBUG)
		printf("BC_ts: %f", BC_ts);
	endif
	% extract from trust matrix
	% the set of the best contributors
	% considering thr threshold BC_ts		 
	%size(tau)
	trust_vector=tau(agi,:);
	ret=find(trust_vector>BC_ts);
	%printf("DEBUG: bestcon(%d)=", agi)
	%ret
	bestCon = ret;
endfunction

% TODO 
function bestRec=BR(agi)
	bestRec = []; 
endfunction 

% COMPUTE THE UNION OF FRIENDS(I) AND ALL THE AGENTS IN GROUPS(I)
function ret=friendsAndGroups(agi)
	ret=union(F(agi), allAgentsFromGroups(agi));
endfunction

% COMPUTE THE UNION OF BC AND BR SETS 
function ret=BCandBR(agi)
	ret=union(BC(agi), BR(agi));
endfunction

% additionalGroup 
% is used to simulate the disadvantage 
% with the addition of an additional cluster 
% of agents 
function r=ratio1Disadvantage(agi, additionalGroup=[])
	global DEBUG 

	if(DEBUG)
		printf("\n agi=%d, ratio1Disadvantage DEBUG..", agi)
	endif

	A=length(setdiff(union(friendsAndGroups(agi), additionalGroup), BCandBR(agi)));
	B=length(union(friendsAndGroups(agi), additionalGroup));

	ret=0;

	if(B>0 && A>0)
		ret = A/B;
	endif

	if(DEBUG)
		printf("\nratio1Disadvantage, lenght(A)=%d, lenght(B)=%d, ret=%f\n", A, B, ret);
	endif
	r=ret;
endfunction

% additionalGroup 
% is used to simulate the disadvantage 
% with the addition of an additional cluster 
% of agents 
function r=ratio2Disadvantage(agi, additionalGroup=[])
	global tau
	global DEBUG

	if(DEBUG)
		printf("\n agi=%d, ratio2Disadvantage DEBUG..", agi)
	endif

	a=setdiff(union(friendsAndGroups(agi), additionalGroup), BCandBR(agi)); % the best alternative 
	if(length(a)>0)
		tv=tau(agi, a); % vector of trust values for all friends and agent groups
		mxt=max(tv); % maximum value of trust (alternative b^0), i.e. trust(agi, b^0)
	else
		mxt=0;
	endif

	% best agents that are not in the set of friends and groups
	s=setdiff(BCandBR(agi), union(friendsAndGroups(agi), additionalGroup));

	ret=0;  % default 
	sm=0; 

	if(length(s)>0)
		trustv=tau(agi,:)-mxt;  		
		%trustv=tau(agi,:)-mxt;  
		size(trustv); 
		sm=sum(trustv(s));
		ret=sm/length(s);
	endif

	if(DEBUG)
		printf("\nratio2Disadvantage, lenght(s)=%d, sm=%d, ret=%f\n", length(s), sm, ret);
	endif

	r=ret; 
endfunction

function ret=disadvantage(agi, additionalGroup=[])
	r=(ratio1Disadvantage(agi, additionalGroup) + ratio2Disadvantage(agi, additionalGroup))/2; 
	assert(r<=1.0);
	assert(r>=0);
	ret=r; 
endfunction

% compute the disadvantage of the nodes 
function computeAllD()
	global N
	global dVector
	for i=1:N
		dVector(i)=disadvantage(i); 
	endfor
endfunction

% compute the performance of all the nodes 
function computePerformance()
	global N
	global dVector
	global dPerformance
	dPerformance=1-dVector;
endfunction

% return a vector of K numbers
% each representing the average mutual 
% trust within a group 
% USED functions: avg_tau() 
function ret=avg_global_tau_comp()

  global K 

  try 

  avgtau = zeros(1,K);

  printf("avg_global_tau_comp()\n");
  
  for i=1:K
    avgtau(1,i) = avg_tau(i);
  endfor

  ret = avgtau;

  catch err
    printf("ERROR (%s): %s - @LINE %d\n", err.identifier, err.message, err.stack.line);
    exit;
  end_try_catch
endfunction

% Mutual avg trust within a group 
function ret=avg_tau(gindex)
  global W
  global groups
  global tau
	global N

  sm=0;
  c=0;
  allAgents = cell2mat(groups(1,gindex));
  for i=1:columns(allAgents)
    ag1 = allAgents(1,i);
    for j=i:columns(allAgents)
      ag2 = allAgents(1,j);
      if(ag1!=ag2)
				sm+=tau(ag1,ag2);
				sm+=tau(ag2,ag1);
				c+=2;
      endif
    endfor
  endfor
  if(c>0)
    ret=sm/c;
  else 
    ret = 0;
  endif
endfunction

% trace disadvantage 
function flush_disadvantage(step, file)
	global dVector; 
  printf("flush_disadvantage(%d, %s)\n", step, file);
  dlmwrite(file, [step dVector], "append", "on", "delimiter", " ", "precision", 3);
endfunction

% trace the performance 
function flush_performance_stats(step, file)
	global dPerformance;
	global N; 

  printf("flush_performance_stats(%d, %s)\n", step, file);
	st=statistics(dPerformance'); 
  dlmwrite(file, [step st(1:5,1)'], "append", "on", "delimiter", " ", "precision", 3);
endfunction

% trace on file the avg global tau within groups
function flush_avg_global_tau(step, file)
  printf("flush_avg_global_tau()", file);
	dlmwrite(file, [step avg_global_tau_comp], "append", "on", "delimiter", " ", "precision", 3);
endfunction

% trace the statistics about the mutual trust 
function flush_avg_global_tau_stat(step, file)
  printf("flush_avg_global_tau_stat(%d, %s)\n", step, file);
  vec = avg_global_tau_comp();
	%vec
  stat=statistics(vec'); % FM trasposto --> vettore colonna
  %minimum, first quartile, median, third quartile, maximum
%  dlmwrite(file, [step stat(1,1:5)], "append", "on", "delimiter", " ", "precision", 3);
   dlmwrite(file, [step stat(1:5,1)'], "append", "on", "delimiter", " ", "precision", 3);
endfunction

%compute the trust between an agent 
%and the agents of a certain group
% agi == agent index 
% gi = group index 
% SUM(tau(agi, j))/S for all j in gi; S=size(gi)
function ret=u2g_tau(agi, gi)
  global W
  global groups
  global tau

  sm=0;c=0;
  allAgents = cell2mat(groups(1,gi));
  if(columns(allAgents)==0)
    ret=-1;
  endif
  for i=1:columns(allAgents)
    if(allAgents(1,i)!=agi)
      sm+=tau(agi,allAgents(1,i));
      c++;
    endif
  endfor
  if(c>0)
    ret=sm/c;
  else
    ret = 0.5;
  endif
endfunction

%compute the trust between a group
% and a single agent 
% agi == agent index 
% gi = group index 
% SUM(tau(j, agi))/S for all j in gi; S=size(gi)
function ret=g2u_tau(gi, agi)
  global W
  global groups
  global tau

  sm=0;
  allAgents  = cell2mat(groups(1,gi));
  if(columns(allAgents)==0)
    ret=0.5;
  endif
  c=0;
  for i=1:columns(allAgents)
    if(allAgents(1,i)!=agi)
      sm+=tau(allAgents(1,i), agi);
      c++;
    endif
  endfor
  if(c>0)
    ret=sm/c;
  else
    ret = 0.5;
  endif 
endfunction


% Sample a feedback for low or high performance
% the feedback is intended for a service provided by j 
function feed = sampleFeedback(j)

  global CH % mapping of the behaviour of the agents, low performance / high performance 
	% mean and standard deviation to simulate the performance of agents 
  global F_HP_NU
  global F_LP_NU  
  global F_HP_SIGMA
  global F_LP_SIGMA
  global lambda_F

  if(CH(1,j)==0) % LP
    feed = min(1,max(0,normrnd(F_LP_NU, F_LP_SIGMA)));
  else
    feed = min(1,max(0,normrnd(F_HP_NU, F_HP_SIGMA)));
  endif

endfunction

% simulate a recommendation (low of high performance) 
% ??? 
function rec = sampleRecc(j)
  global CH
  global F_HP_NU
  global F_LP_NU  
  global F_HP_SIGMA
  global F_LP_SIGMA
  global lambda_F
  global FEED_NU
  global FEED_SIGMA

  if(CH(1,j)==0) % LP
    rec = min(1,max(0,normrnd(F_LP_NU, F_LP_SIGMA)));
  else
    rec = min(1,max(0,normrnd(F_HP_NU, F_HP_SIGMA)));
  endif
endfunction

%MAIN LOOP

% Return a vector containing all the groups of the agent 
% with index agi 
function ret=grps(agi)
  global groups
  global K 

	%groups
  tmp=[];
  c = 0;
  for i=1:length(groups)
    if(any(cell2mat(groups(i))==agi)) % agi in groups
			%printf("\n Agent %d in group %d", agi, i); 
      c=c+1;
      tmp(c)=i;
			%tmp, size(tmp);
    endif
  endfor
  ret=tmp;
endfunction

%%% MERITOCRATIC ACTIVE TASK 

function meritocraticActiveTask(agi)
	% compute th best set BCandBR
	agi 
	best=BCandBR(agi)
	frgr=friendsAndGroups(agi)
	% compute BCandBR - friendsAndGroups
	candidates=setdiff(best,frgr)
	for i=1:length(candidates)
		c=candidates(i)
		ret=askFriendship(agi, c)
		if(ret==1)
			removeWorstFriend(agi); 
		else
			% explore the addition of any cluster 
			% containing b
			% stop with the first cluster
			clustersC = grps(c); % groups of agent c 
			for k=1:length(clustersC)
				tmp=agentsInGroup(clustersC(k)) % agents in the group 
				d1=disadvantage(agi)
				d2=disadvantage(agi, tmp);
				if(d2<=d1) % try to join group clusters(k)
					ret=askJoinGroup(clustersC(k)); 
				endif
			endfor
		endif
	endfor
endfunction 

function ret=askJoinGroup(agi, gi)
	global groups

	% compute all the disadvantage 
	% of all the agent in the group 

	% compute all the disadvantage 
	% of all the agents in the group 
	% with the addition of the agent agi

	% voting of all the agents 
	% is simulated 
	% by the relation 
	% d2<=d1  

	% add or not the agent agi into the group 

	% TODO 
	ret=0;
endfunction


% agent agi asks a friendshipt to 
% agent c
function ret=askFriendship(agi, c)
	%currentD = disadvantage(c);
	best=BCandBR(c)
	% accept friendship request 	
	% will not increase the disadvantage
	% only if the agent agi 
	% is already in the set BCandBR(c)
	r=0; 
	if(any(find(best==agi))) % agi is already in the set best 
		% accept the friendship request 
		addFriend(agi, c);
		removeWorstFriend(c); 
		r=1; 
	endif	
	ret = r; 
endfunction 

function removeWorstFriend(agi)
	global tau; 
	fr=F(agi);
	[mn,i]=min(tau(agi, fr));
	% remove friend i 
	removeFriend(fgi, fr(i));
endfunction 

%%% MERITOCRATIC PASSIVE TASK 

function meritocraticPassiveTask(gi)
	% TODO 
endfunction 

function trace_groups(s, file)
  global K
  global groups
  global CH

  ch0 = 0;
  ch1 = 0;

  printf("trace_groups(%d, %s)..\n", s, file);
  for gi=1:K %K groups
    ch0 = 0;
    ch1 = 0;
    gr = cell2mat(groups(gi));
    for i=1:columns(gr)
      if(CH(1,gr(1,i))==0)
	ch0+=1;
      else
	ch1+=1;
      endif 
    endfor
    dlmwrite(file, [s gi cell2mat(groups(gi)) ch0 ch1], "delimiter", " ", "append", "on");
  endfor
endfunction

function trace_agents_groups(s, file)
  global N
  global groups
  global CH

  printf("trace_agents_groups(%d, %s)..\n", s, file);
  for agi=1:N %K groups
    dlmwrite(file, [s agi CH(agi) numel(grps(agi)) grps(agi)], "delimiter", " ", "append", "on");
  endfor
endfunction

function trace_stat_voting(s, file)
  global SUM_VOTING
  global SUM_RANKING
  global CURRENT_STEP

  printf("trace_stat_voting(%d, %s)..\n", s, file);
  dlmwrite(file, [ CURRENT_STEP SUM_RANKING SUM_VOTING ], "delimiter", " ", "append", "on");
endfunction

function initGroups()
	global groups 
	global W 
	global K 
	global N 

	maxW=min(W, round(N/K))
	printf("***1. Generating %d groups of %d random agents..***\n", K, maxW);
	agi=1;
	for gi=1:K %K groups
  	printf("Group %d <-- ", gi);
  	ic = 1;
  	while(ic<=maxW && agi <=N)
    	gr = cell2mat(groups(gi)); % vector representing group gi
    	if(!(columns(gr)>=maxW)) % group is full!!
      	if(!any(gr==agi)) % agent i NOT  in group gr
					gr=union(gr, agi); 
					groups(gi) = gr;
					assert(any(grps(agi)==gi));
					printf("%d, ", agi);
					ic+=1;
      	endif
    	endif    
    	agi+=1;
  	endwhile
  	% = min(round(rand(1,W)*N), ones(1,W)*N); %sampling of W agents per group
  	printf("\n");
	endfor
endfunction 

function markAgentBehaviours()
	global N 	
	global CH 
	global P_high_perf

	printf("***2. Marking agent behaviours (high perf, low perf, etc)...\n");
	if(exist("ch.txt", "file") && exist("trust.txt", "file"))
  	printf("Reading file ch.txt\n");
  	CH=dlmread("ch.txt", " ");
	else
  	CH = zeros(1,N);
  	for i=1:N
			x = rand(); 
			%printf("\n x(prob high perf)=%f", x); 
    	if(x<=P_high_perf)
      	CH(1,i) = 1; %High Performance
    	else
      	CH(1,i) = 0; %Low Performance
    	endif
  	endfor
  	dlmwrite("ch.txt", CH, "delimiter", " ", "precision", 3);
	endif
endfunction

function initTrustMatrix()
	global tau 
	global N 

	printf("***3. Initializing trust matrix (%d x %d)..\n", N, N);
	if(exist("trust.txt", "file") && exist("ch.txt", "file"))
  	printf("Reading file trust.txt\n");
  	tau=dlmread("trust.txt", " ");
	else
  	for i=1:N
    	if(mod(i,100)==0)
      	printf("%d..\n", i);
    	endif  
    	for j=i+1:N
				tau(i,j) = sampleFeedback(j); % feedback from i to j
				tau(j,i) = sampleFeedback(i); % feedback from j to i
    	endfor
  	endfor
  	dlmwrite("trust.txt", tau, "delimiter", " ", "precision", 3);
	endif
endfunction 

function test()

	%% TEST 

	global N 
	printf("\n *** TEST *** "); 
	%agi = randi([1,N])
	agi=20
		
	printf("\n Test BC for agent %d..", agi);
	bc = BC(agi);
	bc; % OK 

	% test BR 

	printf("\n Test BR for agent %d..", agi)
	br = BR(agi);
	br; % OK

	% test BC

	printf("\n Test BC for agent %d..", agi)
	bc = BC(agi);
	bc; % OK

	% test BCandBR
	printf("\n Test BC and BR for agent %d..", agi)
	bcbr = BCandBR(agi); 
	bcbr;  % OK 

	% test friends 

	printf("\n Test friends for agent %d..", agi)
	fr = F(agi);
	fr; % OK 

	% test frgr

	printf("\n Test friends and groups for agent %d..", agi)
	frgr = friendsAndGroups(agi);
	frgr; % OK 

	% test disadvantage 
	d = disadvantage(agi);
	printf("\n Disadvantage of agent %d", agi); 
	d; 

	printf("\n *** TEST *** "); 

endfunction 


function main()
  global tau 
  global W
  global K
  global groups
  global alpha
  global N
  global CH

	global NI

  global NSTEPS
  global ALGON

  global CURRENT_STEP
  global SUM_VOTING
  global SUM_RANKING

	global P_high_perf

	global ALGO
	global MERIT
	global U2G

	%INIT Friends 
	initFriends(); 

	%INIT GROUPS
	initGroups() 

	% MARK AGENT BEHAVIOURS 
	markAgentBehaviours()
	
	%MARK AGENT BEHAVIOURS 
	initTrustMatrix()

  %collect data at step ZERO
	collect_data(0); 

	% TEST 
	%test(); exit;

	computeAllD();
	computePerformance(); 

	printf("***4. MAIN SIMULATION, %d EPOCHS.. \n",NSTEPS);
	for s=1:NSTEPS
  	SUM_VOTING=0;
  	SUM_RANKING=0;
  	CURRENT_STEP=s;
  	printf("*** s = %d*** \n", s);
  	printf("4.1 Simulating %d random interactions.. \n",NI);
  	for c=1:NI % simulate 20 random interactions
    	index1 = randi([1,N]);
    	index2 = randi([1,N]);
    	if(index1!=index2)
      	%simulate feedback
      	feed = sampleFeedback(index2); %feedback of i on j
      	recc = sampleRecc(index2);
     		printf("[step = %d, c=%d], (%d, %d), Feed: %.2f, recc: %.2f\n",s,c,index1, index2, feed, recc);
      	old_tau=tau(index1,index2);
      	tau(index1,index2) = (old_tau + alpha*feed+(1-alpha)*recc)/2;
    		printf("\t tau(%d, %d)=%.2f<--%.2f\n", index1, index2,old_tau,tau(index1,index2));
    	endif
  	endfor

  	%run the group formation
		% different group formation may be executed 
		% 1. U2G 
		% 2. active and passive task for minimizing the disadvantage	 
  	printf("4.2 Launching %d random execution of the CLUSTERING ALGORITHM...\n", ALGON);
  	rndi = randi([1,N], 1, ALGON); % execution of the ALGON invocation of algoA
  	for j=1:ALGON
   		% printf("Execution %d..\n", j);
			if(ALGO==U2G)
    		algoA(rndi(1, j));
			else 
				if(ALGO==MERIT)
				% TODO 
					meritocraticActiveTask(rndi(1,j));
					exit; % TEST
					%printf("\n Warning, TODO meritocratic active task!"); 
				endif 
			endif
  	endfor

		computeAllD(); 
		computePerformance(); 
   
		collect_data(s); 

endfor % END MAIN SIMULATION 
endfunction

% Computation of 
% -- average trust within all the groups 
% -- statistics of average trust within the groups 
% -- group composition 
% -- statistics about voting 
function collect_data(i)
  printf("Computing average trust \n");
  flush_avg_global_tau(i, "average_tau.txt");
  flush_avg_global_tau_stat(i, "stat_tau.txt");
  
  printf("Tracing group composition\n");
  trace_groups(i, "groups.txt");

  printf("Tracing agents groups \n");
  trace_agents_groups(i, "agents_groups.txt");

  printf("Tracing statistics about voting and ranking.. \n");
  trace_stat_voting(i, "stat_voting.txt");

  printf("Tracing disadvantage values.. \n");
  flush_disadvantage(i, "D.txt");

  printf("Tracing performance stats.. \n");
  flush_performance_stats(i, "AP_stats.txt");
endfunction

% MAIN INVOCATION 
main
