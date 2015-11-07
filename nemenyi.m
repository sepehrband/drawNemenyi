function [p, nemenyi, meanrank, CDa, rankmean] = nemenyi(data, reps, varargin)
% Evaluate models using Freidman/Kruskal-Wallis and Nemenyi tests
%
% Syntax [p, nemenyi] = nemenyi(data, 'Properties', ...)
%
% Inputs:   data = Error distributions. Each column is a different
%                  model. Each row is a different measurement.
%                  
%                  If data is a cell array then ranks are inputted directly. The
%                  structure of the cell array is {ranks,N}, where ranks is 1xK
%                  ranks of K models and N is the number of time series.
%                  In this case reps is not needed.
%           reps = Number of replicates of each measurement. Consult
%                  friedman.m for more details. 
%   
%           For example, assessing 3 models, with 50 initialisations each,
%           across 5 time series requires data to be 5*50x3 and reps = 50.
%
%           Properties
%               'alpha'             significance level for nemenyi test. 
%                                   (Statistical tables are hardcoded)
%                                   0.01 | [0.05] | 0.10
%               'labels'            Optional. Cell array with names for
%                                   each model. It must contain equal
%                                   elements to the number of columns.
%               'ploton'            'line' | 'oline' | 'vline' | ['voline']
%                                   'matrix' | 'detail' | 'odetail' |
%                                   'mcb' | 'omcb' | 'off'
%               'colormap'          colourmap for 'detail' plot. 
%                                   Default: winter
%
% Outputs: p = Friedman/Kruskal-Wallis p-value
%          nemenyi  = nemenyi test result
%                     0 = insignificant difference
%                     1 = significant difference
%          meanrank = meanrank per model
%          CDa      = nemenyi test critical distance for a = [0.01; 0.05; 0.10]
%          rankmean = model ranking
%
% 2012 - Added critical values for up to 100 models. Use qtukey in R to get
% studentised range values e.g. qtukey(0.95,2:100,Inf) and then divide by
% sqrt(2).
%
% Nikolaos Kourentzes 2008 - Updated 2011 - Updated 2012
% http://nikolaos.kourentzes.com
%% Argument Check
if nargin<2 && ~iscell(data)
    error('Too few arguments')
end
%% Input variables check
% Default
cmap = 'winter';
alph = 0.05;
ploton = 'voline';
labels = [];
% Default names
if length(varargin) > 1
      okargs =   {'colormap' 'alpha' 'labels' 'ploton'};
      defaults = { cmap       alph    labels   ploton};
      [~, ~, cmap, alph, labels, ploton] = statgetargs(okargs,defaults,varargin{:});
end
% Find size of data matrix - columns are different models
if ~iscell(data)
    [r, k] = size(data);
    % Test for number of rows vs reps
    r = r/reps;
    if (floor(r) ~= r), 
       error('stats:friedman:BadSize',...
             'The number of rows must be a multiple of REPS.');
    end
else
    k = numel(data{1});
    r = data{2};
    % Rows vs reps is not tested in this case
    meanrank = data{1};
end
% Prepare labels
if isempty(labels)
    labels = [repmat('Model ',k,1) num2str((1:k)')];
    labels = mat2cell(labels,ones(k,1),size(labels,2));
elseif numel(labels)~=k
    error('nemenyi:labels','Incorrect number of labels')
end
%% Test Freidman or Kruskal-Wallis
if ~iscell(data)
    % If results data is provided
    if r>1
        p = friedman(data,reps,'off');
        lblt = 'Friedman ';
    else
        p = kruskalwallis(data,[],'off');
        lblt = 'Kruskal-Wallis ';
    end
else
    % If mean ranks are provided
    % Input is already average ranks for each model, so I test using the
    % following form of Friedman's
    % Koning et al, 2005 IJF
    S = ((12*r)/(k*(k+1)))*sum((meanrank-((k+1)/2)).^2);
    p = 1 - chi2cdf(S, k-1);
    lblt = 'Friedman ';
end
%% Nemenyi Test
% Hardcoded tables for nemenyi test 
% 1st row is for alpha = 0.01, second is for 0.05 and third for 0.10
% For up to 50 different models (columns of data<=50)
% q = [2.576 2.9130000 3.113 3.255 3.364 3.452 3.526 3.59 3.646 3.696 3.741 3.781 3.818 3.853 3.884 3.914 3.941 3.967 3.992 4.015 4.037 4.057 4.077 4.096 4.114 4.132 4.148 4.164 4.179 4.194 4.208 4.222 4.236 4.249 4.261 4.273 4.285 4.296 4.307 4.318 4.329 4.339 4.349 4.359 4.368 4.378 4.387 4.395 4.404; ...
%      1.96009999744911 2.34335187285222 2.56891893605073 2.7280179618177 2.84964032818179 2.9486352775479 3.03065966416554 3.10207744906538 3.16359573902861 3.21875006796116 3.26824754264422 3.31279526985897 3.35380746316779 3.39128412257068 3.42593235484882 3.4584592667834 3.48886485837443 3.51714912962189 3.54401918730698 3.569040161 3.592946027 3.615646276 3.637252631 3.657860551 3.677556303 3.696413427 3.71449839 3.731869175 3.748578108 3.764671858 3.780192852 3.795178566 3.809663649 3.823679212 3.837254248 3.850413505 3.863181025 3.875578729 3.887627121 3.899344587 3.910747391 3.921852503 3.932673359 3.943224099 3.953518159 3.963566147 3.973379375 3.98296845 3.992343271 ; ...
%      1.64473037303991 2.05202387900336 2.29102597104441 2.45931738496681 2.58871792592395 2.69266262275837 2.77963675684432 2.85459007565009 2.91964389951925 2.97762665557655 3.02995255738436 3.07662160494267 3.11975511859505 3.15935309834149 3.19541554418201 3.23006377646015 3.26117647483236 3.29087495964219 3.31915923088965 3.345675735 3.370711558 3.39447671 3.417089277 3.438651085 3.459252641 3.478971727 3.497877641 3.516032608 3.533492489 3.550305367 3.566516497 3.58216477 3.597287662 3.611916995 3.626083879 3.639814478 3.653134249 3.666065818 3.678630398 3.690847789 3.702736375 3.714311713 3.725589359 3.736584163 3.747309558 3.757777567 3.767999503 3.777987386 3.787749702];
% For up to 100 different models (columns of data<=100)
q = [2.57582949100000,2.91349419200000,3.11325044300000,3.25468594200000,3.36374019200000,3.45221268500000,3.52647091800000,3.59033892400000,3.64629157700000,3.69602098200000,3.74073346500000,3.78131856600000,3.81845086500000,3.85265432700000,3.88434331700000,3.91385017600000,3.94144643200000,3.96735694600000,3.99176980800000,4.01484199500000,4.03670927200000,4.05748760500000,4.07727528100000,4.09616068900000,4.11421948900000,4.13151885600000,4.14811818800000,4.16406910300000,4.17941968400000,4.19421235800000,4.20848389400000,4.22226894100000,4.23559861100000,4.24850118800000,4.26100212900000,4.27312476800000,4.28489102400000,4.29631999100000,4.30743005300000,4.31823818000000,4.32875992900000,4.33900944200000,4.34899944700000,4.35874337800000,4.36825184300000,4.37753615500000,4.38660550600000,4.39547050400000,4.40413892600000,4.41261925800000,4.42091857000000,4.42904605500000,4.43700666400000,4.44480746600000,4.45245482500000,4.45995440000000,4.46731113900000,4.47452999200000,4.48161732300000,4.48857596100000,4.49541156200000,4.50212837000000,4.50872921200000,4.51521833000000,4.52159996900000,4.52787695600000,4.53405212000000,4.54012970200000,4.54611182600000,4.55200202500000,4.55780242200000,4.56351513800000,4.56914370800000,4.57469025300000,4.58015689600000,4.58554575700000,4.59085966400000,4.59609932500000,4.60126756900000,4.60636580900000,4.61139687400000,4.61648167800000,4.62126101300000,4.62609833100000,4.63087413000000,4.63559053200000,4.64024683000000,4.64484726700000,4.64939184200000,4.65388197000000,4.65831906500000,4.66270383400000,4.66703769200000,4.67132275900000,4.67555832900000,4.67974652200000,4.68388875400000,4.68798502300000,4.69203674500000;...
     1.95996423300000,2.34370047600000,2.56903207300000,2.72777471700000,2.84970538200000,2.94831990800000,3.03087886700000,3.10173026000000,3.16368342000000,3.21865390100000,3.26800359100000,3.31273870100000,3.35361795900000,3.39123038200000,3.42604124900000,3.45842461900000,3.48868454600000,3.51707276200000,3.54379927700000,3.56904016100000,3.59294602700000,3.61564627600000,3.63725263100000,3.65786055100000,3.67755630300000,3.69641342700000,3.71449839000000,3.73186917500000,3.74857810800000,3.76467185800000,3.78019285200000,3.79517856600000,3.80966364900000,3.82367921200000,3.83725424800000,3.85041350500000,3.86318102500000,3.87557872900000,3.88762712100000,3.89934458700000,3.91074739100000,3.92185250300000,3.93267335900000,3.94322409900000,3.95351815900000,3.96356614700000,3.97337937500000,3.98296845000000,3.99234327100000,4.00151232500000,4.01048480300000,4.01926777600000,4.02786973000000,4.03629702900000,4.04455603600000,4.05265453000000,4.06059675300000,4.06838977700000,4.07603784400000,4.08354731800000,4.09092102800000,4.09816604400000,4.10528448800000,4.11228201600000,4.11916145800000,4.12592705600000,4.13258234500000,4.13913156800000,4.14557613900000,4.15192100800000,4.15816829700000,4.16432083300000,4.17038073800000,4.17635225500000,4.18223679700000,4.18803648700000,4.19375486000000,4.19939262200000,4.20495260300000,4.21043763000000,4.21584841100000,4.22118706700000,4.22645572000000,4.23165649000000,4.23679079300000,4.24185933400000,4.24686494300000,4.25180903400000,4.25669231300000,4.26151619600000,4.26628280200000,4.27099284100000,4.27564843200000,4.28024957500000,4.28479839300000,4.28929488500000,4.29374188000000,4.29813937700000,4.30248879100000;...
     1.64485341000000,2.05229258000000,2.29134134100000,2.45951608200000,2.58852064300000,2.69273191900000,2.77988353700000,2.85460633900000,2.91988855800000,2.97776807700000,3.02969446300000,3.07673332800000,3.11969360000000,3.15919894900000,3.19574364200000,3.22972365800000,3.26146143900000,3.29122427000000,3.31923277000000,3.34567573500000,3.37071155800000,3.39447671000000,3.41708927700000,3.43865108500000,3.45925264100000,3.47897172700000,3.49787764100000,3.51603260800000,3.53349248900000,3.55030536700000,3.56651649700000,3.58216477000000,3.59728766200000,3.61191699500000,3.62608387900000,3.63981447800000,3.65313424900000,3.66606581800000,3.67863039800000,3.69084778900000,3.70273637500000,3.71431171300000,3.72558935900000,3.73658416300000,3.74730955800000,3.75777756700000,3.76799950300000,3.77798738600000,3.78774970200000,3.79729705800000,3.80663793900000,3.81578153700000,3.82473492300000,3.83350516800000,3.84210075800000,3.85052664200000,3.85879059900000,3.86689757900000,3.87485323700000,3.88266323100000,3.89033321900000,3.89786673400000,3.90526872800000,3.91254344300000,3.91969582800000,3.92672941900000,3.93364704500000,3.94045294700000,3.94715137000000,3.95374443300000,3.96023567400000,3.96662862600000,3.97292470500000,3.97912815300000,3.98524038400000,3.99126493400000,3.99720392300000,4.00305876800000,4.00883300100000,4.01452804000000,4.02014671000000,4.02568972000000,4.03115989800000,4.03655865800000,4.04188741500000,4.04714899700000,4.05234481700000,4.05747558400000,4.06254341800000,4.06754973400000,4.07249523900000,4.07738276100000,4.08221300800000,4.08698668600000,4.09170520900000,4.09637070000000,4.10098315700000,4.10554470300000,4.11005533700000;]; 
 if ~iscell(data)
    % Get a matrix of ranks.  For the unusual case of replicated
    % measurements, rank together all replicates in the same row.  This
    % is the advice given by Zar (1996), "Biostatistical Analysis."
    m = nan(reps*r,k);
    for j=1:r
       jrows = reps * (j-1) + (1:reps);
       v = data(jrows,:);
       a = tiedrank(v(:));
       m(jrows,:) = reshape(a, reps, k);
    end
    meanrank=mean(m);
 end
% Find critical distance
    q = q(:,k-1);
    N = r;
    CDa=q*sqrt((k*(k+1))/(6*N));
    switch alph
        case 0.01
            CD = CDa(1);
        case 0.05
            CD = CDa(2);
        case 0.10
            CD = CDa(3);
        otherwise
            error('nemenyi:alpha','Alpha must be either 0.01, 0.05 or 0.10')
    end
% Produce comparison matrix
    % 0 = same, 1 = different
    nemenyi = abs(bsxfun(@minus,meanrank,meanrank'));
    nemenyi = ~(nemenyi<CD);
    nemenyi = single(nemenyi);
    nemenyi(1:k+1:end) = nan(k,1);
%% Find model ranking
if nargout>4 || ~strcmp(ploton,'off')
    rankmean = tiedrank(meanrank);
end
% 
% %% Plots
% % clf
% switch ploton
%     %% Horizontal line based on ranking and statistical difference
%     case 'line'
%         % Start plotting
%         plot([0 k+1],[1 1],'k');
%         hold on
%         % Plot dots for models
%         plot(rankmean,1,'or','markerfacecolor','k');
%         % Find groups (not significantly different) for each model
%             rline = nan(k,k);
%             srankmean = sort(rankmean); % used to take care of ties
%             for i=1:k
%                 idx = rankmean==srankmean(i); % used to take care of ties
%                 idx = find(idx==1); % used to take care of ties
%                 idx = idx(1); % used to take care of ties
%                 rline(i,:) = abs(meanrank - meanrank(idx))<CD;
%             end   
%         % Get rid of exactly equal groups
%             rline = unique(rline,'rows');
%         % Get rid of groups with only one member
%             rline = rline(sum(rline,2)~=1,:);
%             K = size(rline,1);
%         % Set colors
%         if p <= alph
%             friedmanH='Different';
%             clmap = zeros(K,3);
%             clmap(:,1) = 1;
%         else
%             friedmanH='All same';
%             clmap = zeros(K,3);
%             clmap(:,3) = 1;
%         end
%         % Plot remaining groups
%             % Plot dotted line
%                 for i=1:K
%                     plot([1 1]*min(rankmean(rline(i,:)==1)), [1 1-i/20], ':k','linewidth', 1);
%                     plot([1 1]*max(rankmean(rline(i,:)==1)), [1 1-i/20], ':k','linewidth', 1);
%                 end            
%             % Plot group lines
%                 for i=1:K
%                     plot(sort(rankmean(rline(i,:)==1)), ones(1,numel(rankmean(rline(i,:)==1)))*1-i/20, 'color', clmap(i,:),'linewidth', 4);
%                 end
%         % Add names & ranks
%         % Find duplicate names
%             [~, m] = unique(rankmean);
%             clear b
%             dubl = rankmean(setdiff(1:k,m));
%             dubl = unique(dubl);
%             dubl = dubl(:);
%             dloc = cell(numel(dubl),1);
%             for i = 1:numel(dubl)
%                 dubl(i,2) = numel(find(rankmean==dubl(i,1)));
%                 dubl(i,3) = (dubl(i,2)-1)/dubl(i,2);
%                 if mod(dubl(i,2),2)
%                    % disp('odd') 
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                        dloc{i}(floor(dubl(i,2)/2)+1) = dubl(i,1);
%                        dloc{i}(floor(dubl(i,2)/2)+2:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                 else
%                     % disp('even')
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))+dubl(i,3)/2;
%                        dloc{i}(floor(dubl(i,2)/2)+1:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))-dubl(i,3)/2;   
%                        dloc{i} = sort(dloc{i});
%                 end
%             end
%             dubl(:,4) = 1; % counter
%             for i=1:k
%                 if sum(rankmean(i)==dubl(:,1))==1
%                     ll = find(rankmean(i)==dubl(:,1),1,'first');
%                     loc = dloc{ll};
%                     loc = loc(dubl(ll,4));
%                     dubl(ll,4) = dubl(ll,4) + 1;
%                 else
%                     loc = rankmean(i);
%                 end
%                 text(loc,1.1,[sprintf('%.2f',meanrank(i)) ' - ' labels{i}],'Rotation',90,'HorizontalAlignment','left');
%             end
%         ylim([min(0,1-0.1-0.05*size(rline,1)) 2])
%         xlim([0 k+1])
%         set(gca,'xtick',[],'ytick',[]);
%         % Title
%         title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
%         hold off
% %% Ordered Horizontal Line Nemenyi Graph
%     case 'oline'
%         % Create a straight black line as axis
%         plot([0 k+1],[1 1],'k');
%         hold on
%         % Plot dots for models
%         srankmean = sort(rankmean);
%         plot(srankmean,1,'or','markerfacecolor','k');
%         % Find groups (not significantly different) for each model
%         rline = nan(k,k);
%         srankmean = sort(rankmean); % used to take care of ties
%         for i=1:k
%             idx = rankmean==srankmean(i); % used to take care of ties
%             idx = find(idx==1); % used to take care of ties
%             idx = idx(1); % used to take care of ties
%             rline(i,:) = abs(meanrank - meanrank(idx))<CD;
%         end   
%         % Get boundaries of rline
%         rankvector = nan(k,2);
%         for i = 1:k
%             temprank = rankmean(rline(i,:)==1);
%             rankvector(i,:) = [min(temprank) max(temprank)];
%         end
%         rankvector = unique(rankvector,'rows');
%         rankvector = rankvector(~(rankvector(:,1) == rankvector(:,2)),:);
%         K = numel(rankvector(:,1));
%         % Set colors
%         if p <= alph
%             friedmanH='Different';
%             cmap = lower(cmap);
%             if ~strcmp(cmap,'red')
%                 clmap = eval(['colormap(' cmap '(' int2str(K) '))']);
%             else
%                 clmap = repmat([1 0 0],K,1);
%             end
%         else
%             friedmanH='All same';
%             clmap = zeros(K,3);
%             clmap(:,3) = 1;
%         end
%         % Plot lines
%         memran = 0;
%         for i= K:-1:1;
%             if memran ~= rankvector(i,1)
%                 plot([rankvector(i,1) rankvector(i,1)],[1-i/20 1],':k')
%             end
%             plot([rankvector(i,2) rankvector(i,2)],[1-i/20 1],':k')
%             plot(rankvector(i,:), ones(1,2)*1-i/20, 'color', clmap(i,:),'linewidth', 4);
%             memran = rankvector(i,1);
%         end  
%         % Add names & ranks
%         % Find duplicate names
%             [~, m] = unique(rankmean);
%             dubl = rankmean(setdiff(1:k,m));
%             dubl = unique(dubl);
%             dubl = dubl(:);
%             dloc = cell(numel(dubl),1);
%             for i = 1:numel(dubl)
%                 dubl(i,2) = numel(find(rankmean==dubl(i,1)));
%                 dubl(i,3) = (dubl(i,2)-1)/dubl(i,2);
%                 if mod(dubl(i,2),2)
%                    % disp('odd') 
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                        dloc{i}(floor(dubl(i,2)/2)+1) = dubl(i,1);
%                        dloc{i}(floor(dubl(i,2)/2)+2:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                 else
%                     % disp('even')
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))+dubl(i,3)/2;
%                        dloc{i}(floor(dubl(i,2)/2)+1:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))-dubl(i,3)/2;   
%                        dloc{i} = sort(dloc{i});
%                 end
%             end
%             dubl(:,4) = 1; % counter
%             for i=1:k
%                 if sum(rankmean(i)==dubl(:,1))==1
%                     ll = find(rankmean(i)==dubl(:,1),1,'first');
%                     loc = dloc{ll};
%                     loc = loc(dubl(ll,4));
%                     dubl(ll,4) = dubl(ll,4) + 1;
%                 else
%                     loc = rankmean(i);
%                 end
%                 text(loc,1.1,[sprintf('%.2f',meanrank(i)) ' - ' labels{i}],'Rotation',90,'HorizontalAlignment','left');
%             end
%         ylim([min(0,1-0.1-0.05*K) 2])
%         xlim([0 k+1])
%         set(gca,'xtick',[],'ytick',[]);
%         % Title
%         title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
%         hold off       
% %% Vertical line based on ranking and statistical difference
%     case 'vline'
%         % Create a straight black line as axis
%             idx = 0:k+1;
%             plot(ones(k+2,1),idx,'k');
%             hold on
%         % Plot dots for models
%             plot(1,rankmean,'or','markerfacecolor','k');
%         % Find groups (not significantly different) for each model
%             rline = nan(k,k);
%             srankmean = sort(rankmean); % used to take care of ties
%             for i=1:k
%                 idx = rankmean==srankmean(i); % used to take care of ties
%                 idx = find(idx==1); % used to take care of ties
%                 idx = idx(1); % used to take care of ties
%                 rline(i,:) = abs(meanrank - meanrank(idx))<CD;
%             end   
%         % Get rid of exactly equal groups
%             rline = unique(rline,'rows');
%         % Get rid of groups with only one member
%             rline = rline(sum(rline,2)~=1,:);
%             K = size(rline,1);
%         % Set colors
%         if p <= alph
%             friedmanH='Different';
%             clmap = zeros(K,3);
%             clmap(:,1) = 1;
%         else
%             friedmanH='All same';
%             clmap = zeros(K,3);
%             clmap(:,3) = 1;
%         end
%         % Plot remaining groups
%             % Plot dotted line
%                 for i=1:K;
%                     plot([1 1+i/20], [1 1]*min(rankmean(rline(i,:)==1)), ':k','linewidth', 1);
%                     plot([1 1+i/20], [1 1]*max(rankmean(rline(i,:)==1)), ':k','linewidth', 1);
%                 end            
%             % Plot group lines
%                 for i=1:K;
%                     plot(ones(1,numel(rankmean(rline(i,:)==1)))*1+i/20, sort(rankmean(rline(i,:)==1)), 'color', clmap(i,:),'linewidth', 4);
%                 end
%         % Add names & ranks
%             % Find duplicate names
%             [~, m] = unique(rankmean);
%             clear b
%             dubl = rankmean(setdiff(1:k,m));
%             dubl = unique(dubl);
%             dubl = dubl(:);
%             dloc = cell(numel(dubl),1);
%             for i = 1:numel(dubl)
%                 dubl(i,2) = numel(find(rankmean==dubl(i,1)));
%                 dubl(i,3) = (dubl(i,2)-1)/dubl(i,2);
%                 if mod(dubl(i,2),2)
%                    % disp('odd') 
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                        dloc{i}(floor(dubl(i,2)/2)+1) = dubl(i,1);
%                        dloc{i}(floor(dubl(i,2)/2)+2:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                 else
%                     % disp('even')
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))+dubl(i,3)/2;
%                        dloc{i}(floor(dubl(i,2)/2)+1:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))-dubl(i,3)/2;   
%                        dloc{i} = sort(dloc{i});
%                 end
%             end
%             dubl(:,4) = 1; % counter
%             for i=1:k
%                 if sum(rankmean(i)==dubl(:,1))==1
%                     ll = find(rankmean(i)==dubl(:,1),1,'first');
%                     loc = dloc{ll};
%                     loc = loc(dubl(ll,4));
%                     dubl(ll,4) = dubl(ll,4) + 1;
%                 else
%                     loc = rankmean(i);
%                 end
%                 text(0.9, loc,[labels{i} ' - ' sprintf('%.2f',meanrank(i))],'HorizontalAlignment','right');
%             end
%         xlim([min(0,1-0.1-0.05*size(rline,1)) 2])
%         ylim([0 k+1])
%         set(gca,'xtick',[],'ytick',[]);
%         set(gca,'YDir','reverse');
%         % Title
%         title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.2f')])
%         hold off
% %% Ordered Vertical Line Nemenyi Graph
%     case 'voline'
%         % Create a straight black line as axis
%         plot([1 1],[0 k+1],'k');
%         hold on
%         % Plot dots for models
%         srankmean = sort(rankmean,'descend');
%         plot(1,srankmean,'or','markerfacecolor','k');
%         % Find groups (not significantly different) for each model
%         rline = nan(k,k);
%         srankmean = sort(rankmean); % used to take care of ties
%         for i=1:k
%             idx = rankmean==srankmean(i); % used to take care of ties
%             idx = find(idx==1); % used to take care of ties
%             idx = idx(1); % used to take care of ties
%             rline(i,:) = abs(meanrank - meanrank(idx))<CD;
%         end   
%         % Get boundaries of rline
%         rankvector = nan(k,2);
%         for i = 1:k
%             temprank = rankmean(rline(i,:)==1);
%             rankvector(i,:) = [min(temprank) max(temprank)];
%         end
%         rankvector = unique(rankvector,'rows');
%         rankvector = rankvector(~(rankvector(:,1) == rankvector(:,2)),:);
%         K = numel(rankvector(:,1));
%         % Set colors
%         if p <= alph
%             friedmanH='Different';
%             cmap = lower(cmap);
%             if ~strcmp(cmap,'red')
%                 clmap = eval(['colormap(' cmap '(' int2str(K) '))']);
%             else
%                 clmap = repmat([1 0 0],K,1);
%             end
%         else
%             friedmanH='All same';
%             clmap = zeros(K,3);
%             clmap(:,3) = 1;
%         end
%         % Plot lines
%         memran = 0;
%         for i= K:-1:1;
%             if memran ~= rankvector(i,1)
%                 plot([1+i/20 1],[rankvector(i,1) rankvector(i,1)],':k')
%             end
%             plot([1+i/20 1],[rankvector(i,2) rankvector(i,2)],':k')
%             plot(ones(1,2)*1+i/20, rankvector(i,:), 'color', clmap(i,:),'linewidth', 4);
%             memran = rankvector(i,1);
%         end  
%         % Add names & ranks
%         % Find duplicate names
%             [~, m] = unique(rankmean);
%             dubl = rankmean(setdiff(1:k,m));
%             dubl = unique(dubl);
%             dubl = dubl(:);
%             dloc = cell(numel(dubl),1);
%             for i = 1:numel(dubl)
%                 dubl(i,2) = numel(find(rankmean==dubl(i,1)));
%                 dubl(i,3) = (dubl(i,2)-1)/dubl(i,2);
%                 if mod(dubl(i,2),2)
%                    % disp('odd') 
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                        dloc{i}(floor(dubl(i,2)/2)+1) = dubl(i,1);
%                        dloc{i}(floor(dubl(i,2)/2)+2:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3));
%                 else
%                     % disp('even')
%                        dloc{i}(1:floor(dubl(i,2)/2)) = dubl(i,1) - cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))+dubl(i,3)/2;
%                        dloc{i}(floor(dubl(i,2)/2)+1:dubl(i,2)) = dubl(i,1) + cumsum(ones(floor(dubl(i,2)/2),1)*dubl(i,3))-dubl(i,3)/2;   
%                        dloc{i} = sort(dloc{i});
%                 end
%             end
%             dubl(:,4) = 1; % counter
%             for i=1:k
%                 if sum(rankmean(i)==dubl(:,1))==1
%                     ll = find(rankmean(i)==dubl(:,1),1,'first');
%                     loc = dloc{ll};
%                     loc = loc(dubl(ll,4));
%                     dubl(ll,4) = dubl(ll,4) + 1;
%                 else
%                     loc = rankmean(i);
%                 end
%                 text(0.9, loc,[labels{i} ' - ' sprintf('%.2f',meanrank(i))],'HorizontalAlignment','right');
%             end
%         xlim([min(0,1-0.1-0.05*K) 2])
%         ylim([0 k+1])
%         set(gca,'xtick',[],'ytick',[],'ydir','reverse');
%         % Title
%         title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
%         hold off       
% %% Detailed nemenyi graph with meanranks and model groups
%     case 'detail'
%         % plot ranks
%         hp = plot(meanrank,1,'o','markersize',14);
%         xx = xlim;
%         hold on
%         % plot groupings
%         name2 = cell(k,1);
%         rankvector = nan(k,2);
%         for i = 1:k
%             ranktemp = [meanrank(nemenyi(i,:)==0) meanrank(i)];
%             rankvector(i,:) = [min(ranktemp) max(ranktemp)];
%         end
%         [urankvector, ~, uid] = unique(rankvector,'rows');
%         K = numel(urankvector(:,1));
%         cmap = lower(cmap);
%         clmap = eval(['colormap(' cmap '(' int2str(K) '))']);
%         for i = 1:k
%             plot(xx,[i+1 i+1],':','color',[0.8 0.8 0.8])
%             plot(rankvector(i,:),[i+1 i+1], 'color', clmap(uid(i),:),'linewidth', 3)
%             name2{i} = [int2str(i) '-' labels{i}];
%             if mean(clmap(uid(i),:))>0.43
%                 text(meanrank(i),1,int2str(i),'horizontalalignment','center');
%             else
%                 text(meanrank(i),1,['\color[rgb]{1 1 1}' int2str(i)],'horizontalalignment','center');
%             end
%             set(hp(i),'color', clmap(uid(i),:),'markerfacecolor', clmap(uid(i),:),'markeredgecolor','k')
%         end
%         hold off
%         ylim([0 k+2])
%         % Set visual properties
%         set(gca,'YDir', 'reverse', 'ytick', 0:k+2, 'yticklabel', [{''}; {'Models:'}; name2; {''}])
%         set(get(gca,'XLabel'),'String','Average Rank')
%         box('off')
%         % Label
%         if p <= alph
%             friedmanH='Different';
%         else
%             friedmanH='All same ';
%         end
%         title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
% %% Ordered Detailed Nemenyi Graph
%     case 'odetail'
%         % Sort ranks
%         [smeanrank sid] = sort(meanrank);     
%         % plot ranks
%         hp = plot(smeanrank,1,'o','markersize',14);
%         xx = xlim;
%         hold on
%         % plot groupings
%         name2 = cell(k,1);
%         rankvector = nan(k,2);
%         for i = 1:k
%             ranktemp = [meanrank(nemenyi(sid(i),:)==0) meanrank(sid(i))];
%             rankvector(i,:) = [min(ranktemp) max(ranktemp)];
%         end
%         [urankvector, ~, uid] = unique(rankvector,'rows');
%         K = numel(urankvector(:,1));
%         cmap = lower(cmap);
%         clmap = eval(['colormap(' cmap '(' int2str(K) '))']);
%         for i = 1:k
%             % plot(ones(2,1)*rankvector(i,1),[1 i+1],':','color', clmap(uid(i),:))
%             plot(xx,[i+1 i+1],':','color',[0.8 0.8 0.8])
%             plot(rankvector(i,:),[i+1 i+1], 'color', clmap(uid(i),:),'linewidth', 3)
%             plot(ones(2,1)*rankvector(i,2),[1 i+1],':','color', clmap(uid(i),:))
%             name2{i} = [int2str(i) '-' labels{sid(i)}];
%             if mean(clmap(uid(i),:))>0.43
%                 text(meanrank(sid(i)),1,int2str(i),'horizontalalignment','center');
%             else
%                 text(meanrank(sid(i)),1,['\color[rgb]{1 1 1}' int2str(i)],'horizontalalignment','center');
%             end
%             set(hp(i),'color', clmap(uid(i),:),'markerfacecolor', clmap(uid(i),:),'markeredgecolor','k')
%         end
%         hold off
%         ylim([0 k+2])
%         % set visual properties
%         set(gca,'YDir', 'reverse', 'ytick', 0:k+2, 'yticklabel', [{''}; {'Models:'}; name2; {''}])
%         % set(gca,'XDir', 'reverse')
%         set(get(gca,'XLabel'),'String','Average Rank')
%         box('off')
%         % Label
%         if p <= alph
%             friedmanH='Different';
%         else
%             friedmanH='All same ';
%         end
%        title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
% %% Ordered and Grouped Detailed Nemenyi Graph
%     case 'ogdetail'
%         % Sort ranks
%         [smeanrank sid] = sort(meanrank);     
%         % plot ranks
%         hp = plot(smeanrank,1,'o','markersize',14);
%         xx = xlim;
%         hold on
%         % plot groupings
%         name2 = cell(k,1);
%         rankvector = nan(k,2);
%         for i = 1:k
%             ranktemp = [meanrank(nemenyi(sid(i),:)==0) meanrank(sid(i))];
%             rankvector(i,:) = [min(ranktemp) max(ranktemp)];
%         end
%         [urankvector, ~, uid] = unique(rankvector,'rows');
%         K = numel(urankvector(:,1));
%         cmap = lower(cmap);
%         clmap = eval(['colormap(' cmap '(' int2str(K) '))']);
%         for i = 1:k
% %             plot(ones(2,1)*rankvector(i,1),[1 i+1],':','color', clmap(uid(i),:))
%             plot(xx,[i+1 i+1],':','color',clmap(uid(i),:),'linewidth',2)
%             plot(ones(2,1)*rankvector(i,2),[1 i+1],':','color', clmap(uid(i),:))
%             name2{i} = [int2str(i) '-' labels{sid(i)}];
%             if mean(clmap(uid(i),:))>0.43
%                 text(meanrank(sid(i)),1,int2str(i),'horizontalalignment','center');
%             else
%                 text(meanrank(sid(i)),1,['\color[rgb]{1 1 1}' int2str(i)],'horizontalalignment','center');
%             end
%             set(hp(i),'color', clmap(uid(i),:),'markerfacecolor', clmap(uid(i),:),'markeredgecolor','k')
%         end
%         for i = 1:K
%             ranktemp = find(uid==i);
%             patch([urankvector(i,:) urankvector(i,end:-1:1)],[ranktemp(1)+1-0.25 ranktemp(1)+1-0.25 ranktemp(end)+1+0.25 ranktemp(end)+1+0.25],clmap(i,:))
%         end
%         hold off
%         ylim([0 k+2])
%         % set visual properties
%         set(gca,'YDir', 'reverse', 'ytick', 0:k+2, 'yticklabel', [{''}; {'Models:'}; name2; {''}])
%         % set(gca,'XDir', 'reverse')
%         set(get(gca,'XLabel'),'String','Average Rank')
%         box('off')
%         % Label
%         if p <= alph
%             friedmanH='Different';
%         else
%             friedmanH='All same ';
%         end
%        title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
% %% Results as a matrix of statistically different models
%     case 'matrix'
%         % create matrix
%         plnemenyi = nemenyi;
%         plnemenyi(isnan(plnemenyi)) = 0.5;
%         imagesc(plnemenyi) 
%         colormap(gray) % black are the same, white are significantly different
%         % Label
%         if p <= alph
%             friedmanH='Different';
%         else
%             friedmanH='All same';
%         end
%         set(gca,'yticklabel',labels,'ytick',1:k,'xtick',1:k)
%         xticklabel_rotate([],90,labels)
%         title([lblt 'p-value: ' num2str(p) ' \bullet Friedman test: ' friedmanH ' \bullet CritDist: ' num2str(CD)])
%         % set(gca,'yticklabel',labels,'xticklabel',labels,'ytick',1:c,'xtick',1:c)
%         hold off
% %% Results similar to Multiple Comparison with Best as in Koning et al 2005 IJF
%     case 'mcb'
%         % Set best color
%         if p <= alph
%             friedmanH='Different';
%             lc = 'r';
%         else
%             friedmanH='All same';
%             lc = 'b';
%         end
%         % Find maximum label length
%         nn = nan(k,1);
%         for i = 1:k
%             nn(i) = numel(labels{i});
%             % Add space at the end of labels if there is none
%             if ~strcmp(labels{i}(end),' ')
%                 labels{i} = [labels{i} ' '];
%             end
%         end
%         nn = max(nn);
%         % Calculate bounds
%         lb = meanrank - CD;
%         ub = meanrank + CD;
%         [mub, id] = min(ub);
%         % Produce plot
%         plot([1:k;1:k],[lb; ub],'-b')
%         hold on
%         plot(1:k,meanrank,'ob','markerfacecolor','b')
%         plot([(1:k)-0.2;(1:k)+0.2],[ub;ub],'-b')
%         plot([(1:k)-0.2;(1:k)+0.2],[lb;lb],'-b')
%         plot([0 k+1],[mub mub],':k')
%         % Plot best model
%         plot([id;id],[lb(id); ub(id)],'-','color',lc)
%         plot(id,meanrank(id),'o','color',lc,'markerfacecolor',lc)
%         plot([id-0.2;id+0.2],[ub(id);ub(id)],'-','color',lc)
%         plot([id-0.2;id+0.2],[lb(id);lb(id)],'-','color',lc)
%         hold off
%         xlim([0.5 k+0.5])
%         ylabel('Average Rank')
%         set(gca,'xtick',1:k,'xticklabel','')
%         yy = ylim;
%         % Adjust to fit large names
%         pos = get(gca,'position');
%         pos(2) = pos(2) + 0.010*nn;
%         pos(4) = pos(4) - 0.010*nn;
%         set(gca,'position',pos)
%         % Print model names
%         for i = 1:k
%             text(i,yy(1),labels{i},'rotation',90,'horizontalalignment','right')
%         end
%         % Title
%         title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
% %% Ordered version of Multiple Comparison with Best as in Koning et al 2005 IJF
%     case 'omcb'
%         % Set best color
%         if p <= alph
%             friedmanH='Different';
%             lc = 'r';
%         else
%             friedmanH='All same';
%             lc = 'b';
%         end
%         % Find maximum label length
%         nn = nan(k,1);
%         for i = 1:k
%             nn(i) = numel(labels{i});
%             % Add space at the end of labels if there is none
%             if ~strcmp(labels{i}(end),' ')
%                 labels{i} = [labels{i} ' '];
%             end
%         end
%         nn = max(nn);
%         % Sort rank
%         [smeanrank, sid] = sort(meanrank);
%         labels = labels(sid);
%         % Calculate bounds
%         lb = smeanrank - CD;
%         ub = smeanrank + CD;
%         [mub, id] = min(ub);
%         % Produce plot
%         plot([1:k;1:k],[lb; ub],'-b')
%         hold on
%         plot(1:k,smeanrank,'ob','markerfacecolor','b')
%         plot([(1:k)-0.2;(1:k)+0.2],[ub;ub],'-b')
%         plot([(1:k)-0.2;(1:k)+0.2],[lb;lb],'-b')
%         plot([0 k+1],[mub mub],':k')
%         % Plot best model
%         plot([id;id],[lb(id); ub(id)],'-','color',lc)
%         plot(id,smeanrank(id),'o','color',lc,'markerfacecolor',lc)
%         plot([id-0.2;id+0.2],[ub(id);ub(id)],'-','color',lc)
%         plot([id-0.2;id+0.2],[lb(id);lb(id)],'-','color',lc)
%         hold off
%         xlim([0.5 k+0.5])
%         ylabel('Average Rank')
%         set(gca,'xtick',1:k,'xticklabel','')
%         yy = ylim;
%         % Adjust to fit large names
%         pos = get(gca,'position');
%         pos(2) = pos(2) + 0.010*nn;
%         pos(4) = pos(4) - 0.010*nn;
%         set(gca,'position',pos)
%         % Print model names
%         for i = 1:k
%             text(i,yy(1),labels{i},'rotation',90,'horizontalalignment','right')
%         end
%         % Title
%         title([lblt 'p-value: ' num2str(p,'%.3f') ' \bullet ' friedmanH ' \bullet CritDist: ' num2str(CD,'%.1f')])
%% Close ploton switch 
end