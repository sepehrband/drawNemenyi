%% drawNemenyi demo
% Farshid Sepehrband - Jeiran Choupan
% nemenyi.m function is borrowed from: http://kourentzes.com/forecasting/2013/04/19/nemenyi-test/

% Let's say we have a matrix of MxN, as the Nemenyi's input:
% M measurements x N techniques

% First we define a cell of strings, including name of technique. 
% for example: 
% Names = {'tech #1','tech #2',...,'tech #n'};

% Then, we define our output folder and output name and run the code.

% Here is an example, assuming you are comparing 4 techniques:

Names = {'tech #1','tech #2','tech #3','tech #4'}; % names
OutputFolder  = '~/Address/to/output';             % Output folder
Outname = 'NemenyiResults';                        % Output name
Draw_Nemenyi(results, Names, OutputFolder, Outname);

% Now you should see a "NemeyiResults.tif" on the output folder.
% You can change the print options (see %%print section of Draw_Nemenyi.m)
% Enjoy!

