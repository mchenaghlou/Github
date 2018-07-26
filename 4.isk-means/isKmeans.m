function [ V, ClusterIndices, iCVs_ff, iCVs1_ff, added, removed,clustering_snap_shot, ClusterLabelDic, Radiis] = isKmeans( X ,record_video, lambda, lamda_delta, emerging_cluster_capacity,myPathVid)
%==========================================================================
% Author: Milad Chenaghlou (mchenaghlou@student.unimelb.edu.au)
% Created: 20018-07-26
% iskM (Incremental Sequential k-Means) 
%   Input Data:
%   X: dataset with d columns (d dimensions) and n rows (data-points).
%   record_video: record a video of the process
%       record_video = 0 : No video
%       record_video = 1 : Show iCVIs any data dimension
%       record_video = 2 : Show iCVIs and data-points. (Supports 2D only)
%   lambda: forgetting factor for iCVIs
%   lambda_delta: forgetting factor for delta
%   emerging_cluster_capacity: not implemented. (a mechanism to handle more 
%       emerging clusters)
%   myPathVid: The path to create the video file.
%==========================================================================

% Ignore the parameter emerging_cluster_capacity by setting it to 0.
emerging_cluster_capacity = 0;


% Initialize the video settings.
if record_video > 0
    vidObj = VideoWriter(myPathVid);
    vidObj.Quality = 20;
    vidObj.FrameRate = 15;
    open(vidObj);
    f = figure('pos',[100 50 1200 900]);
    p = uipanel('Parent',f,'BorderType','none');
    p.FontSize = 12;
    p.FontWeight = 'bold';
end


% The number of dimensions
D = size(X, 2);

% The number of data-points
n = size(X, 1);

% Start Current-skM with 1 prototype.
k = 1;

% set the radii of the first cluster in Current-skM from the first 10
% data-points
stab_period = 10;

% start Control-skM with 2 prototypes
kp = k + 1;

% initialize the prototypes in Current-skM 
V=X(1:k,1:D);

% Initialize the prototypes in Control-skM
Vp=X(1:kp,1:D);

% iCVI values for Current-skM with forgetting factor
iCVs_ff = zeros(n,1);

% iCVI values for Control-skM with forgetting factor
iCVs1_ff = zeros(n,1);

% iCVI values for Current-skM without forgetting factor
iCVs = zeros(n,1);

% iCVI values for Control-skM without forgetting factor
iCVs1 = zeros(n,1);

% initialize cluster indices for Current-skM
ClusterIndices=zeros(n,1);
ClusterIndices(1:k)=1:1:k;
ClusterIndices(k+1:k+1)=k;

% initialize cluster indices for Control-skM
ClusterIndices1=zeros(n,1);
ClusterIndices1(1:k+1)=1:1:k+1;

% initialize the delta means and variances with zeros.
delta_means = zeros(n,1);
delta_vars = zeros(n,1);

% initialize the means and variances of iCVI values with forgetting factor
% for Current-skM and Control-skM
iCVs_ff_variance = zeros(n,1);
iCVs_ff_means = zeros(n,1);
iCVs1_ff_variance = zeros(n,1);
iCVs1_ff_means = zeros(n,1);

% initialize number of data-points for each prototype for Current-skM and
% Control-skM
NIs = ones(k,1);
NIsp = ones(kp,1);

% initialize the separation values for Current-skM and Control-skM
% calculate_Hq needs documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hq(k+1) = calculate_Hq(X(1:k+1, 1:D), V, 0);
hqp(k+1) = calculate_Hq(Vp, 0, 0);


% initialize iCVI values for Current-skM and Control-skM
iCVs_ff(k+1) = XB_batch(X(1:k+1, 1:D), V(1:k, :));
iCVs1_ff(k+1) = XB_batch(X(1:kp, 1:D), Vp(1:kp, :));

% initialize figures showing the process
if record_video > 0
    err_bar_handle1 = plot(0,0,'.');
    err_bar_handle2 = plot(0,0,'.');

    prototype_handlers = [];
        
    subplot(2,3,1,'Parent',p)
    title('Current-sKM')
    xlabel('Feature 1', 'FontSize',8,'FontWeight','bold')
    ylabel('Feature 2', 'FontSize',8,'FontWeight','bold')
    hold on
    axis equal
    
    if record_video > 1
        plot(X(1:k+1, 1), X(1:k+1, 2), '.');
        for i = 1:k
            subplot(2,3,1,'Parent',p)
            prototype_handlers(i) = plot(V(i, 1), V(i, 2), '*k', 'linewidth',5);
        end
    end
    
    prototype_handlers_p = [];
    subplot(2,3,4,'Parent',p)
    title('Control-sKM')
    xlabel('Feature 1', 'FontSize',8,'FontWeight','bold')
    ylabel('Feature 2', 'FontSize',8,'FontWeight','bold')
    hold on
    axis equal
    if record_video > 1
        plot(X(1:k+1, 1), X(1:k+1, 2), '.');
        for i = 1:kp
            prototype_handlers_p(i) = plot(Vp(i, 1), Vp(i, 2), '*k', 'linewidth',5);
        end
    end

    if record_video > 0
        subplot(2,3,2,'Parent',p)
        hold on
        title('Current-sKM iXB_{\lambda}')
        xlabel('Time', 'FontSize',8,'FontWeight','bold')
        ylabel('iXB_{\lambda}', 'FontSize',8,'FontWeight','bold')        
        axis([0 1 0 1])

        subplot(2,3,5,'Parent',p)
        hold on
        title('Control-sKM iXB_{\lambda}')
        xlabel('Time', 'FontSize',8,'FontWeight','bold')
        ylabel('iXB_{\lambda}', 'FontSize',8,'FontWeight','bold')                
        axis([0 1 0 1])

        subplot(2,3,[3,6])
        title('\delta value with 3\sigma rule')
        set(gca,'XTick',[])
    end
end

%% IDB and IXB Masud Initialization
oCentersp = Vp;
GVecp = zeros(kp, D);
CVecp = zeros(1, kp);
PsiVecp = zeros(1, kp);
oMindistp = 0;

oCenters = V;
GVec = zeros(k, D);
CVec = zeros(1, k);
PsiVec = zeros(1, k);
oMindist = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The added prototype numbers
added = [];

% The removed prototype numbers
removed = [];

% initialize the horizon on the figure for iCVI values
horizon = 30;

% initialie the one sided 3-sigma test value
sigma_count = 1.5;


ClusterLabelDic = [1,1];

% initialize the raduis of the first cluster
Radiis = [10];

% initialize clustering_snap_shot to be used to show the evolution of
% clusters over time
clustering_snap_shot.prototypes = {};
clustering_snap_shot.data = {};
snapshot_time = n/4;


for i = k+1+1 : n
    
    if rem(i, snapshot_time) == 0;
        clustering_snap_shot.prototypes = [clustering_snap_shot.prototypes, {V}];
        clustering_snap_shot.data = [clustering_snap_shot.data, {X(1:i,:)}];
    end
    
    if rem(i, 1000)== 0
        i
    end
    %     if rem(i,200) == 0
    %         milad = 1;
    %         past_time = current;
    %         current = now;
    %         current - past_time
    %     end
    if i >= horizon
        start = i-horizon+1;
    else
        start = 1;
    end
    
    xq1 = X(i, 1:D);
    
    if record_video > 1
        subplot(2,3,1,'Parent',p)
        plot(xq1(:, 1), xq1(:, 2), '.');
        %         drawnow();
        
        subplot(2,3,4,'Parent',p)
        plot(xq1(:, 1), xq1(:, 2), '.');
        %         drawnow();
    end
    
    %% Cluster Update ...
    [ label, V, NIs, ~, Radiis] = KMeansStep(xq1, V, NIs, Radiis);
    ClusterIndices(i) = ClusterLabelDic((ClusterLabelDic(:, 1) == label), 2);
%     ClusterIndices(i) = label;
    if i == 4475
        milad = 1;
    end
    if(ClusterIndices(i) == 0)
        milad = 1;
    end    
    uq = calculate_cluster_membership( xq1, V, 'fuzzy');
    
    [ ClusterIndices1(i), Vp, NIsp, ~, ~] = KMeansStep(xq1, Vp, NIsp, ones(100)); % V^{prime}
    uq1 = calculate_cluster_membership( xq1, Vp, 'fuzzy');
    
    if record_video > 1
        subplot(2,3,1,'Parent',p)
        delete(prototype_handlers(label))
        prototype_handlers(label) = plot(V(label, 1), V(label, 2), '*k', 'linewidth',5);
        
        subplot(2,3,4,'Parent',p)
        delete(prototype_handlers_p(ClusterIndices1(i)))
        prototype_handlers_p(ClusterIndices1(i)) = plot(Vp(ClusterIndices1(i), 1), Vp(ClusterIndices1(i), 2), '*k', 'linewidth',5);
    end
    
    %% Calculate iXB_lambda directly
    [iCVs_ff(i), hq(i), uq] = IXB_ff_milad( xq1, V, hq(i-1), iCVs_ff(i-1), lambda, size(V, 1), 'fuzzy');
    %     [iCVs(i), hq(i)] = IXB_milad( xq1, V, hq(i-1), iCVs(i-1), size(V, 1), i, 'fuzzy');
    
    [iCVs1_ff(i), hqp(i), uq1] = IXB_ff_milad( xq1, Vp, hqp(i-1), iCVs1_ff(i-1), lambda, size(Vp, 1), 'fuzzy');
    %     [iCVs1(i), hqp(i)] = IXB_milad( xq1, Vp, hqp(i-1), iCVs1(i-1), size(Vp, 1), i, 'fuzzy');
    
    %% Calculate iXB_lambda from complex forumula (Fuzzy Within Cluster Dispersion).
    %     [ XB,Mindist,GVec,CVec,PsiVec] = IXB( xq1, uq, oCenters, V, 0.9, GVec,CVec,PsiVec,size(V, 1), i,oMindist);
    %     oMindist = Mindist;
    %     oCenters = V;
    %     iCVs_ff(i) = XB;
    %     [ XB1,Mindistp,GVecp,CVecp,PsiVecp] = IXB( xq1, uq1, oCentersp, Vp, 0.9, GVecp,CVecp,PsiVecp,size(Vp, 1), i,oMindistp);
    %     oMindistp = Mindistp;
    %     oCentersp = Vp;
    %     iCVs1_ff(i) = XB1;
    
    %%  Calculate the Davies-Boulding index
    %     [ DB1,GVec,CVec,PsiVec] = IDB( xq1, uq, oCenters, V, 0.9, GVec,CVec,PsiVec,size(V, 1));
    %     oCenters = V;
    %     iCVs_ff(i) = DB1;
    %     [ DBp1,GVecp,CVecp,PsiVecp] = IDB( xq1, uq1, oCentersp, Vp, 0.9, GVecp,CVecp,PsiVecp,size(Vp, 1));
    %     oCentersp = Vp;
    %     iCVs1_ff(i) = DBp1;
    
    
    %% Calculate Within_Cluster_Dispersion with Masuds code
    [ iXB_ff(i),Mindist,GVec,CVec,PsiVec] = IXB( xq1, uq, oCenters, V, 0.9, GVec,CVec,PsiVec,size(V, 1), i,oMindist);
    oCenters = V;
    oMindist = Mindist;
    cohesions = CVec ./ NIs';
    
    
%     [ iCVs1_ff(i),Mindistp,GVecp,CVecp,PsiVecp] = IXB( xq1, uq1, oCentersp, Vp, 0.9, GVecp,CVecp,PsiVecp,size(Vp, 1), i,oMindistp);
%     oCentersp = Vp;
%     oMindistp = Mindistp;

    
    
    %     if record_video == 1
    %         subplot(2,3,4,'Parent',p)
    %         delete(cohesion_handlers(ClusterIndices1(i)))
    %         cohesion_handlers(ClusterIndices1(i)) = text(Vp(ClusterIndices1(i), 1)-0.1, Vp(ClusterIndices1(i), 2)-0.2, num2str(cohesions(ClusterIndices1(i))));
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if record_video > 0
        subplot(2,3,2,'Parent',p)
        %         axis([start i , 0 max(max(max(iCVs_ff(start:i)), max(iCVs1_ff(start:i))), max(max(iCVs(start:i)), max(iCVs1(start:i))))])
        plot([i-1:i], iCVs_ff(i-1:i), '-b', 'linewidth', 2, 'DisplayName','iXB_{\lambda}')
        if ishandle(err_bar_handle1)
            delete(err_bar_handle1)
        end
        err_bar_handle1 = errorbar(i-1, iCVs_ff_means(i-1), sigma_count * sqrt(iCVs_ff_variance(i-1-emerging_cluster_capacity)),'x', 'linewidth',1.5, 'Color', 'k')
        %         plot([i-1:i], iCVs(i-1:i), 'r--', 'linewidth', 2, 'DisplayName','iXB')
        
        the_max1 = max(max(max(iCVs_ff(start:i)), max(iCVs1_ff(start:i))), max(max(iCVs(start:i)), max(iCVs1(start:i))));
%         the_max2 = max(max(max(iCVs_ff(start:i)), max(iCVs1_ff(start:i))), max(max(iCVs(start:i)), max(iCVs1(start:i))));
        the_max2 = max(iCVs1_ff_means(i-1) + sigma_count * sqrt(iCVs1_ff_variance(i-1-emerging_cluster_capacity)), ...
            iCVs_ff_means(i-1) + sigma_count * sqrt(iCVs_ff_variance(i-1-emerging_cluster_capacity)));
        
        axis([start i , 0 max(the_max1, the_max2)])
        
        subplot(2,3,5,'Parent',p)
        plot([i-1:i], iCVs1_ff(i-1:i), '-b', 'linewidth', 2, 'DisplayName','iXB_{\lambda}')
        if ishandle(err_bar_handle2)
            delete(err_bar_handle2)
        end
        err_bar_handle2 = errorbar(i-1, iCVs1_ff_means(i-1), sigma_count * sqrt(iCVs1_ff_variance(i-1-emerging_cluster_capacity)),'x', 'linewidth',1.5, 'Color', 'k')
        %         plot([i-1:i], iCVs1(i-1:i), 'r--', 'linewidth', 2, 'DisplayName','iXB')
        
        the_max1 = max(max(max(iCVs_ff(start:i)), max(iCVs1_ff(start:i))), max(max(iCVs(start:i)), max(iCVs1(start:i))));
%         the_max2 = max(max(max(iCVs_ff(start:i)), max(iCVs1_ff(start:i))), max(max(iCVs(start:i)), max(iCVs1(start:i))));
        the_max2 = max(iCVs1_ff_means(i-1) + sigma_count * sqrt(iCVs1_ff_variance(i-1-emerging_cluster_capacity)), ...
            iCVs_ff_means(i-1) + sigma_count * sqrt(iCVs_ff_variance(i-1-emerging_cluster_capacity)));
        
        
        
        
        axis([start i , 0 max(the_max1, the_max2)])
    end
    
    if i == 307
        milad = 1;
    end
    if i > emerging_cluster_capacity + 1
        %% Calculate delta
        delta(i) = iCVs1_ff(i) - iCVs_ff(i);
        if record_video > 0
            if i > stab_period
                subplot(2,3,[3,6], 'Parent',p)
                hold on
                cla
                errorbar(delta_means(i-1-emerging_cluster_capacity),sigma_count * sqrt(delta_vars(i-1-emerging_cluster_capacity)),'x', 'linewidth',1.5)
                plot(1, delta(i), 'xk', 'linewidth',5);
                hold off
                
% Horizontal Errorbar                
%                 errorbar(delta_means(i-1-emerging_cluster_capacity),1,sigma_count * sqrt(delta_vars(i-1-emerging_cluster_capacity)),'horizontal', 'x', 'linewidth',1.5)
%                 hold on
%                 plot(delta(i), 1, 'xk', 'linewidth',5);

% Extract subplots into new standalone figure.
% hAx = findobj('type', 'axes');
% for iAx = 1:length(hAx)
% fNew = figure;
% hNew = copyobj(hAx(iAx), fNew);
% % Change the axes position so it fills whole figure
% set(hNew, 'pos', [0.23162 0.2233 0.72058 0.63107])
% end
            end
        end
        
        % For adding prototypes
        if  i > 15 && (delta(i) >= comparing_mean + sigma_count * sqrt(comparing_variance)) && ...
                (iCVs1_ff(i) >= iCVs1_ff_means(i-1) + sigma_count * sqrt(iCVs1_ff_variance(i-1)))
            cen = xq1;
            cluster_center_dists = pdist2(cen, V);
            [min_dist, min_ind] = min(cluster_center_dists );
            SF = squareform(pdist(Vp));
            SF(SF == 0 ) = NaN;
            [v] = min(SF(:));
            if min_dist < v
                % do not add anything.
            else
                added = [added, i];
                Radiis = [Radiis, min_dist];
                V = [V; cen];
                CVec = [CVec, 0];
                GVec = [GVec; zeros(1, D)];
                PsiVec = [PsiVec, 0];
                oCenters = [oCenters; cen];
                cohesions = [cohesions, 0];
                NIs = [NIs; 1];
                
                %%%%%%%% For the label dictionary
                label = max(ClusterLabelDic(:, 2)) + 1;
                ClusterLabelDic = [ClusterLabelDic; [size(V, 1), label]];
                
                % update these for control
                Vp = [Vp; cen];
                CVecp = [CVecp, 0];
                GVecp = [GVecp; zeros(1, D)];
                PsiVecp = [PsiVecp, 0];
                oCentersp = [oCentersp; cen];
                
                NIsp = [NIsp; 1];
                
                
                
                if record_video > 1
                    subplot(2,3,1,'Parent',p)
                    prototype_handlers(length(prototype_handlers) + 1) = plot(cen(:, 1), cen(:, 2), '*k', 'linewidth',5);
                    
                    subplot(2,3,4,'Parent',p)
                    prototype_handlers_p(length(prototype_handlers_p) + 1) = plot(xq1(:, 1), xq1(:, 2), '*k', 'linewidth',5);
                    %                 cohesion_handlers(length(cohesion_handlers) + 1) = text(Vp(ClusterIndices1(i), 1), Vp(ClusterIndices1(i), 2), '0');
                end
                
                if record_video > 0
                    currFrame = getframe(gcf);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                    writeVideo(vidObj,currFrame);
                end
            end
            
        end
        
        %% For removing prototypes
        if  i > 15 && (delta(i) <= comparing_mean - sigma_count * sqrt(comparing_variance)) && ...
                (iCVs_ff(i) >= iCVs_ff_means(i-1) + sigma_count * sqrt(iCVs_ff_variance(i-1)))
            if size(V, 1) > 1
                %           This is on V, because the jump was seen in simple kmeans
                %% calculate the local XB and choose the prototype with the max
                SF = squareform(pdist(V));
                SF(SF == 0 ) = NaN;
                closest_prototype_dists = [];
                for j = 1:size(V, 1)
                    closest_prototype_dists = [closest_prototype_dists, min(SF(j, :))];
                end
                local_iXB = cohesions ./ closest_prototype_dists;
                [max_val, max_ind] = max(local_iXB);
                [min_val, min_ind] = min(SF(max_ind, :));
                to_be_merged_ind = min_ind;
                to_be_removed_ind = max_ind;
                
                %% Choose the last moved prototype
                %                 last_moved = ClusterIndices(i);
                %                 temp_V = V;
                %                 temp_V(ClusterIndices(i), :) = [];
                %                 dists = pdist2(temp_V, V(ClusterIndices(i), :));
                %                 [min_val, min_ind] = min(dists);
                %                 to_be_removed_ind = max(ClusterIndices(i), min_ind);
                %                 to_be_merged_ind = min(ClusterIndices(i), min_ind);
                
                
                %% Choose the two closest prototypes to merger
%                     SF = squareform(pdist(V));
%                     SF(SF == 0 ) = NaN;
%                     [v] = min(SF(:));
%                     [row,~] = find(SF==v);
%                     to_be_removed_ind = max(row);
%                     to_be_merged_ind = min(row);



                new_v = (NIs(to_be_removed_ind) * V(to_be_removed_ind, :) + NIs(to_be_merged_ind) * V(to_be_merged_ind, :)) ./ (NIs(to_be_merged_ind) + NIs(to_be_removed_ind));
                new_vp = (NIsp(to_be_removed_ind+1) * Vp(to_be_removed_ind+1, :) + NIsp(to_be_merged_ind+1) * Vp(to_be_merged_ind+1, :)) ./ (NIsp(to_be_merged_ind+1) + NIsp(to_be_removed_ind+1));
                
                if to_be_removed_ind > 1 && ...
                        (NIs(to_be_merged_ind) < 5000000 && NIs(to_be_removed_ind) < 5000000)
                    
                    removed = [removed, i];
                    NIs(to_be_merged_ind) = NIs(to_be_merged_ind) + NIs(to_be_removed_ind);
                    V(to_be_merged_ind, :) = new_v;
                    cohesions(to_be_merged_ind) = (NIs(to_be_removed_ind) * cohesions(to_be_removed_ind) + NIs(to_be_merged_ind) * cohesions(to_be_merged_ind)) ./ (NIs(to_be_merged_ind) + NIs(to_be_removed_ind));
                    CVec(to_be_merged_ind) = (NIs(to_be_removed_ind) * CVec(to_be_removed_ind) + NIs(to_be_merged_ind) * CVec(to_be_merged_ind)) ./ (NIs(to_be_merged_ind) + NIs(to_be_removed_ind));
                    GVec(to_be_merged_ind,:) = (NIs(to_be_removed_ind) * GVec(to_be_removed_ind) + NIs(to_be_merged_ind) * GVec(to_be_merged_ind)) ./ (NIs(to_be_merged_ind) + NIs(to_be_removed_ind));
                    PsiVec(to_be_merged_ind) = (NIs(to_be_removed_ind) * PsiVec(to_be_removed_ind) + NIs(to_be_merged_ind) * PsiVec(to_be_merged_ind)) ./ (NIs(to_be_merged_ind) + NIs(to_be_removed_ind));
                    
                    NIsp(to_be_merged_ind+1) = NIsp(to_be_merged_ind+1) + NIsp(to_be_removed_ind+1);
                    Vp(to_be_merged_ind+1, :) = new_vp;

                    CVecp(to_be_merged_ind+1) = (NIsp(to_be_removed_ind+1) * CVecp(to_be_removed_ind+1) + NIsp(to_be_merged_ind+1) * CVecp(to_be_merged_ind+1)) ./ (NIsp(to_be_merged_ind+1) + NIsp(to_be_removed_ind+1));
                    GVecp(to_be_merged_ind+1,:) = (NIsp(to_be_removed_ind+1) * GVecp(to_be_removed_ind+1) + NIsp(to_be_merged_ind+1) * GVecp(to_be_merged_ind+1)) ./ (NIsp(to_be_merged_ind+1) + NIsp(to_be_removed_ind+1));
                    PsiVecp(to_be_merged_ind+1) = (NIsp(to_be_removed_ind+1) * PsiVecp(to_be_removed_ind+1) + NIsp(to_be_merged_ind+1) * PsiVecp(to_be_merged_ind+1)) ./ (NIsp(to_be_merged_ind+1) + NIsp(to_be_removed_ind+1));

                    %%%%%%%% For the label dictionary
                    ClusterLabelDic(ClusterLabelDic(:, 1) ==  to_be_removed_ind, 1) = 1000 * to_be_removed_ind;
                    ClusterLabelDic((ClusterLabelDic(:, 1) > to_be_removed_ind) & (ClusterLabelDic(:, 1) < 1000), 1)  = ClusterLabelDic((ClusterLabelDic(:, 1) > to_be_removed_ind) & (ClusterLabelDic(:, 1) < 1000), 1) - 1;
                    Radiis(to_be_removed_ind) = [];
                    NIs(to_be_removed_ind) = [];
                    V(to_be_removed_ind, :) = [];
                    cohesions(to_be_removed_ind) = [];
                    CVec(to_be_removed_ind) = [];
                    GVec(to_be_removed_ind,:) = [];
                    PsiVec(to_be_removed_ind) = [];
                    
                    
                    NIsp(to_be_removed_ind+1) = [];
                    Vp(to_be_removed_ind+1, :) = [];
                    
                    CVecp(to_be_removed_ind+1) = [];
                    GVecp(to_be_removed_ind+1,:) = [];
                    PsiVecp(to_be_removed_ind+1) = [];
                    
                    if record_video > 1
                        
                        delete(prototype_handlers(to_be_removed_ind))
                        prototype_handlers(to_be_removed_ind) = [];
                        delete(prototype_handlers_p(to_be_removed_ind+1))
                        prototype_handlers_p(to_be_removed_ind+1) = [];
                        %                     delete(cohesion_handlers(to_be_removed_ind+1))
                        %                     cohesion_handlers(to_be_removed_ind+1) = [];
                    end
                    if record_video > 0
                        currFrame = getframe(gcf);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                        writeVideo(vidObj,currFrame);
                    end
                end
            end
        end
        
        
        delta_vars(i) = lamda_delta * delta_vars(i-1) + (((1-lamda_delta^2)/2) * ((delta(i) - delta_means(i-1)) ^ 2));
        comparing_variance = delta_vars(i-emerging_cluster_capacity);
        delta_means(i) = (lamda_delta) * delta_means(i-1) + (1-lamda_delta) * delta(i);
        comparing_mean = delta_means(i-emerging_cluster_capacity);
        
        iCVs_ff_variance(i) = lamda_delta * iCVs_ff_variance(i-1) + (((1-lamda_delta^2)/2) * ((iCVs_ff(i) - iCVs_ff_means(i-1)) ^ 2));
        iCVs_ff_means(i) = (lamda_delta) * iCVs_ff_means(i-1) + (1-lamda_delta) * iCVs_ff(i);
        
        iCVs1_ff_variance(i) = lamda_delta * iCVs1_ff_variance(i-1) + (((1-lamda_delta^2)/2) * ((iCVs1_ff(i) - iCVs1_ff_means(i-1)) ^ 2));
        iCVs1_ff_means(i) = (lamda_delta) * iCVs1_ff_means(i-1) + (1-lamda_delta) * iCVs1_ff(i);
        
        %     t11 = var(delta(kp+1:i))
        %     t22 = mean(delta(kp+1:i))
        
        
        %     i = i + 1;
        
        
        if record_video > 0 % && rem(i, 1) == 0
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end
    end
    
end

end

