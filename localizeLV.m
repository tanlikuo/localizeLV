function [bw4d_LV, bw3d_LV] = localizeLV(im4d, res_xy)
%localizeLV Localizes cardiac LV by estimating 3D centrepoint
%
%   Estimates the cardiac left ventricle (LV) centre by roughly 
%   segmenting the LV bloodpool. Only works on MRI short-axis steady state 
%   free precession (SSFP) images.
%
%   [bw3d_LV, bw4d_LV] = localizeLV(im4d, res_xy) 
%    in : im4d is a 4D volume with dimensions y,x,t,z
%    in : res_xy is the x/y resolution in mm/pixel
%    out: bw4d_LV is a 4D binary mask representing the estimated LV bloodpool
%    out: bw3d_LV is adapted from bw4d_LV, but averaged across time
%
%   NOTE: this function requires the MATLAB Image Processing and Statistics
%         toolboxes to be installed

% NOTE: This is based off a cleaned up version of TLK_findLV (dated 20160626).

assert(nargin == 2, 'Function requires two input arguments.');
assert(ndims(im4d) == 4, 'Input volume must be 4D (y,x,t,z)');
assert(isscalar(res_xy), 'x/y resolution should be scalar value');

vol = double(im4d);
[sz_y, sz_x, sz_t, sz_z] = size(vol);

% get the DC (average) volume and the "minimum" (seed) volume
vol_dc = squeeze(mean(vol, 3));
vol_min = squeeze(min(vol, [], 3));

% calculate magnitude and phase of 1st harmonic
[h1_mag, h1_phi] = goCalcH1(vol);

% use 1st harmonic motion to get initial mask
mask_valid = goThresholdMotion(h1_mag);

% get value for threshold segmentation
[mask_4d, thres_hilo, thres_h1] = goThresholdVol(vol, h1_mag, mask_valid);
mask_dc = (vol_dc > thres_hilo) & mask_valid;
mask_min = (vol_min > thres_hilo) & mask_valid;

% localize the lung
[Cxy_lung, bw3d_cvxhull] = goCalcLungCxy(vol_dc, thres_hilo);

% separate and group "core" objects and "full" objects
cc_full = bwconncomp(mask_4d, 8);
[info_full, px2d_full] = goParseBinaryObjs(cc_full);
cc_core = bwconncomp(mask_min, 8);
[info_core, px2d_core] = goParseBinaryObjs(cc_core);
% filter out noise / tiny objects (min 6 px)
arr = [ info_full.area ] < 6;
info_full(arr) = [];  px2d_full(arr) = [];  cc_full.PixelIdxList(arr) = [];  cc_full.NumObjects = sum(~arr);
arr = [ info_core.area ] < 6;
info_core(arr) = [];  px2d_core(arr) = [];  cc_core.PixelIdxList(arr) = [];  cc_core.NumObjects = sum(~arr);

% map out relationships between the "core" and "full" objects
[core2full, full2core, f2c_num, core2core] = mapCore2Full(info_core, info_full, px2d_core, cc_full, sz_t);

% collect the connecting core objects into groups
[core_conn, idx_conn, conn_grp] = getCoreConnections(f2c_num, core2full, core2core, info_core);

% merge a subset of connected core objects together (targetting single
% objects which "break off" small chunks due to noise)
[cc_core_merged, cc_full_merged] = goMergeObjs(cc_core, cc_full, info_core, conn_grp, core2full);

% redo core <-> full object mapping. All remaining core-to-core connections
% should be between "significant" objects
cc_core = cc_core_merged;
cc_full = cc_full_merged;
[info_core, px2d_core] = goParseBinaryObjs(cc_core);
[info_full, px2d_full] = goParseBinaryObjs(cc_full);
[core2full, full2core, f2c_num, core2core] = mapCore2Full(info_core, info_full, px2d_core, cc_full, sz_t);

% redo collection of the connecting core objects into groups
[core_conn, idx_conn, conn_grp] = getCoreConnections(f2c_num, core2full, core2core, info_core);

% split "full" objects in the remaining cases where there are still core
% connections
cc_full_split = goSplitObjs(cc_full, px2d_full, info_core, px2d_core, conn_grp, core2full);

% make final core <-> full object mapping, there should be no more connections / leaks between core objects 
cc_full = cc_full_split;
[info_full, px2d_full] = goParseBinaryObjs(cc_full);
[core2full, full2core, f2c_num, core2core] = mapCore2Full(info_core, info_full, px2d_core, cc_full, sz_t);

% SCORE: get mean area and calculate normalized area range
area_full_mean = zeros(length(info_core), 1);
scores_area = zeros(length(info_core), 1);
for n_core = 1:length(info_core)
    areas = [ info_full(core2full(n_core,:)).area ];
    area_full_mean(n_core) = mean(areas);
    scores_area(n_core) = max(areas) - min(areas);
end
scores_area = scores_area ./ area_full_mean;  % normalize by mean area

% SCORE: get average circularity / eccentricity
cc = struct('Connectivity', 8, 'ImageSize', [sz_y  sz_x], 'NumObjects', length(info_full), 'PixelIdxList', []);
cc.PixelIdxList = px2d_full;
rp = regionprops(cc, 'Eccentricity');
eccen_full = [ rp.Eccentricity ];
scores_circ = zeros(length(info_core), 1);
for n_core = 1:length(info_core)
    scores_circ(n_core) = mean(eccen_full(core2full(n_core,:)));
end
scores_circ = 1 - scores_circ;  % invert so pure circle == 1

% SCORE: get distance from "dark centroid"
Cx = permute(Cxy_lung(:,1), [3 2 1]);  Cy = permute(Cxy_lung(:,2), [3 2 1]);
[xx, yy, ~] = meshgrid(1:sz_x, 1:sz_y, 1:sz_z);
xx = bsxfun(@minus, xx, Cx);  yy = bsxfun(@minus, yy, Cy);
dist_map = sqrt( xx.^2 + yy.^2 );
dist_map = bsxfun(@minus, mean(mean(dist_map,1),2), dist_map);  % invert and normalize to average (possible) distance
dist_map(dist_map < 0) = 0;
% get minimum distance from each "core object" to inverse weighted centroid
scores_dist = zeros(length(info_core), 1);
for n_core = 1:length(info_core)
    scores_dist(n_core) = max(dist_map(cc_core.PixelIdxList{n_core}));
end

% SCORE (filter): check max centrepoint displacement in time
scores_ctrpnt = zeros(length(info_core), 1);
for n_core = 1:length(info_core)
    Cxy_all = vertcat( info_full(core2full(n_core,:)).Cxy );
    Cxy_med = median(Cxy_all, 1);
    Cxy_dist = sqrt( (Cxy_all(:,1) - Cxy_med(1)).^2 + (Cxy_all(:,2) - Cxy_med(2)).^2 );
    scores_ctrpnt(n_core) = max(Cxy_dist);
end
core_filt = false(length(info_core), 1);
core_filt(scores_ctrpnt > (15 / res_xy)) = true;  % threshold at 15 mm

% get combined score, and apply some thresholds
scores = scores_area .* scores_circ .* scores_dist;
min_area = 400 / res_xy / res_xy;  % set a 400 mm^2 minimum area threshold (value determined empirically)
scores(area_full_mean < min_area) = 0;
scores(core_filt) = 0;  % zero scores for objects that have been filtered out

% create "mean 3D object" for chaining calculations
px_mobj = cell(length(info_core), 1);
for n_core = 1:length(info_core)  % find pixels that are "active" at least 50% of cardiac phases
    px2d = vertcat( px2d_full{core2full(n_core, :)} );
    px_range = min(px2d):max(px2d);
    px_cnt = histc(px2d, px_range);
    px2d = px_range(px_cnt >= sz_t / 2);
    px_mobj{n_core} = px2d' + sz_x*sz_y*(info_core(n_core).slice - 1);  % apply z-offset
end
cc_mobj = struct('Connectivity', 8, 'ImageSize', [sz_y  sz_x  sz_z], 'NumObjects', length(info_core), 'PixelIdxList', []);
cc_mobj.PixelIdxList = px_mobj;
[info_mobj, px2d_mobj] = goParseBinaryObjs(cc_mobj);

% GROUPING: chain / group objects together
chains = goFindChains(scores, info_mobj, px2d_mobj, cc_mobj, core_filt, res_xy);
scores_chain = zeros(length(chains), 1);
for n_chn = 1:length(chains)
    scores_chain(n_chn) = sum(scores(chains{n_chn}));
end
idx_chain = find(scores_chain == max(scores_chain), 1, 'first');  % highest scoring chain / group

idx = core2full(chains{idx_chain}, :);  idx = idx(:);  % get matching 4D objects
arr = false(length(info_full), 1);  arr(idx) = true;
cc = cc_full;  cc.PixelIdxList(~arr) = {[]};
bw4d_LV = labelmatrix(cc) > 0;
bw3d_LV = squeeze( mean(bw4d_LV, 3) ) >= 0.5;

return  % convenience for setting breakpoints



function [h1_mag, h1_phi] = goCalcH1(vol)
% perform 1-cycle mini Fourier transform (first harmonic)
% can also use built-in FFT function, but this seemed more straightforward
[~, ~, sz_t, ~] = size(vol);
theta = shiftdim(2*pi*((1:sz_t)'-1) / sz_t, -2);  % sinusoid of frequency 1, along 3rd dimension
h1_cos = squeeze(sum(bsxfun(@times, vol,  cos(theta)), 3));
h1_sin = squeeze(sum(bsxfun(@times, vol, -sin(theta)), 3));
h1_mag = sqrt(h1_cos.^2 + h1_sin.^2);  % magnitude
h1_phi = atan2(-h1_cos, h1_sin);  % phase



function mask = goThresholdMotion(vol)
% use 1st harmonic motion to get initial mask

[sz_y, sz_x, sz_z] = size(vol);
vol = imfilter(vol, fspecial('average'));
pd = fitdist(vol(:), 'Exponential');
vol(vol < icdf(pd, 0.95)) = 0;  % threshold at p=0.95
thres = prctile(vol(vol > 0), 99);
vol(vol > thres) = thres;  % threshold at max 99%

ctr3d_xy = zeros(1, 3);
[xx, yy, zz] = meshgrid(1:sz_x, 1:sz_y, 1:sz_z);
max_iter = 5;  % maximum allowed iterations
for n = 1:max_iter  
    % Get 2D weighted centroid (should be identical to regionprops)
    xx_wt = xx .* vol;  yy_wt = yy .* vol;  zz_wt = zz .* vol;
    ctr_xy = [ squeeze( sum(sum(xx_wt,1),2) )  squeeze( sum(sum(yy_wt,1),2) ) ];
    vol_sum = squeeze( sum(sum(vol,1),2) );
    ctr_xy = ctr_xy ./ [ vol_sum  vol_sum ];
    
    slices = 1:sz_z;
    arr = all(~isnan(ctr_xy), 2);  % check for invalid cases (i.e., no valid pixels)
    
    % NOTE: not very sure if this is the best way... it appears that svd is
    % more appropriate?
    px = polyfit(slices(arr), ctr_xy(arr,1)', 1);  % get x/y line of best fit
    py = polyfit(slices(arr), ctr_xy(arr,2)', 1);

    ctr_xy_fit = [ polyval(px, 1:sz_z)'  polyval(py, 1:sz_z)' ];  % fitted centrepoints
    
    % get 3D weighted centroid and compare distance from previous 3D centroid
    % (should be identical to MATLAB's regionprops)
    ctr3d_new = [ sum(xx_wt(:))  sum(yy_wt(:))  sum(zz_wt(:)) ] / sum(vol(:));
    ctr3d_diff = ctr3d_xy - ctr3d_new;
    ctr3d_diff = sqrt(sum(ctr3d_diff.^2));
    ctr3d_xy = ctr3d_new;
    
    if ctr3d_diff < 1;  % stop iteration when 3D centroid differs by < 1 px
        break;
    elseif n == max_iter
        fprintf('Warning: 3D centroid did not stabilize\n');
    end  
    
    % calculate the distance from every voxel to its 2D centroid
    diff_xx = bsxfun(@minus, xx, shiftdim(ctr_xy_fit(:,1), -2));
    diff_yy = bsxfun(@minus, yy, shiftdim(ctr_xy_fit(:,2), -2));
    vol_dist = sqrt(diff_xx.^2 + diff_yy.^2);
    
    % calculate combined intensity & distance index (higher intensity &
    % shorter distance is better)
    vol_comb = vol_dist .* (max(vol(:)) - vol);  % invert intensity to be compatible w/ distance
    vol_samp = vol_comb((vol > 0) & (vol_comb > 0));  % remove ommited voxels. Also, a few distributions don't allow zero values

    % the distribution is strictly positive, with a long tail to the right
    % a bunch of distributions match reasonably well, in order:
    % LogLogistic, Gamma, Weibull, Nakagami, Rayleigh
    % The first 2 fit better, while the last 3 are more "rounded/damped"
    pd = fitdist(vol_samp, 'LogLogistic');
    vol(vol_comb > pd.icdf(0.9)) = 0;  % threshold at one-tailed p = 0.9
end

% construct final mask
xx = xx(:,:,1);  yy = yy(:,:,1);
mask = vol > 0;
for n = 1:sz_z
    mask2d = mask(:,:,n);
    if ~any(mask2d(:));  continue;  end
    
    Cxy_all = [ mean(xx(mask2d))  mean(yy(mask2d)) ];
    cc = bwconncomp(mask2d);
    if cc.NumObjects == 1;  continue;  end
    
    px_areas = cellfun(@numel, cc.PixelIdxList);
    cc_dist = zeros(cc.NumObjects, 1);
    for n_cc = 1:cc.NumObjects
        Cxy = [ mean(xx(cc.PixelIdxList{n_cc}))  mean(yy(cc.PixelIdxList{n_cc})) ];
        cc_dist(n_cc) = sqrt(sum( (Cxy - Cxy_all).^2 ));
    end
    
    score = ( px_areas' / max(px_areas) + 1 - (cc_dist / max(cc_dist)) ) / 2;
    cc_all = cc;
    cc_all.NumObjects = 1;
    cc_all.PixelIdxList = { vertcat( cc.PixelIdxList{score >= .01} ) };
    rp = regionprops(cc_all, 'BoundingBox', 'ConvexImage', 'ConvexArea');
    mask2d = false(sz_y, sz_x);
    bb = round( rp.BoundingBox );
    bb_x = [ bb(1)  bb(1) + bb(3) - 1 ];
    bb_y = [ bb(2)  bb(2) + bb(4) - 1 ];
    mask2d(bb_y(1):bb_y(2), bb_x(1):bb_x(2)) = rp.ConvexImage;
    
    mask(:,:,n) = mask2d;
end



function [bw4d, thres_hilo, thres_h1] = goThresholdVol(vol4d, vol3d_h1mag, bw3d_valid)
% first get the binary volume of "strong moving objects" from H1 magnitude
data = vol3d_h1mag(bw3d_valid);
data_max = prctile(data, 99);
data(data > data_max) = data_max;  % threshold at 99%
level = graythresh(data / data_max);
thres_h1 = level * data_max;
bw3d_h1mag = (vol3d_h1mag >= thres_h1)  &  bw3d_valid;

% now get the intensity threshold by iteratively scanning the maximum
% intensity projection volume against the binary "strong moving objects"
% volume until we get a match
vol3d_max = squeeze(max(vol4d, [], 3));
low_high = prctile(vol3d_max(bw3d_valid), [1 99]);
levels = linspace(low_high(1), low_high(2), 256);  % keep to 256 levels for speed
data = vol3d_max(bw3d_h1mag);
data_cnt = numel(data);
for thres_hilo = levels(:)'  % iterate until at most 95% coverage
    if (sum(data >= thres_hilo) / data_cnt) <= 0.95;  break;  end
end
sz_t = size(vol4d, 3);
bw4d_valid = repmat(permute(bw3d_valid, [1 2 4 3]), [1 1 sz_t 1]);
bw4d = (vol4d >= thres_hilo)  &  bw4d_valid;



function [core2full, full2core, f2c_num, core2core] = mapCore2Full(info_core, info_full, px2d_core, cc_full, sz_t)
% map out relationships between the "core" and "full" objects
core2full = zeros(length(info_core), sz_t);
full2core = zeros(length(info_full), 30);  % maximum 30 "core" objects in each "full" object
f2c_num = zeros(length(info_full), 1);  % counter for full2core
L4d = labelmatrix(cc_full);
for n_core = 1:length(info_core)
    n_z = info_core(n_core).slice;
    for n_t = 1:sz_t
        % map "core objects" to "full objects"
        L = L4d(:,:,n_t,n_z);
        n_full = L(px2d_core{n_core}(1));  % just one pixel should be enough
        core2full(n_core, n_t) = n_full;
        
        % map "full objects" to "core objects"
        f2c_num(n_full) = f2c_num(n_full) + 1;
        full2core(n_full, f2c_num(n_full)) = n_core;
    end
end

% trim full2core size
max_objs = max( sum(full2core ~= 0, 2) );
full2core = full2core(:, 1:max_objs);

% map out which "core" objects connect or "leak" to each other
core2core = cell(length(info_core), sz_t);
idx_conn = find(f2c_num > 1);  % find "full" objects with >1 "core" counterparts
for n_full = idx_conn(:)'
    n_t = info_full(n_full).phase;
    idx = full2core(n_full, :);  idx = idx(idx > 0);  % get list of connected "core" objects
    for n_core = idx(:)'
        core2core{n_core, n_t} = idx(idx ~= n_core);  % add to list (omit parent object)
    end
end



function [core_conn, idx_conn, conn_grp] = getCoreConnections(f2c_num, core2full, core2core, info_core)
% collect the connecting core objects into groups. A previous
% implementation failed to account for 2nd level connections and above:
% i.e., object 3 connects to 4, but 4 also connects to 3 and 6.
% this implementation uses a more exhaustive search which should solve that
core_conn = f2c_num(core2full);  % number of connections at each phase (size: n_cores x n_t)
idx_conn = find( any(core_conn > 1, 2) );  % find cores with phase leaks / extra connections at any phase
conn_grp = cell(length(idx_conn), 1);  % keep track of connected groups
for n_grp = 1:length(idx_conn)  % first collect all "connected" cores
    n_core = idx_conn(n_grp);
    conn_grp{n_grp} = [ n_core  unique([ core2core{n_core, :} ]) ];  % get list of "connected" cores at any phase
end
do_merge = true;
while do_merge  % repeat until no more valid merges found
    do_merge = false;
    slices = cellfun(@(x) info_core(x(1)).slice, conn_grp);  % get slice index of each group
    for n_grp_A = 1:length(conn_grp)-1
        n_z = slices(n_grp_A);
        for n_grp_B = n_grp_A+1:length(conn_grp)
            if n_z ~= slices(n_grp_B);  continue;  end  % if different slice, members can't be connected, so skip
            if any(ismember(conn_grp{n_grp_A}, conn_grp{n_grp_B}))  % if matching members, merge group
                conn_grp{n_grp_A} = unique([ conn_grp{n_grp_A}  conn_grp{n_grp_B} ]);
                conn_grp(n_grp_B) = {[]};
                do_merge = true;
            end
        end
    end
    sz_grp = cellfun(@numel, conn_grp);  % remove empty groups (from merging)
    conn_grp(sz_grp == 0) = [];
end



function [cc_core_merged, cc_full_merged] = goMergeObjs(cc_core, cc_full, info_core, conn_grp, core2full)
% merge a subset of connected core objects together (targetting single
% objects which "break off" small chunks due to noise)
sz_t = cc_full.ImageSize(3);
cc_core_merged = cc_core;
cc_full_merged = cc_full;
for n_grp = 1:length(conn_grp)
    idx_core = conn_grp{n_grp};
    areas = [ info_core(idx_core).area ];
    [areas, idx] = sort(areas, 'descend');  % sort by largest area
    idx_core = idx_core(idx);
    
    arr = areas > areas(1) * 0.1;  % "significant" objects defined as >10% area
    idx_core_big = idx_core(arr);  idx_core_sm = idx_core(~arr);
    if length(idx_core_big) == 1  % only one primary object, merge all others
        idx_merge = { idx_core };
    else
        idx_merge = num2cell( idx_core_big(:) );
        if ~isempty(idx_core_sm)  % merge all "insignificant" objects to closest "significant" objects
            dist_table = zeros(length(idx_core_big), length(idx_core_sm));
            Cxy_sm = vertcat( info_core(idx_core_sm).Cxy );
            for n = 1:length(idx_core_big)
                Cxy = info_core(idx_core_big(n)).Cxy;
                dist_table(n, :) = sqrt( (Cxy_sm(:,1) - Cxy(1)).^2 + (Cxy_sm(:,2) - Cxy(2)).^2 );
            end
            [~, idx_sm2bg] = min(dist_table, [], 1);  % get closest "significant" object
            for n = 1:length(idx_core_sm)  % update merge table
                idx_merge{idx_sm2bg(n)} = [ idx_merge{idx_sm2bg(n)}  idx_core_sm(n) ];
            end
        end
    end
    
    % apply (partial) merge to "core" and "full" objects
    % NOTE: an earlier version tried to handle the complete "full" object
    % merge-and-split here as well, but there're too many complications
    for n = 1:length(idx_merge)
        idx = idx_merge{n};  % 1st index is primary "significant" object
        if length(idx) <= 1;  continue;  end
        
        px = vertcat(  cc_core_merged.PixelIdxList{idx} );
        cc_core_merged.PixelIdxList{idx(1)} = px;
        cc_core_merged.PixelIdxList(idx(2:end)) = {[]};
        
        idx_full = core2full(idx, :);  % get list of affected "full" objects
        for n_t = 1:sz_t  % merge "full" objects in the same phase
            idx = unique(idx_full(1, n_t), 'stable');
            if all(idx(1) == idx(2:end));  continue;  end
            cc_full_merged.PixelIdxList{idx(1)} = vertcat( cc_full_merged.PixelIdxList{idx} );
            cc_full_merged.PixelIdxList(idx(2:end)) = {[]};
        end
    end
end
arr = cellfun(@isempty, cc_core_merged.PixelIdxList);
cc_core_merged.PixelIdxList(arr) = [];
cc_core_merged.NumObjects = sum(~arr);
arr = cellfun(@isempty, cc_full_merged.PixelIdxList);
cc_full_merged.PixelIdxList(arr) = [];
cc_full_merged.NumObjects = sum(~arr);



function cc_full_split = goSplitObjs(cc_full, px2d_full, info_core, px2d_core, conn_grp, core2full)
% split "full" objects in the remaining cases where there are still core
% connections
sz_x = cc_full.ImageSize(2);  sz_y = cc_full.ImageSize(1);  sz_t = cc_full.ImageSize(3);
cc_full_split = cc_full;
[xx, yy] = meshgrid(1:sz_x, 1:sz_y);
for n_grp = 1:length(conn_grp)
    idx_core = conn_grp{n_grp};
    idx_full = core2full(idx_core, :);
    for n_t = 1:sz_t
        idx_phase = idx_full(:, n_t);
        for idx_f = unique(idx_phase)'
            arr = idx_f == idx_phase;
            idx_c = idx_core(arr);
            if length(idx_c) <= 1;  continue;  end  % no "shared / connecting" cores
            
            % split "full" object by calculating per-pixel distance to
            % "core" object centroids, and also 
            px2d = px2d_full{idx_f};
            dist_table = zeros(length(px2d), length(idx_c));
            for n = 1:length(idx_c)
                Cxy = info_core(idx_c(n)).Cxy;
                dist_table(:,n) = sqrt( (xx(px2d) - Cxy(1)).^2 + (yy(px2d) - Cxy(2)).^2 );
                arr = ismember(px2d, px2d_core{idx_c(n)});
                dist_table(arr,n) = -1;  % make sure "full" pixels which coincide with the "core" object are set to minimum
            end
            [~, dist_table] = min(dist_table, [], 2);
            n_full_num = cc_full_split.NumObjects;
            n_full_new = length(idx_c);
            for n = 1:n_full_new
                cc_full_split.PixelIdxList{n_full_num + n} = cc_full_split.PixelIdxList{idx_f}(dist_table == n);
            end
            cc_full_split.PixelIdxList(idx_f) = {[]};
            cc_full_split.NumObjects = n_full_num + n_full_new;
        end
    end
end
arr = cellfun(@isempty, cc_full_split.PixelIdxList);
cc_full_split.PixelIdxList(arr) = [];
cc_full_split.NumObjects = sum(~arr);



function [Cxy_lung, bw3d_valid] = goCalcLungCxy(vol_dc, thres_hilo)
% now try to detect the lung (actually just the dark areas)
[sz_y, sz_x, sz_z] = size(vol_dc);

% calculate the "valid area" by omitting the background air and limbs
bw3d_high = vol_dc >= thres_hilo;
bw3d_valid = bw3d_high;
[xx, yy] = meshgrid(1:sz_x, 1:sz_y);
for n_z = 1:sz_z
    % get centrepoint of binary "bright" image
    bw = bw3d_valid(:,:,n_z);
    Cxy = [ mean(xx(bw))  mean(yy(bw)) ];
    d_map = sqrt( (xx - Cxy(1)).^2 + (yy - Cxy(2)).^2 );
    
    cc = bwconncomp(bw);
    if cc.NumObjects == 0;  continue;  end
    
    % calculate the maximum convex hull area
    cc_full = cc;
    cc_full.NumObjects = 1;
    cc_full.PixelIdxList = { vertcat(cc_full.PixelIdxList{:}) };
    rp = regionprops(cc_full, 'ConvexHull', 'ConvexArea');
    area_max = rp.ConvexArea;
    bw = roipoly(bw, rp.ConvexHull(:,1), rp.ConvexHull(:,2));
    
    % calculate the distance of each "bright" object from the binary centrepoint
    d_list = zeros(cc.NumObjects, 1);
    for n = 1:cc.NumObjects
        d_list(n) = min( d_map(cc.PixelIdxList{n}) );
    end
    [~, idx] = sort(d_list);
    cc.PixelIdxList = cc.PixelIdxList(idx);
    
    % starting from the closest-to-centre object, keep adding objects until
    % the combined convex hull is 75% of the maximum
    cc_sub = cc_full;
    cc_sub.PixelIdxList = { [] };
    for n = 1:cc.NumObjects
        cc_sub.PixelIdxList{1} = vertcat(cc_sub.PixelIdxList{1}, cc.PixelIdxList{n});
        rp = regionprops(cc_sub, 'ConvexHull', 'ConvexArea');
        if rp.ConvexArea / area_max >= 0.75  % area threshold 75%
            bw = roipoly(bw, rp.ConvexHull(:,1), rp.ConvexHull(:,2));
            break
        end
    end
    bw3d_valid(:,:,n_z) = bw;
end

% calculate the weighted centroid of the dark regions
[xx, yy] = meshgrid(1:sz_x, 1:sz_y);
ctr_x = zeros(sz_z, 1);  ctr_y = zeros(sz_z, 2);
ctr_x(:) = sz_x/2;  ctr_y(:) = sz_y/2;  % default to centre of image
for n_z = 1:sz_z
    im = vol_dc(:,:,n_z);
    bw = bw3d_valid(:,:,n_z)  &  ~bw3d_high(:,:,n_z);
    if ~any(bw(:));  continue;  end
    im = max(im(bw)) - im;
    im(~bw) = 0;
    im = im .^ 2;  % emphasize the darker intensities
    xx_wt = xx .* im;  yy_wt = yy .* im;
    ctr_x(n_z) = sum(xx_wt(:)) / sum(im(:));
    ctr_y(n_z) = sum(yy_wt(:)) / sum(im(:));
end

Cxy_lung = [ ctr_x  ctr_y ];



% parse binary objects for 3d/4d volume
function [info_vol, px2d_vol] = goParseBinaryObjs(cc_vol)
sz = cc_vol.ImageSize;
assert(length(sz) == 3  ||  length(sz) == 4);
if length(sz) == 3;
    sz_y = sz(1);  sz_x = sz(2);  sz_t = 1;      sz_z = sz(3);
else
    sz_y = sz(1);  sz_x = sz(2);  sz_t = sz(3);  sz_z = sz(4);
end

[xx, yy] = meshgrid(1:sz_x, 1:sz_y);
px_areas = cellfun(@numel, cc_vol.PixelIdxList);  % get areas
px_sample = cellfun(@(x) x(1), cc_vol.PixelIdxList);  % just get 1st pixel of every object ...
[~, ~, n_t, n_z] = ind2sub([sz_y, sz_x, sz_t, sz_z], px_sample);  % ... to determine phase and slice
offsets = sz_y*sz_x*sz_t * (n_z - 1) + sz_y*sz_x * (n_t - 1);
info_vol(1:cc_vol.NumObjects,1) = struct('slice', 0, 'phase', 0, 'area', 0, 'Cxy', [0 0]);
px2d_vol = cell(cc_vol.NumObjects, 1);
for n = 1:cc_vol.NumObjects
    info_vol(n).slice = n_z(n);
    info_vol(n).phase = n_t(n);
    info_vol(n).area = px_areas(n);
    px2d = cc_vol.PixelIdxList{n} - offsets(n);
    info_vol(n).Cxy = [ mean(xx(px2d))  mean(yy(px2d)) ];
    px2d_vol{n} = px2d;
end



function chains = goFindChains(scores, info_obj, px2d_obj, cc_obj, obj_filt, res_xy)
sz_z = cc_obj.ImageSize(end);

% create list of unassigned objects, and sort by score, descending
idx_search = 1:length(info_obj);
[~, I] = sort(scores, 'descend');
idx_search = idx_search(I);
L3d_search = labelmatrix(cc_obj);  % labelmatrix for matching objects

% start chaining through all objects, updating the list of unassigned
% objects as we go along
chains = cell(length(info_obj), 1);
n_chn = 1;
% continue until unassigned list is empty, or until none of the remaining
% objects have a score > 0 (i.e., don't search based on score 0 objects)
while ~isempty(idx_search)  &&  any(scores(idx_search) > 0)
    % start new chain search on (latest) highest scoring object
    idx_ref = idx_search(1);
    sl_ref = info_obj(idx_ref).slice;
    arr_chain = false(length(info_obj), 1);  % to keep track of current chain
    arr_chain(idx_ref) = true;
    
    % search up and down the slices
    slices = sl_ref+1 : +1 : sz_z;
    [idx_search, L3d_search, arr_chain] = goFindChains_step(idx_search, L3d_search, arr_chain, idx_ref, slices, info_obj, px2d_obj, obj_filt, res_xy);
    slices = sl_ref-1 : -1 : 1;
    [idx_search, L3d_search, arr_chain] = goFindChains_step(idx_search, L3d_search, arr_chain, idx_ref, slices, info_obj, px2d_obj, obj_filt, res_xy);
    chains(n_chn) = { find(arr_chain) };
    n_chn = n_chn + 1;
end
chains(n_chn:end) = [];



function [idx_search, L3d_search, arr_chain] = goFindChains_step(idx_search, L3d_search, arr_chain, idx_ref, slices, info_obj, px2d_obj, obj_filt, res_xy)
[sz_y, sz_x, ~] = size(L3d_search);
sl_ref = info_obj(idx_ref).slice;
px_ref = px2d_obj{idx_ref};
area_ref = info_obj(idx_ref).area;
area_all = vertcat( info_obj.area );
cxy_ref = info_obj(idx_ref).Cxy;
cxy_all = vertcat( info_obj.Cxy );
binrange = 1:length(info_obj);  % for object-matching histogram

% update list of unassigned objects
idx_search(idx_search == idx_ref) = [];
z_offset = sz_y*sz_x * (sl_ref - 1);
L3d_search(px2d_obj{idx_ref} + z_offset) = 0;

for n_z = slices
    z_offset = sz_y*sz_x * (n_z - 1);
    X = L3d_search(px_ref + z_offset);  % get overlapping pixels
    X = X(X > 0);  % remove zero-valued pixels (no object label)
    if isempty(X);  break;  end  % no overlapping objects, stop search
    
    F = histc(X, binrange);  % get list/labels of overlapping objects
    area_min = min(area_all, area_ref);  % get minimum areas of reference & comparison objects
    F = F(:) ./ area_min;  % normalized area of overlap (force vertical vector or will get error when X is length 1)
    cxy_diff = sqrt( (cxy_all(:,1) - cxy_ref(1)).^2 + (cxy_all(:,2) - cxy_ref(2)).^2 );  % difference in centrepoints
    arr = true(length(F), 1);
    arr = arr  &  F >= 0.6;  % minimum 60% intersection
    arr = arr  &  cxy_diff <= (15 / res_xy);  % maximum 15mm centrepoint difference
    arr = arr  &  ~obj_filt;  % ... and they must not have been filtered out due to other criteria
    if isempty(arr);  break;  end  % no matching objects, stop search
    
    % update list of assigned and unassigned objects
    % allow for multiple objects per slice (WARNING: may be buggy)
    arr_chain(arr) = true;
    px_ref = vertcat( px2d_obj{arr} );
    area_ref = sum([ info_obj(arr).area ]);
    cxy_ref = mean(vertcat( info_obj(arr).Cxy ), 1);
    L3d_search(px_ref + z_offset) = 0;
    for n_obj = find(arr)'
        idx_search(idx_search == n_obj) = [];
    end
end
