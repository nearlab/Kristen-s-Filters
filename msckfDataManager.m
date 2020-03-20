function M = msckfDataManager()
%MSCKFDATAMANAGER Manages data required for implementation of MSCKF
    % Properties
    M.pfPkgList = []; % all the point packages currently being tracked
    M.tagList = []; % tags corresponding to the point packages
    M.poseTimestampList = []; % timestamps corresponding to all the poses in the state vector
    M.nPoses = 0;
    
    % Methods
    M.removePackage = @removePackage;
end


% Old
% methods
%     function obj = addCameraMeasurement(obj, pfTagsObserved, pfTagsNew, cm, I, t_curr, pose)
%         % make new packages
%         nNewTags = length(pfTagsNew)
%         newPkgList = pointPackage.empty;
%         newTagList = zeros(1,nNewTags);
%         for ii = 1:nNewTags
%             tag = pfTagsNew(ii)
%             meas = cm(:,ii)
%             pkg = pointPackage(I,tag,meas,t_curr,pose)
%             newPkgList(ii) = pkg
%             newTagList(ii) = tag
%         end
%         obj.pfPkgList = [obj.pfPkgList, newPkgList];
%         obj.tagList = [obj.tagList, newTagList];
%         % update existing packages
%         oldTagsObserved = setdiff(pfTagsObserved, pfTagsNew)
%         nOldTagsObserved = length(oldTagsObserved)
%         for jj = 1:nOldTagsObserved
%             pkgInd = find(obj.tagList == oldTagsObserved(jj))
%             pkg_curr = obj.pfPkgList(pkgInd)
%             pkg_curr.addObservation(cm(:,jj),t_curr,pose)
%         end    
%     end
%     
% %     function obj = set.pfPkgList(obj,pkgList)
% %         obj.pfPkgList = pkgList;
% %     end
% end
% 
% end
% 
% 
% 

