function jointAngles_lie = getLieModelJAs(mdl_lie)
% "jointAngles_lie" will have following angles, in order:
% back_FB, rshoulder_Elev, rshoulder_Abd, rshjoulder_ExtRot, relbow_flex, 
% lshoulder_Elev, lshoulder_Abd, lshjoulder_ExtRot, lelbow_flex,
% rhip_Flex, rknee_Ext, rankle_Dorsi,
% lhip_Flex, lknee_Ext, lankle_Dorsi,

tr = mdl_lie.transforms;





mdl_lie.transforms(1, 28).frame_in.t(1:3,4);



